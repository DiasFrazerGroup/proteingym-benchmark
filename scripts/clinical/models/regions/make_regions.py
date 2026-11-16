"""Region archives: split one protein's Genebass archive(s) into small archives the unchanged zoo alignment models
can score (eve, evmutation, evmutation2, gemme, rsalor, structuredca).

For a protein (all archives of one UniProt entry, variants pooled as unique mutants):
- Regions: tiles of W residues every W/2 over the wild type (last tile shifted to end at L; one tile [0, L) if
  L <= W). Each variant goes to the tile where it is most central (largest distance to the tile ends; same rule as
  the PLM windows and the structure fragments). Only tiles with variants are written.
- MSA: the archive's a3m is streamed once (never held as Python objects): lowercase insertions dropped, match
  columns [s, e) taken, rows with < `min_cov` non-gap residues in the region dropped, duplicates dropped, query
  first, then a uniform random sample (reservoir, fixed seed) of at most `max_depth` rows. Depths recorded.
- Structure: the AFDB fragment that fully contains the region and in which the region centre is most central
  (fragments.place puts every fragment on the wild type), cropped to the region and renumbered 1..W. If no
  fragment contains the region, the region archive has no structure (rsalor/structuredca cannot score it).
- Assay: the unique mutants assigned to the region, sequences cropped, `mutant` renumbered (region position =
  UniProt position - s). `variants.csv` maps every UniProt mutant to its region and region mutant.

Usage (local env with proteingym-base):
  make_regions.py --uniprot P04114 --window 500 --out $WORK/models/regions/archives
Output: <out>/<UNIPROT>/W<window>/{regions.json, variants.csv, <region>.splits.pgdata}
Documented in docs/clinical/models_long_proteins.md (§3.4).
"""
import argparse
import importlib.util
import io
import json
import random
import shutil
import tempfile
import time
import tomllib
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from biotite.structure.io.pdb import PDBFile

REPO = Path(__file__).resolve().parents[4]
WORK = Path("/no_backup/diasfrazer/egalpern_pgb_clinical")
INDEX = Path("/users/diasfrazer/egalpern/project/pgb_clinical/datasets/clinical/genebass/genebass_fig6b.csv")

# fragments.py of the structure-window models (placement of AFDB fragments on the wild type), reused as is
_spec = importlib.util.spec_from_file_location(
    "fragments", REPO / "models/saprot_window/src/proteingym/models/saprot_window/fragments.py")
fragments = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fragments)

AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H",
       "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V"}
LOWER = bytes(range(ord("a"), ord("z") + 1))
# match states only: keep the 20 aa, X and '-'; any other upper-case letter (B, Z, U, O, J) becomes X
_UPPER_FIX = bytes.maketrans(b"BZUOJ", b"XXXXX")


# ---------------------------------------------------------------------------------------------------------------
# regions
# ---------------------------------------------------------------------------------------------------------------
def tiles(length: int, window: int) -> list[tuple[int, int]]:
    """[a, b) tiles of `window` residues every window//2 over [0, length); one tile if length <= window."""
    if length <= window:
        return [(0, length)]
    step = window // 2
    s = list(range(0, length - window + 1, step))
    if s[-1] != length - window:
        s.append(length - window)
    return [(a, a + window) for a in s]


def assign_tile(positions: list[int], tl: list[tuple[int, int]]) -> int:
    """Index of the tile containing all (0-based) positions in which their centre is most central."""
    lo, hi = min(positions), max(positions)
    c = (lo + hi) / 2
    cand = [(min(c - a, b - 1 - c), -i) for i, (a, b) in enumerate(tl) if a <= lo and hi < b]
    if not cand:
        raise ValueError(f"no tile contains positions {lo + 1}-{hi + 1}")
    return -max(cand)[1]


def parse_mutant(m: str) -> list[tuple[str, int, str]]:
    out = []
    for x in m.split(":"):
        out.append((x[0], int(x[1:-1]), x[-1]))
    return out


# ---------------------------------------------------------------------------------------------------------------
# archive reading
# ---------------------------------------------------------------------------------------------------------------
class Archive:
    """One *.splits.pgdata (outer zip: dataset.pgdata + slices.json), inner dataset extracted to a temp dir."""

    def __init__(self, path: Path, tmp: Path):
        self.path = path
        with zipfile.ZipFile(path) as z:
            z.extract("dataset.pgdata", tmp / path.stem)
        self.inner = zipfile.ZipFile(tmp / path.stem / "dataset.pgdata")
        self.manifest = tomllib.loads(self.inner.read("manifest.lock").decode())
        self.name = self.manifest["name"]
        seq = self.manifest["sequences"][0]
        fa = self.inner.read(seq["path"]).decode().splitlines()
        self.wt = "".join(x.strip() for x in fa if not x.startswith(">"))
        self.uniprot = seq.get("uniprot_id")
        a = self.manifest["assays"][0]
        self.assay = pd.read_csv(io.BytesIO(self.inner.read(a["path"])))
        self.target = a["targets"][0]["name"] if a.get("targets") else self.manifest["assay_targets"][0]["name"]
        self.msa_member = self.manifest["msas"][0]["path"] if self.manifest.get("msas") else None

    def structures(self):
        for s in self.manifest.get("structures", []):
            pdb = PDBFile.read(io.StringIO(self.inner.read(s["path"]).decode()))
            yield s["name"], pdb.get_structure(model=1)

    def a3m_records(self):
        """Yield match-state rows (bytes, upper case, insertions removed) of the MSA, streamed."""
        with self.inner.open(self.msa_member) as fh:
            buf = []
            for line in fh:
                line = line.rstrip(b"\r\n")
                if not line:
                    continue
                if line[:1] == b">":
                    if buf:
                        yield b"".join(buf).translate(_UPPER_FIX, LOWER)
                    buf = []
                else:
                    buf.append(line)
            if buf:
                yield b"".join(buf).translate(_UPPER_FIX, LOWER)


# ---------------------------------------------------------------------------------------------------------------
# MSA per region (one streaming pass for all regions)
# ---------------------------------------------------------------------------------------------------------------
def region_msas(records, wt: str, regions: list[tuple[int, int]], min_cov: float, max_depth: int, seed: int):
    """Per region: query row + reservoir sample of unique rows with >= min_cov non-gap residues in the region."""
    L = len(wt)
    wtb = wt.encode()
    first = next(records)
    if first != wtb:
        raise ValueError("first MSA row (match states) differs from the wild type")
    st = [dict(seen={hash(wtb[a:b])}, res=[], n_unique=0, n_cov=0, rng=random.Random(seed + i),
               min_res=int(np.ceil(min_cov * (b - a))))
          for i, (a, b) in enumerate(regions)]
    n_rows = 1
    bad_len = 0
    for row in records:
        n_rows += 1
        if len(row) != L:
            bad_len += 1
            continue
        for (a, b), s in zip(regions, st):
            seg = row[a:b]
            if (b - a) - seg.count(b"-") < s["min_res"]:
                continue
            s["n_cov"] += 1
            h = hash(seg)
            if h in s["seen"]:
                continue
            s["seen"].add(h)
            s["n_unique"] += 1
            if len(s["res"]) < max_depth:
                s["res"].append(seg)
            else:
                j = s["rng"].randrange(s["n_unique"])
                if j < max_depth:
                    s["res"][j] = seg
    if bad_len:
        raise ValueError(f"{bad_len} MSA rows have a match-state length != {L}")
    out = []
    for (a, b), s in zip(regions, st):
        rows = [wtb[a:b]] + s["res"]
        out.append(dict(rows=rows, depth_total=n_rows, depth_cov=s["n_cov"], depth_unique=s["n_unique"] + 1,
                        depth_kept=len(rows)))
    return out


def gap_excluded_columns(rows: list[bytes], max_gap_fraction: float = 0.5) -> list[int]:
    """EVmutation's rule (evedesign EVcouplings._build_alignment): columns with unweighted gap fraction
    > max_gap_fraction over all rows (query included) are excluded from the model. 1-based region positions."""
    m = np.frombuffer(b"".join(rows), dtype=np.uint8).reshape(len(rows), -1)
    frac = (m == ord("-")).mean(axis=0)
    return [int(i) + 1 for i in np.where(frac > max_gap_fraction)[0]]


# ---------------------------------------------------------------------------------------------------------------
# structures
# ---------------------------------------------------------------------------------------------------------------
def place_structures(arch: Archive, wt: str):
    out = []
    for name, atoms in arch.structures():
        atoms = atoms[~atoms.hetero]
        ca = atoms[atoms.atom_name == "CA"]
        letters = "".join(AA3.get(r, "X") for r in ca.res_name)
        p = fragments.place(name, np.asarray(ca.res_id, dtype=int), letters, wt)
        out.append((p, atoms))
    return out


def region_structure(placed, wt: str, a: int, b: int):
    """Most-central fragment fully containing [a, b), cropped and renumbered 1..(b-a); None if none contains it."""
    c = (a + b - 1) / 2
    cand = [(min(c - p.start, p.end - c), i) for i, (p, _) in enumerate(placed) if p.start <= a and b - 1 <= p.end]
    if not cand:
        return None
    _, i = max(cand)
    p, atoms = placed[i]
    wt_idx = atoms.res_id + p.shift
    crop = atoms[(wt_idx >= a) & (wt_idx < b)].copy()
    crop.res_id = crop.res_id + p.shift - a + 1
    crop.chain_id[:] = "A"
    ca = crop[crop.atom_name == "CA"]
    mism = [f"{wt[a + r - 1]}{r}{AA3.get(n, 'X')}" for r, n in zip(ca.res_id, ca.res_name)
            if AA3.get(n, "X") != wt[a + r - 1]]
    missing = (b - a) - len(set(ca.res_id.tolist()))
    return dict(name=p.name, atoms=crop, mismatches=mism, missing=missing, centrality=float(min(c - p.start, p.end - c)))


# ---------------------------------------------------------------------------------------------------------------
# writing one region archive with proteingym-base
# ---------------------------------------------------------------------------------------------------------------
def toml_str(s: str) -> str:
    return json.dumps(s)


def write_region_archive(d: Path, out: Path, name: str, uniprot: str, wt_region: str, assay: pd.DataFrame,
                         target: str, rows: list[bytes], struct: dict | None, meta: dict) -> Path:
    from proteingym.base import Dataset, Manifest
    from proteingym.base.dataset import Subsets

    (d / "assays").mkdir(parents=True)
    assay.to_csv(d / "assays" / f"{name}.csv", index=False)
    (d / f"{uniprot}.fasta").write_text(f">{uniprot}\n{wt_region}\n")
    with open(d / "msa.a3m", "wb") as f:
        f.write(f">{uniprot}\n".encode() + rows[0] + b"\n")
        for i, r in enumerate(rows[1:], 1):
            f.write(b">s%d\n" % i + r + b"\n")
    non_targets = [c for c in assay.columns if c not in ("sequence", target)]
    m = [f"name = {toml_str(name)}", f"description = {toml_str(meta['description'])}", 'version = "0.1.0"', "",
         "[[assay_targets]]", f'name = "{target}"', "",
         "[[assays]]", f"name = {toml_str(name)}", f'path = "assays/{name}.csv"',
         'sequence_alphabet = "AA"', 'target_phenotype = "other"', f"number_of_variants = {len(assay)}", "",
         "[[assays.targets]]", f'name = "{target}"', "",
         *[x for c in non_targets for x in ("[[assays.non_targets]]", f'name = "{c}"', "")],
         "[[sequences]]", 'type = "wild_type"', 'alphabet = "AA"', f'path = "{uniprot}.fasta"',
         f'uniprot_id = "{uniprot}"', ""]
    if struct is not None:
        pdb = PDBFile()
        pdb.set_structure(struct["atoms"])
        pdb.write(d / "structure.pdb")
        m += ["[[structures]]", 'path = "structure.pdb"', f'name = "{struct["name"]}_region"',
              f"description = {toml_str(meta['structure_description'])}", ""]
    m += ["[[msas]]", 'path = "msa.a3m"', f'name = "{name}_msa"', 'format = "a3m"',
          f'reference_sequence_name = "{uniprot}"', f"description = {toml_str(meta['msa_description'])}", ""]
    (d / "manifest.toml").write_text("\n".join(m))
    ds = Dataset.from_manifest(Manifest.from_path(d / "manifest.toml"))
    Subsets(dataset=ds, slices={}).dump(path=out)
    return out / f"{name}.splits.pgdata"


# ---------------------------------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--uniprot", required=True)
    ap.add_argument("--window", type=int, required=True, help="region size W (>= L: one full-length region)")
    ap.add_argument("--out", type=Path, default=WORK / "models/regions/archives")
    ap.add_argument("--archive-dir", type=Path, default=WORK / "genebass/splits_layers")
    ap.add_argument("--index", type=Path, default=INDEX)
    ap.add_argument("--min-cov", type=float, default=0.5, help="min non-gap fraction of a row in the region")
    ap.add_argument("--max-depth", type=int, default=50000, help="random depth cap (query included)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--tmp", type=Path, default=None)
    a = ap.parse_args()

    t0 = time.time()
    idx = pd.read_csv(a.index)
    names = idx.loc[idx.uniprot_id == a.uniprot, "dataset"].tolist()
    if not names:
        raise SystemExit(f"{a.uniprot} not in {a.index}")
    out = a.out / a.uniprot / f"W{a.window}"
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    with tempfile.TemporaryDirectory(dir=a.tmp) as tmp:
        tmp = Path(tmp)
        archs = [Archive(a.archive_dir / f"{n}.splits.pgdata", tmp) for n in names]
        wt = archs[0].wt
        assert all(x.wt == wt for x in archs), "archives of one protein disagree on the wild type"
        target = archs[0].target
        L = len(wt)
        # unique mutants (pooled over the protein's archives); target value of the first occurrence (unused by
        # zero-shot models, needed by the archive format)
        rows = pd.concat([x.assay[["mutant", "sequence", target]].assign(dataset=x.name) for x in archs])
        mut = rows.drop_duplicates("mutant").reset_index(drop=True)
        tl = tiles(L, a.window)
        recs = []
        for m in mut.mutant:
            sites = parse_mutant(m)
            for w, p, _ in sites:
                assert wt[p - 1] == w, f"{m}: wild type at {p} is {wt[p - 1]}"
            recs.append(assign_tile([p - 1 for _, p, _ in sites], tl))
        mut["tile"] = recs
        used = sorted(mut.tile.unique())
        regions = [tl[i] for i in used]
        print(f"{a.uniprot} L={L} W={a.window}: {len(mut)} unique mutants from {len(archs)} archives, "
              f"{len(tl)} tiles, {len(regions)} with variants", flush=True)

        t1 = time.time()
        msa_arch = next(x for x in archs if x.msa_member)
        msas = region_msas(msa_arch.a3m_records(), wt, regions, a.min_cov, a.max_depth, a.seed)
        t_msa = time.time() - t1
        print(f"MSA streamed in {t_msa:.0f}s: {msas[0]['depth_total']} rows", flush=True)

        placed = place_structures(archs[0], wt)
        info, var = [], []
        for k, ((s, e), ti, ms) in enumerate(zip(regions, used, msas)):
            rid = f"r{s + 1:05d}-{e:05d}"
            name = f"{a.uniprot}_{rid}"
            sub = mut[mut.tile == ti].copy()
            st = region_structure(placed, wt, s, e)

            def renum(m):
                return ":".join(f"{w}{p - s}{t}" for w, p, t in parse_mutant(m))

            assay = pd.DataFrame({"sequence": sub.sequence.str[s:e].to_numpy(),
                                  target: sub[target].to_numpy(),
                                  "mutant": [renum(m) for m in sub.mutant],
                                  "uniprot_mutant": sub.mutant.to_numpy()})
            assert (assay.sequence.str.len() == e - s).all()
            excl = gap_excluded_columns(ms["rows"])
            meta = dict(
                description=f"Region {s + 1}-{e} of {a.uniprot} (L={L}, W={a.window}) for alignment models",
                structure_description=(f"{st['name']} cropped to wild-type residues {s + 1}-{e}, renumbered 1-{e - s}"
                                       if st else ""),
                msa_description=(f"match columns {s + 1}-{e} of {msa_arch.msa_member}; rows with >= {a.min_cov:.0%} "
                                 f"non-gap in the region, deduplicated, random cap {a.max_depth} (seed {a.seed})"))
            path = write_region_archive(tmp / "build" / name, out, name, a.uniprot, wt[s:e], assay, target,
                                        ms["rows"], st, meta)
            info.append(dict(region=rid, archive=path.name, start=s + 1, end=e, n_mutants=len(sub),
                             depth_total=ms["depth_total"], depth_cov=ms["depth_cov"],
                             depth_unique=ms["depth_unique"], depth_kept=ms["depth_kept"],
                             evmutation_gap_excluded=excl,
                             structure=st["name"] if st else None,
                             structure_centrality=st["centrality"] if st else None,
                             structure_mismatches=st["mismatches"] if st else None,
                             structure_missing=st["missing"] if st else None))
            for m, rm in zip(sub.mutant, assay.mutant):
                var.append(dict(mutant=m, region=rid, region_start=s + 1, region_end=e, region_mutant=rm))
            print(f"  {rid}: {len(sub)} mutants, depth {ms['depth_cov']} cov / {ms['depth_unique']} unique -> "
                  f"{ms['depth_kept']} kept, {len(excl)} gap-excluded cols, structure "
                  f"{st['name'] if st else 'NONE'}", flush=True)
        pd.DataFrame(var).to_csv(out / "variants.csv", index=False)
        rows.drop(columns=["sequence"]).to_csv(out / "rows.csv", index=False)
        (out / "regions.json").write_text(json.dumps(dict(
            uniprot=a.uniprot, length=L, window=a.window, min_cov=a.min_cov, max_depth=a.max_depth, seed=a.seed,
            archives=names, target=target, n_mutants=len(mut), t_msa_s=round(t_msa, 1),
            t_total_s=round(time.time() - t0, 1), regions=info), indent=1))
    print(f"done in {time.time() - t0:.0f}s -> {out}", flush=True)


if __name__ == "__main__":
    main()
