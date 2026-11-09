"""Add annotation layers to ProteinGym dataset archives (.splits.pgdata / .pgarchive), for any benchmark.

Layers are extra non-target columns on every assay row, for comparison and for masking an analysis. They are never
model inputs or targets.

  scores  (external predictors, as published, per variant; empty for multi-mutants or when not covered)
    score_alphamissense   AlphaMissense pathogenicity, 0-1, higher = more pathogenic           [dbNSFP]
    score_popeve          popEVE, more negative = more deleterious                             [dbNSFP]
    score_cadd_phred      CADD PHRED-scaled, higher = more deleterious                         [dbNSFP]
    score_primateai3d     PrimateAI-3D, higher = more pathogenic        [Illumina file; only if configured]
  masks   (per row)
    mask_plddt            AlphaFold pLDDT of the mutated residue in the archive's own structure(s) (the most
                          central fragment covering it); multi-mutants: minimum over the mutated sites
    mask_protein_length   length of the full canonical UniProt protein (not the archive's sequence segment)
  helpers
    uniprot_pos           mutated position(s) in canonical UniProt numbering (':'-joined for multi-mutants)
    score_snv_source      how SNVs were found: genomic (the row's own GRCh38 SNV) | protein (all SNVs causing the
                          protein change, via the AlphaMissense hg38 table; scores averaged) | none

Per row, the variant is read from the `mutant` column (e.g. A12V or A12V:K40E, numbered on the archive's wild type),
or derived from the mutated sequence. The archive wild type is placed on the canonical UniProt sequence (exact,
substring, or global alignment), so DMS archives on a domain still get UniProt positions and the full length.
If the row has a GRCh38 SNV column (default `variant_id`, e.g. 2-21041028-G-A), external scores use that exact SNV.

dbNSFP fields are per transcript; the transcript is chosen by UniProt accession and protein position (then by
accession and amino acids, e.g. Ensembl/UniProt numbering offsets), and is recorded only through score_snv_source.

Rows and their order are unchanged, so the archive's splits (slices.json) stay valid; the result is reloaded with
proteingym-base as a check. Sources and versions: sources.toml (also written into each column's description).

Run: python add_layers.py ARCHIVE... --out DIR [--layers scores,masks] [--snv-col variant_id] [--force]
"""
import argparse
import gzip
import io
import re
import subprocess
import sys
import tempfile
import tomllib
import zipfile
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.dataset as pads
import toml
from Bio.Align import PairwiseAligner

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from build_genebass import read_sprot  # noqa: E402

AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H",
       "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V", "MSE": "M"}
SNV_RE = re.compile(r"^(?:chr)?([0-9]{1,2}|X|Y|MT?)[-:_](\d+)[-:_]([ACGT])[-:_>/]([ACGT])$", re.I)
MUT_RE = re.compile(r"^([A-Z])(\d+)([A-Z])$")


# --- sequences and positions -------------------------------------------------------------------------------------
def position_map(query, ref, min_identity=0.9):
    """1-based query position -> 1-based ref position (None if unaligned), via identity, substring or alignment."""
    if query == ref:
        return {i: i for i in range(1, len(query) + 1)}, 1.0
    k = ref.find(query)
    if k >= 0:
        return {i: i + k for i in range(1, len(query) + 1)}, 1.0
    al = PairwiseAligner(mode="local", match_score=2, mismatch_score=-1, open_gap_score=-5, extend_gap_score=-0.5)
    a = al.align(query, ref)[0]
    m, same = {}, 0
    for (qs, qe), (rs, re_) in zip(*a.aligned):
        for d in range(qe - qs):
            m[qs + d + 1] = rs + d + 1
            same += query[qs + d] == ref[rs + d]
    ident = same / max(len(m), 1)
    return (m, ident) if ident >= min_identity else ({}, ident)


def parse_mutant(m):
    """'A12V' or 'A12V:K40E' -> [(wt, pos, alt), ...]; None for synonymous/WT/unparseable."""
    if not isinstance(m, str) or not m or m.upper() in {"WT", "_WT", "WILDTYPE"}:
        return None
    out = []
    for x in m.split(":"):
        g = MUT_RE.match(x.strip())
        if not g:
            return None
        out.append((g[1], int(g[2]), g[3]))
    return out


def mutant_from_sequence(seq, wt):
    if not isinstance(seq, str) or len(seq) != len(wt):
        return None
    return ":".join(f"{a}{i + 1}{b}" for i, (a, b) in enumerate(zip(wt, seq)) if a != b) or None


# --- structures ----------------------------------------------------------------------------------------------------
def read_ca(pdb_bytes):
    seq, plddt = [], []
    for line in pdb_bytes.decode().splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA" and line[16] in " A":
            seq.append(AA3.get(line[17:20], "X"))
            plddt.append(float(line[60:66]))
    return "".join(seq), np.array(plddt)


def site_plddt(structs, wt):  # noqa: C901
    """wt position -> pLDDT from the most central structure covering it (structs: [(seq, plddt, span_start)])."""
    best = {}
    for s, pl, start in structs:
        if start and wt[start - 1:start - 1 + len(s)] == s:
            m = {start + i: i + 1 for i in range(len(s))}          # wt pos -> structure index (1-based)
        else:
            fwd, ident = position_map(s, wt)
            m = {w: q for q, w in fwd.items()}
        if not m:
            continue
        lo, hi = min(m), max(m)
        for w, q in m.items():
            centrality = min(w - lo, hi - w)
            if w not in best or centrality > best[w][0]:
                best[w] = (centrality, pl[q - 1])
    return {w: v for w, (_, v) in best.items()}


# --- external scores ----------------------------------------------------------------------------------------------
class SnvScores:
    """dbNSFP (AlphaMissense, popEVE, CADD) and optional PrimateAI-3D for GRCh38 SNVs, resolved per UniProt protein."""

    def __init__(self, cfg):
        self.cfg = cfg
        self.tabix = cfg["tabix"]["path"]
        self.db = cfg["dbnsfp"]["path"]
        with gzip.open(self.db, "rt") as f:
            self.hdr = {c: i for i, c in enumerate(f.readline().rstrip("\n").split("\t"))}
        need = ["#chr", "pos(1-based)", "ref", "alt", "aaref", "aaalt", "aapos", "Uniprot_acc",
                "AlphaMissense_score", "popEVE_score", "CADD_phred"]
        missing = [c for c in need if c not in self.hdr]
        assert not missing, f"dbNSFP columns missing: {missing}"
        self.pai = cfg["primateai3d"] if cfg["primateai3d"].get("path") else None

    def dbnsfp(self, snvs):
        """{(chrom,pos,ref,alt)} -> DataFrame of raw dbNSFP rows (strings)."""
        if not snvs:
            return pd.DataFrame()
        with tempfile.NamedTemporaryFile("w", suffix=".tsv") as reg:
            for c, p in sorted({(s[0], s[1]) for s in snvs}, key=lambda x: (x[0], x[1])):
                reg.write(f"{c}\t{p}\t{p}\n")
            reg.flush()
            out = subprocess.run([self.tabix, "-R", reg.name, self.db], capture_output=True, text=True, check=True)
        h = self.hdr
        rows = [r.split("\t") for r in out.stdout.splitlines()]
        rows = [r for r in rows if (r[h["#chr"]], int(r[h["pos(1-based)"]]), r[h["ref"]], r[h["alt"]]) in snvs]
        cols = ["#chr", "pos(1-based)", "ref", "alt", "aaref", "aaalt", "aapos", "Uniprot_acc",
                "AlphaMissense_score", "popEVE_score", "CADD_phred"]
        d = pd.DataFrame([[r[h[c]] for c in cols] for r in rows], columns=cols).drop_duplicates()
        d["key"] = list(zip(d["#chr"], d["pos(1-based)"].astype(int), d.ref, d.alt))
        return d

    @staticmethod
    def resolve(row, uniprot, upos, wt_aa, alt_aa):
        """Pick the transcript entry of a dbNSFP row for this protein; returns (am, popeve, cadd) floats."""
        def split(x):
            return x.split(";")
        acc, apos = split(row.Uniprot_acc), split(row.aapos)
        n = max(len(acc), len(apos))
        get = (lambda v, i: v[i] if len(v) == n else v[0])
        aa_ok = row.aaref == wt_aa and row.aaalt == alt_aa
        idx = [i for i in range(n) if get(acc, i).split("-")[0] == uniprot and get(apos, i) == str(upos)]
        if not idx and aa_ok:
            idx = [i for i in range(n) if get(acc, i).split("-")[0] == uniprot] or \
                  [i for i in range(n) if get(apos, i) == str(upos)]
        if not idx:
            return np.nan, np.nan, np.nan

        def val(field):
            v = split(getattr(row, field))
            x = get(v, idx[0])
            return float(x) if x not in {".", ""} else np.nan

        return val("AlphaMissense_score"), val("popEVE_score"), \
            (float(row.CADD_phred) if row.CADD_phred not in {".", ""} else np.nan)

    def primateai3d(self, snvs):
        if not self.pai or not snvs:
            return {}
        c = self.pai
        want = pd.DataFrame(list(snvs), columns=["chrom", "pos", "ref", "alt"])
        hits = []
        for ch in pd.read_csv(c["path"], sep=None, engine="python", chunksize=2_000_000,
                              usecols=[c["chrom"], c["pos"], c["ref"], c["alt"], c["score"]]):
            ch = ch.rename(columns={c["chrom"]: "chrom", c["pos"]: "pos", c["ref"]: "ref", c["alt"]: "alt",
                                    c["score"]: "s"})
            ch["chrom"] = ch.chrom.astype(str).str.removeprefix("chr")
            hits.append(ch.merge(want, on=["chrom", "pos", "ref", "alt"]))
        h = pd.concat(hits)
        return {(r.chrom, r.pos, r.ref, r.alt): r.s for r in h.itertuples()}


@lru_cache(maxsize=64)
def snv_index_for(path, uniprot):
    t = pads.dataset(path).to_table(filter=pads.field("uniprot_id") == uniprot).to_pandas()
    return t.groupby("protein_variant")[["chrom", "pos", "ref", "alt"]].apply(
        lambda g: [tuple(x) for x in g.itertuples(index=False)]).to_dict()


# --- archive I/O -------------------------------------------------------------------------------------------------
def read_archive(path):
    with zipfile.ZipFile(path) as z:
        outer = {i.filename: (z.read(i), i.compress_type) for i in z.infolist()}
    with zipfile.ZipFile(io.BytesIO(outer["dataset.pgdata"][0])) as z:
        inner = {i.filename: (z.read(i), i.compress_type) for i in z.infolist()}
    return outer, inner


def zip_bytes(members):
    b = io.BytesIO()
    with zipfile.ZipFile(b, "w") as z:
        for name, (data, ct) in members.items():
            z.writestr(name, data, compress_type=ct or zipfile.ZIP_DEFLATED)
    return b.getvalue()


def fmt(x):
    return "" if x is None or (isinstance(x, float) and np.isnan(x)) else (
        str(int(x)) if isinstance(x, (int, np.integer)) else f"{x:.6g}" if isinstance(x, float) else str(x))


# --- main ---------------------------------------------------------------------------------------------------------
def annotate(path, out_dir, cfg, sprot, scorer, layers, snv_col, mutant_col, force):
    outer, inner = read_archive(path)
    man = toml.loads(inner["manifest.lock"][0].decode())
    wt_sec = next(s for s in man["sequences"] if s.get("type", "wild_type") == "wild_type")
    wt = "".join(inner[wt_sec["path"]][0].decode().split("\n", 1)[1].split())
    uniprot = (wt_sec.get("uniprot_id") or "").split("-")[0]
    canon = sprot.get(uniprot, {}).get("sequence")
    pmap, ident = position_map(wt, canon) if canon else ({}, 0.0)
    plddt, plddt_src = {}, "none"
    if "masks" in layers:
        # proteingym-base rewrites structures when it builds an archive and zeroes the B-factor (pLDDT) column, so
        # archive structures are used only if they still carry pLDDT; otherwise the AlphaFold DB v6 models of the
        # UniProt entry (all fragments), placed on the archive wild type by alignment.
        structs = [(*read_ca(inner[s["path"]][0]), int(s.get("metadata", {}).get("span_start", 0) or 0))
                   for s in man.get("structures", [])]
        if structs and any(pl.max() > 0 for _, pl, _ in structs):
            plddt, plddt_src = site_plddt(structs, wt), "archive"
        elif uniprot:
            files = sorted(Path(cfg["plddt"]["afdb_dir"]).glob(f"AF-{uniprot}-F*-model_v6.pdb"))
            st = [(*read_ca(f.read_bytes()), 0) for f in files]
            if st:
                plddt, plddt_src = site_plddt(st, wt), "afdb_v6"
    v = cfg["dbnsfp"]["version"]
    desc = {
        "uniprot_pos": f"mutated position(s) on the canonical UniProt sequence ({cfg['uniprot']['version']})",
        "score_alphamissense": f"AlphaMissense pathogenicity (0-1, higher = more pathogenic); {v}",
        "score_popeve": f"popEVE score (more negative = more deleterious); {v}",
        "score_cadd_phred": f"CADD PHRED (higher = more deleterious); {v}",
        "score_primateai3d": f"PrimateAI-3D score (higher = more pathogenic); {cfg['primateai3d']['version']}",
        "score_snv_source": "SNVs used for the scores: genomic (row's own GRCh38 SNV), protein (all SNVs causing "
                            "the protein change, AlphaMissense hg38 table; scores averaged), none",
        "mask_plddt": "AlphaFold pLDDT of the mutated residue in this archive's structure(s), most central "
                      "fragment; multi-mutants: minimum over sites",
        "mask_plddt_source": "pLDDT source: archive (its own structures) | afdb_v6 (AlphaFold DB v6 models of the "
                             "UniProt entry; archive structures had no pLDDT) | none",
        "mask_protein_length": f"length of the full canonical UniProt protein ({cfg['uniprot']['version']}), "
                               "not of the archive sequence",
    }
    summary = []
    for asec in man["assays"]:
        raw = inner[asec["path"]][0].decode()
        df = pd.read_csv(io.StringIO(raw), dtype=str, keep_default_na=False)
        seq_col = asec.get("sequence_alias") or "sequence"
        muts = (df[mutant_col] if mutant_col in df else
                pd.Series([mutant_from_sequence(s, wt) for s in df[seq_col]], index=df.index))
        parsed = [parse_mutant(m) for m in muts]
        new = {"uniprot_pos": [":".join(str(pmap.get(p, "")) for _, p, _ in m) if m else "" for m in parsed]}
        if "masks" in layers:
            new["mask_plddt"] = [min((plddt.get(p, np.nan) for _, p, _ in m), default=np.nan) if m else np.nan
                                 for m in parsed]
            new["mask_plddt_source"] = [plddt_src if m else "none" for m in parsed]
            new["mask_protein_length"] = [len(canon) if canon else np.nan] * len(df)
        if "scores" in layers:
            single = [m[0] if m and len(m) == 1 else None for m in parsed]
            snvs_of, src = [], []
            idx = snv_index_for(cfg["snv_index"]["path"], uniprot) if uniprot else {}
            for i, s in enumerate(single):
                g = SNV_RE.match(df[snv_col].iat[i]) if snv_col in df and s else None
                if g:
                    snvs_of.append([(g[1].upper().removeprefix("CHR"), int(g[2]), g[3].upper(), g[4].upper())])
                    src.append("genomic")
                elif s and pmap.get(s[1]):
                    hits = idx.get(f"{s[0]}{pmap[s[1]]}{s[2]}", [])
                    snvs_of.append(hits)
                    src.append("protein" if hits else "none")
                else:
                    snvs_of.append([])
                    src.append("none")
            allsnv = {x for xs in snvs_of for x in xs}
            d = scorer.dbnsfp(allsnv)
            byk = {k: g for k, g in d.groupby("key")} if len(d) else {}
            pai = scorer.primateai3d(allsnv)
            am, pe, cadd, p3 = [], [], [], []
            for s, xs in zip(single, snvs_of):
                vals = [SnvScores.resolve(r, uniprot, pmap.get(s[1]), s[0], s[2])
                        for x in xs for r in byk.get(x, pd.DataFrame()).itertuples()] if s else []
                arr = np.array(vals, dtype=float).reshape(-1, 3)
                with np.errstate(all="ignore"), __import__("warnings").catch_warnings():
                    __import__("warnings").simplefilter("ignore", RuntimeWarning)   # all-NaN columns -> NaN
                    m = np.nanmean(arr, axis=0) if len(arr) else np.full(3, np.nan)
                am.append(m[0]); pe.append(m[1]); cadd.append(m[2])
                pv = [pai[x] for x in xs if x in pai]
                p3.append(float(np.mean(pv)) if pv else np.nan)
            new.update(score_alphamissense=am, score_popeve=pe, score_cadd_phred=cadd, score_snv_source=src)
            if scorer.pai:
                new["score_primateai3d"] = p3
        clash = [c for c in new if c in df.columns]
        if clash and not force:
            raise SystemExit(f"{path.name}: columns already present {clash} (use --force to replace)")
        for c, vals in new.items():
            df[c] = [fmt(x) for x in vals]
        inner[asec["path"]] = (df.to_csv(index=False).encode(), inner[asec["path"]][1])
        nt = [x for x in asec.get("non_targets", []) if x["name"] not in new]
        asec["non_targets"] = nt + [{"name": c, "description": desc[c]} for c in new]
        cov = {c: float(np.mean([x != "" for x in df[c]])) for c in new}
        summary.append(dict(archive=path.name, assay=asec["name"], rows=len(df), uniprot=uniprot,
                            wt_to_uniprot_identity=round(ident, 4), **{f"cov_{c}": round(x, 4) for c, x in cov.items()}))
    inner["manifest.lock"] = (toml.dumps(man).encode(), inner["manifest.lock"][1])
    outer["dataset.pgdata"] = (zip_bytes(inner), outer["dataset.pgdata"][1])
    out = out_dir / path.name
    out.write_bytes(zip_bytes(outer))
    # check: proteingym-base reads it, same rows, new columns present
    from proteingym.base import Subsets
    ds = Subsets.from_path(out).dataset
    for a_, s_ in zip(ds.assays, summary):
        names = {f.name for f in a_.fields}
        assert len(a_.records) == s_["rows"], (out, len(a_.records), s_["rows"])
        assert {"uniprot_pos"} <= names, names
    return summary


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("archives", nargs="+", type=Path)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--sources", type=Path, default=Path(__file__).with_name("sources.toml"))
    ap.add_argument("--layers", default="scores,masks")
    ap.add_argument("--snv-col", default="variant_id", help="GRCh38 SNV column, if the assay has one")
    ap.add_argument("--mutant-col", default="mutant")
    ap.add_argument("--force", action="store_true", help="replace layer columns already present")
    a = ap.parse_args()
    cfg = tomllib.loads(a.sources.read_text())
    layers = set(a.layers.split(","))
    a.out.mkdir(parents=True, exist_ok=True)
    assert all(p.resolve().parent != a.out.resolve() for p in a.archives), "--out must differ from the input folder"
    sprot, _ = read_sprot(cfg["uniprot"]["path"])
    scorer = SnvScores(cfg) if "scores" in layers else None
    rows = []
    for p in a.archives:
        rows += annotate(p, a.out, cfg, sprot, scorer, layers, a.snv_col, a.mutant_col, a.force)
        print("\t".join(f"{k}={v}" for k, v in rows[-1].items()), flush=True)
    pd.DataFrame(rows).to_csv(a.out / "layers_summary.tsv", sep="\t", index=False)
    print(f"done: {len(a.archives)} archives -> {a.out}")


if __name__ == "__main__":
    main()
