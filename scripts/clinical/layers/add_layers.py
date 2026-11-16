"""Add annotation layers to ProteinGym dataset archives (.splits.pgdata / .pgarchive), for any benchmark.

Layers are extra non-target columns on every assay row, for comparison and for masking an analysis. They are never
model inputs or targets.

  scores  (external predictors, as published, per variant; empty for multi-mutants or when not covered)
    score_alphamissense   AlphaMissense pathogenicity, 0-1, higher = more pathogenic  [AlphaMissense hg38 table]
    score_alphamissense_source  snv (UniProt + the row's GRCh38 SNV) | protein (UniProt + protein change) | none
    score_popeve          popEVE, more negative = more deleterious               [popEVE v1.1 zip, per RefSeq protein]
    score_cadd_phred      CADD PHRED-scaled, higher = more deleterious                         [dbNSFP]
    score_primateai3d     PrimateAI-3D, higher = more pathogenic        [Illumina file; only if configured]
  masks   (per row)
    mask_plddt            AlphaFold pLDDT of the mutated residue in the archive's own structure(s) (the most
                          central fragment covering it); multi-mutants: minimum over the mutated sites
    mask_protein_length   length of the full canonical UniProt protein (not the archive's sequence segment)
  helpers
    uniprot_pos           mutated position(s) in canonical UniProt numbering (':'-joined for multi-mutants)
    score_snv_source      SNVs used for CADD / PrimateAI-3D: genomic (the row's own GRCh38 SNV) | protein (all SNVs
                          causing the protein change, via the AlphaMissense hg38 table; scores averaged) | none

Per row, the variant is read from the `mutant` column (e.g. A12V or A12V:K40E, numbered on the archive's wild type),
or derived from the mutated sequence. The archive wild type is placed on the canonical UniProt sequence (exact,
substring, or global alignment), so DMS archives on a domain still get UniProt positions and the full length.
If the row has a GRCh38 SNV column (default `variant_id`, e.g. 2-21041028-G-A), SNV-level scores use that exact SNV.

Each score comes from its own complete release (paths and versions in sources.toml):
  AlphaMissense  the official hg38 table (canonical isoforms), as a parquet sorted by UniProt (build_snv_index.py);
                 looked up by (UniProt, SNV), else by (UniProt, protein change). Its sequence is rebuilt from the
                 wild-type letters of protein_variant and placed on the current UniProt sequence by alignment if it
                 differs.
  popEVE         the v1.1 zip, one CSV per RefSeq protein, read in place. UniProt -> RefSeq from Swiss-Prot DR
                 lines (displayed isoform first); the RefSeq sequence (from `mutant`) is checked against UniProt and
                 positions mapped by alignment if it differs. The RefSeq used is in the column description and in
                 layers_summary.tsv.
  CADD           dbNSFP (SNV-level).

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
def sequence_from_mutants(mutants):
    """Wild-type sequence implied by a list of 'A12V' substitutions ('X' at positions without any)."""
    wt = {}
    for m in mutants:
        g = MUT_RE.match(m)
        if g:
            wt[int(g[2])] = g[1]
    return "".join(wt.get(i, "X") for i in range(1, max(wt, default=0) + 1))


def place(canon, other):
    """canonical UniProt position -> position on `other` (a source's own sequence); plus a short note."""
    if not canon or not other:
        return {}, "no sequence"
    if len(canon) == len(other) and all(a == b or b == "X" for a, b in zip(canon, other)):
        return {i: i for i in range(1, len(canon) + 1)}, "identical to the UniProt canonical sequence"
    m, ident = position_map(canon, other)
    diff = sum(canon[i - 1] != other[j - 1] for i, j in m.items())
    return m, (f"differs from the UniProt canonical sequence (lengths {len(canon)}/{len(other)}); positions mapped "
               f"by alignment: {len(m)}/{len(canon)} canonical positions placed, {diff} mismatched, identity "
               f"{ident:.3f}")


class AlphaMissense:
    """AlphaMissense hg38 table (parquet from build_snv_index.py), per UniProt accession: score by (UniProt, GRCh38
    SNV), else by protein change (mean over its SNVs; AlphaMissense is protein-level so they agree). Also the SNV
    index (protein change -> causing SNVs) for SNV-level scores."""

    def __init__(self, cfg):
        self.path = cfg["alphamissense"]["table"]
        self.cache = {}

    def protein(self, uniprot, canon):
        if uniprot not in self.cache:
            t = pads.dataset(self.path).to_table(filter=pads.field("uniprot_id") == uniprot).to_pandas()
            seq = sequence_from_mutants(t.protein_variant.unique())
            pos_map, note = place(canon, seq)
            snv = dict(zip(zip(t.chrom, t.pos, t.ref, t.alt), zip(t.protein_variant, t.am_pathogenicity)))
            by_pv = {pv: (float(g.am_pathogenicity.mean()), list(zip(g.chrom, g.pos, g.ref, g.alt)))
                     for pv, g in t.groupby("protein_variant")}
            self.cache = {uniprot: dict(seq=seq, map=pos_map, note=note, snv=snv, pv=by_pv)}   # one protein at a time
        return self.cache[uniprot]

    @staticmethod
    def key(p, wt_aa, upos, alt_aa):
        """protein_variant on AlphaMissense numbering for a change at canonical position upos (None if not placed)."""
        q = p["map"].get(upos) if upos else None
        return f"{wt_aa}{q}{alt_aa}" if q and p["seq"][q - 1] == wt_aa else None


class PopEVE:
    """popEVE per RefSeq protein, read straight from the release zip; UniProt -> RefSeq via Swiss-Prot DR lines."""

    def __init__(self, cfg, refseq):
        c = cfg["popeve"]
        self.zip = zipfile.ZipFile(c["zip"])
        self.score = c["score"]
        self.members = {Path(n).stem.split(".")[0]: n for n in self.zip.namelist() if n.endswith(".csv")}
        self.refseq = refseq
        self.cache = {}

    def protein(self, uniprot, canon):
        """-> dict(np, seq, map, note, scores{mutant on RefSeq numbering: score}) or None if no RefSeq in popEVE."""
        if uniprot in self.cache:
            return self.cache[uniprot]
        best = None
        for np_acc, displayed in self.refseq.get(uniprot, []):
            name = self.members.get(np_acc.split(".")[0])
            if not name:
                continue
            d = pd.read_csv(self.zip.open(name), usecols=["mutant", self.score])
            seq = sequence_from_mutants(d.mutant)
            pos_map, note = place(canon, seq)
            rank = (note.startswith("identical"), displayed, len(pos_map))
            if best is None or rank > best[0]:
                used = Path(name).stem
                if used != np_acc:
                    note += f"; Swiss-Prot lists {np_acc}, popEVE has {used}"
                best = (rank, dict(np=used, seq=seq, map=pos_map, note=note,
                                   scores=dict(zip(d.mutant, d[self.score].astype(float)))))
        self.cache = {uniprot: best[1] if best and best[1]["map"] else None}               # one protein at a time
        return self.cache[uniprot]

    @staticmethod
    def get(p, wt_aa, upos, alt_aa):
        q = p["map"].get(upos) if p and upos else None
        return p["scores"].get(f"{wt_aa}{q}{alt_aa}", np.nan) if q and p["seq"][q - 1] == wt_aa else np.nan


def read_refseq(path):
    """Swiss-Prot human DR RefSeq lines: accession -> [(NP_xxx.y, is_displayed_isoform)], displayed isoforms first.
    Lines without an isoform marker belong to entries with a single (displayed) isoform."""
    out, acc, disp, rows = {}, None, None, []
    with gzip.open(path, "rt") as f:
        for line in f:
            tag = line[:2]
            if tag == "AC" and acc is None:
                acc = line[5:].split(";")[0].strip()
            elif tag == "CC" and (g := re.search(r"IsoId=([^;,]+)[^;]*; Sequence=Displayed", line)):
                disp = g[1]
            elif tag == "DR" and line[5:].startswith("RefSeq;"):
                g = re.match(r"RefSeq; (NP_\d+\.\d+);[^\[]*(?:\[([^\]]+)\])?", line[5:])
                if g:
                    rows.append((g[1], g[2]))
            elif tag == "//":
                out[acc] = sorted({(n, iso is None or iso == disp) for n, iso in rows}, key=lambda x: not x[1])
                acc, disp, rows = None, None, []
    return out


class Dbnsfp:
    """dbNSFP for GRCh38 SNVs: CADD PHRED (SNV-level, no transcript choice), and optional PrimateAI-3D."""

    def __init__(self, cfg):
        self.tabix = cfg["tabix"]["path"]
        self.db = cfg["dbnsfp"]["path"]
        with gzip.open(self.db, "rt") as f:
            self.hdr = {c: i for i, c in enumerate(f.readline().rstrip("\n").split("\t"))}
        need = ["#chr", "pos(1-based)", "ref", "alt", "CADD_phred"]
        missing = [c for c in need if c not in self.hdr]
        assert not missing, f"dbNSFP columns missing: {missing}"
        self.pai = cfg["primateai3d"] if cfg["primateai3d"].get("path") else None

    def cadd(self, snvs):
        """{(chrom,pos,ref,alt)} -> {snv: CADD PHRED}."""
        if not snvs:
            return {}
        with tempfile.NamedTemporaryFile("w", suffix=".tsv") as reg:
            for c, p in sorted({(s[0], s[1]) for s in snvs}):
                reg.write(f"{c}\t{p}\t{p}\n")
            reg.flush()
            out = subprocess.run([self.tabix, "-R", reg.name, self.db], capture_output=True, text=True, check=True)
        h, res = self.hdr, {}
        for line in out.stdout.splitlines():
            r = line.split("\t")
            k = (r[h["#chr"]], int(r[h["pos(1-based)"]]), r[h["ref"]], r[h["alt"]])
            if k in snvs and r[h["CADD_phred"]] not in {".", ""}:
                res[k] = float(r[h["CADD_phred"]])
        return res

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
    am_p = pe_p = None
    if "scores" in layers and uniprot and canon:
        am_p = scorer["am"].protein(uniprot, canon)
        pe_p = scorer["popeve"].protein(uniprot, canon)
    am_info = (f"UniProt {uniprot}: AlphaMissense sequence {am_p['note']}" if am_p and am_p["map"] else
               f"UniProt {uniprot} not in AlphaMissense")
    pe_info = (f"RefSeq {pe_p['np']} for UniProt {uniprot}: {pe_p['note']}" if pe_p else
               f"no popEVE RefSeq protein for UniProt {uniprot}")
    desc = {
        "uniprot_pos": f"mutated position(s) on the canonical UniProt sequence ({cfg['uniprot']['version']})",
        "score_alphamissense": "AlphaMissense pathogenicity (0-1, higher = more pathogenic), by UniProt + GRCh38 SNV, "
                               f"else by protein change; {cfg['alphamissense']['version']}; {am_info}",
        "score_alphamissense_source": "how the AlphaMissense value was found: snv (UniProt + the row's GRCh38 SNV) | "
                                      "protein (UniProt + protein change; mean over its SNVs) | none",
        "score_popeve": f"popEVE score (more negative = more deleterious), by protein change; "
                        f"{cfg['popeve']['version']}; {pe_info}",
        "score_cadd_phred": f"CADD PHRED (higher = more deleterious); {cfg['dbnsfp']['version']}",
        "score_primateai3d": f"PrimateAI-3D score (higher = more pathogenic); {cfg['primateai3d']['version']}",
        "score_snv_source": "SNVs used for the SNV-level scores (CADD, PrimateAI-3D): genomic (row's own GRCh38 SNV), "
                            "protein (all SNVs causing the protein change, AlphaMissense hg38 table; scores "
                            "averaged), none",
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
            snvs_of, src, am, am_src, pe = [], [], [], [], []
            for i, s in enumerate(single):
                upos = pmap.get(s[1]) if s else None
                key = AlphaMissense.key(am_p, s[0], upos, s[2]) if s and am_p else None
                g = SNV_RE.match(df[snv_col].iat[i]) if snv_col in df and s else None
                snv = (g[1].upper().removeprefix("CHR"), int(g[2]), g[3].upper(), g[4].upper()) if g else None
                if snv:
                    snvs_of.append([snv])
                    src.append("genomic")
                elif key and key in am_p["pv"]:
                    snvs_of.append(am_p["pv"][key][1])
                    src.append("protein")
                else:
                    snvs_of.append([])
                    src.append("none")
                if snv and am_p and snv in am_p["snv"]:
                    am.append(float(am_p["snv"][snv][1]))
                    am_src.append("snv")
                    scorer["am_snv_vs_protein"].append((path.name, f"{s[0]}{s[1]}{s[2]}", am_p["snv"][snv][0], key))
                elif key and key in am_p["pv"]:
                    am.append(am_p["pv"][key][0])
                    am_src.append("protein")
                else:
                    am.append(np.nan)
                    am_src.append("none")
                pe.append(PopEVE.get(pe_p, s[0], upos, s[2]) if s else np.nan)
            allsnv = {x for xs in snvs_of for x in xs}
            cadd_of = scorer["dbnsfp"].cadd(allsnv)
            pai = scorer["dbnsfp"].primateai3d(allsnv)
            cadd, p3 = [], []
            for xs in snvs_of:
                cv = [cadd_of[x] for x in xs if x in cadd_of]
                cadd.append(float(np.mean(cv)) if cv else np.nan)
                pv = [pai[x] for x in xs if x in pai]
                p3.append(float(np.mean(pv)) if pv else np.nan)
            new.update(score_alphamissense=am, score_alphamissense_source=am_src, score_popeve=pe,
                       score_cadd_phred=cadd, score_snv_source=src)
            if scorer["dbnsfp"].pai:
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
        extra = dict(am_sequence=am_p["note"] if am_p else "", popeve_refseq=pe_p["np"] if pe_p else "",
                     popeve_sequence=pe_p["note"] if pe_p else "") if "scores" in layers else {}
        summary.append(dict(archive=path.name, assay=asec["name"], rows=len(df), uniprot=uniprot,
                            wt_to_uniprot_identity=round(ident, 4), **{f"cov_{c}": round(x, 4) for c, x in cov.items()},
                            **extra))
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
    scorer = dict(am=AlphaMissense(cfg), popeve=PopEVE(cfg, read_refseq(cfg["uniprot"]["path"])), dbnsfp=Dbnsfp(cfg),
                  am_snv_vs_protein=[]) if "scores" in layers else None
    rows = []
    for p in a.archives:
        rows += annotate(p, a.out, cfg, sprot, scorer, layers, a.snv_col, a.mutant_col, a.force)
        print("\t".join(f"{k}={v}" for k, v in rows[-1].items()), flush=True)
    pd.DataFrame(rows).to_csv(a.out / "layers_summary.tsv", sep="\t", index=False)
    if scorer:   # AlphaMissense found by SNV: does the SNV's protein change agree with the row's (UniProt numbering)?
        chk = pd.DataFrame(scorer["am_snv_vs_protein"], columns=["archive", "variant", "am_snv_protein_variant",
                                                                 "expected_protein_variant"])
        bad = chk[chk.am_snv_protein_variant != chk.expected_protein_variant]
        bad.to_csv(a.out / "layers_am_snv_mismatch.tsv", sep="\t", index=False)
        print(f"AlphaMissense by SNV: {len(chk)} rows, {len(bad)} whose SNV gives another protein change "
              f"(layers_am_snv_mismatch.tsv)")
    print(f"done: {len(a.archives)} archives -> {a.out}")


if __name__ == "__main__":
    main()
