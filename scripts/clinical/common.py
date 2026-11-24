"""Shared helpers for the clinical benchmark scripts: Swiss-Prot, SNV -> protein mapping, a3m and AFDB checks."""
import gzip
import re

import numpy as np
import polars as pl

MAX_LEN = 2700  # longest protein with a single full-length AlphaFold DB model
AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H",
       "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V"}


def read_sprot(path):
    """Swiss-Prot human flat file: accession -> dict(entry_name, gene, pfam_ids, sequence)."""
    out, rec = {}, None
    with gzip.open(path, "rt") as f:
        for line in f:
            tag = line[:2]
            if tag == "ID":
                rec = dict(entry_name=line.split()[1], acc=None, gene=None, pfam_ids=[], seq=[], in_seq=False)
            elif tag == "AC" and rec["acc"] is None:
                rec["acc"] = line[5:].split(";")[0].strip()
            elif tag == "GN" and rec["gene"] is None and (g := re.search(r"Name=([^;{ ]+)", line)):
                rec["gene"] = g[1]
            elif tag == "DR" and line[5:].startswith("Pfam;"):
                rec["pfam_ids"].append(line[5:].split(";")[1].strip())
            elif tag == "SQ":
                rec["in_seq"] = True
            elif tag == "//":
                out[rec["acc"]] = dict(entry_name=rec["entry_name"], gene=rec["gene"], pfam_ids=rec["pfam_ids"],
                                       sequence="".join(rec["seq"]))
            elif rec is not None and rec["in_seq"]:
                rec["seq"].append(line.replace(" ", "").strip())
    return out


def snv_to_protein(index, variant_ids, sp):
    """SNV ids (chrom-pos-ref-alt, no 'chr') -> DataFrame(variant_id, uniprot_id, mutant): the protein change each
    SNV causes on a Swiss-Prot canonical protein, from the AlphaMissense SNV index (build_snv_index.py). Kept only
    when the wild-type residue is the Swiss-Prot one."""
    v = pl.DataFrame({"variant_id": sorted(set(variant_ids.dropna()))}).with_columns(
        pl.col("variant_id").str.split_exact("-", 3).struct.rename_fields(["chrom", "pos", "ref", "alt"]).alias("s")
    ).unnest("s").with_columns(pl.col("pos").cast(pl.Int64, strict=False))
    t = pl.scan_parquet(index).join(v.lazy(), on=["chrom", "pos", "ref", "alt"]).select(
        "variant_id", "uniprot_id", mutant="protein_variant").collect().to_pandas()
    t = t[t.uniprot_id.isin(sp)]
    ok = [sp[u]["sequence"][int(m[1:-1]) - 1:int(m[1:-1])] == m[0] for u, m in zip(t.uniprot_id, t.mutant)]
    return t[ok]


def mutate(wt, mutant):
    p = int(mutant[1:-1])
    assert wt[p - 1] == mutant[0] != mutant[-1], mutant
    return wt[:p - 1] + mutant[-1] + wt[p:]


def a3m_depth(path, wt):
    """Number of sequences of an a3m whose first record must be the wild type; None if it is not."""
    n, first = 0, []
    with open(path, "rb") as f:
        for line in f:
            if line.startswith(b">"):
                n += 1
            elif n == 1:
                first.append(line.strip().decode())
    return n if "".join(first) == wt else None


def afdb_models(u, wt, afdb):
    """AlphaFold DB v6 models of `u` as [(path, start, end, mean pLDDT)], checked against the wild type.
    Proteins up to 2,700 aa have one full-length model; longer ones come as 1,400-aa fragments, one every 200 aa."""
    files = sorted(afdb.glob(f"AF-{u}-F*-model_v6.pdb"), key=lambda p: int(re.search(r"-F(\d+)-", p.name)[1]))
    out = []
    for p in files:
        ca = [(AA3.get(x[17:20], "X"), float(x[60:66])) for x in open(p) if x.startswith("ATOM") and x[12:16] == " CA "]
        seq = "".join(a for a, _ in ca)
        start = 200 * (int(re.search(r"-F(\d+)-", p.name)[1]) - 1) + 1
        ref = wt[start - 1:start - 1 + len(seq)]
        if len(ref) != len(seq) or sum(a != b for a, b in zip(ref, seq)) > 0.01 * len(seq):
            return []                                              # model of another sequence version
        out.append((p, start, start + len(seq) - 1, round(float(np.mean([b for _, b in ca])), 2)))
    return out
