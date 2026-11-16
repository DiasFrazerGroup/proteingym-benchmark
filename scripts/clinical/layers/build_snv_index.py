"""One-time: the AlphaMissense tables as parquets sorted by UniProt accession ([alphamissense] in sources.toml).

table          from AlphaMissense_hg38.tsv.gz: uniprot_id, protein_variant (UniProt canonical numbering as used by
               AlphaMissense, e.g. V2L), chrom (no 'chr'), pos, ref, alt, am_pathogenicity. One row per possible
               missense SNV of the canonical transcript of each UniProt canonical isoform (71.7M rows, 19,117
               proteins). add_layers.py uses it for AlphaMissense by (UniProt, GRCh38 SNV), and as the SNV index
               (protein change -> causing SNVs) for SNV-level scores (CADD, PrimateAI-3D) of variants given only as
               protein changes.
protein_table  from AlphaMissense_aa_substitutions.tsv.gz: uniprot_id, protein_variant, am_pathogenicity. Every
               single substitution of every UniProt canonical isoform (216M rows, 20,516 proteins), also those
               without a hg38 canonical transcript (e.g. LPA) and substitutions no SNV can cause. add_layers.py uses
               it for AlphaMissense by (UniProt, protein change).
Row groups of ~1M rows, so reading one protein touches one or two row groups.
Run (cluster, ~10 min, ~60 GB): python build_snv_index.py [--sources sources.toml] [--only table|protein_table]
"""
import argparse
import gzip
import tomllib
from pathlib import Path

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.csv as pacsv
import pyarrow.parquet as pq

KEEP = {"CHROM": "chrom", "POS": "pos", "REF": "ref", "ALT": "alt", "uniprot_id": "uniprot_id",
        "protein_variant": "protein_variant", "am_pathogenicity": "am_pathogenicity"}


def build(source, out, columns):
    with gzip.open(source, "rt") as f:                           # 3 licence comment lines, then the header
        lines = [f.readline() for _ in range(4)]
    names = lines[3].lstrip("#").rstrip("\n").split("\t")
    use = [c for c in names if c in KEEP]
    t = pacsv.read_csv(
        source,
        read_options=pacsv.ReadOptions(skip_rows=4, column_names=names, block_size=1 << 26),
        parse_options=pacsv.ParseOptions(delimiter="\t"),
        convert_options=pacsv.ConvertOptions(
            include_columns=use,
            column_types={"POS": pa.int64(), "am_pathogenicity": pa.float32(), "CHROM": pa.string()}))
    t = t.rename_columns([KEEP[c] for c in t.column_names])
    if "chrom" in t.column_names:
        t = t.set_column(t.column_names.index("chrom"), "chrom", pc.replace_substring_regex(t["chrom"], "^chr", ""))
    t = t.select(columns)
    t = t.sort_by([("uniprot_id", "ascending")] + ([("pos", "ascending")] if "pos" in columns else []))
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(".tmp.parquet")
    pq.write_table(t, tmp, compression="zstd", row_group_size=1_000_000)
    tmp.rename(out)
    print(f"wrote {out} ({t.num_rows:,} rows, {pc.count_distinct(t['uniprot_id']).as_py():,} proteins)", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sources", type=Path, default=Path(__file__).with_name("sources.toml"))
    ap.add_argument("--only", choices=["table", "protein_table"])
    a = ap.parse_args()
    cfg = tomllib.loads(a.sources.read_text())["alphamissense"]
    if a.only in (None, "table"):
        build(cfg["source"], Path(cfg["table"]),
              ["uniprot_id", "protein_variant", "chrom", "pos", "ref", "alt", "am_pathogenicity"])
    if a.only in (None, "protein_table"):
        build(cfg["protein_source"], Path(cfg["protein_table"]), ["uniprot_id", "protein_variant", "am_pathogenicity"])


if __name__ == "__main__":
    main()
