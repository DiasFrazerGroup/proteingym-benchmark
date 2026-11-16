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
isoform_table  from AlphaMissense_isoforms_hg38.tsv.gz (optional): chrom, pos, ref, alt, transcript_id,
               protein_variant, am_pathogenicity for non-canonical GENCODE v32 transcripts (144.6M rows), kept in
               genomic order. add_layers.py uses it only for proteins whose AlphaMissense canonical sequence is not
               the current UniProt one (e.g. LPA): the row's SNV, on a transcript giving the same protein change.
Row groups of ~1M rows (sorted by UniProt) or 0.5M rows (genomic order), so a lookup touches few row groups.
Run (cluster, ~15 min, ~60 GB): python build_snv_index.py [--sources sources.toml] [--only table|protein_table|isoform_table]
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
        "transcript_id": "transcript_id", "protein_variant": "protein_variant", "am_pathogenicity": "am_pathogenicity"}


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


def build_stream(source, out, columns):
    """Large genomic table, already sorted by position: streamed (low memory), row groups keep chrom/pos statistics."""
    with gzip.open(source, "rt") as f:
        lines = [f.readline() for _ in range(4)]
    names = lines[3].lstrip("#").rstrip("\n").split("\t")
    rd = pacsv.open_csv(
        source,
        read_options=pacsv.ReadOptions(skip_rows=4, column_names=names, block_size=1 << 26),
        parse_options=pacsv.ParseOptions(delimiter="\t"),
        convert_options=pacsv.ConvertOptions(
            include_columns=[c for c in names if c in KEEP],
            column_types={"POS": pa.int64(), "am_pathogenicity": pa.float32(), "CHROM": pa.string()}))
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp, writer, n = out.with_suffix(".tmp.parquet"), None, 0
    for b in rd:
        t = pa.Table.from_batches([b])
        t = t.rename_columns([KEEP[c] for c in t.column_names])
        t = t.set_column(t.column_names.index("chrom"), "chrom", pc.replace_substring_regex(t["chrom"], "^chr", ""))
        t = t.select(columns)
        writer = writer or pq.ParquetWriter(tmp, t.schema, compression="zstd")
        writer.write_table(t, row_group_size=500_000)
        n += t.num_rows
    writer.close()
    tmp.rename(out)
    print(f"wrote {out} ({n:,} rows)", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sources", type=Path, default=Path(__file__).with_name("sources.toml"))
    ap.add_argument("--only", choices=["table", "protein_table", "isoform_table"])
    a = ap.parse_args()
    cfg = tomllib.loads(a.sources.read_text())["alphamissense"]
    if a.only in (None, "table"):
        build(cfg["source"], Path(cfg["table"]),
              ["uniprot_id", "protein_variant", "chrom", "pos", "ref", "alt", "am_pathogenicity"])
    if a.only in (None, "protein_table"):
        build(cfg["protein_source"], Path(cfg["protein_table"]), ["uniprot_id", "protein_variant", "am_pathogenicity"])
    if a.only in (None, "isoform_table") and cfg.get("isoform_table"):
        build_stream(cfg["isoform_source"], Path(cfg["isoform_table"]),
                     ["chrom", "pos", "ref", "alt", "transcript_id", "protein_variant", "am_pathogenicity"])


if __name__ == "__main__":
    main()
