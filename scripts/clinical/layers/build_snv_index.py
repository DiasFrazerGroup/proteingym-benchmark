"""One-time: the AlphaMissense hg38 table as a parquet sorted by UniProt accession (score source + SNV index).

Columns: uniprot_id, protein_variant (e.g. V2L, UniProt canonical numbering as used by AlphaMissense), chrom (no
'chr'), pos, ref, alt, am_pathogenicity. One row per possible missense SNV of the canonical transcript of each
UniProt canonical isoform (71.7M rows). add_layers.py uses it
  - as the AlphaMissense score source (by UniProt accession + GRCh38 SNV, or by protein change), and
  - as the SNV index: protein change -> causing SNVs, for SNV-level scores (CADD via dbNSFP, PrimateAI-3D) of
    variants given only as protein changes.
Sorted by uniprot_id with ~1M-row groups, so reading one protein touches one or two row groups.
Run (cluster, ~5 min, ~20 GB): python build_snv_index.py [--sources sources.toml]
"""
import argparse
import gzip
import tomllib
from pathlib import Path

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.csv as pacsv
import pyarrow.parquet as pq


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sources", type=Path, default=Path(__file__).with_name("sources.toml"))
    a = ap.parse_args()
    cfg = tomllib.loads(a.sources.read_text())["alphamissense"]
    out = Path(cfg["table"])
    out.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(cfg["source"], "rt") as f:                    # 3 licence comment lines, then the header
        lines = [f.readline() for _ in range(4)]
    names = lines[3].lstrip("#").rstrip("\n").split("\t")
    keep = {"CHROM": "chrom", "POS": "pos", "REF": "ref", "ALT": "alt", "uniprot_id": "uniprot_id",
            "protein_variant": "protein_variant", "am_pathogenicity": "am_pathogenicity"}
    t = pacsv.read_csv(
        cfg["source"],
        read_options=pacsv.ReadOptions(skip_rows=4, column_names=names, block_size=1 << 26),
        parse_options=pacsv.ParseOptions(delimiter="\t"),
        convert_options=pacsv.ConvertOptions(
            include_columns=list(keep),
            column_types={"POS": pa.int64(), "am_pathogenicity": pa.float32(), "CHROM": pa.string()}))
    t = t.rename_columns([keep[c] for c in t.column_names])
    t = t.set_column(t.column_names.index("chrom"), "chrom", pc.replace_substring_regex(t["chrom"], "^chr", ""))
    t = t.select(["uniprot_id", "protein_variant", "chrom", "pos", "ref", "alt", "am_pathogenicity"])
    t = t.sort_by([("uniprot_id", "ascending"), ("pos", "ascending")])
    tmp = out.with_suffix(".tmp.parquet")
    pq.write_table(t, tmp, compression="zstd", row_group_size=1_000_000)
    tmp.rename(out)
    print(f"wrote {out} ({t.num_rows:,} SNVs, {pc.count_distinct(t['uniprot_id']).as_py():,} proteins)")


if __name__ == "__main__":
    main()
