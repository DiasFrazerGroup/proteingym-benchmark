"""One-time: AlphaMissense's SNV table as a parquet index, used to place genomic variants on UniProt.

AlphaMissense_hg38.tsv.gz / AlphaMissense_hg19.tsv.gz (Zenodo 10.5281/zenodo.8208688) list every possible missense
SNV of the canonical transcript of each UniProt canonical isoform. The index keeps
uniprot_id, protein_variant (UniProt numbering, e.g. V2L), chrom (no 'chr'), pos, ref, alt (71.7M rows per build).

Run (~10 min, ~40 GB): python build_snv_index.py AlphaMissense_hg38.tsv.gz am_hg38_snv_index.parquet
"""
import sys

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.csv as pacsv
import pyarrow.parquet as pq

source, out = sys.argv[1:3]
t = pacsv.read_csv(
    source,
    read_options=pacsv.ReadOptions(skip_rows=3, block_size=1 << 26),          # 3 licence lines, then the header
    parse_options=pacsv.ParseOptions(delimiter="\t"),
    convert_options=pacsv.ConvertOptions(
        include_columns=["uniprot_id", "protein_variant", "#CHROM", "POS", "REF", "ALT"],
        column_types={"POS": pa.int64(), "#CHROM": pa.string()}))
t = t.rename_columns(["uniprot_id", "protein_variant", "chrom", "pos", "ref", "alt"])
t = t.set_column(2, "chrom", pc.replace_substring_regex(t["chrom"], "^chr", ""))
pq.write_table(t.sort_by([("uniprot_id", "ascending"), ("pos", "ascending")]), out, compression="zstd",
               row_group_size=1_000_000)
print(f"wrote {out} ({t.num_rows:,} rows, {pc.count_distinct(t['uniprot_id']).as_py():,} proteins)")
