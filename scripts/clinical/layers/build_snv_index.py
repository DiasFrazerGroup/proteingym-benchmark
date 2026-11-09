"""One-time: protein change -> causing SNVs (GRCh38), from AlphaMissense_hg38.tsv.gz, as parquet.

Columns: uniprot_id, protein_variant (e.g. V2L, UniProt canonical numbering as used by AlphaMissense), chrom (no
'chr'), pos, ref, alt. Used by add_layers.py to look up SNV-level scores (dbNSFP, PrimateAI-3D) for variants given
only as protein changes.
Run: python build_snv_index.py [--sources sources.toml]
"""
import argparse
import tomllib
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sources", type=Path, default=Path(__file__).with_name("sources.toml"))
    a = ap.parse_args()
    cfg = tomllib.loads(a.sources.read_text())["snv_index"]
    out = Path(cfg["path"])
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(".tmp.parquet")
    writer, n = None, 0
    reader = pd.read_csv(cfg["source"], sep="\t", skiprows=3, chunksize=5_000_000, dtype=str,
                         usecols=["#CHROM", "POS", "REF", "ALT", "uniprot_id", "protein_variant"])
    for ch in reader:
        ch = ch.rename(columns={"#CHROM": "chrom", "POS": "pos", "REF": "ref", "ALT": "alt"})
        ch["chrom"] = ch.chrom.str.removeprefix("chr")
        ch["pos"] = ch.pos.astype("int64")
        t = pa.Table.from_pandas(ch[["uniprot_id", "protein_variant", "chrom", "pos", "ref", "alt"]],
                                 preserve_index=False)
        writer = writer or pq.ParquetWriter(tmp, t.schema, compression="zstd")
        writer.write_table(t)
        n += len(ch)
        print(f"{n:,} rows", flush=True)
    writer.close()
    tmp.rename(out)
    print(f"wrote {out} ({n:,} SNVs)")


if __name__ == "__main__":
    main()
