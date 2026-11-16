"""One-time: the protein sequence of every RefSeq protein in the popEVE release zip ([popeve] sequences).

The sequence is rebuilt from the wild-type letters of each CSV's `mutant` column ('X' where a position has none).
add_layers.py uses it to find the popEVE protein of a UniProt entry by sequence when none of the RefSeq accessions
listed in Swiss-Prot's DR lines is in the zip (e.g. the RefSeq accession changed). Columns: np (NP_xxx.y), length,
sequence.
Run (cluster, ~10 min): python build_popeve_index.py [--sources sources.toml] [--workers 8]
"""
import argparse
import tomllib
import zipfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd

from add_layers import sequence_from_mutants


def seqs(args):
    path, names = args
    z = zipfile.ZipFile(path)
    return [(Path(n).stem, sequence_from_mutants(pd.read_csv(z.open(n), usecols=["mutant"]).mutant)) for n in names]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sources", type=Path, default=Path(__file__).with_name("sources.toml"))
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    cfg = tomllib.loads(a.sources.read_text())["popeve"]
    names = [n for n in zipfile.ZipFile(cfg["zip"]).namelist() if n.endswith(".csv")]
    chunks = [(cfg["zip"], names[i::a.workers * 4]) for i in range(a.workers * 4)]
    with ProcessPoolExecutor(a.workers) as ex:
        rows = [r for part in ex.map(seqs, chunks) for r in part]
    d = pd.DataFrame(rows, columns=["np", "sequence"]).sort_values("np")
    d.insert(1, "length", d.sequence.str.len())
    out = Path(cfg["sequences"])
    out.parent.mkdir(parents=True, exist_ok=True)
    d.to_parquet(out, index=False)
    print(f"wrote {out} ({len(d):,} proteins)")


if __name__ == "__main__":
    main()
