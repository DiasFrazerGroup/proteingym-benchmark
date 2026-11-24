"""FASTA (to stdout) of the benchmark proteins that have no MSA yet, shortest first.

Run: python queries.py <work dir of prepare.py> <MSA folder> <uniprot_sprot_human.dat.gz> > queries.fasta
"""
import sys
from pathlib import Path

import pandas as pd

from common import read_sprot

work, msas, sprot = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3]
sp = read_sprot(sprot)
ids = set()
for f in work.glob("*/variants.tsv.gz"):
    ids |= set(pd.read_csv(f, sep="\t", usecols=["uniprot_id"]).uniprot_id)
todo = sorted((u for u in ids if not (msas / f"{u}.a3m").exists()), key=lambda u: (len(sp[u]["sequence"]), u))
for u in todo:
    print(f">{u}\n{sp[u]['sequence']}")
print(f"{len(ids)} proteins, {len(todo)} without MSA", file=sys.stderr)
