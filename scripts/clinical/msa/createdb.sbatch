#!/usr/bin/env bash
# UniRef100 (UniProt 2026_01, cluster copy in /db/uniprot) -> MMseqs2 sequence DB. No precomputed index: the search
# builds its k-mer index on the fly, split to fit memory (an index would add several hundred GB on disk).
set -euo pipefail
FASTA=/db/uniprot/2026_01/uniref/uniref100/blast/db/uniref100.fasta
if [ -f "$DB/uniref100.dbtype" ]; then echo "exists: $DB/uniref100"; exit 0; fi
date
"$MMSEQS" createdb "$FASTA" "$DB/uniref100.tmp" --threads "$SLURM_CPUS_PER_TASK"
for f in "$DB"/uniref100.tmp*; do mv "$f" "${f/uniref100.tmp/uniref100}"; done   # atomic-ish: only complete DBs get the final name
du -sh "$DB"; date
echo "CREATEDB OK"
