#!/usr/bin/env bash
# One chunk of card proteins vs UniRef100 2026_01, PG2 MSA protocol (scripts/clinical/msa/colabfold_search_uniref100.py).
# QUERIES = chunk FASTA; writes <ID>.a3m into $WORK/msas.
set -euo pipefail
source scripts/clinical/cards/env.sh
date
"$PY" -u scripts/clinical/msa/colabfold_search_uniref100.py "$QUERIES" "$UNIREF100_DB" "$WORK/msas" \
  --mmseqs "$MMSEQS" --threads "$SLURM_CPUS_PER_TASK" --max-seqs 1000000 --split-memory-limit 300G \
  --work "$WORK/search_work/$(basename "$QUERIES" .fasta)"
ls "$WORK/msas"/*.a3m | wc -l; du -sh "$WORK/msas"; date
echo "SEARCH OK"
