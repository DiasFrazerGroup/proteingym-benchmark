#!/usr/bin/env bash
# 61 Genebass wild types vs UniRef100 2026_01: ColabFold monomer search, filter 0 (PG2 MSA protocol).
# QUERIES / SEARCH_WORK override the query FASTA and work dir (e.g. a later batch of proteins); outputs are one
# <ID>.a3m each in $WORK/genebass_msas, so existing MSAs of other proteins are left untouched.
set -euo pipefail
QUERIES=${QUERIES:-$WORK/genebass_queries.fasta}
SEARCH_WORK=${SEARCH_WORK:-$WORK/search_work}
date
"$PY" scripts/clinical/msa/colabfold_search_uniref100.py "$QUERIES" "$DB/uniref100" \
  "$WORK/genebass_msas" --mmseqs "$MMSEQS" --threads "$SLURM_CPUS_PER_TASK" --max-seqs 1000000 \
  --split-memory-limit 300G --work "$SEARCH_WORK"
ls "$WORK/genebass_msas"/*.a3m | wc -l
du -sh "$WORK/genebass_msas"; date
echo "SEARCH OK"
