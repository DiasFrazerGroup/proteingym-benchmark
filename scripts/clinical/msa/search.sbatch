#!/usr/bin/env bash
# 61 Genebass wild types vs UniRef100 2026_01: ColabFold monomer search, filter 0 (PG2 MSA protocol).
set -euo pipefail
date
"$PY" scripts/clinical/msa/colabfold_search_uniref100.py "$WORK/genebass_queries.fasta" "$DB/uniref100" \
  "$WORK/genebass_msas" --mmseqs "$MMSEQS" --threads "$SLURM_CPUS_PER_TASK" --max-seqs 1000000 \
  --split-memory-limit 300G --work "$WORK/search_work"
ls "$WORK/genebass_msas"/*.a3m | wc -l
du -sh "$WORK/genebass_msas"; date
echo "SEARCH OK"
