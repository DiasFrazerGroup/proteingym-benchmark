#!/bin/bash
#SBATCH --job-name=msa_headers
#SBATCH --partition=genoa64
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=/no_backup/diasfrazer/egalpern_pgb_clinical/logs/msa_headers_%j.log
# Restore the ColabFold a3m headers (bit scores) in the final archives: splits_layers -> splits_layers_hdr.
# Streams and checks every sequence; ~60 GB of a3m in total, deflate-bound (~10-30 min on 16 cores).
set -euo pipefail
WORK=/no_backup/diasfrazer/egalpern_pgb_clinical
PY=/users/diasfrazer/egalpern/project/proteingym_clinical/env/bin/python
cd /users/diasfrazer/egalpern/project/pgb_clinical
$PY scripts/clinical/msa/restore_msa_headers.py "$WORK"/genebass/splits_layers/*.splits.pgdata \
  --msa-dir "$WORK/genebass_msas" --out "$WORK/genebass/splits_layers_hdr"
n_in=$(ls "$WORK"/genebass/splits_layers/*.splits.pgdata | wc -l); n_out=$(ls "$WORK"/genebass/splits_layers_hdr/*.splits.pgdata | wc -l)
du -sh "$WORK/genebass/splits_layers" "$WORK/genebass/splits_layers_hdr"
[ "$n_in" = "$n_out" ] && echo "HEADERS OK $n_out archives" || echo "HEADERS MISMATCH $n_in vs $n_out"
