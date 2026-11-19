#!/bin/bash
#SBATCH --job-name=filter_msas
#SBATCH --partition=genoa64
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=02:00:00
#SBATCH --output=/no_backup/diasfrazer/egalpern_pgb_clinical/logs/filter_msas_%j.log
# Final archives (<= 2,700 aa, EVE-protocol filtered alignments): splits_layers -> final_le2700 (filter_msas.py).
# Largest a3m in the set: ABCA1 2.6 GB (~10 GB peak per worker); ~20 GB of a3m in total.
set -euo pipefail
WORK=/no_backup/diasfrazer/egalpern_pgb_clinical
PY=/users/diasfrazer/egalpern/project/proteingym_clinical/env/bin/python
cd /users/diasfrazer/egalpern/project/pgb_clinical
D=datasets/clinical/genebass
$PY scripts/clinical/msa/filter_msas.py --index $D/genebass_fig6b.csv --archives "$WORK/genebass/splits_layers" \
  --msa-dir "$WORK/genebass_msas" --out "$WORK/genebass/final_le2700" \
  --out-index $D/genebass_le2700.csv --out-filter $D/genebass_le2700_msa_filter.tsv --max-length 2700
du -sh "$WORK/genebass/splits_layers" "$WORK/genebass/final_le2700"
echo "FILTER MSAS DONE"
