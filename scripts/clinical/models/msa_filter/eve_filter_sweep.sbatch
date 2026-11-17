#!/bin/bash
#SBATCH --job-name=eve_sweep
#SBATCH --partition=genoa64
#SBATCH --cpus-per-task=34
#SBATCH --mem=96G
#SBATCH --time=02:00:00
#SBATCH --output=/no_backup/diasfrazer/egalpern_pgb_clinical/logs/eve_sweep_%j.log
# 68 a3m files, 40 GB, streamed once each (one process per file; the largest are 3.9 GB / ~1M rows).
set -euo pipefail
W=/no_backup/diasfrazer/egalpern_pgb_clinical
mkdir -p $W/models/msa_filter
/users/diasfrazer/egalpern/miniforge3/envs/frustratometer_galpern2025/bin/python /users/diasfrazer/egalpern/project/pgb_clinical_models/scripts/clinical/models/msa_filter/eve_filter_sweep.py $W/genebass_msas $W/models/msa_filter/eve_sweep.tsv
