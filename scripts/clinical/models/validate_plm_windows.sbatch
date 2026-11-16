#!/bin/bash
# Windowed PLM validation + scores (validate_plm_windows.py) in the zoo base image, window package mounted.
# Usage: sbatch [--gres ...] scripts/clinical/models/validate_plm_windows.sbatch esm2|esmc [--windows 2048,0] [--limit N]
#SBATCH --job-name=plm_windows
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
set -euo pipefail
M=$1; shift
REPO=${REPO:-$PWD}
P=/users/diasfrazer/egalpern/project/proteingym_clinical/containers   # prototype images: zoo esm2/esmc, evedesign 0.0.5
ARCH=${ARCH:-/no_backup/diasfrazer/egalpern_pgb_clinical/genebass_api/build/splits}
OUT=${OUT:-/no_backup/diasfrazer/egalpern_pgb_clinical/models/plm_windows/$M}
mkdir -p "$OUT"
apptainer exec --nv --cleanenv --env HF_HOME=/hf,TORCH_HOME=/hf,PYTHONPATH=/src \
  --bind "$P/hf:/hf" --bind "$REPO/models/${M}_window/src:/src:ro" --bind "$REPO/models/${M}_window/README.md:/card.md:ro" --bind "$REPO/scripts/clinical/models:/scripts:ro" \
  --bind "$ARCH:/data:ro" --bind "$OUT:/out" \
  "$P/$M.sif" bash -c "cd /opt/program && .venv/bin/python /scripts/validate_plm_windows.py $M /data /card.md --out /out $*" 2>&1 \
  | grep --line-buffered -v "Loading weights"
echo "PLM WINDOWS DONE"
