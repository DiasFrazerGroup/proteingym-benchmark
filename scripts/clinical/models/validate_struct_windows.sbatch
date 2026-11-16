#!/bin/bash
# Structure models for long proteins (validate_struct_windows.py) in the zoo base image, *_window package mounted.
# Usage: sbatch [--gres ...] scripts/clinical/models/validate_struct_windows.sbatch saprot|esmif1 [--limit N]
#SBATCH --job-name=struct_windows
#SBATCH --partition=gpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
set -euo pipefail
M=$1; shift
REPO=${REPO:-$PWD}
P=/users/diasfrazer/egalpern/project/proteingym_clinical/containers
case $M in saprot) F=saprot_window; SIF=saprot.sif;; esmif1) F=esm-if1_window; SIF=esm-if1.sif;; esac
ARCH=${ARCH:-/no_backup/diasfrazer/egalpern_pgb_clinical/genebass_api/build/splits}
OUT=${OUT:-/no_backup/diasfrazer/egalpern_pgb_clinical/models/struct_windows/$M}
mkdir -p "$OUT"
apptainer exec --nv --cleanenv --env HF_HOME=/hf,TORCH_HOME=/hf,PYTHONPATH=/src \
  --bind "$P/hf:/hf" --bind "$REPO/models/$F/src:/src:ro" --bind "$REPO/models/$F/README.md:/card.md:ro" \
  --bind "$REPO/scripts/clinical/models:/scripts:ro" --bind "$ARCH:/data:ro" --bind "$OUT:/out" \
  "$P/$SIF" bash -c "cd /opt/program && .venv/bin/python /scripts/validate_struct_windows.py $M /data /card.md --out /out $*" 2>&1 \
  | grep -v "Loading weights\|UserWarning\|labeled_dataset_from_df"
echo "STRUCT WINDOWS DONE"
