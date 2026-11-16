#!/bin/bash
# Score region archives with an unchanged zoo model in its own image (run_regions.py, one container start).
# Usage: sbatch [--partition/--gres/--mem/--time ...] scripts/clinical/models/regions/run_regions.sbatch \
#          <model> <tag> <region archive glob or dir> [<more>...]
#   model: image name in $IMG (eve evmutation evmutation2 gemme rsalor structuredca); tag: output subfolder
#   HP='{"key": value}' overrides card hyper_parameters (e.g. HP='{"max_msa_depth": 50000}').
#SBATCH --job-name=run_regions
#SBATCH --partition=genoa64
#SBATCH --qos=short
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
set -uo pipefail
M=$1; TAG=$2; shift 2
REPO=${REPO:-$PWD}
IMG=${IMG:-/no_backup/diasfrazer/egalpern_pgb_clinical/containers}
ARCH=${ARCH:-/no_backup/diasfrazer/egalpern_pgb_clinical/models/regions/archives}
OUT=${OUT:-/no_backup/diasfrazer/egalpern_pgb_clinical/models/regions/preds/$M/$TAG}
mkdir -p "$OUT"
NV=""; [ -n "${CUDA_VISIBLE_DEVICES:-}" ] && NV="--nv"
FILES=()
for x in "$@"; do
  if [ -d "$ARCH/$x" ]; then FILES+=("$ARCH/$x"/*.splits.pgdata); else FILES+=($ARCH/$x); fi
done
IN=("${FILES[@]/#$ARCH//arch}")
echo "host $(hostname) job ${SLURM_JOB_ID:-none} model $M tag $TAG ${#FILES[@]} regions cpus ${SLURM_CPUS_PER_TASK:-?} gpu ${CUDA_VISIBLE_DEVICES:-none}"
export APPTAINERENV_REGIONS_HP="${HP:-{\}}"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8} NUMBA_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
apptainer exec $NV --cleanenv --env OMP_NUM_THREADS=$OMP_NUM_THREADS,NUMBA_NUM_THREADS=$NUMBA_NUM_THREADS \
  --bind "$REPO/scripts/clinical/models/regions:/scripts:ro" --bind "$ARCH:/arch:ro" --bind "$OUT:/out" \
  "$IMG/$M.sif" bash -c "cd \$(mktemp -d) && /opt/program/.venv/bin/python -u /scripts/run_regions.py /opt/program/README.md /out ${IN[*]}" 2>&1 \
  | grep --line-buffered -v "UserWarning\|warnings.warn\|it/s\]$"
echo "RUN REGIONS DONE"
