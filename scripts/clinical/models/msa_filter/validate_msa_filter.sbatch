#!/bin/bash
# One validate_msa_filter.py run (build + score on one archive) in a zoo model's image, with the *_filtered package
# mounted on PYTHONPATH. Usage:
#   CARD=zoo|filtered [TIMEOUT=2h] sbatch [--partition/--gres/--mem/--time] validate_msa_filter.sbatch \
#        <image> <filtered folder> <archive file name> <tag> [--hp JSON] [--class X] [--singles all|N]
#   image: eve evmutation rsalor structuredca (in $IMG); CARD=zoo uses the image's baked zoo card and class.
#   The original ColabFold a3m files ($A3M) are mounted at /a3m (use --hp {"a3m_headers_dir": "/a3m"}).
#SBATCH --job-name=msa_filter_val
#SBATCH --partition=genoa64
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
set -uo pipefail
M=$1; F=$2; A=$3; TAG=$4; shift 4
REPO=/users/diasfrazer/egalpern/project/pgb_clinical_models
IMG=${IMG:-/no_backup/diasfrazer/egalpern_pgb_clinical/containers}
ARCH=${ARCH:-/no_backup/diasfrazer/egalpern_pgb_clinical/genebass/splits_layers}
OUT=${OUT:-/no_backup/diasfrazer/egalpern_pgb_clinical/models/msa_filter/val}
mkdir -p "$OUT"
CARDP=/card.md; [ "${CARD:-filtered}" = zoo ] && CARDP=/opt/program/README.md
NV=""; [ -n "${CUDA_VISIBLE_DEVICES:-}" ] && NV="--nv"
echo "host $(hostname) job ${SLURM_JOB_ID:-none} image $M pkg $F archive $A tag $TAG card $CARDP timeout ${TIMEOUT:-none} args $*"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8} NUMBA_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
ARGS=()
for x in "$@"; do ARGS+=("$x"); done
${TIMEOUT:+timeout $TIMEOUT} apptainer exec $NV --cleanenv \
  --env PYTHONPATH=/src,OMP_NUM_THREADS=$OMP_NUM_THREADS,NUMBA_NUM_THREADS=$NUMBA_NUM_THREADS \
  --bind "$REPO/models/$F/src:/src:ro" --bind "$REPO/models/$F/README.md:/card.md:ro" \
  --bind "$REPO/scripts/clinical/models/msa_filter:/scripts:ro" --bind "$ARCH:/arch:ro" --bind "$OUT:/out" --bind "${A3M:-/no_backup/diasfrazer/egalpern_pgb_clinical/genebass_msas}:/a3m:ro" \
  "$IMG/$M.sif" bash -c 'cd $(mktemp -d) && exec /opt/program/.venv/bin/python -u /scripts/validate_msa_filter.py run "$@"' _ \
  "$CARDP" "/arch/$A" /out "$TAG" "${ARGS[@]}" 2>&1 \
  | grep --line-buffered -v "UserWarning\|warnings.warn\|it/s\]$\|s/it\]$\| | INFO \|SyntaxWarning\|invalid escape"
rc=${PIPESTATUS[0]}
echo "EXIT $rc (124 = timeout)"
[ $rc -eq 124 ] && echo "{\"tag\": \"$TAG\", \"timeout\": \"${TIMEOUT}\"}" > "$OUT/$TAG.timeout.json"
echo "MSA FILTER VAL DONE"
