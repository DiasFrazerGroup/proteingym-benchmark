#!/bin/bash
# Build Apptainer images of zoo model folders (Singularity.def, --fakeroot) into $IMG.
# Usage: sbatch scripts/clinical/models/build_images.sbatch eve evmutation ...
#SBATCH --job-name=build_images
#SBATCH --partition=genoa64
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
set -uo pipefail
REPO=${REPO:-$PWD}
IMG=${IMG:-/no_backup/diasfrazer/egalpern_pgb_clinical/containers}
export APPTAINER_CACHEDIR=$IMG/cache APPTAINER_TMPDIR=$IMG/tmp
mkdir -p "$IMG/cache" "$IMG/tmp" "$IMG/logs"
echo "host $(hostname) job ${SLURM_JOB_ID:-none}"
for m in "$@"; do
  [ -f "$IMG/$m.sif" ] && { echo "$m: exists"; continue; }
  t=$(date +%s)
  if (cd "$REPO/models/$m" && apptainer build --fakeroot "$IMG/$m.sif" Singularity.def) > "$IMG/logs/build_$m.log" 2>&1; then
    echo "$m: BUILT in $(( ($(date +%s) - t) / 60 )) min"
  else
    echo "$m: BUILD FAILED ($IMG/logs/build_$m.log)"; tail -5 "$IMG/logs/build_$m.log"
  fi
done
echo "BUILD IMAGES DONE"
