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
  SRC="$REPO/models/$m"
  # evmutation2's Singularity.def copies .python-version, which the folder does not have (upstream bug): build
  # from a copy with the version of its Dockerfile base image
  if grep -q "\.python-version" "$SRC/Singularity.def" && [ ! -f "$SRC/.python-version" ]; then
    rm -rf "$IMG/tmp/src_$m"; cp -r "$SRC" "$IMG/tmp/src_$m"; SRC="$IMG/tmp/src_$m"
    sed -n 's/^FROM python:\([0-9.]*\).*/\1/p' "$SRC/Dockerfile" | head -1 > "$SRC/.python-version"
    echo "$m: building from a copy with .python-version $(cat "$SRC/.python-version")"
  fi
  # GEMME/JET2: lcqb.upmc.fr is down (2026-11). Local copies of the same tarballs ($IMG/tarballs: GEMME.tgz from the
  # Wayback Machine snapshot 20231001225718, JET2.tgz from the authors' image elodielaine/gemme@sha256:a35099f7...);
  # the Singularity.def's own sha256 pins still check them. Build from a copy whose %setup puts them in /opt.
  if grep -q "lcqb.upmc.fr/GEMME/package" "$SRC/Singularity.def" && [ -f "$IMG/tarballs/GEMME.tgz" ] \
      && [ -f "$IMG/tarballs/JET2.tgz" ]; then
    rm -rf "$IMG/tmp/src_$m"; cp -r "$SRC" "$IMG/tmp/src_$m"; SRC="$IMG/tmp/src_$m"
    cp "$IMG/tarballs/GEMME.tgz" "$IMG/tarballs/JET2.tgz" "$SRC/"
    sed -i -e '/curl -fsSLk https:\/\/www.lcqb.upmc.fr/d' \
      -e 's#^    mkdir -p ${SINGULARITY_ROOTFS}/opt/program$#&\n    cp GEMME.tgz JET2.tgz ${SINGULARITY_ROOTFS}/opt/#' \
      -e 's#/tmp/\(GEMME\|JET2\)\.tgz#/opt/\1.tgz#g' "$SRC/Singularity.def"
    echo "$m: building from a copy with local GEMME/JET2 tarballs (sha256 still checked)"
  fi
  if (cd "$SRC" && apptainer build --fakeroot "$IMG/$m.sif" Singularity.def) > "$IMG/logs/build_$m.log" 2>&1; then
    echo "$m: BUILT in $(( ($(date +%s) - t) / 60 )) min"
  else
    echo "$m: BUILD FAILED ($IMG/logs/build_$m.log)"; tail -5 "$IMG/logs/build_$m.log"
  fi
done
echo "BUILD IMAGES DONE"
