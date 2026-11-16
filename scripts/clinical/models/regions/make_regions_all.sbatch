#!/bin/bash
# Region archives for many proteins at one region size (make_regions.py), e.g. every Genebass protein.
# Usage: sbatch scripts/clinical/models/regions/make_regions_all.sbatch <W> <UNIPROT> [<UNIPROT> ...]
#   OUT (default $WORK/models/regions/full/archives), EXTRA (more make_regions.py args)
#SBATCH --job-name=make_regions_all
#SBATCH --partition=genoa64
#SBATCH --qos=short
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
set -uo pipefail
REPO=${REPO:-$PWD}
PY=/users/diasfrazer/egalpern/project/proteingym_clinical/env/bin/python
OUT=${OUT:-/no_backup/diasfrazer/egalpern_pgb_clinical/models/regions/full/archives}
W=$1; shift
echo "host $(hostname) job ${SLURM_JOB_ID:-none} W $W, $# proteins"
for U in "$@"; do
  [ -f "$OUT/$U/W$W/regions.json" ] && { echo "$U: exists"; continue; }
  $PY -u "$REPO/scripts/clinical/models/regions/make_regions.py" --uniprot "$U" --window "$W" --out "$OUT" \
    --tmp "${TMPDIR:-/tmp}" ${EXTRA:-} 2>&1 | grep --line-buffered -v "^  r"
  [ "${PIPESTATUS[0]}" -eq 0 ] || echo "$U: FAILED"
done
echo "MAKE REGIONS ALL DONE"
