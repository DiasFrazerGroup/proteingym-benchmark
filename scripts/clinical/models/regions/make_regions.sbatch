#!/bin/bash
# Region archives for one protein and several region sizes (make_regions.py), local proteingym-base env.
# Usage: sbatch scripts/clinical/models/regions/make_regions.sbatch <UNIPROT> <W> [<W> ...]   (extra args: EXTRA=...)
#SBATCH --job-name=make_regions
#SBATCH --partition=genoa64
#SBATCH --qos=short
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
set -uo pipefail
REPO=${REPO:-$PWD}
PY=/users/diasfrazer/egalpern/project/proteingym_clinical/env/bin/python
OUT=${OUT:-/no_backup/diasfrazer/egalpern_pgb_clinical/models/regions/archives}
U=$1; shift
echo "host $(hostname) job ${SLURM_JOB_ID:-none} $U windows $*"
for W in "$@"; do
  /usr/bin/time -v $PY -u "$REPO/scripts/clinical/models/regions/make_regions.py" --uniprot "$U" --window "$W" \
    --out "$OUT" --tmp "${TMPDIR:-/tmp}" ${EXTRA:-} 2>&1 | grep --line-buffered -v "^\s*\(Page\|Voluntary\|Involuntary\|Swaps\|File system\|Socket\|Signals\|Average\|Minor\|Major\|Exit\|Percent\|Command\|System time\|User time\)"
done
echo "MAKE REGIONS DONE"
