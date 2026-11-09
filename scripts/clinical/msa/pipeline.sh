#!/usr/bin/env bash
# Submit the Genebass MSA chain on the CRG cluster (all jobs independent of any interactive session):
#   smoke     colabfold_search_uniref100.py on 2 proteins vs a Swiss-Prot-human mini DB (recipe + a3m check)
#   createdb  UniRef100 2026_01 FASTA (/db/uniprot) -> MMseqs2 DB                       (parallel with smoke)
#   search    61 Genebass wild types vs UniRef100, PG2 protocol                          (afterok smoke+createdb)
#   rebuild   archives with MSAs + Fig. 6b reproduction check                            (afterok search)
# Usage: bash scripts/clinical/msa/pipeline.sh            (from the repo root)
set -euo pipefail
REPO=$(cd "$(dirname "$0")/../../.." && pwd)
WORK=/no_backup/diasfrazer/egalpern_pgb_clinical
DB=/no_backup/diasfrazer/egalpern_uniref100_2026_01_mmseqs
export REPO WORK DB
# MMseqs2 18 (CRG EasyBuild 18-git20261027-gompi-2025b-int64) copied with its 6 GCCcore/zlib/bzip2 libs behind a
# wrapper: /software differs between genoa64 nodes (genoa64-06 only has v17), so the module path is not portable.
export MMSEQS=$WORK/mmseqs18/bin/mmseqs
export PY=/users/diasfrazer/egalpern/project/proteingym_clinical/env/bin/python
S="$REPO/scripts/clinical/msa"
common=(--partition=genoa64 --export=ALL --chdir="$REPO")

j0=$(sbatch --parsable "${common[@]}" --job-name=msa_smoke --qos=vshort --time=00:45:00 --cpus-per-task=8 \
     --mem=32G --output="$WORK/logs/smoke_%j.log" "$S/smoke.sbatch")
j1=$(sbatch --parsable "${common[@]}" --job-name=u100_createdb --qos=normal --time=10:00:00 --cpus-per-task=4 \
     --mem=32G --output="$WORK/logs/createdb_%j.log" "$S/createdb.sbatch")
j2=$(sbatch --parsable "${common[@]}" --job-name=msa_search --qos=vlong --time=2-00:00:00 --cpus-per-task=64 \
     --mem=400G --dependency=afterok:$j0:$j1 --output="$WORK/logs/search_%j.log" "$S/search.sbatch")
j3=$(sbatch --parsable "${common[@]}" --job-name=gb_rebuild --qos=short --time=06:00:00 --cpus-per-task=8 \
     --mem=64G --dependency=afterok:$j2 --output="$WORK/logs/rebuild_%j.log" "$S/rebuild.sbatch")
echo "smoke $j0 | createdb $j1 | search $j2 (after $j0,$j1) | rebuild $j3 (after $j2)" | tee -a "$WORK/logs/jobs.txt"
