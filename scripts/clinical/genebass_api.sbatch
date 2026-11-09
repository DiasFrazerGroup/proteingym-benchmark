#!/usr/bin/env bash
# Genebass benchmark from Genebass itself (API), same Fig. 6b rules; then archives (no MSA), Fig. 6b check, layers.
# Submit from the repo root:
#   sbatch -p genoa64 --qos=short -t 04:00:00 -c 4 --mem=48G --chdir=<repo> scripts/clinical/genebass_api.sbatch
set -euo pipefail
PY=${PY:-/users/diasfrazer/egalpern/project/proteingym_clinical/env/bin/python}
WORK=${WORK:-/no_backup/diasfrazer/egalpern_pgb_clinical}
O=$WORK/genebass_api
F=/no_backup/diasfrazer/db/vesm_djzn25_figdata/extracted/data/Fig6_data
date
"$PY" scripts/clinical/genebass_api.py --fetch "$O/fetched" --fig6 "$F" --out "$O"
"$PY" scripts/clinical/build_genebass.py --fig6 "$F" --hits "$O/genebass_api_hits.csv" \
  --betas "$O/genebass_api_betas.csv.gz" \
  --sprot /no_backup/diasfrazer/db/uniprot_sprot_human/uniprot_sprot_human.dat.gz \
  --afdb /users/diasfrazer/egalpern/project/human_proteome_de/data/structures --out "$O/build"
"$PY" scripts/clinical/check_genebass.py --fig6 "$F" --out "$O/build" || echo "CHECK: differs from paper (see above)"
"$PY" scripts/clinical/layers/add_layers.py "$O"/build/splits/*.splits.pgdata --out "$O/build/splits_layers"
date
echo "GENEBASS_API JOB OK"
