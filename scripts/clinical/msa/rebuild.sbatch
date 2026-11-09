#!/usr/bin/env bash
# Rebuild the Genebass archives with the ColabFold/UniRef100 MSAs, then rerun the Fig. 6b reproduction check.
# Archives go to $WORK/genebass (large: MSAs are copied into every archive); the index CSV is copied into the repo.
set -euo pipefail
F=/no_backup/diasfrazer/db/vesm_djzn25_figdata/extracted/data/Fig6_data
"$PY" scripts/clinical/build_genebass.py --fig6 "$F" \
  --sprot /no_backup/diasfrazer/db/uniprot_sprot_human/uniprot_sprot_human.dat.gz \
  --afdb /users/diasfrazer/egalpern/project/human_proteome_de/data/structures \
  --msa-dir "$WORK/genebass_msas" --out "$WORK/genebass"
"$PY" scripts/clinical/check_genebass.py --fig6 "$F" --out "$WORK/genebass"
cp "$WORK/genebass/genebass_fig6b.csv" datasets/clinical/genebass/genebass_fig6b.csv
du -sh "$WORK/genebass/splits"
echo "REBUILD OK"
