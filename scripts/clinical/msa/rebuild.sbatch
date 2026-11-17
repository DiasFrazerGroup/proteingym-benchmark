#!/usr/bin/env bash
# Rebuild the Genebass archives with the ColabFold/UniRef100 MSAs, then rerun the Fig. 6b reproduction check,
# then add the annotation layers (scripts/clinical/layers).
# Archives go to $WORK/genebass (large: MSAs are copied into every archive); the index CSV is copied into the repo.
set -euo pipefail
F=/no_backup/diasfrazer/db/vesm_djzn25_figdata/extracted/data/Fig6_data
# association data from Genebass itself (scripts/clinical/genebass_api.sbatch must have run)
"$PY" scripts/clinical/build_genebass.py --fig6 "$F" \
  --hits "$WORK/genebass_api/genebass_api_hits.csv" --betas "$WORK/genebass_api/genebass_api_betas.csv.gz" \
  --sprot /no_backup/diasfrazer/db/uniprot_sprot_human/uniprot_sprot_human.dat.gz \
  --afdb /users/diasfrazer/egalpern/project/human_proteome_de/data/structures \
  --msa-dir "$WORK/genebass_msas" --out "$WORK/genebass"
# all rows: MISMATCH is expected (Genebass has more pairs and rows than the paper)
"$PY" scripts/clinical/check_genebass.py --fig6 "$F" --out "$WORK/genebass" || echo "CHECK: differs from paper"
# the paper's exact set, recovered from the in_vesm_paper flags alone: must print REPRODUCED
"$PY" scripts/clinical/check_genebass.py --fig6 "$F" --out "$WORK/genebass" --paper-set
cp "$WORK/genebass/genebass_fig6b.csv" datasets/clinical/genebass/genebass_fig6b.csv
# annotation layers (external scores + masks): final archives in $WORK/genebass/splits_layers
"$PY" scripts/clinical/layers/add_layers.py "$WORK"/genebass/splits/*.splits.pgdata --out "$WORK/genebass/splits_layers"
# ColabFold a3m headers (bit scores) back into the final archives (proteingym-base MSA.dump drops them)
"$PY" scripts/clinical/msa/restore_msa_headers.py "$WORK"/genebass/splits_layers/*.splits.pgdata \
  --msa-dir "$WORK/genebass_msas" --out "$WORK/genebass/splits_layers_hdr"
rm -rf "$WORK/genebass/splits_layers" && mv "$WORK/genebass/splits_layers_hdr" "$WORK/genebass/splits_layers"
du -sh "$WORK/genebass/splits" "$WORK/genebass/splits_layers"
echo "REBUILD OK"
