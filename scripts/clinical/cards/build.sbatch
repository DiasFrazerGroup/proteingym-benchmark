#!/usr/bin/env bash
# Archives for the given cards (build_cards.py; args passed through, e.g. --cards clinvar_v0 [--genes 3]).
# No layers or model scores: the archives hold benchmark data only.
set -euo pipefail
source scripts/clinical/cards/env.sh
date
"$PY" -u scripts/clinical/cards/build_cards.py --prepared "$WORK/cards" --msa-dir $MSA_DIRS "$@"
du -sh "$WORK"/cards/*/splits; date
echo "BUILD OK"
