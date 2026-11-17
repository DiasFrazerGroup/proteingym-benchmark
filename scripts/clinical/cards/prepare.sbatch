#!/usr/bin/env bash
# Card variant sets -> UniProt-mapped rows, protein list, MSA queries (prepare_cards.py).
set -euo pipefail
source scripts/clinical/cards/env.sh
date
"$PY" -u scripts/clinical/cards/prepare_cards.py --curation "$CURATION" --out "$WORK/cards" --msa-dirs $MSA_DIRS
date; echo "PREPARE OK"
