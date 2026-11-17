#!/usr/bin/env bash
# Check every card's archives against clinical-curation's benchmark_cards/summary.tsv (check_cards.py).
set -euo pipefail
source scripts/clinical/cards/env.sh
"$PY" -u scripts/clinical/cards/check_cards.py --prepared "$WORK/cards" --curation "$CURATION"
