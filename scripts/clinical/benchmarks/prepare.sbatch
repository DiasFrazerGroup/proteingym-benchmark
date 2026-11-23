#!/usr/bin/env bash
# Variant selection and UniProt mapping for the three benchmarks, protein list, MSA queries (prepare.py).
set -euo pipefail
source scripts/clinical/benchmarks/env.sh
date
"$PY" -u $S/prepare.py --curation "$CURATION" --out "$OUT" --msa-dirs $MSA_DIRS
date; echo "PREPARE OK"
