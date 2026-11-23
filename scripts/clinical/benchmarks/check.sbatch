#!/usr/bin/env bash
# Check the benchmark archives and report the reference models' results (check.py).
set -euo pipefail
source scripts/clinical/benchmarks/env.sh
"$PY" -u $S/check.py --prepared "$OUT" --curation "$CURATION" --results "$OUT/results.tsv"
