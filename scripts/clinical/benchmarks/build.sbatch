#!/usr/bin/env bash
# Archives for the given benchmarks (build.py; args passed through, e.g. --benchmarks clinvar [--genes 3]).
set -euo pipefail
source scripts/clinical/benchmarks/env.sh
date
"$PY" -u $S/build.py --prepared "$OUT" --msa-dir $MSA_DIRS "$@"
du -sh "$OUT"/*/splits; date
echo "BUILD OK"
