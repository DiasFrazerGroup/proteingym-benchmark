#!/usr/bin/env bash
# Submit the ClinVar / de novo / gnomAD benchmark jobs on genoa64 (from the repo root):
#   bash scripts/clinical/benchmarks/submit.sh prepare        variant selection + UniProt mapping + MSA queries
#   bash scripts/clinical/benchmarks/submit.sh msa [N]        MSA search for proteins without one, in chunks of N
#                                                             (default 150); then build + check (afterok)
#   bash scripts/clinical/benchmarks/submit.sh build [args]   build now (one job per benchmark) + check; extra args go
#                                                             to build.py (--genes N: smoke run, no check)
#   bash scripts/clinical/benchmarks/submit.sh check          check only
# ONLY="clinvar denovo" restricts the build to those benchmarks (the check always covers all three).
set -euo pipefail
source scripts/clinical/benchmarks/env.sh
mkdir -p "$WORK/logs" "$WORK/msas" "$WORK/queries"
common=(--partition=genoa64 --export=ALL --chdir="$(pwd)")
BENCH=(${ONLY:-clinvar denovo gnomad})

submit_check() {  # $1 = dependency spec (may be empty)
  local d=()
  [ -n "$1" ] && d=(--dependency="$1")
  j=$(sbatch --parsable "${common[@]}" "${d[@]}" --job-name=bench_check --qos=short --time=03:00:00 \
      --cpus-per-task=2 --mem=64G --output="$WORK/logs/check_%j.log" $S/check.sbatch)
  echo "check $j $1" | tee -a "$WORK/logs/jobs.txt"
}

submit_build() {  # $1 = dependency spec (may be empty); rest = build.py args
  local dep=$1; shift
  local d=() ids=()
  [ -n "$dep" ] && d=(--dependency="$dep")
  for b in "${BENCH[@]}"; do
    j=$(sbatch --parsable "${common[@]}" "${d[@]}" --job-name=bench_build_$b --qos=normal --time=24:00:00 \
        --cpus-per-task=2 --mem=64G --output="$WORK/logs/build_${b}_%j.log" $S/build.sbatch --benchmarks "$b" "$@")
    ids+=("$j"); echo "build $b $j $*" | tee -a "$WORK/logs/jobs.txt"
  done
  [[ " $* " != *" --genes "* ]] && submit_check "afterok:$(IFS=:; echo "${ids[*]}")"
  return 0
}

case "$1" in
  prepare)
    j=$(sbatch --parsable "${common[@]}" --job-name=bench_prepare --qos=short --time=03:00:00 --cpus-per-task=4 \
        --mem=96G --output="$WORK/logs/prepare_%j.log" $S/prepare.sbatch)
    echo "prepare $j" | tee -a "$WORK/logs/jobs.txt" ;;
  msa)
    n=${2:-150}
    rm -f "$WORK"/queries/chunk_*.fasta
    awk -v n="$n" -v d="$WORK/queries" '/^>/{if(c%n==0){f=sprintf("%s/chunk_%03d.fasta",d,c/n)} c++} {print > f}' \
      "$OUT/queries.fasta"
    deps=()
    for q in "$WORK"/queries/chunk_*.fasta; do
      j=$(QUERIES=$q sbatch --parsable "${common[@]}" --job-name=bench_msa_$(basename $q .fasta) --qos=vlong \
          --time=2-00:00:00 --cpus-per-task=64 --mem=400G --output="$WORK/logs/msa_%j.log" $S/msa_search.sbatch)
      deps+=("$j"); echo "msa $(basename $q) $j" | tee -a "$WORK/logs/jobs.txt"
    done
    submit_build "afterok:$(IFS=:; echo "${deps[*]}")" ;;
  build)
    shift
    submit_build "" "$@" ;;
  check)
    submit_check "" ;;
esac
