#!/usr/bin/env bash
# Submit the card archive chain on genoa64 (from the repo root):
#   bash scripts/clinical/cards/submit.sh prepare        prepare_cards.py (card rows -> UniProt-mapped rows, MSA queries)
#   bash scripts/clinical/cards/submit.sh msa [N]        split queries.fasta into chunks of N proteins (default 150), one
#                                                        MSA search job each; then build + check (afterok all chunks)
#   bash scripts/clinical/cards/submit.sh build [args]   build now (one job per card) + check; extra args go to
#                                                        build_cards.py (--genes N: smoke run, no check)
set -euo pipefail
source scripts/clinical/cards/env.sh
mkdir -p "$WORK/logs" "$WORK/msas" "$WORK/queries"
common=(--partition=genoa64 --export=ALL --chdir="$(pwd)")
S=scripts/clinical/cards
CARDS=(clinvar_v0 dnm_aptest_ge1 dnm_aptest_ge2 dnm_aptest_random300 gnomad_af_clinvar_dnm_ge1
       gnomad_af_clinvar_dnm_ge2 gnomad_af_random300)

submit_build() {  # $1 = dependency spec (may be empty); rest = build_cards.py args
  local dep=$1; shift
  local d=() ids=()
  [ -n "$dep" ] && d=(--dependency="$dep")
  for c in "${CARDS[@]}"; do
    j=$(sbatch --parsable "${common[@]}" "${d[@]}" --job-name=cards_build_$c --qos=normal --time=24:00:00 \
        --cpus-per-task=2 --mem=64G --output="$WORK/logs/build_${c}_%j.log" $S/build.sbatch --cards "$c" "$@")
    ids+=("$j"); echo "build $c $j $*" | tee -a "$WORK/logs/jobs.txt"
  done
  if [[ " $* " != *" --genes "* ]]; then
    j=$(sbatch --parsable "${common[@]}" --dependency=afterok:$(IFS=:; echo "${ids[*]}") --job-name=cards_check \
        --qos=short --time=03:00:00 --cpus-per-task=2 --mem=64G --output="$WORK/logs/check_%j.log" $S/check.sbatch)
    echo "check $j (after ${ids[*]})" | tee -a "$WORK/logs/jobs.txt"
  fi
}

case "$1" in
  prepare)
    j=$(sbatch --parsable "${common[@]}" --job-name=cards_prepare --qos=short --time=03:00:00 --cpus-per-task=4 \
        --mem=96G --output="$WORK/logs/prepare_%j.log" $S/prepare.sbatch)
    echo "prepare $j" | tee -a "$WORK/logs/jobs.txt" ;;
  msa)
    n=${2:-150}
    rm -f "$WORK"/queries/chunk_*.fasta
    awk -v n="$n" -v d="$WORK/queries" '/^>/{if(c%n==0){f=sprintf("%s/chunk_%03d.fasta",d,c/n)} c++} {print > f}' \
      "$WORK/cards/queries.fasta"
    deps=()
    for q in "$WORK"/queries/chunk_*.fasta; do
      j=$(QUERIES=$q sbatch --parsable "${common[@]}" --job-name=cards_msa_$(basename $q .fasta) --qos=vlong \
          --time=2-00:00:00 --cpus-per-task=64 --mem=400G --output="$WORK/logs/msa_%j.log" $S/msa_search.sbatch)
      deps+=("$j"); echo "msa $(basename $q) $j" | tee -a "$WORK/logs/jobs.txt"
    done
    submit_build "afterok:$(IFS=:; echo "${deps[*]}")" ;;
  build)
    shift
    submit_build "" "$@" ;;
esac
