#!/usr/bin/env bash
# Rebuild the clinical benchmarks (clinvar, denovo, gnomad, genebass) from the repo root:
#
#   bash scripts/clinical/run.sh prepare    select the variants, place them on UniProt, list the proteins without MSA
#   bash scripts/clinical/run.sh msa        MSAs for those proteins (hours; GPU=1 for the GPU search)
#   bash scripts/clinical/run.sh build      archives + index CSVs (datasets/clinical/<benchmark>/<benchmark>.csv)
#   bash scripts/clinical/run.sh check      verify the archives, report the reference models
#
# Paths and tools: scripts/clinical/env.sh. Benchmarks to run: BENCHMARKS="clinvar denovo" (default: all four).
# One-time inputs: build_snv_index.py (AlphaMissense SNV tables), genebass_api.py (Genebass tables), and the
# UniRef100 MMseqs2 database (see msa.sh).
set -euo pipefail
S=scripts/clinical
source $S/env.sh
B=${BENCHMARKS:-clinvar denovo gnomad genebass}

case $1 in
  index)
    mkdir -p "$SNV_INDEX"
    for g in hg38 hg19; do
      $PY $S/build_snv_index.py "$AM_RAW/AlphaMissense_$g.tsv.gz" "$SNV_INDEX/am_$g.parquet"
    done ;;
  prepare)
    $PY $S/prepare.py $B --out "$WORK" --curation "$CURATION" --genebass "$GENEBASS" --fig6 "$FIG6" \
      --sprot "$SPROT" --snv-index "$SNV_INDEX"
    $PY $S/queries.py "$WORK" "$MSAS" "$SPROT" > "$WORK/queries.fasta" ;;
  msa)
    bash $S/msa.sh "$WORK/queries.fasta" "$MSAS" ;;
  build)
    for b in $B; do
      $PY $S/build.py $b --variants "$WORK/$b/variants.tsv.gz" --out "$WORK/$b/splits" \
        --index datasets/clinical/$b/$b.csv --sprot "$SPROT" --afdb "$AFDB" --msas "$MSAS" --jobs "$JOBS"
    done ;;
  check)
    $PY $S/check.py $B --work "$WORK" --index datasets/clinical --curation "$CURATION" ;;
esac
