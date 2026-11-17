#!/bin/bash
# Genebass runs of a *_filtered model: one validate_msa_filter.py run per protein (one container start each, so the
# per-protein time and peak RSS are clean), scoring the union of the protein's Genebass mutants over all its archives.
# Usage: sbatch [--mem/--time/--cpus-per-task ...] run_genebass.sbatch <image> <filtered folder> <UNIPROT> [<UNIPROT>...]
#   proteins and archives from $OUT/proteins.tsv (uniprot_id, gene, L, msa, n_mut, archives); outputs
#   $OUT/<folder>__<UNIPROT>__genebass.{csv,json|failed.json}; finished proteins are skipped.
#SBATCH --job-name=genebass_filtered
#SBATCH --partition=genoa64
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
set -uo pipefail
M=$1; F=$2; shift 2
REPO=/users/diasfrazer/egalpern/project/pgb_clinical_models
export OUT=${OUT:-/no_backup/diasfrazer/egalpern_pgb_clinical/models/msa_filter/genebass}
for U in "$@"; do
  TAG="${F}__${U}__genebass"
  if [ -f "$OUT/$TAG.json" ] || [ -f "$OUT/$TAG.failed.json" ]; then echo "$U: done, skipped"; continue; fi
  ARCHS=$(awk -F'\t' -v u="$U" '$1 == u {print $6}' "$OUT/proteins.tsv")
  [ -z "$ARCHS" ] && { echo "$U: not in proteins.tsv"; continue; }
  FIRST=${ARCHS%%,*}
  echo "=== $U $(date +%T) archives $ARCHS"
  bash "$REPO/scripts/clinical/models/msa_filter/validate_msa_filter.sbatch" "$M" "$F" "$FIRST" "$TAG" \
    --singles none --mutant-archives "$ARCHS"
done
echo "GENEBASS FILTERED DONE"
