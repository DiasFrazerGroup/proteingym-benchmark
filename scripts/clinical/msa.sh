#!/usr/bin/env bash
# MSAs of the PG2 protocol: ColabFold search on UniRef100 (2026_04) with no post-alignment filtering.
#
#   bash msa.sh queries.fasta out_dir      ->  out_dir/<ID>.a3m per query (query first, header = its ID)
#
# The monomer search of colabfold_search (colabfold 1.5.5, mmseqs_search_monomer) with filter=0, run straight on
# UniRef100: no expandaln (UniRef100 has no clusters to expand), no environmental database, and --max-seqs
# 1,000,000 instead of 10,000 (ColabFold gets its depth from the cluster expansion).
#
# Environment: MMSEQS (MMseqs2 >= 18), UNIREF100 (MMseqs2 database), THREADS, GPU=1 for the GPU prefilter, and
# for the CPU search SPLIT_MEMORY (memory for the k-mer index, e.g. 300G; default: all of the machine).
# One-time database (current UniRef100 release; here 2026_04, 227M sequences):
#   mmseqs databases UniRef100 uniref100 tmp --gpu 1    # --gpu 1 writes the GPU-ready (padded) format; drop it for CPU
set -euo pipefail
Q=$1 OUT=$2
M=${MMSEQS:-mmseqs} DB=$UNIREF100 T=${THREADS:-$(nproc)}
W=${MSA_TMP:-$OUT/tmp_$(basename "$Q" .fasta)}
if [ "${GPU:-0}" = 1 ]; then PREFILTER="--gpu 1 --prefilter-mode 1"; else PREFILTER="-s 8 --split-memory-limit ${SPLIT_MEMORY:-0}"; fi
rm -rf "$W"; mkdir -p "$OUT" "$W"

$M createdb "$Q" "$W/qdb" --shuffle 0
$M search "$W/qdb" "$DB" "$W/res" "$W/tmp" --threads "$T" --num-iterations 3 --db-load-mode 0 -a -e 0.1 \
  --max-seqs 1000000 $PREFILTER
$M mvdb "$W/tmp/latest/profile_1" "$W/prof_res"
$M lndb "$W/qdb_h" "$W/prof_res_h"
$M align "$W/prof_res" "$DB" "$W/res" "$W/res_realign" --threads "$T" --db-load-mode 0 -e 10 \
  --max-accept 1000000 --alt-ali 10 -a
$M filterresult "$W/qdb" "$DB" "$W/res_realign" "$W/res_filt" --threads "$T" --db-load-mode 0 --qid 0 --qsc -20 \
  --diff 0 --max-seq-id 1.0 --filter-min-enable 100
$M result2msa "$W/qdb" "$DB" "$W/res_filt" "$W/msa" --threads "$T" --db-load-mode 0 --msa-format-mode 6 \
  --filter-msa 0 --filter-min-enable 1000 --diff 3000 --qid 0.0,0.2,0.4,0.6,0.8,1.0 --qsc 0 --max-seq-id 0.95
$M unpackdb "$W/msa" "$W/a3m" --unpack-name-mode 0 --unpack-suffix .a3m
$M rmdb "$W/msa"

i=0                                    # unpackdb numbers the files in query order
for id in $(sed -n 's/^>\([^ ]*\).*/\1/p' "$Q"); do
  tr -d '\0' < "$W/a3m/$i.a3m" | sed "1s/.*/>$id/" > "$OUT/$id.a3m.part"
  mv "$OUT/$id.a3m.part" "$OUT/$id.a3m"
  rm "$W/a3m/$i.a3m"
  i=$((i + 1))
done
rm -rf "$W"
echo "MSA OK: $i a3m in $OUT"
