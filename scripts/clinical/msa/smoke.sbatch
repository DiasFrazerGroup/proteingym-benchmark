#!/usr/bin/env bash
# Smoke test of the search recipe: Swiss-Prot human as a mini target DB, 2 queries, then the builder's a3m check.
set -euo pipefail
T="$WORK/smoke"; rm -rf "$T"; mkdir -p "$T"
"$PY" - "$T" <<'EOF'
import sys
sys.path.insert(0, "scripts/clinical")
from build_genebass import read_sprot
t = sys.argv[1]
sp, _ = read_sprot("/no_backup/diasfrazer/db/uniprot_sprot_human/uniprot_sprot_human.dat.gz")
open(f"{t}/sprot.fasta", "w").write("".join(f">{k}\n{v['sequence']}\n" for k, v in sp.items()))
open(f"{t}/q.fasta", "w").write("".join(f">{k}\n{sp[k]['sequence']}\n" for k in ["P35504", "P02647"]))
EOF
"$MMSEQS" createdb "$T/sprot.fasta" "$T/sprotdb" > /dev/null
"$PY" scripts/clinical/msa/colabfold_search_uniref100.py "$T/q.fasta" "$T/sprotdb" "$T/out" --mmseqs "$MMSEQS" \
  --threads "$SLURM_CPUS_PER_TASK" --max-seqs 1000
"$PY" - "$T" <<'EOF'
import sys
sys.path.insert(0, "scripts/clinical")
from build_genebass import a3m_check, read_sprot
t = sys.argv[1]
seqs = {l[1:].strip(): s.strip() for l, s in zip(*[iter(open(f"{t}/q.fasta"))] * 2)}
for u, wt in seqs.items():
    print(u, a3m_check(f"{t}/out/{u}.a3m", wt))
print("SMOKE OK")
EOF
