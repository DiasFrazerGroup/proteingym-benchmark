#!/usr/bin/env bash
# Add the annotation layers to every archive (*.pgdata, *.pgarchive) in IN_DIR, writing OUT_DIR. Cluster job; submit
# from the repo root (or pass --chdir=<repo>). Builds the AlphaMissense parquets and the popEVE sequence index first if
# they do not exist yet.
# Usage: sbatch -p genoa64 --qos=short -t 04:00:00 -c 4 --mem=48G \
#          scripts/clinical/layers/layers.sbatch IN_DIR OUT_DIR [extra add_layers.py args]
set -euo pipefail
shopt -s nullglob
IN=$1; OUT=$2; shift 2
PY=${PY:-/users/diasfrazer/egalpern/project/proteingym_clinical/env/bin/python}
D=scripts/clinical/layers
IDX=$("$PY" -c "import tomllib;c=tomllib.load(open('$D/sources.toml','rb'))['alphamissense'];print(c['table'],c['protein_table'],c.get('isoform_table',''))")
for f in $IDX; do [ -f "$f" ] || { "$PY" "$D/build_snv_index.py" --sources "$D/sources.toml"; break; }; done
PSEQ=$("$PY" -c "import tomllib;print(tomllib.load(open('$D/sources.toml','rb'))['popeve'].get('sequences',''))")
[ -z "$PSEQ" ] || [ -f "$PSEQ" ] || (cd "$D" && "$PY" build_popeve_index.py --workers "${SLURM_CPUS_PER_TASK:-4}")
files=("$IN"/*.pgdata "$IN"/*.pgarchive)
[ ${#files[@]} -gt 0 ] || { echo "no archives in $IN"; exit 1; }
"$PY" "$D/add_layers.py" "${files[@]}" --out "$OUT" --sources "$D/sources.toml" "$@"
echo "LAYERS OK"
