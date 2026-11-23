# Shared settings for the ClinVar / de novo / gnomAD benchmark jobs (sourced by the .sbatch files and submit.sh).
export WORK=${WORK:-/no_backup/diasfrazer/cpugh/pgb_clinical_benchmarks}
export CURATION=${CURATION:-/users/diasfrazer/cpugh/projects/pg2_clinical/clinical-curation}
export PY=${PY:-/users/diasfrazer/egalpern/project/proteingym_clinical/env/bin/python}
# MMseqs2 18 + UniRef100 2026_01 MMseqs2 DB built for the Genebass MSAs (read-only, group access)
export MMSEQS=${MMSEQS:-/no_backup/diasfrazer/egalpern_pgb_clinical/mmseqs18/bin/mmseqs}
# (not $DB: user profiles on the cluster set DB to the shared db folder)
export UNIREF100_DB=${UNIREF100_DB:-/no_backup/diasfrazer/egalpern_uniref100_2026_01_mmseqs/uniref100}
# MSAs: ours first, then the Genebass ones (reused when their first record is our wild type)
export MSA_DIRS="${MSA_DIRS:-$WORK/msas /no_backup/diasfrazer/egalpern_pgb_clinical/genebass_msas}"
S=scripts/clinical/benchmarks
OUT=$WORK/benchmarks
