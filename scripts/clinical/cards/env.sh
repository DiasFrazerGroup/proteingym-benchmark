# Shared settings for the card archive jobs (sourced by the .sbatch files and submit.sh).
export WORK=${WORK:-/no_backup/diasfrazer/cpugh/pgb_clinical_cards}
export CURATION=${CURATION:-/users/diasfrazer/cpugh/projects/pg2_clinical/clinical-curation}
export PY=${PY:-/users/diasfrazer/egalpern/project/proteingym_clinical/env/bin/python}
# MMseqs2 18 + UniRef100 2026_01 MMseqs2 DB built for the Genebass MSAs (read-only, group access)
export MMSEQS=${MMSEQS:-/no_backup/diasfrazer/egalpern_pgb_clinical/mmseqs18/bin/mmseqs}
export DB=${DB:-/no_backup/diasfrazer/egalpern_uniref100_2026_01_mmseqs}
# MSAs: ours first, then the Genebass ones (reused when their first record is our wild type)
export MSA_DIRS="${MSA_DIRS:-$WORK/msas /no_backup/diasfrazer/egalpern_pgb_clinical/genebass_msas}"
