# Paths for the clinical benchmark build (sourced by run.sh; override any of them in the environment).
export WORK=${WORK:-/no_backup/diasfrazer/cpugh/pgb_clinical}          # prepared tables, MSAs, archives
export PY=${PY:-/users/diasfrazer/egalpern/project/proteingym_clinical/env/bin/python}   # pandas, pyarrow, proteingym-base
export JOBS=${JOBS:-16}                                                 # parallel archive builds
DB=/no_backup/diasfrazer/db
# inputs
export CURATION=${CURATION:-/users/diasfrazer/cpugh/projects/pg2_clinical/clinical-curation}  # ClinVar / de novo / gnomAD tables
export FIG6=${FIG6:-$DB/vesm_djzn25_figdata/extracted/data/Fig6_data}   # VESM paper Fig. 6 data (figshare 31162471)
export GENEBASS=${GENEBASS:-/no_backup/diasfrazer/egalpern_pgb_clinical/genebass_api}      # genebass_api.py output
export SPROT=${SPROT:-$DB/uniprot_sprot_human/uniprot_sprot_human.dat.gz}                  # Swiss-Prot human 2026_03
export AFDB=${AFDB:-/users/diasfrazer/egalpern/project/human_proteome_de/data/structures}  # AF-<ID>-F<k>-model_v6.pdb
export AM_RAW=${AM_RAW:-$DB/alphamissense/raw}                          # AlphaMissense_hg38.tsv.gz, AlphaMissense_hg19.tsv.gz
export SNV_INDEX=${SNV_INDEX:-$WORK/snv_index}                          # am_hg38 / am_hg19 parquet (build_snv_index.py)
# MSAs: one <UNIPROT>.a3m per protein
export MSAS=${MSAS:-$WORK/msas}
export MMSEQS=${MMSEQS:-mmseqs}                                         # MMseqs2 >= 18
export UNIREF100=${UNIREF100:-$WORK/uniref100_2026_04/uniref100}   # MMseqs2 DB, UniRef100 2026_04 (see msa.sh)
export GPU=${GPU:-0}                                                    # 1: GPU search (database made with --gpu 1)
