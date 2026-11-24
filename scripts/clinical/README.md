# Clinical benchmark scripts

These build the four clinical benchmarks (`clinvar`, `denovo`, `gnomad`, `genebass`) as PG2 archives. What the
benchmarks are: [docs/clinical/benchmarks.md](../../docs/clinical/benchmarks.md).

## Run

From the repo root, with the paths in [env.sh](env.sh) (override any of them in the environment):

```shell
bash scripts/clinical/run.sh prepare
bash scripts/clinical/run.sh msa
bash scripts/clinical/run.sh build
bash scripts/clinical/run.sh check
```

`BENCHMARKS="clinvar genebass"` limits a step to some benchmarks. `JOBS` sets the number of parallel archive
builds.

| Step | Script | Reads | Writes |
|---|---|---|---|
| `prepare` | [prepare.py](prepare.py), [queries.py](queries.py) | clinical-curation tables, Genebass tables, VESM Fig. 6 data, Swiss-Prot, SNV index | `$WORK/<benchmark>/variants.tsv.gz`, `$WORK/queries.fasta` |
| `msa` | [msa.sh](msa.sh) | `$WORK/queries.fasta`, UniRef100 MMseqs2 database | `$MSAS/<UNIPROT>.a3m` |
| `build` | [build.py](build.py) | variants, Swiss-Prot, AlphaFold DB models, MSAs | `$WORK/<benchmark>/splits/*.splits.pgdata`, `datasets/clinical/<benchmark>/<benchmark>.csv` |
| `check` | [check.py](check.py) | variants, archives, index | pass / fail, `$WORK/reference_results.tsv` |

The steps are plain commands. On a cluster, submit them however you like, for example
`sbatch -c 64 --wrap "JOBS=60 bash scripts/clinical/run.sh build"`.

## One-time inputs

- **SNV index:** `bash scripts/clinical/run.sh index` runs [build_snv_index.py](build_snv_index.py) on the
  AlphaMissense hg38 and hg19 tables (10 minutes, 40 GB of memory).
- **Genebass tables:** [genebass_api.py](genebass_api.py) downloads them from Genebass.
- **UniRef100 database:** `mmseqs createdb uniref100.fasta uniref100` (UniProt 2026_01). For the GPU search,
  also `mmseqs makepaddedseqdb uniref100 uniref100_gpu`, and point `UNIREF100` at `uniref100_gpu`.

## MSA search on a GPU

`GPU=1 bash scripts/clinical/msa.sh queries.fasta out_dir` uses MMseqs2's GPU prefilter. It needs MMseqs2 18
built with CUDA and the padded database. MMseqs2 spreads the database over the GPUs in `CUDA_VISIBLE_DEVICES`.
The GPU prefilter always runs at full sensitivity, so `-s` isn't passed.

To split the work, cut `queries.fasta` into chunks and run `msa.sh` on each; every chunk writes its own
`<UNIPROT>.a3m` files.
