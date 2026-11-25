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
- **UniRef100 database:** `mmseqs databases UniRef100 uniref100 tmp --gpu 1` (downloads the current release; the
  archives use 2026_04). `--gpu 1` writes the GPU-ready format; drop it for a CPU search.

## MSA search on a GPU

`GPU=1 bash scripts/clinical/msa.sh queries.fasta out_dir` uses MMseqs2's GPU prefilter. It needs MMseqs2 18
built with CUDA and the GPU-ready database. MMseqs2 spreads the database over the GPUs in `CUDA_VISIBLE_DEVICES`.
The GPU prefilter always runs at full sensitivity, so `-s` isn't passed.

After the GPU prefilter, the alignment of one protein runs on one CPU thread, and its time grows with the square of
the protein length: about 8 minutes per pass at 1,000 aa, four passes. Two things keep it fast:

- **Put the database on local disk.** The alignment reads target sequences at random, which is slow from a
  parallel filesystem.
- **Use many proteins per run and send long proteins to the CPU search.** The archives used the GPU up to
  1,250 aa and the CPU search (`GPU=0`, `SPLIT_MEMORY=300G`) above.

To split the work, cut `queries.fasta` into chunks and run `msa.sh` on each; every chunk writes its own
`<UNIPROT>.a3m` files.
