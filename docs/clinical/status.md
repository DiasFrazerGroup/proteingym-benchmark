# Clinical benchmark: status

Snapshot: 2026-11-11.

## Repository

| | |
|---|---|
| Fork | [DiasFrazerGroup/proteingym-benchmark](https://github.com/DiasFrazerGroup/proteingym-benchmark) (fork of ProteinGym/proteingym-benchmark) |
| Base branch | `docs/clinical-benchmark-design`: design notes `docs/clinical_benchmark.md`, by the model-API collaborator |
| Our branch | `clinical/genebass-data`, local clone `project/pgb_clinical`. **Not pushed yet.** |
| Env | `proteingym_clinical/env` (Python 3.12, proteingym-base `main` ada3c78); cluster jobs on genoa64 |
| Work dir | `/no_backup/diasfrazer/egalpern_pgb_clinical` (logs and `logs/jobs.txt`, MSAs, archives, MMseqs2 v18 copy, tabix) |

Commits on our branch, oldest first:

- Genebass datasets (VESM Fig. 6b) — builder, check, index.
- ColabFold/UniRef100 MSAs as a SLURM chain.
- Annotation layers: scores and masks.
- Genebass tables straight from the Genebass API, plus the fixes found while running it (phenotype filter, pair
  rule, UGT1A-locus annotation, SNV-based mapping, check robustness).
- These reports, in `docs/clinical/`.

The earlier prototype (parallel folder `project/proteingym_clinical/`) is kept aside for reference: windowed zoo
models, structure-model runs, correlation variants and the Spearman/per-gene analyses.

## Done

- [x] **Genebass archives, Fig. 6b setting:** builder, Fig. 6b reproduction check and index.
  - Paper-sourced version: 103 archives, reproduces the paper's printed numbers exactly.
  - Genebass-sourced version: 109 archives (see [genebass.md](genebass.md)).
- [x] **Association data straight from Genebass:** all of the paper's rows are identical in Genebass. Pair rule
  recovered: portal top pLoF list on the paper's 27 phenotypes.
- [x] **MSAs, PG2 protocol** (ColabFold on UniRef100 2026_01, no post-alignment filtering) for the 61 paper
  proteins: 4 h 24 min on 64 CPUs, 38 GB.
- [x] **Annotation layers** for any archive: AlphaMissense, popEVE and CADD (dbNSFP 5.3.1a); pLDDT, full length and
  UniProt positions. Validated against the paper, AlphaMissense's own tables and raw AFDB files.
- [x] **Reports:** [genebass.md](genebass.md), [long_proteins.md](long_proteins.md), this file. The project overview is in
  [project.md](project.md).

## Running

| Job | ID | What | Expected |
|---|---|---|---|
| `gb_rebuild` | 31809884 (genoa64, 8 CPU, 200 GB, ≤12 h) | Final Genebass archives: Genebass tables, MSAs, layers → `$WORK/genebass/splits_layers/` | Several hours (archives with multi-GB MSAs). The log ends `REBUILD OK`. |
| `msa_new5` | 31810198 (genoa64, 64 CPU, 400 GB, ≤6 h) | MSAs for CST3, GPLD1, GPT, SHBG, SLC22A12 (`QUERIES=$WORK/genebass_queries_new5.fasta`, wild types taken from the built archives) → `$WORK/genebass_msas/` | 1–2 h (database scan). The log ends `SEARCH OK`, 66 `.a3m`. |
| `gb_rebuild` | 31810199 (same as 31809884; after `msa_new5` OK and after 31809884 ends) | Rebuild with all 66 MSAs; supersedes 31809884's archives | Several hours; then run the `--paper-rows-only` check. |

## Pending

1. **MSAs for the 5 new proteins** (CST3, GPLD1, GPT, SHBG, SLC22A12). A short search job, then rebuild those
   archives.
2. **Structures, PG2 protocol:** AF2 with default settings. Needs the exact AF2 setup and a rule for the proteins
   over 2,700 aa (APOB, LRP2, PKD1). The archives currently use AFDB v6.
3. **PrimateAI-3D:** needs the Illumina academic licence and the file. Then set `[primateai3d] path` in
   `scripts/clinical/layers/sources.toml` and rerun the layers.
4. **Clinical evaluation stage in the repo:** per-pair Pearson and −log10 P, overall and per group, plus Spearman
   and per-gene variants. It needs model predictions, which in turn need the model API for long proteins.
5. **Sharing:** push the branch; host the archives (Hugging Face vs cluster only); licence handling for the score
   layer.
6. **Later tasks:** ClinVar (needs proteome scale), BRaVa (feasibility analysed; variant-level files not downloaded
   yet), de novo variants.
7. **proteingym-base issue to report:** archive structures lose pLDDT (B-factors set to 0).
8. **Open from earlier work:** the dark_energy ESM-IF context-chain bug (`esmif_core.py:214-244`), verified but not
   fixed.

## Decisions taken

| Decision | Where from |
|---|---|
| Work in the group fork, on our own branch; never push without asking | User, 2026-11-06 |
| Genebass first, exactly the VESM Fig. 6b setting, with the 3 biomarker groups as the paper's panel (the "Other" pairs only in the overall mean) | User; paper notebook |
| One archive per gene–phenotype pair, PG2 DMS style, with a trivial split | Ours, following the PG2 conventions |
| Association data from Genebass itself, not the paper's extraction; the paper's tables only for groups and comparison | User, 2026-11-09 |
| Keep the 6 Genebass-only pairs, flagged `in_vesm_paper`; whether they enter the headline metric is still open | Ours (pending user) |
| Rows that can't be mapped to a UniProt wild type are dropped and listed | Ours (an archive needs a wild-type position) |
| MSA: ColabFold on UniRef100 Jan 2026, no post-alignment filtering; ignore PG2's jackhmmer MSA | User, 2026-11-09 |
| ColabFold recipe adapted to UniRef100: no cluster expansion, no environmental DB, `--max-seqs 1e6` | Ours (inferred from PG2 depths), to confirm with the PG2 builders |
| Structures: PG2 protocol, AF2 default settings (not done yet) | User |
| Layers: external scores (AlphaMissense, popEVE, CADD, PrimateAI-3D) and masks (pLDDT, protein length) on every row | User, 2026-11-09 |
| pLDDT from the original AFDB files | Ours (proteingym-base zeroes B-factors) |
| No PrimateAI-3D scraping; licensed download only | Ours (licence terms) |
| Model scope for now: masked marginal, 1,022 context; alignment models deferred, now to be included | User, 2026-11-06 → 2026-11-09 |
| The work split is not fixed; data and model changes go on separate branches, never coupled (see [project.md](project.md#5-organisation-of-the-work)) | User, 2026-11-11 |
| All compute as sbatch jobs on genoa64, not in the interactive session | User, 2026-11-09 |

## Picking this up

Everything below runs on the CRG cluster as sbatch jobs on genoa64; don't compute in an interactive session.
Python: `/users/diasfrazer/egalpern/project/proteingym_clinical/env/bin/python` (proteingym-base and its
dependencies). `WORK=/no_backup/diasfrazer/egalpern_pgb_clinical`. Job IDs are appended to `$WORK/logs/jobs.txt`.

1. **Check the final rebuild (job 31809884).**
   - `squeue -u $USER`, then `tail $WORK/logs/rebuild_31809884.log`.
   - Success ends with `REBUILD OK`, with 109 archives in `$WORK/genebass/splits/` and the layered versions in
     `$WORK/genebass/splits_layers/`.
   - The log's own check prints **MISMATCH. That is expected:** the Genebass set has more pairs and rows than the
     paper. The real test is the exact one below, which must print `REPRODUCED`. It did for the archives without
     MSAs (`$WORK/genebass_api/build`, job 31810196):

     ```shell
     sbatch -p genoa64 --qos=vshort -t 00:30:00 -c 2 --mem=32G --chdir=<repo> --wrap \
       "$PY scripts/clinical/check_genebass.py --fig6 /no_backup/diasfrazer/db/vesm_djzn25_figdata/extracted/data/Fig6_data \
        --out $WORK/genebass --paper-rows-only"
     ```
   - If the job ran out of memory or time: raise `--mem` or `--time`. MSAs are up to 3.9 GB each and are parsed by
     proteingym-base. Then resubmit `scripts/clinical/msa/rebuild.sbatch` with the exports used in
     `scripts/clinical/msa/pipeline.sh` (REPO, WORK, DB, MMSEQS, PY).
2. **MSAs for the 5 proteins without one** (CST3 P01034, GPLD1 P80108, GPT P24298, SHBG P04278, SLC22A12 Q96S37):
   - write their wild types to a FASTA;
   - run `scripts/clinical/msa/search.sbatch` with that FASTA in place of `genebass_queries.fasta` (same DB, output
     to `$WORK/genebass_msas`; minutes to an hour for 5 short proteins, but the database scan dominates);
   - then rerun the rebuild.
3. **Copy the index into the repo.** After a successful rebuild, `datasets/clinical/genebass/genebass_fig6b.csv` in
   the repo is overwritten with the Genebass-sourced index (109 rows). Commit it on `clinical/genebass-data`.
4. **Open decisions** (ask the user): whether the 6 Genebass-only pairs count in the headline metric; the AF2
   structure setup; pushing the branch; hosting and licences; PrimateAI-3D.
5. **Branch discipline** (see [project.md](project.md#5-organisation-of-the-work)): data work only on
   `clinical/*-data` branches (`datasets/clinical/`, `scripts/clinical/`, `docs/clinical/`); model or evaluation
   changes go on their own branches. Don't push without the user's OK.
