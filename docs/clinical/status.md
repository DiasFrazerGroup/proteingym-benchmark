# Clinical benchmark: status

Snapshot: 2026-11-16.

## Repository

| | |
|---|---|
| Fork | [DiasFrazerGroup/proteingym-benchmark](https://github.com/DiasFrazerGroup/proteingym-benchmark) (fork of ProteinGym/proteingym-benchmark) |
| Base branch | `docs/clinical-benchmark-design`: design notes `docs/clinical_benchmark.md`, by the model-API collaborator |
| Our branch | `clinical/genebass-data`, local clone `project/pgb_clinical`. Pushed to the fork (2026-11-16). |
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
- [x] **Annotation layers** for any archive: AlphaMissense (its hg38 table), popEVE (v1.1 release zip) and CADD (dbNSFP 5.3.1a); pLDDT, full length and
  UniProt positions. Validated against the paper, AlphaMissense's own tables and raw AFDB files.
- [x] **Final benchmark set (2026-11-19):** proteins <= 2,700 aa (100 pairs, 63 proteins, 6,066 rows), archives
  with the EVE-protocol filtered alignment (model side's `msa_filter.py`, vendored verbatim): `$WORK/genebass/final_le2700/`
  (868 MB), index `datasets/clinical/genebass/genebass_le2700.csv`, filter table `genebass_le2700_msa_filter.tsv`.
  Built by `scripts/clinical/msa/filter_msas.py` (job 32028266), checked by `check_final.py` (job 32028275,
  `CHECK FINAL OK`). See [genebass.md](genebass.md#final-dataset--2700-aa-filtered-alignments).
- [x] **Reports:** [genebass.md](genebass.md), [long_proteins.md](long_proteins.md), this file. The project overview is in
  [project.md](project.md).

## Running

Nothing on the data side (final set built 2026-11-19, see Done). Model work (long proteins) is on branch `clinical/model-api`, documented in
`docs/clinical/models_long_proteins.md` there.

Finished 2026-11-16:
- Rebuild 31919891 with all 66 MSAs and the `in_vesm_paper` flags: `paper set (in_vesm_paper): 103 of 109 archives,
  6,977 rows (paper: 103 pairs, 6,977 rows)`, all-rows check MISMATCH (expected), `--paper-set` check **REPRODUCED**,
  `REBUILD OK` (2 h 06 min).
- Layers with the new sources (branch `clinical/layers-sources`, merged here): job 31922391 on the final archives,
  `LAYERS OK` (42 min). AlphaMissense 100% of rows, popEVE 73.6%, CADD 99.9% (see
  [layers_sources_progress.md](layers_sources_progress.md)).
- `$WORK/genebass/splits/` (63 GB, uncompressed archives without layers) deleted after the checks, as agreed.
  The final archives are `$WORK/genebass/splits_layers/`; to rebuild them, rerun `scripts/clinical/msa/rebuild.sbatch`.

Earlier jobs: rebuild 31809884 and 31810199, MSAs for the 5 new proteins 31810198 (`SEARCH OK`, 66 `.a3m`).

## Pending

1. **Hosting:** files stay local on the cluster for now (user, 2026-11-16). Licence handling for the score layer
   (AlphaMissense CC BY-NC-SA 4.0, CADD non-commercial) before any public release.
2. **Clinical evaluation stage in the repo:** per-pair Pearson and −log10 P, overall and per group, plus Spearman
   and per-gene variants; headline on all archives and rows, with a switch for the paper's set (`in_vesm_paper`).
   It needs model predictions, which in turn need the model API for long proteins.
4. **Later tasks:** ClinVar (needs proteome scale), BRaVa (feasibility analysed; variant-level files not downloaded
   yet), de novo variants.
5. **proteingym-base issue to report:** archive structures lose pLDDT (B-factors set to 0).
6. **Open from earlier work:** the dark_energy ESM-IF context-chain bug (`esmif_core.py:214-244`), verified but not
   fixed.

## Decisions taken

| Decision | Where from |
|---|---|
| Work in the group fork, on our own branch; never push without asking | User, 2026-11-06 |
| Genebass first, exactly the VESM Fig. 6b setting, with the 3 biomarker groups as the paper's panel (the "Other" pairs only in the overall mean) | User; paper notebook |
| One archive per gene–phenotype pair, PG2 DMS style, with a trivial split | Ours, following the PG2 conventions |
| Association data from Genebass itself, not the paper's extraction; the paper's tables only for groups and comparison | User, 2026-11-09 |
| Headline metric on all 109 pairs and all rows; the paper's exact set recoverable from `in_vesm_paper` flags on pairs (index) and rows (assay) | User, 2026-11-16 |
| Rows that can't be mapped to a UniProt wild type are dropped and listed | Ours (an archive needs a wild-type position) |
| MSA: ColabFold on UniRef100 Jan 2026; ignore PG2's jackhmmer MSA | User, 2026-11-09 |
| Final set capped at 2,700 aa (one full-length AFDB structure per protein; APOB, LRP2, PKD1 dropped) | User, 2026-11-19 |
| Final archives store the EVE-protocol filtered alignment (same code as the models' filter, vendored verbatim), not the full ColabFold search | User, 2026-11-19 (replaces "no post-alignment filtering" for the final set) |
| Index `archive` column points at the layered archives (`splits_layers/`); the final index at `final_le2700/` | Ours, 2026-11-19 (bug fix) |
| ColabFold recipe adapted to UniRef100: no cluster expansion, no environmental DB, `--max-seqs 1e6` | Ours (inferred from PG2 depths), to confirm with the PG2 builders |
| Archives keep the ColabFold a3m headers (bit scores for the EVE-protocol filter); `restore_msa_headers.py` after the layer step | User, 2026-11-17 |
| Structures: AFDB v6, AFDB fragments for proteins > 2,700 aa; no AF2 re-folding | User, 2026-11-16 (replaces "PG2 protocol, AF2") |
| PrimateAI-3D: not for now | User, 2026-11-16 |
| Push `clinical/genebass-data` to the fork | User, 2026-11-16 |
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

0. **Final set (current):** `$WORK/genebass/final_le2700/` + `datasets/clinical/genebass/genebass_le2700.csv`. To
   rebuild it from `splits_layers/` alone: `sbatch scripts/clinical/msa/filter_msas.sbatch` (1.5 min), then the
   `check_final.py` command in [genebass.md](genebass.md#final-dataset--2700-aa-filtered-alignments) (must print
   `CHECK FINAL OK`). A full rebuild (`rebuild.sbatch`) now ends with these two steps. The steps below are the history
   of the full 109-pair build.
1. **Check the rebuild (job 31919891).**
   - `squeue -u $USER`, then `grep -E "MISMATCH|REPRODUCED|REBUILD OK|paper set" $WORK/logs/rebuild_31919891.log`.
   - Success: 109 archives in `$WORK/genebass/splits/` and `splits_layers/`, the line `paper set (in_vesm_paper):
     103 of 109 archives, 6,977 rows`, then **MISMATCH** (all rows; expected, Genebass has more pairs and rows than
     the paper), then **REPRODUCED** (the `--paper-set` check, using only the archives' flags), then `REBUILD OK`.
   - To rerun the exact check alone:

     ```shell
     sbatch -p genoa64 --qos=short -t 01:00:00 -c 2 --mem=64G --chdir=<repo> --wrap \
       "$PY scripts/clinical/check_genebass.py --fig6 /no_backup/diasfrazer/db/vesm_djzn25_figdata/extracted/data/Fig6_data \
        --out $WORK/genebass --paper-set"
     ```
   - If the job ran out of memory or time: raise `--mem` or `--time` and resubmit `scripts/clinical/msa/rebuild.sbatch`
     with REPO, WORK and PY exported (as in `scripts/clinical/msa/pipeline.sh`).
2. **Commit and push the index.** The job copies the new index into `datasets/clinical/genebass/genebass_fig6b.csv`;
   commit it on `clinical/genebass-data` and push to the fork.
3. **Free space:** delete `$WORK/genebass/splits/` (63 GB) once step 1 passed.
4. **Open decision** (ask the user): hosting of the archives and licence handling for the score layer.
5. **Branch discipline** (see [project.md](project.md#5-organisation-of-the-work)): data work only on
   `clinical/*-data` branches (`datasets/clinical/`, `scripts/clinical/`, `docs/clinical/`); model or evaluation
   changes go on their own branches. Don't push without the user's OK.
