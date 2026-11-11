# ProteinGym clinical benchmark: project overview

Status: 2026-11-11. Fork: [DiasFrazerGroup/proteingym-benchmark](https://github.com/DiasFrazerGroup/proteingym-benchmark).
Related reports: [genebass.md](genebass.md), [long_proteins.md](long_proteins.md), [status.md](status.md), and the
design notes [clinical_benchmark.md](../clinical_benchmark.md).

## 1. Motivation

ProteinGym2 benchmarks a zoo of variant-effect models on deep mutational scans (DMS). Clinical genetics asks a
related but different question: do model scores separate, or track, variants by their effect on human disease and
traits? Answering it with the same models needs:

- data from biobanks and clinical databases rather than DMS;
- full-length human proteins, often much longer than DMS constructs;
- sparse variants, a handful to a few hundred per gene;
- metrics that work across genes, not only within one assay.

The aim is to add clinical tasks to ProteinGym2 so that the existing zoo runs **unchanged**, in the **PG2 data
format**, and the result could later be adopted by the ProteinGym2 maintainers.

## 2. Goals

1. **Clinical tasks:**
   - **Biobank association statistics:** continuous per-variant effects on traits. Genebass first, BRaVa later.
   - **ClinVar:** pathogenic vs benign.
   - **De novo variants:** severe developmental disorder cases vs controls.
2. **Same zoo, same format:**
   - standard PG2 archives (`.splits.pgdata`);
   - wild types, structures and MSAs built with PG2's own protocols;
   - no clinical-specific model inputs.
3. **The zoo working on human clinical proteins:** long proteins (windows for language models, fragments for
   structures, MSA policies) and sparse variants (score only the positions needed).
4. **Clinical evaluation:**
   - metrics pooled across genes and per gene;
   - each task's reference metric reproduced exactly (e.g. VESM Fig. 6b for Genebass);
   - masks for low-confidence sites (pLDDT), short proteins and similar.
5. **Comparison layers:** published predictors (AlphaMissense, popEVE, CADD, PrimateAI-3D) as reference columns on
   every variant, without running those models.
6. **Reproducible and shareable:**
   - builders that read primary sources (e.g. Genebass itself, not a paper's extract);
   - validation scripts;
   - cluster jobs that do not depend on an interactive session;
   - hosting on Hugging Face like the PG2 DMS collection, with licences checked per data layer.

## 3. Tasks and their data sources

| Task | Label per variant | Source | Evaluation unit | Status |
|---|---|---|---|---|
| Genebass | single-variant `beta` (continuous) | Genebass API (UK Biobank, 394,841 exomes); VESM Fig. 6b rules | gene–phenotype pair; mean −log10 P (Pearson), overall and per biomarker group | **data built** ([genebass.md](genebass.md)) |
| BRaVa | single-variant effect (continuous) | BRaVa consortium release (multi-biobank, about 1.2M people) | as Genebass | feasibility analysed; variant files not downloaded |
| ClinVar | 1 pathogenic / 0 benign | ClinVar (review-status filtered) | pooled AUROC/AUPRC; per gene where both classes exist | not started (proteome scale) |
| De novo variants | 1 case / 0 control | SDD cohorts vs controls | pooled AUROC with a gene-matched comparison | not started |

The main pitfalls per task are listed in [clinical_benchmark.md](../clinical_benchmark.md#clinical-sub-tasks):

- noisy single-variant statistics for rare variants;
- ClinVar circularity with predictors trained on it;
- gene-level confounding in ClinVar and de novo variants.

## 4. Design

**Archives.** One `.splits.pgdata` per evaluation unit:

- per gene–phenotype pair for biobank tasks, like PG2's one archive per (protein, condition);
- per gene for ClinVar and de novo variants.

Each archive holds the full-length wild type, PG2-protocol inputs, and a trivial split (zero-shot models ignore it).
Clinical archives live apart from the DMS ones (`datasets/clinical/`), so the DMS games never pick them up.

**Inputs, following PG2's protocols:**

- **Wild type:** canonical UniProt.
- **MSA:** ColabFold on UniRef100 (January 2026), no post-alignment filtering. PG2's second, jackhmmer MSA is left
  out.
- **Structure:** AlphaFold2 with default settings, with a fragment rule for proteins too long to fold in one piece.
  The current archives use AlphaFold DB v6 until this is done.

**Variant identity.** Every row carries:

- `mutant` (UniProt numbering);
- the genomic variant where one exists (GRCh38);
- `hgvsp`;
- `mapping`, which records how the variant was placed on the protein.

These let clinical labels, DMS scores and model predictions be joined per variant across benchmarks.

**Layers.** Extra columns on every row, never model inputs or targets, added by one script for any archive
(`scripts/clinical/layers`):

- **Scores:** AlphaMissense, popEVE, CADD and PrimateAI-3D, with source and version recorded.
- **Masks:** pLDDT of the site and the full canonical protein length, used to restrict an analysis (e.g. pLDDT > 70,
  length ≥ 100).

**Evaluation.** A clinical stage beside the existing DVC games:

- it scores each (protein, model) once;
- it evaluates every target from those predictions;
- it computes pooled and per-unit metrics, dropping missing values per target.

## 5. Organisation of the work

The split of tasks between people is not fixed. What is fixed is that **data and model changes never share a
branch**, so each can be reviewed, tested and merged on its own.

| Kind of change | Touches | Branch (proposed) |
|---|---|---|
| Datasets: builders, layers, MSA/structure pipelines, dataset docs | `datasets/clinical/`, `scripts/clinical/`, `docs/clinical/` | `clinical/genebass-data` (now); one `clinical/<task>-data` per new task |
| Model API: windows, fragments, sparse scoring, MSA region policy | `models/`, evedesign subclasses | `clinical/model-api` (or one per model) |
| Evaluation stage: clinical DVC game, pooled metrics | `benchmark/clinical/`, metric code | `clinical/evaluation` |
| Design notes | `docs/clinical_benchmark.md` | `docs/clinical-benchmark-design` (base) |

Rules:

1. A data branch never edits `models/` or `benchmark/`, and a model branch never edits `datasets/` or
   `scripts/clinical/`. Today `clinical/genebass-data` touches only `datasets/clinical/`, `scripts/clinical/` and
   `docs/clinical/`.
2. Branches meet only through the archive format. Models read standard archives, and data builders don't assume any
   model.
3. Model changes are validated against the zoo on short proteins (identical scores) before any clinical run (see
   [long_proteins.md](long_proteins.md)).
4. Archives are generated, not committed (Git LFS / Hugging Face). Branches carry code, index files and docs.

## 6. Milestones

1. **Genebass data (done, final rebuild running):** 109 archives from Genebass, reproducing VESM Fig. 6b on the
   paper's rows; MSAs; layers.
2. **Genebass inputs complete:** MSAs for the 5 added proteins; AF2 (PG2 protocol) structures with a long-protein
   fragment rule.
3. **Model API for long proteins and sparse variants:** all zoo models scoring the Genebass archives, with
   validations.
4. **Clinical evaluation stage:** Fig. 6b metrics in the repo (Pearson; plus Spearman and per-gene variants), first
   Genebass leaderboard.
5. **BRaVa:** same design on the BRaVa release.
6. **ClinVar and de novo variants:** proteome-scale data builders, pooled metrics, gene-level controls.
7. **Release:** Hugging Face hosting and a dataset card; a proposal to ProteinGym2.

## 7. Open decisions

- **Genebass:** whether the 6 Genebass-only pairs enter the headline metric or stay a flagged extension.
- **ClinVar:** curation (review stars, conflicting interpretations) and avoiding circularity with predictors trained
  on ClinVar.
- **De novo variants:** which case and control cohorts; handling of variants seen in both.
- **Metrics:** Pearson (paper) vs Spearman; per pair vs per gene; weighted alternatives (inverse-variance).
- **Structures:** the exact AF2 setup and the fragment rule above about 2,700 aa.
- **MSAs:** confirming the ColabFold-on-UniRef100 recipe with the PG2 builders; region policy for long proteins.
- **Hosting and licences:** private vs public Hugging Face; how to ship score layers with non-commercial or
  no-redistribution terms (PrimateAI-3D).
