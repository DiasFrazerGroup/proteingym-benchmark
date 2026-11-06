# Clinical benchmark: design notes

Status: **draft / for discussion**. This document surveys how a clinical
variant benchmark (benign vs pathogenic labels, plus correlations with
per-variant summary statistics) fits into the existing ProteinGym2 benchmark
and dataset format. Nothing here has been implemented or tested end to end
yet.

Code references point at this repo, at
[proteingym-base](https://github.com/ProteinGym/proteingym-base) (`main`), and
at [evedesign](https://github.com/evedesignbio/evedesign). The models lock
evedesign 0.0.5; the evedesign observations below were read from the 0.0.6
source, so details should be confirmed against 0.0.5.

## Contents

1. [Goals](#goals)
2. [How models score today](#how-models-score-today)
3. [Long proteins](#long-proteins)
4. [Sparse scoring](#sparse-scoring)
5. [Metrics](#metrics)
6. [Do we need proteingym-base changes?](#do-we-need-proteingym-base-changes)
7. [Datasets: the `.pgdata` format](#datasets-the-pgdata-format)
8. [Where clinical datasets live](#where-clinical-datasets-live)
9. [Multiple assays per gene and multiple benchmarks per protein](#multiple-assays-per-gene-and-multiple-benchmarks-per-protein)
10. [Open questions](#open-questions)
11. [Proposed next steps](#proposed-next-steps)

## Goals

- Score clinically observed missense variants with the existing zero-shot
  model zoo.
- Evaluate against:
  - **binary labels**: pathogenic / likely pathogenic vs benign / likely
    benign (e.g. ClinVar);
  - **summary statistics**: continuous per-variant values (e.g. biobank
    effect sizes or odds ratios). The exact statistic is still to be decided,
    see [Open questions](#open-questions).
- Follow proteingym-base standards (dataset archives, model cards, runner)
  so that every existing model runs unchanged.

Two properties of clinical data differ from DMS data:

- **Long proteins.** Many disease genes are far longer than the 1022-residue
  context of ESM-2-style models.
- **Sparse variants.** Each gene usually has a handful to a few hundred
  labelled variants rather than a near-complete saturation scan.

## How models score today

Every model runs through the shared evedesign runner (`__main__.py`, see
[models/README.md](../models/README.md)). For zero-shot models it:

1. loads a `.splits.pgdata` with `proteingym.base.Subsets`;
2. converts it with `evedesign.proteingym.dataset_to_evedesign` into
   - a single-protein `System` whose `rep` is the wild-type sequence, whose
     MSA is the first MSA in the archive and whose structures are all
     structures in the archive (`first_index` is hardcoded to `1`);
   - a test set of **full-length mutated sequences**, one per assay row
     (splits and `--test-fold` are ignored for zero-shot models);
3. calls `build(system)` and then `score(instances)`;
4. writes `<dataset>_predictions.pgdata`.

Nothing in the runner or the conversion handles length. Each model class
decides for itself.

## Long proteins

Behaviour above ~1022 residues, per model:

| Model | Behaviour on long proteins |
|---|---|
| `esm2` (`evedesign.models.esm2:ESM2`) | **Fails.** `max_seq_length = 1022` is checked in `can_model` and again per instance in `_validate_instances_and_max_length`; both raise. |
| `kermut`, `pkermut` | **Fail.** Both check `_ESM_MAX_LENGTH = 1022`. They are supervised, so are probably out of scope for clinical. |
| `esmc` | Runs. `MAX_SEQ_LENGTH = None` deliberately (rotary embeddings). Only limited by memory. |
| `saprot` | No length check found, so it attempts the full sequence. ESM-style rotary embeddings mean it won't crash, but it was trained on ≤1024 tokens, so quality beyond that is untested. Also needs a structure (Foldseek 3Di tokens). |
| `esm-if1`, `pottsmpnn`, `structuredca` | No cap, but they condition on the whole structure, so memory and time grow with length. They depend on structure coverage and quality (long AlphaFold models often have low-pLDDT regions). |
| `evmutation`, `evmutation2`, `eve`, `gemme`, `rsalor` | No hard cap, but they need an MSA whose match states align to the full `rep` (evcouplings raises on a length mismatch). Cost grows with L, steeply for coupling models (EVmutation's couplings are L²·q²). Real MSAs for long proteins often cover only a domain, and `system_from_dataset` has an open TODO for "MSA target region ≠ full sequence". |
| `blosum62` | Fine at any length. |

There is **no sliding-window scoring** (as used by ProteinGym v1 for ESM-1v,
Tranception and others). The models either refuse long sequences or run
them in one pass.

### What is needed

- **Windowed scoring** for fixed-context PLMs (ESM-2 and similar): a window
  of ≤1022 residues centred on the mutated position, scored with masked
  marginals. The cleanest place for this is an evedesign subclass or wrapper,
  following the `esmc` pattern of overriding `max_seq_length` and validation.
- **A policy for MSA focus regions**: decide whether variants outside the
  MSA-covered region are dropped or reported as missing (NaN) predictions.
  This also requires fixing the `first_index=1` assumption.
- **A policy for structures** of long or disordered proteins (e.g. a pLDDT
  threshold, or masking low-confidence regions for SaProt).
- **A decision on which models must support long proteins** for the clinical
  leaderboard, and how missing predictions are handled in the metrics.

## Sparse scoring

### Models: works today, and is cheap

Scoring operates on whatever list of sequences the dataset contains, so a
few variants per protein is fine. Most models only do work for the
positions that are actually mutated:

- `ESM2.score_mutants` / `esmc` and `saprot.score`: one masked forward pass
  per **unique mutated position**, not per residue of the protein.
- `gemme`, `rsalor`, `blosum62`, `structuredca`: compute a full
  profile or scan once at `build`, then look scores up. Build cost does not
  depend on how many variants there are.
- `ESM2.score` (plain pseudo-log-likelihood) and `esm-if1`: one full-length
  forward pass per variant. That's fine for few variants, but each pass is
  still full length, so long-protein limits apply.
- `kermut` / `pkermut`: precompute masked marginals at *every* position
  (L forward passes). That is wasteful for sparse data, but these are supervised.

### Metrics: the real gap

Metrics are computed **per dataset**
(`proteingym-base/src/proteingym/base/metrics.py`), and the zero-shot game
only runs `spearman` (`benchmark/zero_shot/default.yaml`). A per-protein
Spearman over a few binary labels is meaningless. See [Metrics](#metrics).

## Metrics

| | Benign vs pathogenic | Summary statistics |
|---|---|---|
| Label | binary 1/0 | continuous (effect size, OR, AF, …) |
| Metric | AUROC, AUPRC (maybe balanced accuracy) | Spearman (maybe Pearson) |
| Per protein? | Usually not: few variants, often only one class | Yes, where genes have enough variants |
| Pooled across proteins? | **Yes, the headline number** | Probably also, with per-gene numbers as a breakdown |

Points to note:

- **Cross-protein aggregation does not exist yet.** The DVC pipeline runs
  one container per (dataset, model), and `calculate_metric` /
  `aggregate_metrics` produce one row per dataset. Pooled metrics need a new
  stage that reads all prediction files of a model, joins them to the
  labels, and computes pooled AUROC/AUPRC and Spearman. It would sit
  alongside `aggregate_metrics` and `generate_metrics_csv`.
- Per-gene averages should only include genes that have both classes (for
  AUROC) or enough variants (for Spearman). Model scores are not calibrated
  across proteins, so pooled and per-gene-averaged numbers can differ a lot.
  Report both.
- **Missing labels.** A variant may have a binary label but no summary
  statistic, or the reverse. `metric_spearman` currently does not drop nulls
  before calling `spearmanr`, so the clinical metric code must filter nulls
  per target.
- **One prediction run for both targets.** The zero-shot DVC matrix is keyed
  on (dataset, target), so it runs a separate model container per target,
  even though zero-shot predictions are identical across targets. With long
  proteins and two targets that doubles the cost. The clinical game should
  score each (protein, model) once and evaluate all targets from those
  predictions.

## Do we need proteingym-base changes?

**Not necessarily.** A minimal route needs none:

1. Build one `.pgdata` per protein with numeric targets (binary encoded as
   `1`/`0`, summary statistic as a float). The archive format already
   supports this (see below).
2. Add a `benchmark/clinical/` game that reuses the model containers to
   write predictions.
3. Add a benchmark-side script that pools predictions across datasets and
   computes AUROC/AUPRC/Spearman.

Where proteingym-base *would* be involved:

- **Per-dataset metrics through the existing CLI.** `proteingym-base
  evaluate-zero-shot --selected-metrics …` only discovers `metric_*`
  functions defined inside `proteingym/base/metrics.py`
  (`_discover_metric_functions`), so a per-dataset `auroc` used through the
  CLI must be added there. [metric.md](metric.md) also documents this as the
  convention.
- **Pooled metrics** do not fit the `DatasetScoringContext` API (one dataset
  per call) anyway. They could live in the benchmark first. If the team wants
  proteingym-base to stay the single source of metric definitions, they can
  be moved there later, possibly as a new "pooled" evaluation mode.

## Datasets: the `.pgdata` format

A `.pgdata` is ProteinGym's standard **dataset archive**: one dataset, in
practice one protein (`proteingym-base/docs/dataset_archive.md`). It is a zip
built from a TOML manifest (`proteingym-base build`, or
`Dataset.from_manifest(...).dump()`):

```text
manifest.lock   # metadata: name, targets, variables, assay descriptions, provenance
assays/         # CSV(s): one row per variant → full mutated sequence + target column(s)
assays_raw/     # optional raw data
sequences/      # wild-type FASTA (becomes the model's reference `rep`)
structures/     # PDB(s), e.g. AlphaFold model (SaProt, ESM-IF1, PottsMPNN, …)
msas/           # a3m MSA + weights (EVE, EVmutation, GEMME, RSALOR, …)
```

A `.splits.pgdata` wraps a `dataset.pgdata` plus a `slices.json` of
cross-validation folds. The benchmark pipelines currently consume
**`.splits.pgdata` only**, even zero-shot ones, which ignore the folds.

### A clinical dataset

One archive per gene/protein. Example manifest:

```toml
name = "BRCA1_HUMAN_clinical"
description = "ClinVar + biobank variants for BRCA1 (P38398)"
version = "1.0.0"

[[ assay_targets ]]
name = "pathogenic"
description = "1 = pathogenic/likely pathogenic, 0 = benign/likely benign"

[[ assay_targets ]]
name = "effect_size"          # the summary statistic
unit = "beta"

[[ assays ]]
name = "BRCA1_clinical"
path = "./BRCA1/variants.csv"
sequence_alias = "mutated_sequence"
sequence_alphabet = "AA"
target_phenotype = "other"

[[ assays.targets ]]
name = "pathogenic"
alias = "clinvar_label"

[[ assays.targets ]]
name = "effect_size"
alias = "beta"

[[ assays.non_targets ]]       # carried along, not predicted
name = "hgvs_p"
[[ assays.non_targets ]]
name = "review_stars"

[[ sequences ]]
path = "./BRCA1/P38398.fasta"
type = "wild_type"
uniprot_id = "P38398"

[[ structures ]]
path = "./BRCA1/AF-P38398.pdb"
name = "P38398"

[[ msas ]]
path = "./BRCA1/msa.a3m"
format = "a3m"
```

Variants CSV:

```csv
mutated_sequence,clinvar_label,beta,hgvs_p,review_stars
MDLSALRVEEV...,1,,p.Cys61Gly,3
MDLSALRVEEV...,0,0.02,p.Ser1613Gly,2
MDLSALRVEEV...,,0.41,p.Arg1699Trp,
```

How each part is used:

- **Assay CSV**: the sparse variant list. Only the variants we have go in.
  A variant with only one label leaves the other blank.
- **Two `assay_targets`**: binary label and summary statistic side by side.
  Zero-shot models never see them; they are only used at evaluation.
- **`non_targets`**: metadata such as HGVS, review stars, source and ancestry,
  for filtering (e.g. ≥2 stars) and stratified reporting in the pooled
  evaluation.
- **Wild-type sequence, structure and MSA**: model inputs. This is where the
  [long-protein](#long-proteins) issues bite: full-length structures and MSAs
  whose match states cover the full sequence are needed.

The DMS-specific assay fields (`readout`, `growth_coupled`,
`library_construction_method`, `transformation`, …) are all optional, so
clinical data fits the schema without changing proteingym-base.

Constraints to respect:

- One wild-type and one chain per archive (`update_structure` asserts a
  single chain).
- Predictions are joined to labels on the **full sequence** (plus assay
  variables). Variants that give the same protein sequence (e.g. different
  nucleotide changes, same amino-acid change) must be collapsed to one row
  when building the archive.
- `Dataset.to_df()` groups rows with the same (sequence, assay variables) and
  averages numeric targets by default, so duplicate rows are silently merged.

## Where clinical datasets live

There is currently **one shared dataset pool** for both existing games:

- `scripts/generate_datasets_list.py` globs
  `datasets/splits/*.splits.pgdata` and writes the **same** `datasets.json`
  into `benchmark/zero_shot/` and `benchmark/supervised/`, with one entry per
  (dataset, split, target).
- The repo itself only holds a couple of example archives as Git LFS
  pointers. The full collection lives on Hugging Face
  ([ProteinGym/ProteinGym2.0](https://huggingface.co/datasets/ProteinGym/ProteinGym2.0)).

Implications:

1. **Keep clinical archives out of `datasets/splits/`.** Otherwise the
   generator puts them into both DMS games. Use e.g. `datasets/clinical/`,
   and have the clinical game use its own list (a generator flag or a small
   separate script).
2. **Splits.** The runner expects a `.splits.pgdata` with at least one split.
   Options:
   - wrap each clinical dataset with one trivial split, so the runner and
     model containers need no changes (recommended to start);
   - or let the clinical pipeline accept plain `.pgdata`.
3. **Hosting.** Publish clinical archives on Hugging Face like the DMS ones
   (possibly as a separate subset). Keep only a couple of small examples in
   the repo for CI.

## Multiple assays per gene and multiple benchmarks per protein

The convention is **one `.pgdata` per (protein, assay)**. The design work is
in sharing the per-protein inputs (wild-type sequence, structure, MSA) across
archives, not in bundling variants together.

### Two DMS on the same gene

This already happens: e.g. `P53_HUMAN` has 4 datasets and `HSP82_YEAST` has 3
(`HSP82_YEAST_Cote-Hammarlof_2020_growth-H2O2`, `HSP82_YEAST_Flynn_2019`,
`HSP82_YEAST_Mishra_2016`). Each is its own archive, named
`<UNIPROT>_<Author>_<Year>[_condition]`, with its own row in `metrics.csv`.

A manifest *can* hold several `[[assays]]`, but `Dataset.to_df()` groups
rows by (sequence, assay variables) and **averages** each target by default.
Two DMS with the same target name and the same variables would be silently
averaged. A combined archive only makes sense when the assays differ by an
assay variable (like pH/temperature in the NEIME example) or use different
target names. For independent studies, keep separate archives.

### Two benchmarks on the same proteins (e.g. DMS and clinical for BRCA1)

Keep **separate archives in separate benchmark folders** (`datasets/splits/`
vs `datasets/clinical/`). Variant sets, labels, metrics, splits and release
cycles all differ, and each pipeline picks up its datasets by folder.

What must be shared is the per-protein inputs, so results are comparable:

- **A single per-protein input store**, keyed by UniProt ID and isoform:
  wild-type FASTA, AlphaFold structure, MSA and weights. Both benchmarks'
  manifests reference the same files. `build` copies them into each archive,
  so data is stored twice but has one source of truth.
- **Mismatched reference sequences.** DMS archives often use a construct or
  domain fragment as the reference, while clinical needs the full-length
  canonical isoform. With `first_index=1` hardcoded in
  `evedesign/proteingym.py` and MSAs that must match the reference sequence,
  BRCA1-DMS and BRCA1-clinical may legitimately have different reference
  sequences and MSAs, and model scores for the "same" variant can differ.
  For cross-benchmark comparison, score DMS variants in full-length context
  too, or at least record the offset.
- **Store variant identity explicitly.** Prediction files are matched to
  labels by full sequence, which has no meaning outside that archive. Store
  `uniprot_id`, `hgvs_p` / `mutant` (e.g. `R1699W`) and the canonical-isoform
  position as `non_targets` in both benchmarks, so DMS scores, clinical
  labels and model scores can be joined per variant.
- **Recomputation.** Each pipeline re-runs every model on its own archive.
  Zero-shot scores depend only on (model, inputs, variant sequence), so they
  could be cached by a hash of those. The pipeline does not do this today. It
  would be an optimisation for later, mostly useful for long clinical
  proteins.

| Case | Archive layout | Shared |
|---|---|---|
| 2 DMS, same gene | 2 archives, same folder, `<UNIPROT>_<study>` names | wild-type, structure, MSA (if same construct) |
| DMS + clinical, same gene | separate archives in separate benchmark folders | per-protein input store; variant IDs in `non_targets` for cross-joins |

## Open questions

- **Which summary statistic?**
  - Per-variant effect sizes or odds ratios from biobanks: continuous, noisy,
    sparse per gene. Pooled Spearman is natural; maybe weight by standard
    error. This fits the per-protein, two-target archive design above.
  - Gene-level or burden statistics (one number per gene): not a per-variant
    target. Model scores would need aggregating to gene level first, which is
    a different evaluation shape that needs its own design.
  - Allele frequency or constraint: a proxy for selection, but confounded by
    mutation rate.
- Which models must support long proteins for the clinical leaderboard, and
  how are missing predictions (model refused, variant outside MSA) counted?
- Label curation: ClinVar review-star threshold, handling of conflicting
  interpretations, and avoiding circularity with models trained or tuned on
  clinical labels.
- Should pooled metrics eventually move into proteingym-base?

## Proposed next steps

1. Build a single example gene archive with `proteingym-base build`
   (wild-type, AlphaFold structure, MSA, variants CSV with both targets),
   wrap it with a trivial split, and run `blosum62` and `esmc` on it to
   check the whole pipeline end to end.
2. Add `benchmark/clinical/` (DVC pipeline keyed on (dataset, model) rather
   than (dataset, target), and its own `datasets.json` sourced from
   `datasets/clinical/`).
3. Add the pooled evaluation script (AUROC/AUPRC/Spearman, pooled and
   per-gene averaged, filtering nulls per target).
4. Prototype windowed ESM-2 scoring for proteins over 1022 residues.
5. Write a builder that produces clinical manifests from the shared
   per-protein input store plus a variant table.
