# Zoo models on long proteins: problem, solutions, tests

Living document for the model side of the clinical benchmark. Started 2026-11-16, when the model side became ours
(it was the model-API collaborator's in the [design notes](../clinical_benchmark.md)). Every change, test and result
goes here, newest results in the [experiment log](#9-experiment-log).

Branch: `clinical/model-api` (from `docs/clinical-benchmark-design`), worktree `project/pgb_clinical_models`.
Model code only; dataset code stays on `clinical/*-data` branches (see [project.md](project.md#5-organisation-of-the-work)
on the data branch).

## Contents

1. [The problem in numbers](#1-the-problem-in-numbers)
2. [Zoo behaviour today](#2-zoo-behaviour-today)
3. [Solutions per model family](#3-solutions-per-model-family)
4. [Validation protocol](#4-validation-protocol)
5. [Implementation in the fork](#5-implementation-in-the-fork)
6. [Running on the CRG cluster](#6-running-on-the-crg-cluster)
7. [What the prototype already showed](#7-what-the-prototype-already-showed)
8. [Open decisions](#8-open-decisions)
9. [Experiment log](#9-experiment-log)

## 1. The problem in numbers

**Genebass benchmark (final build, 109 gene–phenotype pairs, 66 proteins, 7,791 variant rows):**

| Length | Proteins | Share of variant rows |
|---|---|---|
| > 1,022 aa (ESM-2 / ESM-1v context) | 18 | 46% |
| > 2,048 aa (ESMC context) | 6 | 33% |
| > 2,700 aa (AFDB fragments only) | 3 (LRP2 4,655; APOB 4,563; PKD1 4,303) | 22% |

Variant-weighted mean length: 1,768 aa. Next longest: STAB1 2,570, SPTA1 2,419, ABCA1 2,261, LPA 2,040.

**The wider Genebass set (332 pairs, 185 proteins) and ClinVar/BRaVa** are worse: 52 of 185 proteins are over
1,022 aa and hold 57% of the variants; 8 are over 2,700 aa (up to DST, 7,570 aa). The methods below must scale to
the proteome, not only to these 66 proteins.

**Inputs in the archives:** full-length UniProt wild type; AFDB v6 structures (one model up to 2,700 aa, otherwise
AFDB fragments F1..Fn of 1,400 aa starting every 200 aa, each with `span_start`/`span_end` in the manifest);
ColabFold MSA on UniRef100 2026_01 with the full wild type as query, up to 1.16M sequences and 3.9 GB (PKD1).
Structures stay AFDB v6 with fragments (user decision 2026-11-16, no AF2 re-folding).

## 2. Zoo behaviour today

From the code (zoo in `models/`, evedesign ≥ 0.0.6). Runner: every model goes through the shared `__main__.py`,
which converts the archive with `evedesign.proteingym.dataset_to_evedesign` (single-protein `System`, first MSA,
all structures, `first_index = 1`) and calls `build(system)` then `score(instances)` on full-length mutant sequences.
Nothing in the runner handles length.

| Model | Today on long proteins | Root cause |
|---|---|---|
| `esm2` | **Fails** above 1,022 aa | `max_seq_length = 1022` checked in `can_model` and per instance |
| `esm2` (any length) | Scores the unmasked pseudo-log-likelihood of the whole mutant, one full pass per mutant | zoo default `score`; masked marginal (`score_mutants`) exists in evedesign but is unused |
| `esmc` | Runs; L² memory; risk of silent truncation | cap removed, but `score_mutants` tokenises with `truncation=True` [inferred] |
| `saprot` | Runs beyond its 1,024-token training length, quality unknown | no cap |
| `esm-if1` | One full autoregressive pass per mutant, encoder re-run each time | zoo `score` |
| `pottsmpnn` | Works (one encoder pass; ~21k aa on 24 GB measured elsewhere) | — |
| structure models with fragments | evedesign asserts one chain, `first_index = 1`; several fragments are passed as separate structures | runner |
| all structure models | **pLDDT is 0 inside every archive** | proteingym-base zeroes B-factors when it writes structures |
| `eve`, `evmutation`, `evmutation2`, `gemme`, `rsalor`, `structuredca` | See [§3.4](#34-alignment-models) (audit in progress) | MSA depth and L-scaling |
| `blosum62` | Fine | — |

## 3. Solutions per model family

### 3.1 Protein language models (ESM-2, ESMC, VESM-style)

**Masked marginal with one window per mutated position.** For each unique mutated position, a window of `window`
residues centred on it (shifted to fit at the protein ends) is passed with the position masked; score =
Σ over substitutions of log p(mut) − log p(wt). For L ≤ window the window is the whole protein, so the result must
equal the zoo's `ESM2.score_mutants` exactly. Cost: one forward pass of ≤ window residues per unique position
(positions barely repeat in clinical data: 79,752 variants at 62,453 positions in the 332-pair set).

**Special tokens at window edges (to decide, both implemented as an option `special_tokens`):**
- `termini`: `<cls>` only if the window starts at residue 1, `<eos>` only if it ends at the C-terminus. This is
  the dark_energy convention (`dark_energy/scripts/esm1v/esm1v_core.py:128-181`) and tells the model the truth
  about where the protein starts and ends.
- `every_window`: `<cls>`/`<eos>` around every window, as the tokenizer does by default (ProteinGym v1 ESM scripts
  and the 2026-11-06 prototype).

Note: `long_proteins.md` (data branch) says the prototype puts the tokens only at the true termini. The prototype
code (`proteingym_clinical/models/clinical_window.py`) actually adds them to every window. This explains the small
residual difference against the dark-energy scores (per-variant r 0.9993 for L ≤ 1,022 vs 0.99997 above it is the
other way round from what one would expect, so the short-protein residual must come from something else, e.g.
batching/padding; to check in test T2).

**Window size:** 1,022 for every PLM for parity (ESM-2 context). ESMC (2,048 trained context) also gets a
native-context card. Test T3 measures the window effect on proteins where the full length fits.

**Scoring mode:** masked marginal is the default (ProteinGym standard). A wild-type-marginal / unmasked option is
cheap (one pass per window, ~400× fewer residues) and on Genebass VESM-650M unmasked beats masked (2.41 vs 2.10,
Fig. 6b), while ESM-2 and ESM-1v prefer masked. Offer it as a separate card, not a replacement.

### 3.2 Structure-conditioned sequence models (SaProt, ESM-IF1)

- **SaProt:** same per-position windows on the paired (residue, 3Di) tokens; 3Di from the structure used for that
  position (see fragments below).
- **ESM-IF1:** tiles of `window` residues every `window/2`; each mutant scored in the tile where it is most
  central, as Σ log p(mutant | tile backbone) − Σ log p(wt | tile backbone), batched. For L ≤ window this is
  L·(zoo score(mut) − zoo score(wt)): a per-protein linear transform of the zoo score.
- **Fragments (proteins > 2,700 aa):** for each position use the fragment where it is most central (never the
  dark_energy unweighted mean over fragments: PottsMPNN scores disagree between fragments by up to 16.6, mean 1.56).
  Fragment spans are in the manifest (`[structures.metadata] span_start/span_end`).
- **pLDDT:** never read from archive structures (zeroed). If a model needs it, pass it from the AFDB files
  (`layers` already does this for `mask_plddt`), and report the proteingym-base bug upstream.

### 3.3 Structure-based Potts models (PottsMPNN, StructureDCA)

PottsMPNN: unchanged on full length up to ~21k aa (one encoder pass). For > 2,700 aa proteins, either a stitched
model per protein or most-central-fragment scoring; the second is consistent with the other models.
StructureDCA: see §3.4 (it also needs the MSA).

### 3.4 Alignment models

Pending the code audit (EVE, EVmutation/plmc, EVmutation2, GEMME, RSALOR, StructureDCA): MSA filters and which
positions they drop, depth caps and weighting cost, L²·q² memory, behaviour on positions without a score, and where
a region/window rebuild fits. Candidate solutions to test:
- depth cap by subsampling (EVE already caps at 200k);
- region models: one model per window of a few hundred residues around the variants (as ProteinGym v1 did for
  EVE with focus regions), most-central region per position;
- positions the model cannot score → missing prediction (NaN), reported as coverage, never silently 0.

## 4. Validation protocol

Every change is validated before any benchmark run:

| Test | What | Pass criterion |
|---|---|---|
| T1 | New model on proteins with L ≤ window vs the unchanged zoo model | identical scores (max abs diff ~1e-6) or a documented exact transform |
| T2 | Windowed PLM vs the dark-energy masked scores (independent implementation) | per-variant r ≥ 0.999; explain any residual |
| T3 | Window vs full length where the model allows it (ESMC ≤ 2,048, SaProt, PottsMPNN) | report r and benchmark difference |
| T4 | Fragment choice: most-central vs other fragments (structure models, > 2,700 aa) | report spread |
| T5 | Alignment models on APOB, LRP2, PKD1: memory, time, coverage | fits the cluster; coverage reported |
| T6 | Genebass benchmark (109 pairs, all rows; paper set via `in_vesm_paper`) | numbers logged in §9 |

## 5. Implementation in the fork

Zoo convention: self-contained model folders, model selected by the card's `model_class`, `hyper_parameters`
forwarded to `__init__`. New behaviour = a subclass + a new card (this is how the zoo's `esmc` differs from `esm2`).

Planned folders (copies of the zoo folder, same Dockerfile/Singularity.def/lockfile, plus the subclass):

| Folder | Class | Based on |
|---|---|---|
| `models/esm2_window` | `ESM2Window` | `evedesign.models.esm2:ESM2` |
| `models/esmc_window` | `ESMCWindow` | zoo `esmc` |
| `models/saprot_window` | `SaProtWindow` | zoo `saprot` |
| `models/esm-if1_window` | `ESMIF1Window` | zoo `esm-if1` |

Shared logic (window bounds, most-central tile/fragment) is duplicated per folder as a small `window.py`, to keep
folders self-contained. Upstream route later: windowed `score_mutants` in evedesign's `ESM2` would fix `esm2` and
`esmc` at once.

## 6. Running on the CRG cluster

- Containers: Apptainer, built from each folder's `Singularity.def` with `apptainer build --fakeroot` on the GPU
  node; the prototype's images are in `project/proteingym_clinical/containers/`.
- GPU: partition `gpu`, `--gres=gpu:7g.80gb:1` (full H100). All compute as sbatch, never in the interactive session.
- One container start per model for many archives (`proteingym_clinical/genebass/scripts/batch_train.py`: identical
  outputs, 1.5 s instead of 37 s overhead per protein).
- Archives: `/no_backup/diasfrazer/egalpern_pgb_clinical/genebass/splits_layers/` (final), or the no-MSA copies in
  `.../genebass_api/build/splits_layers/` for sequence models.

## 7. What the prototype already showed

Prototype (2026-11-06, `project/proteingym_clinical/`, subclasses mounted into the unchanged zoo containers;
Genebass paper set, 101 pairs at the time):

| Model | Adaptation | Validation | Fig. 6b mean −log10 P |
|---|---|---|---|
| esm2_window | 1,022 window per position, masked marginal | = zoo `score_mutants` on short proteins (2.4e-7); r 0.9993 / 0.99997 vs dark energy | 2.15 |
| esmc_window | same | = zoo (4e-7); APOB window vs full length r = 0.921 (48 variants) | 1.99 |
| saprot_clinical | window on residue + 3Di tokens | = zoo exactly | 1.51 |
| pottsmpnn_clinical | unchanged, fragments most-central | — | 1.28 |
| esm_if1_clinical | tiles 1,022 / stride 511 | = L·(zoo − zoo_wt), r = 1.000000 | 1.13 |
| blosum62 | — | — | 0.78 |

The same predictions on the final Genebass-sourced set (2026-11-16, 109 pairs, all rows): esm2_window 2.35,
esmc_window 2.16, SaProt 1.57, PottsMPNN 1.33, ESM-IF1 1.22, BLOSUM62 0.79; reference VESM++ 2.67, AlphaMissense
2.65 (`proteingym_clinical/genebass/scripts/16_eval_genebass_source.py`).

## 8. Open decisions

1. Special tokens at window edges: `termini` or `every_window` (default proposed: `termini`, after T2).
2. ESMC: 1,022 for parity only, or also a native 2,048 card.
3. Unmasked/wt-marginal cards for the PLMs in addition to masked.
4. Alignment models: region policy and depth caps (after the audit and T5).
5. Missing predictions in the metric: drop the variant for that model (current evaluation) or penalise coverage.

## 9. Experiment log

| Date | What | Result | Where |
|---|---|---|---|
| 2026-11-16 | Doc started; prototype results re-evaluated on the final Genebass set | §7 | `proteingym_clinical/genebass/results/genebass_source_*` |

## 10. Picking this up (handoff, 2026-11-16)

State when the interactive session ended:
- **Model side (this branch):** this doc only; no model code ported yet. Next: (1) finish §3.4 from a code audit of
  the alignment models (the audit was running in a background agent that died with the session; rerun it: read
  each `models/<m>/` card + evedesign `models/{eve,evcouplings,evmutation2}.py` + `proteingym.py` for MSA filters,
  dropped positions, depth caps, L-scaling, first_index assumptions); (2) port the prototype subclasses
  (`proteingym_clinical/models/clinical_window.py`, `clinical_struct.py`) into the folders of §5 with the
  `special_tokens` option; (3) run T1–T2 on the GPU (§6).
- **Layer coverage fix (data side):** worktree `project/pgb_clinical_layers`, branch `clinical/layers-sources`.
  AlphaMissense from `/no_backup/diasfrazer/db/alphamissense/raw/AlphaMissense_hg38.tsv.gz`, popEVE from the local
  release zip via Swiss-Prot RefSeq cross-references, CADD stays dbNSFP. Was being done by a background agent; check
  `git log`/`git status` in that worktree and `$WORK/logs/jobs.txt` for what it finished.
- **Data rebuild:** SLURM job 31919891 (in_vesm_paper flags), log `$WORK/logs/rebuild_31919891.log`; then commit
  and push the index on `clinical/genebass-data` (see `docs/clinical/status.md` there).
