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
| `eve`, `evmutation`, `evmutation2`, `gemme`, `rsalor`, `structuredca` | **Not feasible** at L ≈ 4,500 with N ≈ 1M (memory 100+ GB, O(N²) weighting); some score unscored positions as WT | MSA depth, L-scaling, no region support; see [§3.4](#34-alignment-models) |
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

Code audit 2026-11-16 (wrappers in `models/`, evedesign 0.0.5 for five models, 0.0.6 for `eve`; upstream packages
not inspected are marked [inferred]). N = MSA depth (up to 1.16M), L = length, q = 21.

**Shared runner and conversion.**
- Exceptions are not caught: one failure kills that dataset's run, no predictions written. Only `eve` pre-filters
  invalid instances; in the other five one invalid instance raises (`evedesign/model.py:200-216`).
- `dataset_to_evedesign` builds one Python `Sequence` per a3m row (`proteingym.py:69-75`) and `remove_inserts()`
  copies it again (`sequence.py:274-290`): several GB of Python strings for the deepest MSAs before any model runs.
- `first_index = 1` hard-coded (`proteingym.py:139`), first MSA only (`:61`), `msa_weights` from the dataset are used
  if present (`:63-67`), single-chain structures asserted (`:94`); "MSA target region != full sequence" is an open TODO
  (`:147-148`).

| | eve | evmutation (plmc) | evmutation2 | gemme | rsalor | structuredca |
|---|---|---|---|---|---|---|
| query == WT check | focus columns only | yes, raises | length only | yes | yes (+ raises if trimmed) | yes |
| depth cap | 200k random (`max_msa_depth`) | none | 2,048 for the encoder, whole MSA parsed first | JET2: first 20k rows; R: full MSA | none [min_seqid 0.35, inferred] | none [min_seqid 0.25, inferred] |
| filters | fragments > 50% gaps dropped; all columns kept | columns > 50% gaps (unweighted) excluded; no fragment filter | — | none | upstream | none |
| weighting | dense numpy O(N²·20L) | numba O(N²L) on the full MSA | — | JET2 | C++ O(N²L) [inferred] | C++ O(N²L) [inferred] |
| unscored positions | none | **excluded columns: mutant = WT (0), not NaN** | none | NA cells add 0 | **no-RSA positions: mutant = WT** | upstream |
| needs structure | no | no | no | dummy | **one full-length chain** | **exactly one structure, one chain** |
| hardware | GPU (weights, one-hot on CPU) | CPU/OpenMP | CPU | CPU, JVM 4 GB | CPU | CPU |

**Why they fail on long proteins (L ≈ 4,300–4,700, N up to 1.16M):**
- **EVE:** one-hot training matrix float64 N×L×20 = ~138 GB at the 200k cap and L = 4,300; weighting is one gemv
  over it per sequence; VAE ~0.9B parameters (~15 GB GPU with Adam) [inferred]. The 50% fragment filter drops most
  domain-local ColabFold hits on long proteins [inferred], so Neff collapses.
- **EVmutation:** alignment as N×L U1 array ~21 GB for LRP2, weights via `np.vectorize` (5×10⁹ Python calls), N²
  weighting ~6×10¹⁵ comparisons; plmc couplings L(L−1)/2·q² ≈ 4×10⁹ floats ≈ 16 GB per vector, L-BFGS keeps several
  (> 100 GB) [inferred]. Columns excluded for gaps silently score 0.
- **EVmutation2:** whole MSA parsed; 16 encoder samples each with an L×L×c_z pair representation (~11 GB each at
  L = 4,655, c_z ≈ 128) [inferred]; ~O(L³) triangle updates on CPU [inferred].
- **GEMME:** R (`computePred.R`) receives the full 1.16M×L FASTA [inferred: far beyond 100 GB]; JET2 gets the top
  20k rows in a3m order, not a random sample.
- **RSALOR, StructureDCA:** need one full-length single-chain structure, which AFDB does not provide above 2,700 aa;
  positions without RSA score as WT.

**Solution to implement (region models):**
1. A region step between `dataset_to_evedesign` and `build`: crop the wild type, the MSA match columns (after
   `remove_inserts`), any weights and the structure (renumbered) to a region [s, e) around the variants; set
   `first_index = s+1`; crop the instances. The query still equals the cropped wild type, so all query checks pass,
   and every model indexes instances by array position, so cropping instances is enough.
2. Inside each region: drop rows that are mostly gaps in the region, then cap depth by random sampling, then compute
   (or supply) sequence weights once per region. Supplying `msa_weights` skips in-model weighting for EVE and
   EVmutation.
3. Regions: windows of a few hundred residues (size per model, set by memory: L²·q² for EVmutation), most-central
   region per variant (same rule as the PLM windows); for structure-needing models, the region must lie inside one
   AFDB fragment.
4. Positions a model cannot score → **NaN**, never 0/WT: patch EVmutation (excluded columns), RSALOR (no RSA) and
   GEMME (NA cells) wrappers to return NaN, and catch per-protein failures in the runner so one protein does not
   lose the dataset.
5. Test T5 first on APOB, LRP2, PKD1 and one mid-size protein: memory, time, coverage, and region vs full length on
   a protein where full length fits.

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

Shared logic is duplicated per folder to keep folders self-contained: `window.py` (PLM windows, `esm2_window`,
`esmc_window`), `fragments.py` (structure placement and most-central assignment, `saprot_window`, `esm-if1_window`).
Validation scripts and job wrappers: `scripts/clinical/models/validate_{plm,struct}_windows.{py,sbatch}` (run in the
prototype's zoo images with the new package mounted on `PYTHONPATH`; proper images still to build). Outputs in
`$WORK/models/`, job IDs in `$WORK/logs/jobs.txt`. Upstream route later: windowed `score_mutants` in evedesign's `ESM2` would fix `esm2` and
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
| 2026-11-16 | Alignment-model code audit | §3.4: none feasible on L ≈ 4,500 / N ≈ 1M; some score unscored positions as WT | — |
| 2026-11-16 | `esm2_window`, `esmc_window` ported (shared `window.py`, `special_tokens` option) | — | `models/esm*_window/` |
| 2026-11-16 | T1 esm2_window, 3 short proteins (job 31921230) | max abs diff vs zoo `score_mutants` 4.8e-7: **PASS** | `$WORK/logs/plm_windows_esm2_smoke_31921230.log` |
| 2026-11-16 | `saprot_window`, `esm-if1_window` ported with `fragments.py` (structure placed on the WT by sequence, most-central structure per position) | APOB: 17 fragments placed at 1, 201, …, 3,201 (identity ≥ 0.9993, Y1422C), no position uncovered (local check) | `models/{saprot,esm-if1}_window/` |
| 2026-11-16 | T1 saprot_window, 3 short proteins (job 31921437) | max diff vs zoo 6.7e-6: **PASS** | `$WORK/logs/struct_windows_saprot_smoke_31921437.log` |
| 2026-11-16 | esm2_window full run on 1g.20gb with the card's batch 64 (job 31921334) | CUDA OOM (64 × 1,022-residue windows); rerun with batch 8 (31921955) | — |
| 2026-11-16 | esm-if1_window smoke (job 31921438) | ImportError: fair-esm imported before the zoo module patches biotite (same bug as the prototype); fixed import order, rerun 31921956 | — |
| 2026-11-16 | esm2_window, 66 proteins, 5,840 unique mutants (job 31923533, 3g.40gb, batch 16) | T1 on all 48 proteins ≤ 1,022 aa: 1.6e-5 **PASS**; termini vs every_window on L > 1,022: r 0.99998 (max 0.18) | `$WORK/models/plm_windows/esm2/` |
| 2026-11-16 | esmc_window windows 1,022 / 2,048 / full (job 31921335, 26 min on 3g.40gb) | T1 48 proteins 2.7e-5 **PASS**; termini vs every_window r 0.9985 (max 1.73): ESMC is more sensitive to edge tokens than ESM2 | `$WORK/models/plm_windows/esmc/` |
| 2026-11-16 | saprot_window, every archive (job 31923536) | T1 1.1e-5 **PASS**; APOB/LRP2/PKD1 use all 16–18 fragments, no position uncovered; same variant in different archives: spread 1e-5 | `$WORK/models/struct_windows/saprot/` |
| 2026-11-16 | esm-if1_window, every archive (job 31923537, 3g.40gb) | T1 1.8e-4 **PASS** (window = L·(zoo − zoo_wt)); vs prototype r 1.00000 on all 5,840 variants (max 0.03 on long proteins); same variant across archives: spread 4e-3 (batch numerics); one caching-allocator OOM warning on APOB, recovered | `$WORK/models/struct_windows/esmif1/` |
| 2026-11-16 | Alignment-model images (job 31923582, genoa64 CPU, `build_images.sbatch`) | rsalor built (3 min); gemme failed: download server gave an empty reply (`curl: (52)`), retry; others in progress | `/no_backup/diasfrazer/egalpern_pgb_clinical/containers/` |
| 2026-11-16 | Per-variant comparison with earlier implementations (`scripts/clinical/models/export_scores.py`) | see "T2/T3 results" below | — |

### T2/T3 results (2026-11-16)

Per variant, Pearson r and max abs difference; L ≤ 1,022: 2,933 variants, L > 1,022: 2,621 variants.

| New | Reference | L ≤ 1,022 | L > 1,022 |
|---|---|---|---|
| esm2_window every_window | prototype esm2_window (every window) | 1.00000 / 1.4e-5 | **1.00000 / 3e-6** |
| esm2_window termini | prototype esm2_window | 1.00000 | 0.99998 / 0.18 |
| esm2_window termini | dark-energy ESM2 650M masked (termini convention, fair-esm code) | 1.00000 / 1.3e-5 | 0.99997 / 1.3 |
| esm2_window every_window | dark-energy ESM2 650M masked | 1.00000 | 0.99995 / 1.28 |
| esmc_window every_window | prototype esmc_window | 1.00000 / 2.7e-5 | **1.00000 / 7e-5** |
| esmc_window termini 1,022 | prototype | 1.00000 | 0.99853 / 1.73 |
| esmc_window termini 2,048 | prototype | 1.00000 | 0.99663 / 1.94 |
| esmc_window full length | prototype | 1.00000 | 0.99159 / 3.6 |
| saprot_window (termini) | prototype saprot_clinical (every window, fragment datasets) | 1.00000 / 1.9e-5 | 0.99920 / 1.46 |

- The ports reproduce the prototype exactly in the same token mode, which confirms the prototype put `<cls>`/`<eos>`
  on every window.
- The residual against the dark-energy scores on long proteins (max 1.3) is the same in both token modes, so it is
  not the edge tokens; open (candidates: the window placement at the ends, fair-esm vs transformers numerics).
- **Benchmark effect of the ESMC choices is small** (Genebass, mean −log10 P, all 109 pairs): 2,048 window 2.200,
  every_window 1,022 2.200, full length 2.192, termini 1,022 2.181 (prototype on its 105 pairs: 2.159).

### Genebass benchmark with the ported models (2026-11-16)

Mean −log10 P of per-pair Pearson r(β, score) (`proteingym_clinical/genebass/scripts/16_eval_genebass_source.py`,
predictions exported by `scripts/clinical/models/export_scores.py` as `v2__*`; layers from their own releases).

| Model | All 109 pairs (coverage) | Paper set, 103 pairs |
|---|---|---|
| AlphaMissense (layer) | 2.695 (100%) | 2.363 (paper: 2.37) |
| VESM++ (paper's scores, reference) | 2.671 (105 pairs) | 2.468 |
| esm2_window, termini 1,022 | 2.403 (100%) | 2.153 (paper ESM2-650M, unmasked: 2.10) |
| popEVE (layer) | 2.372 (97 pairs, 74% of rows) | 2.053 (91 pairs) |
| esmc_window, termini 1,022 | 2.181 (100%) | 1.973 |
| CADD (layer) | 1.912 (99.9%) | 1.682 |
| saprot_window | 1.640 (100%) | 1.529 |
| Allele frequency | 1.629 | 1.351 |
| esm-if1_window | 1.273 (100%) | 1.140 |

All four ported models now cover every pair (the prototype missed CST3, GPLD1 and UGT1A5). PottsMPNN and BLOSUM62
are not ported yet (PottsMPNN needs the fragment logic; BLOSUM62 needs nothing).

## 10. Picking this up

State on 2026-11-16 (evening):
- **Done, validated (T1–T3), benchmarked:** `esm2_window`, `esmc_window`, `saprot_window`, `esm-if1_window`
  (§9). They run in the prototype's zoo images (`project/proteingym_clinical/containers/*.sif`) with the new package
  mounted (`scripts/clinical/models/validate_{plm,struct}_windows.sbatch`). Their own images are not built yet
  (`build_images.sbatch esm2_window ...`).
- **Data side:** final archives `$WORK/genebass/splits_layers/` (MSAs, layers from AlphaMissense/popEVE releases,
  `in_vesm_paper` flags); the model runs used the no-MSA copies in `$WORK/genebass_api/build/splits/` (same rows).
- **Next, in order:**
  1. Alignment-model images: job 31923582 (`$WORK/logs/build_images_*.log`); rebuild `gemme` (download failed).
  2. T5 on APOB, LRP2, PKD1 + one mid-size protein: region archives (§3.4: crop WT, MSA match columns, structure
     fragment; re-filter rows by coverage in the region; random depth cap; weights once per region), run the
     unchanged zoo models on them, map predictions back. Then NaN instead of WT for unscored positions (EVmutation,
     RSALOR, GEMME) and per-protein failure handling.
  3. Port PottsMPNN with `fragments.py` (most-central fragment) and BLOSUM62 unchanged; decide `special_tokens`
     (§8) — ESMC is the only model where it matters (benchmark 2.181 termini vs 2.200 every_window).
  4. Open: residual vs dark-energy ESM2 scores on long proteins (max 1.3, r 0.99997).
