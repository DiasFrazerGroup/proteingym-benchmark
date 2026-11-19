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

**Special tokens at window edges** (option `special_tokens`; **decided 2026-11-16: `termini`, real ends only**, the default in every card):
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

**User decision 2026-11-16: no derived datasets in the final design.** The region archives below are the T5 /
first-benchmark vehicle only (numbers are identical either way). The final design moves the region logic inside
the alignment-model wrappers, like `esm2_window`/`saprot_window`: each model reads the same archive as every other
model and its card declares the policy (`region: 1000`, `min_coverage: 0.5`, `max_depth: 50000`, `seed`). Planned
changes per model: [§10](#10-picking-this-up).

**What was built: region archives (2026-11-16, `scripts/clinical/models/regions/`).** The six models stay
unchanged; a preprocessing step turns one protein into small archives they can score, and the scores are mapped
back.

1. **`make_regions.py`** (local env with proteingym-base; `make_regions.sbatch <UNIPROT> <W>...`): all archives of
   one UniProt entry are pooled (variants as unique mutants, so several phenotype archives of a protein share one
   model fit). Output `$WORK/models/regions/archives/<UNIPROT>/W<W>/`: one `<UNIPROT>_r<s>-<e>.splits.pgdata` per
   region, `regions.json` (spans, MSA depths, structure, EVmutation-excluded columns), `variants.csv` (UniProt
   mutant -> region, region mutant), `rows.csv` (original archive rows).
   - **Regions:** tiles of W residues every W/2 (last tile shifted to end at L; W >= L gives one full-length
     region); each mutant goes to the tile where it is most central (`tiles`/most-central rule of `esm-if1_window`);
     only tiles with variants are written.
   - **MSA:** the archive's a3m is streamed once for all regions (bytes, no Python objects per sequence): lowercase
     insertions dropped (`bytes.translate`), match columns [s, e), rows with < `--min-cov` (0.5) non-gap residues
     in the region dropped, duplicates dropped (hash set), query first, then a uniform reservoir sample of
     `--max-depth` (50,000, query included) rows with a fixed seed. Depth total / covering / unique / kept recorded.
     B/Z/U/O/J become X.
   - **Structure:** every AFDB fragment placed on the wild type with `fragments.place` (the `saprot_window` module,
     loaded by path); the region gets the fragment that fully contains it and in which the region centre is most
     central, cropped to [s, e) and renumbered 1..W (chain A). Mismatches to the wild type and missing residues are
     recorded. No containing fragment -> no structure in the region archive (cannot happen for W <= 1,200 with
     AFDB's 1,400/200 fragments).
   - **Assay:** cropped sequences, `mutant` renumbered (position - s), `uniprot_mutant` kept; target column = the
     first archive's value (unused by zero-shot models, needed by the format). Written with proteingym-base
     (`Dataset.from_manifest` -> `Subsets(slices={}).dump`).
2. **`run_regions.py`** (inside the model's image, `run_regions.sbatch <model> <tag> <UNIPROT>/W<W> ...`, one
   container start for many regions): the zoo runner's own steps (`Subsets.from_path` -> `dataset_to_evedesign`
   -> card `model_class(**hyper_parameters)` -> `build` -> `score`), plus
   - the region wild type is scored with the variants and **pred = score(variant) - score(WT)**: eve (ELBO),
     evmutation (Hamiltonian), evmutation2 (log-likelihood) and structuredca (-energy) return absolute sequence
     scores whose offset differs between regions; gemme and rsalor are WT-relative already (score(WT) = 0);
   - **unscored flags** (become NaN, never WT/0): a mutated position outside `model.positions()` (evmutation:
     columns with unweighted gap fraction > 0.5 over the region MSA; eve: non-focus columns, none with the card's
     `threshold_focus_cols_frac_gaps = 1.0`), or a NaN cell of the model's per-position matrix `_pssm` at
     (position, mutant aa) (rsalor: no RSA at the position; gemme: NA cell from pred.R). The wrappers would add 0
     for all of these;
   - per-region `try/except`: a failed region is logged (`status.jsonl`, `<region>.error.log`) and the others are
     scored; per-region wall time and peak RSS (self and child processes: plmc, R, java) are recorded.
   Card hyper-parameters can be overridden with `HP='{...}'` (only used where stated).
3. **`collect_regions.py`** (local env): maps region predictions back to UniProt mutants (NaN for unscored or
   failed regions, with the reason), checks the evmutation flag against the gap rule recomputed from the region
   MSA, writes the T5 summary (`summary`), region-vs-full comparisons (`compare`) and the Genebass export
   (`export` -> `proteingym_clinical/genebass/predictions/v2__<model>/<UNIPROT>.csv`, `mutant,pred`, NaN rows
   dropped = missing for that model).

Plumbing check (local, APOA1): one full-length region with `--min-cov 0` scored by rsalor equals the unchanged
model on the original archive (max abs diff 0.0; rsalor drops duplicates itself).

### 3.5 What EVE / popEVE did for the human proteome (checked 2026-11-17)

User decision 2026-11-17: no short MSA segments; learn from the EVE/popEVE proteome run and apply it. The region
work of §3.4 is stopped (the wrapper refactor agent was stopped before changing anything; region scripts stay as
the T5 record).

**Sources:** popEVE paper (Orenbuch et al., Nat. Genet. 2025, Methods "Multiple sequence alignments", Europe PMC
PMC12695638); EVE preprint (Frazer et al., bioRxiv 10.1101/2020.12.21.423785, Methods); the popEVE v1.1 release in
`/no_backup/diasfrazer/db/popeve/raw/popEVE_ukbb_20250312.zip` (18,331 RefSeq proteins).

**Alignment protocol (EVE 2021, reused by popEVE):** EVcouplings pipeline, jackhmmer, 5 iterations, UniRef100
(March 2022 for popEVE). Bit-score thresholds explored around 0.3 bits/residue; the chosen alignment has coverage
L_cov >= 0.8 L and 10 L <= N <= 100,000, relaxed step by step when no threshold satisfies it. Sequences kept only if
they align to >= 50% of the target; columns kept with >= 70% occupancy. Sequence weights with identity threshold
T = 0.2. Canonical UniProt isoform. Ours: ColabFold/MMseqs2 on UniRef100, no filtering, up to 1.16M rows (user's data
decision 2026-11-09, PG2 protocol).

**Long proteins: they were not solved, they were left out.**
- EVE 2021: "avoided extremely large genes such as TTN"; no explicit cutoff stated.
- popEVE paper: no long-protein handling described (no segments, windows or domains).
- popEVE release: the longest protein is 2,302 aa (NP_057425.3); 131 proteins > 2,000 aa, 41 > 2,200 aa, none above
  2,302. EVE covers 100% of positions in every file checked, i.e. full-length models, no segments. ESM-1v is missing
  for several proteins of 1,800–2,230 aa (popEVE falls back to EVE there).
- Consistent with our layer coverage: popEVE lacks APOB (4,563), LRP2 (4,655), PKD1 (4,303), SPTA1 (2,419) and STAB1
  (2,570), all above 2,302 aa.

**What this means for us:**
1. Up to ~2,300 aa, the EVE recipe is full length with a filtered, capped alignment (N <= 100k, sequences >= 50%
   coverage, columns >= 70% occupancy). This filtering can be done inside the model wrapper on our ColabFold MSA
   (model-side, no derived data). Re-searching with jackhmmer would be a data change (user decision needed).
2. Above ~2,300 aa there is no precedent in EVE/popEVE. Full-length options: EVE with the one-hot built per batch
   (int8 N x L instead of float64 N x L x 20: ~0.5 GB at N = 100k, L = 4,600; VAE input 92k, first layer ~0.2B
   parameters, fits a 7g.80gb slice [estimate]); sequence weights on GPU in blocks or by MMseqs2 clustering.
3. Validation available locally: our EVE on proteins <= 2,300 aa vs the release's `EVE` column, per protein.

**EVE filter on our ColabFold MSAs (2026-11-17, job 31937015, 1.3 min; `scripts/clinical/models/msa_filter/`).**
User decision: keep ColabFold as the base alignment, apply the EVE filter, give every MSA model the same filtered
alignment; compaction (one-hot per batch) is a second step. The a3m headers carry the MMseqs2 bit score, identity,
e-value and aligned ranges, so EVE's bit-score sweep is reproducible (bits/residue = bitscore / L; MMseqs2 scale, not
HMMER's, hence a sweep 0-1.5). Per protein and threshold: N after the >= 50% coverage filter, L_cov (columns >= 70%
occupied). Table: `$WORK/models/msa_filter/eve_sweep.tsv`.
- The >= 50% coverage filter alone takes every MSA below 100k (1M rows -> 2.5k-100k). The bit-score threshold then
  mostly serves L_cov.
- Strict EVE criteria (L_cov >= 0.8 L, 10 L <= N <= 100k) met by 44 / 66 proteins.
- Long multidomain proteins lose most rows to the coverage filter (domain-local hits): APOB 2,533 of 21,926 rows
  (N/L 0.55), PKD1 4,186 of 883,818 (N/L 1.0), LRP2 18,958 of 768,451 (L_cov 0.66 at 0; 10,474 at 0.05 with L_cov 0.88),
  LPA 311 (kringle repeats). Up to 2,300 aa most proteins keep N/L >= 5 (ABCA1 2,261 aa: 18,775 at 0.5, L_cov 0.83).

### 3.6 EVE-protocol MSA filter inside the alignment-model wrappers (2026-11-17)

**Rules (user, 2026-11-17).** This is a benchmark of existing models: innovations minimal, the starting point is the
EVE/popEVE implementation, every model treated the same and not improved, no benchmark pairs lost to model limits (a
model that cannot score a protein gets NaN there and is reported with its coverage). Where EVE (2021 preprint/code)
and popEVE differ, take the one with higher coverage ("popEVE choice"); where both are silent, take the simplest
option ("ours"). Earlier decisions kept: ColabFold MSA as the base alignment, the same filtered alignment for every
MSA model, filter model-side (card hyper-parameters, no derived data), no MSA segments, compaction later (§10).

**Sources checked.** EVE preprint Methods (bioRxiv 10.1101/2020.12.21.423785): hmmer, 5 iterations, UniRef100;
"we explore a range of bit scores, using 0.3 bits per residue as a reference"; keep sequences aligned to >= 50% of the
target and columns with >= 70% residue occupancy; "prioritize alignments with coverage L_cov >= 0.8L and 100,000 >=
N >= 10L ... if these requirements cannot be met, we sequentially relax them down to L_cov >= 0.7L and N <= 200,000",
then "the best available" alignment; weights with Hamming cutoff T = 0.2. EVE code (`utils/data_utils.py`,
OATML-Markslab/EVE commit 460d70e, vendored verbatim in evedesign 0.0.6 `models/evemodel/`): fragments = gap fraction
> 0.5 (`threshold_sequence_frac_gaps`), focus columns = gap fraction <= 0.3 (`threshold_focus_cols_frac_gaps`),
sequences with an indeterminate residue (B, J, X, Z, ...) in a focus column removed, theta 0.2. popEVE paper (Nat.
Genet. 2025, Europe PMC full text): EVcouplings/jackhmmer, UniRef100 March 2022, T = 0.2; no alignment-selection or
column rules beyond EVE's. popEVE release (`popEVE_ukbb_20250312.zip`, 203 random proteins + APOA1, PCSK9, ABCA1):
the `EVE` column is present at 99.5% of positions (one protein all-NaN), including positions whose `gap frequency`
is up to 1.0 (in 80% of proteins the scored maximum is > 0.9); on average 24% of positions have gap frequency > 0.3,
so EVE's code default (0.3) would leave them unscored. The release therefore models all columns, i.e.
`threshold_focus_cols_frac_gaps = 1.0`, which is the zoo `eve` card's value.

**Filter (`msa_filter.py`, identical copy in every `*_filtered` folder).**

| Step | Rule | Source |
|---|---|---|
| Insertions | lowercase / '.' removed, match states only; rows of the wrong length dropped (malformed) | EVE (a2m match states) |
| Fragments | gap fraction > 0.5 dropped (`min_seq_coverage` 0.5) | EVE code + paper (50% of the target) |
| Indeterminate residues | rows with any letter other than the 20 amino acids and '-' dropped | EVE code (all columns are focus columns, see Columns) |
| Query | row 0 always kept, first | EVE / zoo wrappers |
| Bit score | bits/residue = MMseqs2 bit score (a3m header) / L; N(t), L_cov(t) for each t in `bits_thresholds` (0–1.5) | EVE paper (sweep); MMseqs2 scale is a deviation |
| L_cov | columns with >= 70% residue occupancy among the rows kept at t (`col_occupancy`) | EVE paper |
| Step 0 "strict" | L_cov >= 0.8 L and 10 L <= N <= 100,000 | EVE paper |
| Step 1 "relaxed" | L_cov >= 0.7 L and 10 L <= N <= 200,000 | EVE paper (relaxation limits) |
| Step 2 "low_depth" | L_cov >= 0.7 L, 2 <= N <= 200,000 (10 L floor dropped) | ours ("best available": coverage first, as the paper prioritises coverage) |
| Step 3 "low_lcov" | 2 <= N <= 200,000 | ours (deepest alignment) |
| Step 4 "subsample" | no t gives N <= 200,000: highest t, then uniform random subsample of 200,000 rows (query first, `seed`) | ours (last resort; EVE raises the threshold, steps 0–3 do that) |
| Which t | the lowest qualifying t (deepest alignment) within the first step that has one | ours (paper silent) |
| Depth > 100k | handled by the threshold (steps 0/1), not by subsampling | EVE |
| Columns | all kept and passed to the model; L_cov only selects t | popEVE choice (release scores 100% of positions) |
| Column use inside a model | the zoo model's own: EVE card `threshold_focus_cols_frac_gaps: 1.0` (popEVE-like); EVmutation excludes > 50%-gap columns itself | zoo cards unchanged |
| Weights | each model's own (EVE theta 0.2 as in the paper; RSALOR 0.8, StructureDCA 0.9, EVmutation 0.8) | zoo cards unchanged |
| Unscored variants | NaN (EVmutation gap-excluded columns, RSALOR no-RSA, GEMME NA, EVE non-focus) instead of the wrappers' WT/0 | user rule |

Per protein the wrapper logs and appends to `<output>/msa_filter_stats.jsonl`: rows in, rows dropped (bad length,
fragments, indeterminate), rows after the row filters, bit-score source, chosen threshold, N, N/L, L_cov (fraction),
relaxation step and rule, subsampled, quality flag (`ok`, `low_depth` = N < 10 L, `low_lcov` = L_cov < 0.7 L), the
whole sweep (N, L_cov per t), time and peak RSS.

**Deviations from EVE/popEVE (all forced by the data or the zoo, none an improvement):**

| Deviation | Reason |
|---|---|
| ColabFold/MMseqs2 search (UniRef100 + environmental DBs, up to 1.16M rows) instead of EVcouplings/jackhmmer, 5 iterations, UniRef100 2022 | user data decision (2026-11-09, PG2 protocol); re-searching is a data change |
| Bit scores on the MMseqs2 scale (bits/residue = header bit score / L), thresholds swept 0–1.5; EVE's 0.3 bits/residue reference does not transfer (an exact self-hit of APOA1 scores 0.53 bits/residue) | MMseqs2 instead of HMMER |
| Thresholds are post-hoc cuts of one search, not separate jackhmmer runs at each bit score (EVE re-searches, so each threshold also changes the profile) | one ColabFold search per protein |
| (resolved) proteingym-base `MSA.dump` writes `seq_<i>` headers, so the first archives had no bit scores (rule `no_bitscores`: row filters + 200k cap only) | **user decision 2026-11-17: archives keep the ColabFold headers** (data branch `scripts/clinical/msa/restore_msa_headers.py`, job 31939179 -> `genebass/splits_layers_hdr/`, sequences verified identical and in order) |
| `a3m_headers_dir` (original ColabFold a3m, matched by query, every row verified) | test fallback only (used for the "bits" validation runs before the header-restored archives existed); never in benchmark runs |
| Steps 2–4 and "lowest qualifying t" | paper silent ("best available") |
| Columns all kept (EVE code default 0.3 would drop ~24% of positions) | popEVE choice |
| Per-model weighting and column handling unchanged | fairness: models not modified |

**Bit scores in the wrapper.** proteingym-base keeps the raw a3m bytes (`MSA.file`), but
`evedesign.proteingym.msa_to_sequences` builds `Sequence(id=None)`. Each `*_filtered` package's `__init__.py` calls
`install_header_hook()`, which wraps that function (the shared runner imports the model package before converting the
dataset) and attaches the per-row bit scores (`pg_bitscores`, parsed from the raw a3m bytes `MSA.file` of the
archive; one header per loaded sequence, same order, count checked) to the `Sequences`; ids stay None, so the wrapped
models see exactly what they saw before. evedesign's `dataset_to_evedesign` does not keep headers or pass `MSA.file`
(checked, evedesign 0.0.5/0.0.6 `proteingym.py`), hence the hook. Archives with `seq_<i>` headers give all-NaN scores,
which the filter treats as "no bit scores". `a3m_headers_dir` (default null) is a test fallback only.

**NaN in the pipeline.** Checked on APOA1: a NaN prediction survives `predictions_delta` -> `dump` -> reload as a
float NaN (not null), and proteingym-base's `metric_spearman` then returns NaN for that dataset (scipy propagates).
The clinical export (`export_scores.py` -> Genebass evaluation) drops NaN rows (missing for that model, coverage
reported). PG2 metrics would need a NaN policy before benchmark runs (open decision 5, §8).

### 3.7 Inverse-folding models: ProteinMPNN and PottsMPNN (2026-11-19)

User request 2026-11-19: inverse-folding models beyond ESM-IF1, checked up to 2,700 aa and beyond. Rule as for the
other structure models: one full-length AFDB structure (every protein <= 2,700 aa) -> the zoo model unchanged; AFDB
fragments (> 2,700 aa: APOB, LRP2, PKD1) -> most-central fragment (`fragments.py`, copied verbatim, §3.2).

- **ProteinMPNN (`models/proteinmpnn`, new zoo folder).** `evedesign.models.mpnn:LigandMPNN` with ProteinMPNN
  weights `v_48_020` (48 neighbours, 0.20 A noise; upstream default of both ProteinMPNN and LigandMPNN), from
  dauparas/ProteinMPNN commit 8907e66 (byte-identical to files.ipd.uw.edu/pub/ligandmpnn/proteinmpnn_v_48_020.pt),
  sha256 pinned, downloaded at image build. Score = mean log-likelihood of the sequence given the backbone, decoded
  in one random order (`score(..., use_sequence=True)`), the quantity of ProteinGym's ProteinMPNN baseline
  (`compute_fitness.py`, which draws a new random order per mutant). The wrapper's `fix_full_decoding_order` draws
  the order with an unseeded `torch.randn`; the only addition (`model.py`, class `ProteinMPNN`) seeds torch before
  build (`seed: 0`), so scores are reproducible. Same proteingym-base commit as the eve/pottsmpnn/esm-if1 images
  (pinned in `pyproject.toml`).
- **`proteinmpnn_window`:** one structure -> zoo unchanged. Fragments: each fragment renumbered onto the wild type
  and the zoo model built on it (own backbone, own decoding order from `seed`); mutant scored in its most-central
  fragment as n_fragment x (mean_ll(mutant) - mean_ll(wt)) (summed log-likelihood difference, as esm-if1_window).
- **`pottsmpnn_window`:** one structure -> zoo unchanged. Fragments: the zoo's Potts energy table per fragment
  (same `get_etab` call), mutant scored as E(wt) - E(mutant) with the table of its most-central fragment
  (energies are additive, no rescaling). Reference sequence = entity wild type over the fragment (APOB's AFDB
  fragment carries C at 1422, wild type Y).
- Mutants not covered by any fragment: NaN (none in Genebass).
- **Length:** both are kNN graphs (48 neighbours): memory and time linear in L; PottsMPNN measured elsewhere up to
  ~21k aa on 24 GB. Both run on CPU.
- Other inverse-folding / structure models: ESM-IF1 (`esm-if1_window`, §3.2) and SaProt (`saprot_window`) done.
  ProteinGym's other structure-based baselines (MIF-ST, ProSST, ...) have no zoo folder or evedesign wrapper; not
  added (would be new zoo models, not long-protein handling).

**Tests and Genebass runs (2026-11-19).** Images `proteinmpnn.sif`, `pottsmpnn.sif` (`build_images.sbatch`, job
32024839, 6 min each); window packages mounted in the zoo images (`validate_msa_filter.sbatch`, `run_genebass.sbatch`,
CPU, genoa64). Outputs `$WORK/models/mpnn/{val,genebass}/`.
- T1 (APOA1, 328 variants, two separate processes): window class == zoo class, max |diff| 0 for both models (also
  shows the seeded ProteinMPNN order is reproducible). Jobs 32025090/91, 32025227/28.
- Genebass, all 66 proteins / 5,840 unique mutants, 0 unscored (APOB, LRP2, PKD1 through fragments). Jobs
  32025130–32 (proteinmpnn_window, 142 CPU-min in total), 32025254 (pottsmpnn_window, 29 min). Predictions
  `proteingym_clinical/genebass/predictions/v2__{proteinmpnn,pottsmpnn}_window/` (+ `status.tsv`).

| Protein (L) | ProteinMPNN score time / peak RSS | PottsMPNN build+score / peak RSS |
|---|---|---|
| ABCA1 (2,261), 279 mutants | 12 min / 10.1 GB | 9 s / 9.3 GB |
| STAB1 (2,570), 199 mutants, one structure | 11 min / 8.8 GB | 8 s / 8.4 GB |
| PKD1 (4,303), fragments | 12 min / 13.3 GB | 24 s / 16.3 GB |
| APOB (4,563), fragments, 507 mutants | 15 min / 2.3 GB | 31 s / 6.0 GB |

Peak RSS includes the archive's MSA, which the runner loads although these models do not use it. ProteinMPNN cost
is one forward pass per mutant on CPU: ~1.7 s at 1,400 residues (APOB fragments), 2.5–3.3 s at 2,261–2,570.
- PottsMPNN vs the prototype `pottsmpnn_clinical` (zoo image, fragments most-central): Pearson 1.0, max |diff| 1e-4,
  63 proteins in common (the prototype missed CST3, GPLD1, UGT1A5).
- ProteinMPNN vs PottsMPNN per protein: Spearman median 0.85 (0.71–0.96), on APOB/LRP2/PKD1 (fragments) 0.87–0.88,
  STAB1 0.90. ProteinMPNN vs ESM-IF1 median 0.71.
- Conclusion: both models are fine at any length in this set: full length up to 2,700 aa (zoo unchanged), AFDB
  fragments above.

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
| `models/rsalor_filtered` | `RSALORFiltered` (`filtered_model.py`) | zoo `rsalor` (`model.py` copied verbatim) |
| `models/structuredca_filtered` | `StructureDCAFiltered` | zoo `structuredca` (`model.py` copied verbatim) |
| `models/evmutation_filtered` | `EVmutationFiltered` | `evedesign.models.evcouplings:EVcouplingsPLM` |
| `models/eve_filtered` | `EVEFiltered` | `evedesign.models.eve:EVE` (evedesign 0.0.6) |
| `models/evmutation2_filtered` | `EVmutation2Filtered` | `evedesign.models.evmutation2:EVmutation2` (not tested: ~1,000 aa GPU limit) |
| `models/gemme_filtered` | `GEMMEFiltered` | zoo `gemme` (`model.py` copied verbatim; not tested: image does not build) |
| `models/proteinmpnn` | `ProteinMPNN` (seeds the decoding order) | `evedesign.models.mpnn:LigandMPNN`, weights v_48_020 (new zoo folder, §3.7) |
| `models/proteinmpnn_window` | `ProteinMPNNWindow` | `models/proteinmpnn` (fragments, §3.7) |
| `models/pottsmpnn_window` | `PottsMPNNWindow` | zoo `pottsmpnn` (fragments, §3.7) |

Shared logic is duplicated per folder to keep folders self-contained: `window.py` (PLM windows, `esm2_window`,
`esmc_window`), `fragments.py` (structure placement and most-central assignment, `saprot_window`, `esm-if1_window`),
`msa_filter.py` (EVE-protocol MSA filter, NaN for unscored variants, header hook; all six `*_filtered` folders, §3.6).
The `*_filtered` folders are copies of the zoo folder (same Dockerfile/Singularity.def/lockfile, console script and
package renamed), plus `__init__.py` (installs the header hook), `msa_filter.py`, `filtered_model.py`
(`class XFiltered(MSAFilterMixin, X)`: filter keys popped from the card's `hyper_parameters`, the rest go to the zoo
class unchanged) and `tests/test_msa_filter.py` (synthetic alignments: row filters, threshold choice, relaxation
order, subsampling, missing bit scores). Validation: `scripts/clinical/models/msa_filter/validate_msa_filter.{py,sbatch}`
(one build + score per job in the zoo image with the package mounted) and `compare_msa_filter.py` (local summary,
popEVE comparison); outputs in `$WORK/models/msa_filter/val/`.
Validation scripts and job wrappers: `scripts/clinical/models/validate_{plm,struct}_windows.{py,sbatch}` (run in the
prototype's zoo images with the new package mounted on `PYTHONPATH`; proper images still to build). Outputs in
`$WORK/models/`, job IDs in `$WORK/logs/jobs.txt`. Upstream route later: windowed `score_mutants` in evedesign's `ESM2` would fix `esm2` and
`esmc` at once.

## 6. Running on the CRG cluster

- Containers: Apptainer, built from each folder's `Singularity.def` with `apptainer build --fakeroot` on the GPU
  node; the prototype's images are in `project/proteingym_clinical/containers/`.
- GPU: partition `gpu`, `--gres=gpu:7g.80gb:1` (full H100). All compute as sbatch, never in the interactive session.
  GPU jobs are limited to 12 h (`--time` above 12:00:00 is rejected for every QoS) and by a group GRES limit
  (`AssocGrpGRES`: a third 3g.40gb job waits). EVE trains 400k steps at ~79 it/s (L 267) and ~38 it/s (L 692) on
  3g.40gb.
  GPU jobs are limited to 12 h (`--time` above 12:00:00 is rejected for every QoS) and by a group GRES limit
  (`AssocGrpGRES`: a third 3g.40gb job waits). EVE trains 400k steps at ~79 it/s (L 267) and ~38 it/s (L 692) on
  3g.40gb.
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

1. ~~Special tokens at window edges~~ **Decided 2026-11-16 (user): `termini`**, `<cls>`/`<eos>` only at the real protein ends, so the model is never told a false start or end. `every_window` stays available for comparability with ProteinGym v1 / the prototype (benchmark effect: ESM2 none, ESMC 2.181 vs 2.200).
2. ESMC: 1,022 for parity only, or also a native 2,048 card.
3. Unmasked/wt-marginal cards for the PLMs in addition to masked.
4. Alignment models: region policy and depth caps (after the audit and T5).
5. Missing predictions in the metric: drop the variant for that model (current evaluation) or penalise coverage.
6. ~~Region archives as derived datasets~~ **decided 2026-11-16 (user):** no derived datasets; region policy inside
   the alignment-model wrappers, declared on the card (§3.4, §10). Region archives stay only for T5/first runs.
7. ~~Bit scores for the EVE alignment selection~~ **decided 2026-11-17 (user):** archives keep the ColabFold
   headers (data branch, `genebass/splits_layers_hdr/`, replaces `splits_layers` once verified); `*_filtered` reads
   the bit scores from the archive's raw a3m (§3.6).
8. EVE columns: popEVE-like all columns (`threshold_focus_cols_frac_gaps: 1.0`, zoo card) kept; EVE 2021 code
   default 0.3 would leave ~24% of positions unscored (§3.6).
9. ~~Inverse-folding models beyond ESM-IF1~~ **done 2026-11-19 (user request):** ProteinMPNN (new zoo folder, v_48_020,
   seeded decoding order) and PottsMPNN, zoo unchanged up to 2,700 aa, most-central AFDB fragment above (§3.7).

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
| 2026-11-17 | `*_filtered` alignment models (EVE-protocol MSA filter, §3.6) + validation (rsalor, structuredca, evmutation, eve) | T1 PASS on 4 models; filter stats, correlations, time/memory, EVE vs popEVE: "MSA filter results" below; archives carry no bit scores (open decision) | `$WORK/models/msa_filter/val/` |
| 2026-11-19 | ProteinMPNN (new zoo folder) + `proteinmpnn_window`, `pottsmpnn_window` (§3.7): T1 on APOA1, Genebass 66 proteins | T1 max diff 0 (both); 66/66 proteins, 0 unscored; PottsMPNN = prototype (r 1.0) | `$WORK/models/mpnn/` |

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

### T5 results: alignment models on region archives (2026-11-17)

`collect_regions.py summary --models rsalor structuredca evmutation evmutation2 eve --tag t5`
(full table: `$WORK/models/regions/t5_summary.tsv`). Region archives from `make_regions.py` (§3.4); unchanged zoo
models; CPU jobs on genoa64, evmutation2 on a GPU slice. Coverage = fraction of the protein's Genebass mutants with a
prediction. Time = sum over regions; RSS = peak per region.

| Model | Protein (L) | W | Regions ok | Coverage | Time total (max region) | Peak RSS |
|---|---|---|---|---|---|---|
| rsalor | APOB (4,563) | 500 / 1,000 | 18/18, 9/9 | 1.00 / 1.00 | 0.4 / 0.3 min | 0.4 GB |
| rsalor | LRP2 (4,655) | 500 / 1,000 | 18/18, 9/9 | 1.00 / 1.00 | 1.7 / 1.5 min | 0.8 GB |
| rsalor | PKD1 (4,303) | 500 / 1,000 | 17/17, 8/8 | 1.00 / 1.00 | 1.0 / 0.7 min | 0.8 GB |
| structuredca | APOB | 500 / 1,000 | all | 1.00 / 1.00 | 3.0 / 4.1 min | 0.6 GB |
| structuredca | LRP2 | 500 / 1,000 | all | 1.00 / 1.00 | 15.7 / 18.0 min (5.0) | 0.9 GB |
| structuredca | PKD1 | 500 / 1,000 | all | 1.00 / 1.00 | 4.4 / 4.6 min | 0.7 GB |
| evmutation | APOB | 500 / 1,000 | all | 0.93 / 0.97 | 31 / 59 min (16) | 3.8 / 19.3 GB |
| evmutation | LRP2 | 500 / 1,000 | all | 0.94 / 0.92 | 225 / 373 min (102) | 4.1 / 15.0 GB |
| evmutation | PKD1 | 500 / 1,000 | all | 0.93 / 0.87 | 128 / 168 min (57) | 4.0 / 14.9 GB |
| evmutation2 (GPU) | APOB, LRP2, PKD1 | 500 | all | 1.00 | 22–25 min (1.6) | 12.9 GB host |
| evmutation2 (GPU) | APOB, LRP2, PKD1 | 1,000 | 1/9, 0/8, 0/9 | 0.13 / 0 / 0 | — | CUDA OOM |
| eve (GPU) | PCSK9 (692) | 250 | 1/5 run | 0.32 | 92 min per region | 4.2 GB host |

- **Every model fits the cluster at W = 500** with coverage 1.00 (rsalor, structuredca, evmutation2) or ≥ 0.93
  (evmutation: its gap rule excludes columns, those mutants are NaN, cross-checked against the rule). Depth after
  re-filtering: median 1.8k (APOB), 10–14k (PKD1), 29–35k (LRP2) sequences per region.
- **evmutation2 at W = 1,000 does not fit a GPU slice** (CUDA OOM, tries to allocate 14.9 GB+); W = 500 does.
- **eve** is too slow for routine use: 92 min for one PCSK9 region of 250 on 3g.40gb; not run on long proteins.
- **gemme**: not tested (image build fails, lcqb.upmc.fr download returns empty replies).

Region size vs full length on PCSK9 (692 aa, fits whole: W = 10,000 is the full-length run), Pearson on the 78
Genebass mutants:

| Model | W500 vs full | W250 vs full | W250 vs W500 |
|---|---|---|---|
| rsalor | 0.970 | 0.884 | 0.888 |
| structuredca | 0.966 | 0.864 | 0.878 |
| evmutation | 0.935 (n 45) | 0.906 (n 53) | 0.932 (n 49) |
| evmutation2 | 0.930 | 0.808 | 0.813 |

W = 500 reproduces the full-length scores at r ≈ 0.93–0.97; W = 250 loses noticeably (0.81–0.91). Default for the
wrapper refactor: region 500 (evmutation2 needs ≤ 500 on a GPU slice anyway); 1,000 for CPU models if wanted.

Full runs over all 66 proteins were dropped by user decision (2026-11-17: a few tests are enough). The one finished
full run (rsalor, W = 1,000) stays exported as `predictions/v2__rsalor`.

### MSA filter results (2026-11-17)

`*_filtered` models (§3.6) in the zoo images with the package mounted (`validate_msa_filter.sbatch`, one build + score
per job, jobs 31937458–31937558 in `$WORK/logs/jobs.txt`), summary by `compare_msa_filter.py`
(`$WORK/models/msa_filter/val/{runs,pairs}.tsv`). Scored: all single substitutions (StructureDCA: a seeded sample of
2,000 + the archive's Genebass mutants), pred = score − score(WT). Variants: **zoo** = unchanged zoo model on the
unfiltered MSA (2 h cap); **off** = `*_filtered` with `msa_filter: false`; **filtered** = the card on the first
archives, whose headers were `seq_<i>` (now `genebass/splits_layers_noheaders/`): row filters only, rule
`no_bitscores`; **bits** = the full EVE protocol with bit scores from the original ColabFold a3m
(`a3m_headers_dir`, test fallback); **hdr** = the card on the header-restored archives (`genebass/splits_layers/`
since 2026-11-17 18:34), bit scores read from the archive itself. All runs except `hdr` (and the ABCA1 EVE scoring
rerun, if needed) used the old archives; sequences and order are identical in both.

**Filter per protein** (rows in -> fragments / indeterminate dropped -> rows after row filters; then the threshold):

| Protein (L) | Variant | Rows in | Fragments / indeterminate | After row filters | t (bits/res) | N | N/L | L_cov | Rule | Quality |
|---|---|---|---|---|---|---|---|---|---|---|
| APOA1 (267) | filtered = bits | 5,882 | 268 / 55 | 5,559 | 0.0 | 5,559 | 20.8 | 0.95 | strict | ok |
| PCSK9 (692) | filtered | 531,387 | 430,931 / 175 | 100,281 | — | 100,281 | 144.9 | 0.48 | (no bit scores) | low_lcov |
| PCSK9 (692) | bits | | | | 0.45 | 1,410 | 2.0 | 0.87 | low_depth | low_depth |
| LDLR (860) | filtered | 656,624 | 563,104 / 831 | 92,689 | — | 92,689 | 107.8 | 0.78 | (no bit scores) | ok |
| LDLR (860) | bits | | | | 0.5 | 17,000 | 19.8 | 0.88 | strict | ok |
| IGF1R (1,367) | bits | 1,035,677 | 1,011,733 / 520 | 23,424 | 0.3 | 12,356 | 9.0 | 0.90 | low_depth | low_depth |
| ABCA1 (2,261) | filtered | 1,100,277 | 1,066,651 / 911 | 32,715 | — | 32,715 | 14.5 | 0.66 | (no bit scores) | low_lcov |
| ABCA1 (2,261) | bits | | | | 0.4 | 21,904 | 9.7 | 0.71 | low_depth | low_depth |

(In `runs.tsv` the `filtered` rows show rule `low_lcov`/`relaxed` at t = 0: those runs predate the fix that reports
all-NaN header scores as `no_bitscores`; the selected rows are the same.) The fragment filter alone removes 82–98% of
the rows of the deep ColabFold MSAs. With bit scores the protocol mostly lands in "relaxed"/"low_depth": the
10 L depth floor is the criterion that fails (PCSK9: coverage >= 0.7 L only at t >= 0.45 with 1,410 rows).

**T1 (filter disabled = zoo):**

| Model | Protein | Max abs diff off vs zoo | Note |
|---|---|---|---|
| rsalor | APOA1 | 0 (5,073 singles) | **PASS** |
| structuredca | APOA1 | 1.2e-4 (2,017) | **PASS** (solver numerics) |
| evmutation | APOA1 | 0 on the 5,016 scored; 57 NaN | **PASS**: the 57 NaN are mutations in gap-excluded columns, scored as WT by the zoo |
| eve | APOA1 (2,000 training steps, 50 samples) | 0 (5,073) | **PASS** |

**Filtered vs unfiltered** (all scored singles; Spearman on the archive's Genebass mutants in brackets):

| Model | Protein | filtered vs zoo | bits vs zoo | filtered vs bits |
|---|---|---|---|---|
| rsalor | APOA1 | 0.994 (0.993) | = filtered | 1.000 |
| rsalor | PCSK9 | 0.987 (0.972) | 0.918 (0.902) | 0.943 (0.953) |
| rsalor | LDLR | 0.998 (0.997) | 0.990 (0.992) | 0.994 (0.992) |
| rsalor | ABCA1 | 0.999 (0.998) | 0.999 (0.997) | 1.000 (0.999) |
| structuredca | APOA1 | 0.993 (0.980) | = filtered | — |
| structuredca | PCSK9 | 0.966 (0.965) | 0.791 (0.714) | 0.840 (0.761) |
| structuredca | LDLR | 0.982 (0.972) | 0.925 (0.910) | 0.936 (0.917) |
| structuredca | ABCA1 | zoo timed out (> 2 h) | — | 0.977 (0.969) |
| evmutation | APOA1 | 0.9996 (0.997) | = filtered | — |

**Time and memory** (build = filter + model fit; score = all singles or 2,000 for StructureDCA; peak RSS of the job
from sacct, dominated by loading the archive's MSA: 1.1M rows ~10 GB for ABCA1):

| Model | Protein (L) | Variant | N | Build (min) | Score (min) | Job peak RSS (GB) | Unscored (NaN) |
|---|---|---|---|---|---|---|---|
| rsalor | PCSK9 (692) | zoo / filtered / bits | 531k / 100k / 1.4k | 1.0 / 0.2 / 0.1 | 0.2 | 3.4 / 2.9 / 2.8 | 0 |
| rsalor | LDLR (860) | zoo / filtered / bits | 657k / 93k / 17k | 3.9 / 0.3 / 0.3 | 0.3 | 4.6 / 3.7 / 3.7 | 0 |
| rsalor | ABCA1 (2,261) | zoo / filtered / bits | 1.1M / 33k / 22k | 10.2 / 0.8 / 0.8 | 2.3 | 16.5 / 13.1 / 13.0 | 0 |
| structuredca | PCSK9 | zoo / filtered / bits | | 35.2 / 7.2 / 0.4 | 0.7 | 3.7 / 2.8 / 2.8 | 0 |
| structuredca | LDLR | zoo / filtered / bits | | 35.6 / 7.8 / 2.7 | 1.1 | 4.7 / 3.8 / 3.7 | 0 |
| structuredca | ABCA1 | zoo / filtered / bits | | > 120 (timeout) / 15.6 / 12.0 | 8.0 | 21.8 / 13.2 / 13.0 | 0 |
| evmutation | LDLR (860) | bits | 17,000 | 11.5 | 0.6 | 14.5 | 589 of 16,340 (3.6%) |
| evmutation | IGF1R (1,367) | bits | 12,356 | 37.5 | 4.2 | 34.3 | 760 of 25,973 (2.9%) |
| eve | APOA1 (267) | filtered (3g.40gb) | 5,559 | 93.3 (400k steps, 79 it/s) | 6.6 | 3.4 | 0 |

**EVE vs the popEVE release `EVE` column** (same sequence as the release's RefSeq protein; the release score is
higher = more damaging, ours is an ELBO difference, so the expected sign is negative):

| Protein | Variant | N | Spearman all singles | Spearman Genebass mutants | Time |
|---|---|---|---|---|---|
| APOA1 (267) | filtered (= bits) | 5,559 | −0.854 (5,073) | −0.735 (30) | 1 h 41 min |
| PCSK9 (692) | bits | 1,410 (Neff 112) | −0.881 (13,148) | −0.887 (61) | 3 h 42 min (build 185 min at 38 it/s, score 36 min; 8.7 GB) |
| ABCA1 (2,261) | bits | 21,904 (N/L 9.7) | −0.902 (42,959) | −0.893 (189) | ~18 h: training 31937558 timed out at 12 h after saving the checkpoint, scoring 31940397 6.1 h (43 GB host RAM) |

Reading:
- The plumbing is exact (T1 on all four models) and every position stays scorable (no columns dropped; NaN only
  where the model itself does not model a column: EVmutation 3–4% of singles).
- Without bit scores (the archives today) the filter keeps all non-fragment rows and changes little
  (Spearman vs zoo 0.97–0.999) while cutting StructureDCA time 4–5x and making ABCA1 feasible (zoo > 2 h).
- With bit scores the EVE protocol can cut the alignment hard where coverage is patchy (PCSK9: 100k -> 1.4k rows,
  Spearman vs zoo 0.79–0.92). Whether that helps is a benchmark question (Genebass), not answered here.
- **Header-restored archives (`hdr`, jobs 31940394/31940395, rsalor on PCSK9 and LDLR):** `bitscore_source =
  dataset a3m headers`, same threshold and N as `bits`, predictions identical (max abs diff 0). The card needs no
  `a3m_headers_dir` on these archives.
- EVE with the protocol alignment agrees with the popEVE release at Spearman 0.85–0.90 (APOA1, PCSK9, ABCA1; one seed,
  our ColabFold alignment vs their jackhmmer alignment), even for PCSK9 where the protocol keeps only 1,410 rows.
  That is the alignment effect, not the code: see "EVE/EVmutation reproduce the released scores" below.
- Memory is dominated by loading the archive's MSA into Python (proteingym-base + evedesign objects), not by the
  models; filtered models on ABCA1 need ~13 GB, IGF1R EVmutation 34 GB (plmc on 1,367 columns).

**EVE / EVmutation reproduce the released scores (2026-11-18).** Goal: separate the implementation from the alignment.
The zoo cards, unchanged, were run on the exact ProteinGym v1 alignments behind ProteinGym's released scores
(`DMS_msa_files.zip`, ProteinGym v1.3 download; the a2m read as EVE reads it with `threshold_focus_cols_frac_gaps=1`:
every column upper-cased, `.` -> `-`, fixed width; put into the ProteinGym2 draft archive in place of its MSA). ProteinGym's EVE
settings match the zoo card (focus-gap threshold 1, default parameters, 20,000 samples, theta 0.2); it averages
seeds 0–4 (`EVE_ensemble`), the zoo runs one seed. Jobs 31967282–85, 31974871; outputs, swap and compare scripts in
`$WORK/models/eve_repro/` (`val/compare_pgv1.tsv`). Spearman on the assay's single mutants:

| | GRB2 (217 aa, 33k seqs) | P53 (393 aa, 4.1k seqs, Neff ~144) |
|---|---|---|
| ours seed 42 vs released EVE_single / EVE_ensemble | 0.990 / **0.997** | 0.887 / 0.888 |
| ours seed 0 vs released EVE_single / EVE_ensemble | | 0.940 / **0.951** |
| ours seed 42 vs ours seed 0 (1,048 assay mutants / all 7,467 singles) | | 0.934 / 0.982 |
| released EVE_single vs EVE_ensemble (noise bar) | 0.995 | 0.960 |
| vs DMS: ours (s42, s0) / released single / ensemble | 0.579 / 0.586 / 0.590 | 0.522, 0.577 / 0.555 / 0.576 |
| ours EVmutation vs released EVmutation | 0.996 (638 mutants at a2m upper-case columns) | 0.972 (Pearson 0.987) |

- EVE reproduces within seed noise. On a deep alignment one seed matches the ensemble (0.997). On a shallow one,
  single seeds scatter by 0.89–0.95 among themselves, ours and ProteinGym's alike.
- EVmutation reproduces on the columns EVcouplings models. GRB2 positions 185–214 are lower case (non-focus) in
  ProteinGym's a2m, yet the released EVmutation has real scores there (Spearman 0.95 with its EVE). They come from
  some other run, which we cannot reproduce from this file (ours vs released 0.05 there). P53: our input keeps the 54
  non-focus columns that EVcouplings drops.
- Released single seed vs 5-seed ensemble over 96 human ProteinGym assays: median 0.987, min 0.876.
- Release vs release, as a yardstick for the alignment effect: ProteinGym `EVE_ensemble` vs the popEVE release `EVE`
  on 44 shared human proteins (full-length alignments; `$WORK/models/eve_repro/proteingym_vs_popeve_eve.tsv`): median
  Spearman 0.937, 10th percentile 0.847, 16 of 44 below 0.90. Our 0.85–0.90 vs popEVE falls in this range.
- Implication: a single EVE seed is noisy on shallow alignments (several Genebass proteins have N/L ~2–10). ProteinGym
  averaged 5 seeds (seeds 0-4); doing so too would cost 5x the GPU time.

## 10. Picking this up

**Progress note (inverse-folding models, 2026-11-19):** ProteinMPNN and PottsMPNN done for every Genebass protein
(§3.7): images built, T1 passed, predictions in `predictions/v2__{proteinmpnn,pottsmpnn}_window/`. Next: add them to
`17_eval_common_variants.py` (they score every variant, so they do not shrink the common set).

**Progress note (MSA filter, 2026-11-17):** `*_filtered` folders built and validated on a few proteins (§3.6, §9
"MSA filter results"). **Next:**
1. **Bit scores in the archives: done.** Header-restored archives are `genebass/splits_layers/` (old ones in
   `splits_layers_noheaders/`); checked (jobs 31940394/31940395): bit scores read from the archive, predictions
   identical to the `bits` runs. Upstream: `MSA.dump` in proteingym-base still drops headers (report it).
   EVE ABCA1 done (−0.90 vs popEVE). EVE and EVmutation reproduce ProteinGym's released scores on its own
   alignments (§9 "EVE / EVmutation reproduce the released scores"). Seeds: the zoo card's single seed is kept (user
   2026-11-19: not our call; the number of seeds belongs to the model card, not the benchmark).
2. **Compaction for proteins > 2,300 aa (EVE, EVmutation):** EVE builds the one-hot training set as float64
   N x L x 20 and, for scoring, the whole mutant set at once (43k x 2,261 x 20 float64 = 15.5 GB for ABCA1); build it
   per batch from an int8 N x L matrix (training) and per scoring batch; weights in blocks. EVmutation: plmc memory
   grows with L^2 q^2 (couplings), measure on one > 2,300 aa protein after compaction. Not done (user: later step).
3. **RSALOR / StructureDCA above 2,700 aa:** they need one full-length single-chain structure, AFDB gives fragments
   only. Decide between most-central fragment (consistent with the structure models, §3.2) and NaN with coverage
   reported (benchmark rule: no new tricks); the region tiling of §3.4 is not used (user decision).
4. EVE: full EVE runs on all proteins <= 2,300 aa only after (1); compare with the popEVE release per protein.
5. NaN policy for the PG2 metrics (§3.6 "NaN in the pipeline", §8 item 5).
6. Images for the `*_filtered` folders (`build_images.sbatch rsalor_filtered ...`); validation ran in the zoo
   images with the package mounted. GEMME still has no image (download server), `gemme_filtered` untested;
   `evmutation2_filtered` untested (~1,000 aa GPU limit, T5).

**Progress note (T5, alignment models, 2026-11-17):** T5 done (§9, T5 results). Region code in
`scripts/clinical/models/regions/`, archives in `$WORK/models/regions/archives/`, predictions in
`$WORK/models/regions/preds/<model>/t5/`. No full runs (user decision). **Next:** move the region logic inside the
alignment-model wrappers (card declares `region`, `min_coverage`, `max_depth`, `seed`; no derived datasets), default
region 500.

State on 2026-11-16 (evening):
- **Done, validated (T1–T3), benchmarked:** `esm2_window`, `esmc_window`, `saprot_window`, `esm-if1_window`
  (§9). They run in the prototype's zoo images (`project/proteingym_clinical/containers/*.sif`) with the new package
  mounted (`scripts/clinical/models/validate_{plm,struct}_windows.sbatch`). Their own images are not built yet
  (`build_images.sbatch esm2_window ...`).
- **Data side:** final archives `$WORK/genebass/splits_layers/` (MSAs, layers from AlphaMissense/popEVE releases,
  `in_vesm_paper` flags); the model runs used the no-MSA copies in `$WORK/genebass_api/build/splits/` (same rows).
- **Next, in order:**
  1. Alignment-model images: job 31923582 (`$WORK/logs/build_images_*.log`); rebuild `gemme` (download failed).
  2. (T5 DONE, §9) T5 on APOB, LRP2, PKD1 + one mid-size protein: region archives (§3.4: crop WT, MSA match columns, structure
     fragment; re-filter rows by coverage in the region; random depth cap; weights once per region), run the
     unchanged zoo models on them, map predictions back. Then NaN instead of WT for unscored positions (EVmutation,
     RSALOR, GEMME) and per-protein failure handling.
  3. Port PottsMPNN with `fragments.py` (most-central fragment) and BLOSUM62 unchanged; decide `special_tokens`
     (§8) — ESMC is the only model where it matters (benchmark 2.181 termini vs 2.200 every_window).
  4. Open: residual vs dark-energy ESM2 scores on long proteins (max 1.3, r 0.99997).
