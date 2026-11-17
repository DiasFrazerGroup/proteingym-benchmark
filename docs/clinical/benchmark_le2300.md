# Genebass benchmark restricted to proteins ≤ 2,300 aa (alternative set "le2300")

Status 2026-11-17. Companion to [models_long_proteins.md](models_long_proteins.md) (§3.5: what EVE/popEVE did).
Principle (user, 2026-11-17): this benchmarks existing models; innovations minimal; the starting point is the
EVE/popEVE implementation; be fair without losing the benchmark to model limitations.

## 1. Why 2,300 aa

- **Precedent, not our choice:** the popEVE v1.1 release (18,331 proteins) has no protein longer than 2,302 aa; EVE
  2021 "avoided extremely large genes such as TTN". 2,300 aa is the range in which EVE was actually run on the
  human proteome, at full length, with no segments.
- **Structures:** AFDB gives one full-length model up to 2,700 aa. All 61 proteins ≤ 2,300 aa have exactly one
  structure in their archive, so the fragment placement (used by `saprot_window`, `esm-if1_window`) is never
  triggered in this set. User decision 2026-11-17: keep the AFDB fragment placement for the longer proteins.
- **popEVE coverage in this set:** 97 / 97 pairs, 99.8% of rows (5,734 / 5,744). popEVE's score on all 109 pairs
  (2.372) is identical to its score here, because the 12 pairs it misses are exactly the proteins > 2,300 aa.

## 2. What the set contains

| | All | ≤ 2,300 aa | Removed |
|---|---|---|---|
| Proteins | 66 | 61 | 5 |
| Pairs | 109 | 97 (89%) | 12 |
| Rows (variants) | 7,791 | 5,744 (73.7%) | 2,047 |
| Paper set pairs (`in_vesm_paper`) | 103 | 91 | 12 |
| Lipid Metabolism / Liver / Renal / Other pairs | 39 / 38 / 15 / 17 | 33 / 37 / 11 / 16 | 6 / 1 / 4 / 1 |

Removed: APOB (4,563 aa, 5 pairs, lipid), LRP2 (4,655, 3, renal), PKD1 (4,303, 1, renal), SPTA1 (2,419, 1, other),
STAB1 (2,570, 2, lipid/liver). Renal Function loses 4 of 15 pairs: the group means are less stable there.

## 3. Does the ranking survive the restriction?

Existing per-pair results (`proteingym_clinical/genebass/results/genebass_source_{headline_all,le2300}_*.tsv`, 36 models
and baselines with ≥ 93 pairs): Spearman of model means, all pairs vs ≤ 2,300 aa = **0.963**.

| Model | All 109 | ≤ 2,300 aa (97) |
|---|---|---|
| AlphaMissense (layer) | 2.695 | 2.715 |
| VESM++ (paper scores, 93 pairs) | 2.671 | 2.652 |
| esm2_window, termini | 2.403 | 2.356 |
| popEVE (layer) | 2.372 (97 pairs) | 2.372 |
| rsalor (T5 region run, W 1,000; to be replaced by the filtered full-length model) | 2.194 | 2.233 |
| esmc_window, termini | 2.181 | 2.226 |
| CADD | 1.912 | 1.969 |
| Allele frequency | 1.629 | 1.746 |
| saprot_window | 1.640 | 1.625 |
| esm-if1_window | 1.273 | 1.280 |

The top tier and the bottom are unchanged. Mid-table models within ~0.1 swap places (ESM2 loses ~0.05 without the
long proteins, ESMC gains ~0.04): differences of that size should not be read as rankings in either set.

## 4. Can every model see every protein ≤ 2,300 aa?

Alignment input: ColabFold MSA filtered with the EVE protocol (≥ 50% sequence coverage, bit-score sweep, ≤ 100k
sequences; implementation in progress, §3.5–3.6 of models_long_proteins.md). Filtered depth over the 61 proteins:
median 19k, max 99k sequences.

| Model | ≤ 1,000 aa (47 proteins, 78 pairs) | 1,000–2,300 aa (14 proteins, 19 pairs) | Evidence / cost |
|---|---|---|---|
| esm2_window, esmc_window | yes | yes | done, benchmarked |
| saprot_window, esm-if1_window | yes | yes (one AFDB structure each) | done, benchmarked |
| PottsMPNN, BLOSUM62 | yes | yes (one structure) | not ported yet; no length issue ≤ 2,700 |
| RSALOR, StructureDCA | yes | yes (one structure) | T5: minutes, < 1 GB per 1,000-residue region; filtered full-length runs being measured |
| EVE | yes | yes, costly | memory ~7 GB at 20k × 2,300 (37 GB at 100k) [estimate]; GPU time per protein being measured on one > 1,500 aa protein |
| EVmutation | yes | yes, slow | ~50–100 GB, hours per protein at 2,300 aa [estimate from T5]; being measured |
| GEMME | probably | probably | untested: image cannot be built (download server) |
| **EVmutation2** | yes | **no** | CUDA OOM already at 1,000 residues on a GPU slice (memory ∝ L²); would need a code change |

So "all models on ≤ 2,300 aa" holds for every model except EVmutation2, which is limited to ≈ 1,000 aa. Options
consistent with the principle: report EVmutation2 with NaN above its limit (coverage 78 / 97 pairs) and show it on the
≤ 1,000 aa subset; no code changes to it.

Weak alignments inside the set (flag, keep scores): LPA (2,040 aa, kringle repeats, 105 sequences after the
filter) and EXOC3L4 (722 aa, 473 sequences).

## 5. Implementation (simple, no data change)

- The set is defined by `seq_length <= 2300` from the index (`genebass_fig6b.csv`), computed at evaluation time:
  no new archives, no flags written to data. Added to the prototype evaluation
  (`proteingym_clinical/genebass/scripts/16_eval_genebass_source.py`, set `le2300`, job 31937497).
- For the fork: a `max_length: 2300` option of the Genebass evaluation (or a named subset next to `in_vesm_paper`
  on the data branch, if the benchmark should ship it as a fixed list). To decide when the clinical pipeline is
  wired.
- Reporting rule for every set: each model with its pair and row coverage, plus the common subset of pairs scored
  by all compared models. A model's mean over its own pairs alone is not a ranking.

## 6. Expanding beyond 2,300 aa (to evaluate, not innovations to the models)

| Strategy | Brings back | Condition to accept |
|---|---|---|
| Alignment compaction (one-hot built per batch, int8 storage) for EVE; same for EVmutation's weighting | EVE / EVmutation on SPTA1, STAB1 (≤ 2,700) and APOB, LRP2, PKD1 | identical scores to the unchanged code on a protein where both run; otherwise those proteins stay NaN as in popEVE |
| Cut-off at 2,700 aa (AFDB single-structure limit) instead of 2,300 | +2 proteins, +3 pairs (SPTA1, STAB1) | EVE / EVmutation must run at 2,570 aa (same compute class as 2,300) |
| AFDB fragment placement for RSALOR / StructureDCA > 2,700 aa | APOB, LRP2, PKD1 for these two | would be our addition (they need one structure): only if minimal and documented as a deviation; default NaN |
| EVmutation2 > 1,000 aa | 19 pairs | needs upstream code change: not planned |

Note: with the EVE filter, the long multidomain proteins keep few sequences (APOB 2,533, N/L 0.55; PKD1 4,047, N/L
0.94), so even when compaction makes them run, their alignment-model scores carry a weak-alignment flag.
