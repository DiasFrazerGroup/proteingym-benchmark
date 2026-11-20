# eve_filtered (zoo `eve`, evedesign 0.0.6)

Status 2026-11-20: implementation validated. **Not run on Genebass** (user decision 2026-11-19: GPU cost, about
1.5–18 h per protein). Not in the benchmark tables; popEVE (layer) is the EVE-family reference.

## Changes

- Only the shared MSA filter ([msa_filter_shared.md](msa_filter_shared.md)) and NaN for non-focus columns.
- The card is the zoo's: theta 0.2, 20,000 samples, `threshold_focus_cols_frac_gaps: 1.0` (all columns, as the
  popEVE release), one seed (`random_seed: 42`).
- Seeds: ProteinGym averages seeds 0–4 (`EVE_ensemble`). The user decided on 2026-11-19 that the number of seeds
  belongs to the model card, not the benchmark ("not our call"), so the single seed is kept.

## Validation

- **T1** (APOA1, filter off vs zoo): max diff 0.
- **Reproduces ProteinGym's released EVE on ProteinGym's own alignments** (v1.3 `DMS_msa_files`, zoo card unchanged):

  | Protein | Ours vs EVE_ensemble | Noise bar: released single vs ensemble |
  |---|---|---|
  | GRB2 (deep MSA) | Spearman 0.997 | 0.995 |
  | P53 (shallow MSA) | 0.951 (seed 0), 0.888 (seed 42) | 0.960 |

  On P53, our two seeds agree with each other at 0.934, so the spread is seed noise, not the code.
- **vs the popEVE release** (our ColabFold-filtered alignment vs their jackhmmer one): Spearman −0.85 (APOA1),
  −0.88 (PCSK9) and −0.90 (ABCA1, 2,261 aa). Release vs release (ProteinGym vs popEVE, 44 proteins) has a median
  of 0.937, so this is within the alignment effect.
- **Feasibility at the cap (STAB1, 2,570 aa, N 10,291, 7g.80gb):**
  - training (400k steps) took 5.8 h at 18.5 it/s with 15 GB host RAM;
  - scoring runs at about 2,000 sequences per 8 min (51 GB host RAM for 48,831 singles);
  - the job was stopped once feasibility was shown;
  - ABCA1 took about 18 h in total (two jobs because of the 12 h GPU limit).

## Above 2,700 aa

- **Memory:** the zoo builds the one-hot training set as float64 N × L × 20, plus the whole scoring set at once
  (ABCA1 scoring: 43 GB host). Above about 3,000 aa this needs compaction: int8 N × L, one-hot built per batch.
  That is a code change, not done.
- **Time:** GPU time grows with L (18.5 it/s at 2,570 aa).
- **No precedent:** popEVE stops at 2,302 aa, so there is nothing to validate against.
- **Weak alignments:** APOB (N/L 0.55) and PKD1 (N/L ~1) would carry the flag.

**Expectation:** feasible with compaction plus two chained jobs per protein, with weak alignments.
