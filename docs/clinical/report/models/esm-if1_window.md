# esm-if1_window (structure-conditioned; zoo `esm-if1`)

Status 2026-11-20: final, on all 63 proteins of the capped set (≤ 2,700 aa). Validated (T1).
Predictions: `predictions/v2__esm-if1_window/`.

## Zoo behaviour

The zoo model scores the mean autoregressive log-likelihood of the full mutant sequence given one backbone. It does
one unbatched pass per mutant and accepts one structure numbered like the wild type.

ESM-IF1 was trained on chains of at most 500 residues (Hsu et al. 2022, appendix A.1). Its positions are sinusoidal,
so nothing stops it from running on longer chains. Beyond residue 500, however, its predictions collapse:

| Region | Full-length WT recovery | 500-residue tiles | 500 + context |
|---|---|---|---|
| 1–500 | 0.485 | 0.475 | 0.486 |
| 501–1000 | 0.123 | 0.504 | 0.519 |
| > 1000 | 0.086–0.094 | 0.45–0.52 | 0.44–0.50 |

The table shows wild-type recovery (argmax = wild type), teacher-forced, on 12 proteins of 656–2,570 aa. Mean
log p(wild type) drops from −1.6 to −3.1 past residue 500 at full length.

## Changes

`models/esm-if1_window/`, class `ESMIF1Window` (`window_model.py`); the zoo `model.py` is copied verbatim.

1. **Structure placement** (`fragments.py`, shared with `saprot_window`). Every structure is placed on the wild type
   by sequence. Each mutant goes to the structure where its position is most central. In the ≤ 2,700 aa set every
   archive has one full-length AFDB structure, so this step never changes anything there.
2. **Tiles.** Tiles are `window` residues long and start every `window/2`. Each mutant is scored in the tile where it
   is most central, as Σ_tile log p(mutant | backbone) − Σ_tile log p(wt | backbone). Scoring is batched
   (`batch_size`).
3. **Context chains (`context: true`): the dark-energy method.** This follows the dark-energy project script
   `dark_energy_mane/scripts/awsem/esm/esm_if_logit_longprot.py`. The target tile is decoded and the rest of the
   protein is given to the encoder as extra chains. We use ESM's own multichain layout
   (`multichain_util._concatenate_coords`): target first, then each context chain after 10 NaN padding residues.
   Only the target is decoded, as in ESM's `score_sequence_in_complex`.
4. **Bugs and suboptimal handling in the dark-energy script, fixed here:**
   - Context chains were contiguous slices of list indices. They repeated target residues and skipped others. Here
     they are the rest of the frame, split into pieces of at most `window` residues, with no overlap with the target
     and no gaps.
   - Only the first 5 context chains were kept, in N→C order, so distant C-terminal parts were dropped. Here all of
     the protein is used.
   - The decoder ran over every context chain. Here only the target is decoded.
   - Log-probabilities were merged across overlapping chunks, scaled by an asymmetric overlap weight on first
     coverage. Here there is no merging: each mutant is scored once, in its most central tile.
5. **Card:** `window: 500`, `context: true`, `batch_size: 1`. Batch size 1 fits an 11 GB GPU at 2,570 aa; batch 4
   runs out of memory above about 1,300 aa.
6. Residues that differ between the structure and the wild type are allowed (the zoo raises). The wild type is the
   reference.

## Validation

- **T1** (L ≤ 500, one tile, no context): equal to L·(zoo(mut) − zoo(wt)), max diff 1.8e-4 on 3 proteins (job
  32042243). This is a per-protein linear transform of the zoo score, so per-protein ranks are unchanged.
- **Genebass run:** 63 proteins, 4,653 unique mutants, 0 unscored (job 32042243, RTX 2080 Ti, 42 min). The same
  variant in several archives differs by at most 1e-3 (batch numerics).
- **Shareable evidence (not part of this report, outside the repo):** the plots in
  `$WORK/models/esmif1_dms/plots/` show WT recovery against position, ProteinGym DMS Spearman per assay, and
  agreement with ProteinMPNN/PottsMPNN. On DMS of long proteins (positions > 500), windows give Spearman 0.30–0.62
  against 0.05–0.06 at full length.

| Genebass (100 pairs, same 6,066 rows) | mean −log10 P |
|---|---|
| full length (`window: 2700`, superseded) | 1.27 |
| 1,022 tiles, no context (first port, superseded) | 1.25 |
| **500 + context (final)** | **1.66** |
| ProteinMPNN / PottsMPNN (no windows needed) | 1.27 / 1.27 |

Agreement with ProteinMPNN on Genebass proteins > 500 aa is the median per-protein Spearman over 46 proteins. It is
0.65 at full length and 0.80 with 500 + context. Superseded predictions are kept as `v2__esm-if1_window__full2700`
and `__tiles1022`.

## Above 2,700 aa

AFDB has only 1,400-residue fragments here. The mutant's most central fragment is used, with 500-residue tiles and
context from the same fragment. The code handles this and ran earlier on APOB/LRP2/PKD1 (tiles of 1,022). Not in the
capped benchmark. Cost grows linearly in tiles, so it is feasible at any length.
