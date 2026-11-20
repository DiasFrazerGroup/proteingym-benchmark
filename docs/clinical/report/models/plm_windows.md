# esm2_window, esmc_window (protein language models; zoo `esm2`, `esmc`)

Status 2026-11-20: final, all proteins (no length limit). Validated (T1–T3).
Predictions: `predictions/v2__esm2_window__termini_1022/`, `v2__esmc_window__termini_1022/`.

## Zoo behaviour

- `esm2` fails above 1,022 aa (`max_seq_length`). It scores an unmasked pseudo-log-likelihood with one full pass per
  mutant. The masked marginal (`score_mutants`) exists in evedesign but is unused.
- `esmc` runs at any length, with L² memory, and risks silent truncation.

## Changes

`models/esm2_window/` (class `ESM2Window`) and `models/esmc_window/` (`ESMCWindow`), with shared `window.py`:

1. **Masked marginal**, the ProteinGym standard: Σ over substitutions of log p(mut) − log p(wt), with the position
   masked.
2. **One window per mutated position.** The window has `window` residues centred on the position, shifted to fit at
   the ends. When L ≤ window it is the whole protein, which equals the zoo's `score_mutants`.
3. **Special tokens `termini`** (user decision 2026-11-16): `<cls>`/`<eos>` only at the real protein ends, as in
   dark energy. `every_window` is kept as an option.
4. **Window 1,022 for both models**, for parity with ESM-2's context. ESMC 2,048 and full length were tested.
5. Cards: `window: 1022`, `special_tokens: termini`. Batch size fits the GPU slice (8–16 on 3g.40gb).

## Validation

| Test | Result |
|---|---|
| T1, 48 proteins ≤ 1,022 aa vs zoo `score_mutants` | ESM2 1.6e-5, ESMC 2.7e-5: **PASS** |
| T2 vs dark-energy ESM2 650M masked (independent code) | r 1.00000 (L ≤ 1,022), 0.99997 (L > 1,022) |
| T3 ESMC window 1,022 / 2,048 / full length | Genebass 2.181 / 2.200 / 2.192: no meaningful effect |
| termini vs every_window | ESM2 r 0.99998; ESMC r 0.9985 (Genebass 2.181 vs 2.200) |

Genebass: 66 proteins, 5,840 unique mutants, 0 unscored. On the capped set (100 pairs, 6,066 rows): ESM2 2.295,
ESMC 2.166.

## Above 2,700 aa

No limit. Windows are 1,022 residues and cost one pass per mutated position. They already score APOB, LRP2 and
PKD1; on all 109 pairs ESM2 scores 2.403 and ESMC 2.181.
