# saprot_window (structure-aware PLM; zoo `saprot`)

Status 2026-11-20: final, all proteins. Validated (T1, T3). Predictions: `predictions/v2__saprot_window/`.

## Zoo behaviour

The zoo runs beyond SaProt's 1,024-token training length (no cap), which is out of distribution. It needs one
structure numbered like the wild type.

## Changes

`models/saprot_window/`, class `SaProtWindow` (`window_model.py`, with `fragments.py` shared with `esm-if1_window`):

1. The PLM per-position window (1,022 residues, `termini` tokens, masked marginal) is applied to the paired
   (residue, 3Di) tokens.
2. The 3Di tokens come from the structure used for that position. Structures are placed on the wild type by
   sequence, and the most central one is used. In the ≤ 2,700 aa set this is always the single full-length AFDB
   model.
3. `allow_missing_structure: true`: positions that no structure covers get the masked `#` 3Di token, so they are
   scored from sequence only. No such position occurs in Genebass.

## Validation

- T1 (L ≤ 1,022 vs zoo): max diff 1.1e-5, **PASS**.
- vs the prototype (every_window, fragment datasets): r 1.00000 (L ≤ 1,022), 0.99920 (L > 1,022).
- Genebass: 66 proteins, 0 unscored. APOB/LRP2/PKD1 use all 16–18 fragments.
- Capped set (100 pairs): 1.594.

pLDDT is not used (proteingym-base zeroes B-factors in archive structures; reported as an upstream issue).

## Above 2,700 aa

This already works with AFDB fragments (most central fragment per position). On all 109 pairs it scores 1.640.
