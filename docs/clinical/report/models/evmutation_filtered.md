# evmutation_filtered (zoo `evmutation`: `evedesign.models.evcouplings:EVcouplingsPLM`, plmc)

Status 2026-11-20: final, 63/63 proteins of the capped set. Predictions: `predictions/v2__evmutation_filtered/`.

## Changes

- Shared MSA filter ([msa_filter_shared.md](msa_filter_shared.md)).
- **NaN for mutants in columns the model excludes** (more than 50% gaps, unweighted). The zoo scores them as the wild
  type (0). That leaves 264 of 4,653 Genebass mutants unscored, and they define the benchmark's "variants scored by every model" (5,750 of 6,066 rows).
- The card is the zoo's: theta 0.8, λ_h = λ_J = 0.01, 100 iterations.

## Validation

- **T1** (APOA1): max diff 0 on the 5,016 scored singles. The 57 NaN are gap-excluded columns that the zoo scored as
  the wild type.
- **Reproduces ProteinGym's released EVmutation** on its own alignment: GRB2 Spearman 0.996 (focus columns), P53
  0.972 (Pearson 0.987; our input keeps 54 columns EVcouplings drops).
- **Genebass:** 63 proteins, CPU (genoa64), 1,939 CPU-min in total. The largest are SPTA1 (215 min, 74 GB),
  ABCA1 (100 min, 69 GB) and STAB1 (32 min, 74 GB).
- Capped set: 2.208 (its scored rows are the variants scored by every model).

## Above 2,700 aa

- **Fails:** plmc segfaults on loading more than about 3,100 modelled sites (APOB 4,314, PKD1 4,156, LRP2 4,152).
  Its 32-bit parameter count overflows: L² q² / 2 > 2³¹. This is a model limitation.
- **Fix would be upstream:** a 64-bit index in plmc. Memory would then be about 16 GB per coupling vector, over 100 GB
  for L-BFGS [estimate].
- **Alternatives:** NaN (current), or region models (T5: W = 500 reproduces full length at r 0.93), which the user
  rejected (no MSA segments).
