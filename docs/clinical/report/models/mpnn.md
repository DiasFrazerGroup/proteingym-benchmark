# proteinmpnn, proteinmpnn_window, pottsmpnn_window (inverse folding)

Status 2026-11-20: final, all proteins. Validated (T1). Predictions: `predictions/v2__proteinmpnn_window/`,
`v2__pottsmpnn_window/`. Details: `models_long_proteins.md` §3.7.

## ProteinMPNN: new zoo folder `models/proteinmpnn`

ProteinGym has a ProteinMPNN baseline, but the zoo had no folder for it.
- The model is `evedesign.models.mpnn:LigandMPNN` with ProteinMPNN weights `v_48_020`, the upstream default: 48
  neighbours, 0.20 Å noise. The weights come from dauparas/ProteinMPNN 8907e66 and are sha256-pinned.
- Score: mean log-likelihood of the sequence given the backbone, in one random decoding order. This is the quantity of
  ProteinGym's `compute_fitness.py`.
- **Only change:** torch is seeded before build (`seed: 0`). The wrapper draws the decoding order with an unseeded
  `torch.randn`, so unseeded scores were not reproducible.

## Window variants (only used above 2,700 aa)

- **`proteinmpnn_window`:** with one structure, the zoo model runs unchanged. With fragments, each fragment is
  renumbered onto the wild type and gets its own model. A mutant is scored in its most central fragment as
  n_fragment × (mean_ll(mut) − mean_ll(wt)).
- **`pottsmpnn_window`:** with one structure, the zoo runs unchanged. With fragments, the zoo's Potts energy table is
  computed per fragment, and a mutant is scored as E(wt) − E(mut) with its most central fragment's table. The
  reference sequence is the entity's wild type.

## Validation and runs

- T1 on APOA1 (window class vs zoo class, two processes): max diff 0 for both.
- Genebass: 66/66 proteins, 5,840 mutants, 0 unscored.
  - ProteinMPNN: 142 CPU-min, at most 12 min and 10 GB per protein.
  - PottsMPNN: 29 min in total.
- PottsMPNN vs the prototype `pottsmpnn_clinical`: Pearson 1.0.
- ProteinMPNN vs PottsMPNN per protein: Spearman median 0.85.
- Capped set (100 pairs): ProteinMPNN 1.268, PottsMPNN 1.266.

## Above 2,700 aa

These models need no windows. Both are k-nearest-neighbour graphs, so cost is linear in L; PottsMPNN has run up to
~21k aa on 24 GB. Above 2,700 aa they use AFDB fragments (most central). On all 109 pairs: ProteinMPNN 1.320,
PottsMPNN 1.356.
