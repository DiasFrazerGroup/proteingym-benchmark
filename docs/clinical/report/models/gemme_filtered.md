# gemme_filtered (zoo `gemme`: JET2 + computePred.R)

Status 2026-11-20: final, 63/63 proteins of the capped set, **two outputs kept** (user decision 2026-11-19).
Predictions: `predictions/v2__gemme_filtered__evolEpi/` (zoo default) and `v2__gemme_filtered__evolCombi/`
(ProteinGym/paper).

## Changes

- Shared MSA filter ([msa_filter_shared.md](msa_filter_shared.md)). NaN for NA cells of the prediction matrix (the
  zoo adds 0).
- **`output` option** (`model.py`, commit 8e12abf):
  - `evolEpi`, the default, is what the zoo reads (`_normPred_evolEpi.txt`);
  - `evolCombi` is GEMME's combined model, which ProteinGym reports and the GEMME paper recommends.
  - Both files come out of the same `computePred.R` run; nothing else changes.
- **Image build** (`scripts/clinical/models/build_images.sbatch`, commit 36f5e05):
  - The GEMME/JET2 download server (lcqb.upmc.fr) is down, so the image is built from local copies of the same
    tarballs.
  - `GEMME.tgz` comes from the Wayback Machine (20231001225718). `JET2.tgz` comes from the authors' image
    `elodielaine/gemme@sha256:a35099f7…`.
  - The `Singularity.def`'s own sha256 pins still check both, so the code is byte-identical to the zoo's.

## Validation

- **Reproduces ProteinGym's released GEMME** on ProteinGym's alignments:

  | Protein | evolCombi vs released (Spearman) | evolEpi vs released (Spearman) |
  |---|---|---|
  | GRB2 | 0.998 | 0.94 |
  | P53 | 0.998 | 0.86 |

  ProteinGym reports evolCombi.
- JET2's iJET resampling is unseeded upstream, so repeated runs differ slightly (Spearman about 0.998).
- **Genebass:** 63 proteins, CPU, 250 min in total, at most 10 min and 11 GB per protein (SPTA1).
- evolEpi vs evolCombi per protein: median Spearman 0.95 (0.83–0.98).
- Capped set (100 pairs, own rows): evolCombi 2.352, evolEpi 2.341.

## Above 2,700 aa

- Not run (cap).
- **Time:** 10 min at 2,419 aa, so time is not a problem.
- **Inputs:** JET2 takes the first 20k rows and R the whole filtered alignment. Both are small after the filter
  (APOB 2.5k rows).
- **Expected to run.** Structure is a dummy input, so no fragments are needed. Weak alignments (APOB, PKD1) would
  carry the flag. Untested.
