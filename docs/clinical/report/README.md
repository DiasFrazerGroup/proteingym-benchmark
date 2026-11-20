# Clinical Genebass benchmark: zoo models on long proteins, final report

Status 2026-11-20. Branches:
- model side: `clinical/model-api` (this branch);
- data side: `clinical/genebass-data` (22c1d1b).

The detailed record, tests and experiment log are in [models_long_proteins.md](../models_long_proteins.md). This report
is the summary. There is one page per model in [models/](models/).

## Principles (user)

- Benchmark existing models. Innovations are minimal, and the starting point is the published implementation
  (EVE/popEVE for alignments, the zoo for everything).
- Be fair, without losing benchmark pairs to model limitations. A model that cannot score a variant gets NaN, and its
  coverage is reported. No imputation.
- Data changes and model changes live on separate branches.

## Decisions

| Date | Decision |
|---|---|
| 2026-11-16 | Structures stay AFDB v6 (fragments above 2,700 aa). PLM windows use `termini` special tokens. No derived (region) datasets. |
| 2026-11-17 | Alignment models get the EVE/popEVE protocol: ColabFold MSA filtered model-side, the same alignment for every model, all columns kept, no MSA segments. Archives keep the a3m headers (bit scores). No imputation. |
| 2026-11-19 | **Benchmark capped at 2,700 aa** (63 proteins, 100 pairs). Archives carry only the filtered alignment rows. |
| 2026-11-19 | EVE: keep the zoo's single seed (the number of seeds belongs to the model card). Do not run EVE on Genebass. |
| 2026-11-19 | GEMME: keep both outputs (evolEpi = zoo, evolCombi = ProteinGym/paper). |
| 2026-11-19 | ESM-IF1: the dark-energy method (500-residue tiles + context chains), bugs fixed. |
| 2026-11-19 | ProteinMPNN and PottsMPNN added (inverse folding beyond ESM-IF1). |

## Why 2,700 aa

- **Structures:** every protein up to 2,700 aa has exactly one full-length AFDB structure (checked on all 100
  archives), so no structure model needs fragments.
- **Alignment models:** every alignment model that was run fits up to 2,700 aa (EVmutation2 only with the chunking
  below and an 80 GB GPU for proteins over ~1,800 aa).
- **Precedent:** popEVE itself stops at 2,302 aa, so SPTA1 and STAB1 have no popEVE.
- **What is lost:**
  - APOB, LRP2 and PKD1: 9 of 109 pairs and 1,725 of 7,791 rows;
  - the renal group loses 4 of its 15 pairs;
  - the model ranking is stable: Spearman of model means, all pairs vs capped, was 0.963 at the 2,300 cut
    ([benchmark_le2300.md](../benchmark_le2300.md)).

## Data and benchmark modifications (data branch `clinical/genebass-data`, 22c1d1b)

Full description: `docs/clinical/genebass.md`, section "Final dataset (<= 2,700 aa, filtered alignments)", on that
branch.

| | |
|---|---|
| Index | `datasets/clinical/genebass/genebass_le2700.csv`: 100 pairs, 63 proteins, 6,066 rows; paper set 94 pairs |
| Alignment filter per protein | `datasets/clinical/genebass/genebass_le2700_msa_filter.tsv` (threshold, N, N/L, L_cov, rule, quality) |
| Archives | `$WORK/genebass/final_le2700/` (100 archives, 868 MB instead of 7.1 GB): MSA = only the rows the models use (EVE-protocol filter, exact original a3m records with headers, original order); everything else byte-identical to the full archives |
| Build | `scripts/clinical/msa/filter_msas.py` (+ `.sbatch`; uses `msa_filter.py` vendored byte for byte from this branch, refuses to run if it differs), `check_final.py` (`CHECK FINAL OK`), `rebuild.sbatch` (whole chain) |
| Weak alignments (N < 10 L, flag `msa_quality = low_depth`, scores kept) | LPA, PCSK9, GPLD1, CGNL1, PAPPA2, STAB1, B4GALNT3, MROH2A, MAP3K15, ANGPTL3, ALB, NPC1L1, RBM47, IGF1R, ABCA1 |
| Layers | unchanged: AlphaMissense and CADD cover all 63 proteins; popEVE lacks SPTA1 and STAB1 (97 of 100 pairs) |

The earlier archive steps are documented on the data branch: Genebass source tables, layers, `in_vesm_paper` flags,
header restore.

## Models

| Model (folder) | Type | What we changed | ≤ 2,700 aa | > 2,700 aa (future) | Page |
|---|---|---|---|---|---|
| `esm2_window`, `esmc_window` | PLM | masked marginal in a 1,022-residue window per position, `termini` tokens | 63/63 | works (windows): already scored | [plm_windows](models/plm_windows.md) |
| `saprot_window` | PLM + 3Di | same windows on (residue, 3Di) tokens, structure placement | 63/63 | works (AFDB fragments): already scored | [saprot_window](models/saprot_window.md) |
| `esm-if1_window` | inverse folding | 500-residue tiles + context chains (dark-energy method, fixed) | 63/63 | works (fragments + tiles) | [esm-if1_window](models/esm-if1_window.md) |
| `proteinmpnn` (new), `proteinmpnn_window`, `pottsmpnn_window` | inverse folding | ProteinMPNN zoo folder (seeded order); most-central fragment above 2,700 aa | 63/63, zoo unchanged | works (fragments): already scored | [mpnn](models/mpnn.md) |
| all `*_filtered` | alignment | EVE-protocol MSA filter, NaN for unscored variants | — | — | [msa_filter_shared](models/msa_filter_shared.md) |
| `evmutation_filtered` | alignment (plmc) | filter + NaN | 63/63 (264 of 4,653 mutants NaN: gap columns) | **fails**: plmc 32-bit limit at ~3,100 sites | [evmutation_filtered](models/evmutation_filtered.md) |
| `gemme_filtered` | alignment | filter + NaN, `output` option, image from pinned local tarballs | 63/63 | expected to run (untested) | [gemme_filtered](models/gemme_filtered.md) |
| `rsalor_filtered`, `structuredca_filtered` | alignment + structure | filter + NaN | 63/63 | **fail**: need one full-length structure | [rsalor_structuredca](models/rsalor_structuredca_filtered.md) |
| `eve_filtered` | alignment (VAE) | filter + NaN | validated, **not run** (decision) | needs compaction (float64 one-hot) | [eve_filtered](models/eve_filtered.md) |
| `evmutation2_filtered` | alignment (MSA transformer) | filter + NaN, exact chunked encoder | 63/63 (0 NaN); > 1,800 aa needs the 80 GB GPU (STAB1 67 GB) | **fails**: GPU memory ~L², about 190–220 GB for APOB/PKD1/LRP2 | [evmutation2_filtered](models/evmutation2_filtered.md) |
| `blosum62` | substitution matrix | none (zoo works at any length) | not run in this round | works | — |

## Results (Genebass, mean −log10 P of r(β, score), VESM Fig. 6b metric)

Final evaluation: `results/genebass_common_variants_final.tsv` (prototype repo, job 32052952). Columns:
- **common** = the 5,750 rows of the 100 pairs that every model scores (EVmutation's 264 gap-column mutants and
  a few reference NaNs removed; popEVE has 97 pairs);
- **own** = each model on all the rows it scores within the cap;
- **all 109** = including the 9 pairs above 2,700 aa, only for models that score them (not comparable).

| Model | Type | common (100 pairs) | own (≤ 2,700 aa) | rows own | all 109 pairs |
|---|---|---|---|---|---|
| AlphaMissense (layer) | reference | 2.441 | 2.643 | 6,066 | 2.695 |
| `structuredca_filtered` | alignment + structure | **2.281** | **2.471** | 6,066 | — |
| `evmutation_filtered` | alignment | 2.208 | 2.208 | 5,750 | — |
| popEVE (layer) | reference | 2.241 (97 pairs) | 2.372 | 5,734 | — |
| `gemme_filtered` evolCombi | alignment | 2.178 | 2.352 | 6,066 | — |
| `gemme_filtered` evolEpi | alignment | 2.155 | 2.341 | 6,066 | — |
| `esm2_window` | PLM | 2.132 | 2.295 | 6,066 | 2.403 |
| `esmc_window` | PLM | 2.037 | 2.166 | 6,066 | 2.181 |
| `rsalor_filtered` | alignment + structure | 1.998 | 2.151 | 6,066 | — |
| `evmutation2_filtered` | alignment | 1.937 | 2.092 | 6,066 | — |
| Allele frequency | reference | 1.801 | 1.710 | 6,063 | 1.629 |
| CADD (layer) | reference | 1.778 | 1.944 | 6,063 | 1.912 |
| `esm-if1_window` (dark energy) | inverse folding | 1.551 | 1.663 | 6,066 | — |
| `saprot_window` | PLM + 3Di | 1.489 | 1.594 | 6,066 | 1.640 |
| `pottsmpnn_window` | inverse folding | 1.186 | 1.266 | 6,066 | 1.356 |
| `proteinmpnn_window` | inverse folding | 1.172 | 1.268 | 6,066 | 1.320 |

- Rows are in the TSV's order (by the 2,300 aa column).
- **Common vs own rows.** Two models move:
  - EVmutation drops from 2nd to 5th among the models. Its own rows are the common rows, because it is the model
    whose 264 gap-column NaN define them.
  - CADD moves above allele frequency.
  - The rest keep their order.
- **ESM-IF1 versions** (own rows): the dark-energy method gives 1.663. Full length gives 1.270, and the old
  1,022-residue tiles give 1.250 (`results/genebass_esmif1_versions.tsv`). Details:
  [esm-if1_window](models/esm-if1_window.md).
- **At 2,300 aa** (97 pairs, column `le2300` in the TSV): every model is 0.02–0.07 higher on common rows.
  - The exception is popEVE: it has the same 97 pairs in both, so its value does not change.
  - As a result popEVE falls behind EVmutation at 2,300 aa. The order of everything else is the same.
- `esm-if1_window` has no "all 109" value: the dark-energy run was made for the capped set only.

## Reproducibility

| Step | Script |
|---|---|
| Model images | `scripts/clinical/models/build_images.sbatch <folder>...` |
| PLM windows | `scripts/clinical/models/validate_plm_windows.{py,sbatch}` |
| Structure windows (SaProt, ESM-IF1) | `scripts/clinical/models/validate_struct_windows.{py,sbatch}` |
| Alignment models (one protein; all of Genebass) | `scripts/clinical/models/msa_filter/validate_msa_filter.{py,sbatch}`; `run_genebass.sbatch` |
| Exports | `scripts/clinical/models/export_scores.py` (windows) and `msa_filter/export_genebass.py` (alignment models, `status.tsv`) |
| Evaluation | `proteingym_clinical/genebass/scripts/17_eval_common_variants.py` (prototype repo) |

Job IDs: `$WORK/logs/jobs.txt`. Run outputs: `$WORK/models/`.

## Open items

1. **proteingym-base, upstream:**
   - `MSA.dump` drops a3m headers (we restored them on the data side);
   - structures are written with B-factors zeroed, so pLDDT is lost in every archive.
   Both should be reported upstream.
2. **evedesign `dataset_to_evedesign`:** does not pass headers to models (we use a hook), hard-codes
   `first_index = 1`, asserts one chain.
3. **NaN policy in the PG2 metrics:** proteingym-base `metric_spearman` returns NaN for a dataset with any NaN
   prediction. The clinical evaluation drops NaN rows per model and reports coverage.
4. **Beyond 2,700 aa:** see the "> 2,700 aa" column above and each model page.
   - Needs code changes: EVmutation (plmc 64-bit), RSALOR/StructureDCA (structure fragments, against the rule),
     EVE (compaction).
   - EVmutation2: GPU memory, not code. It needs about 2.5–3× the 80 GB GPU for APOB/PKD1/LRP2; more chunking in the
     vendored AF3 encoder or a larger GPU.
5. **EVE on Genebass:** validated and feasible to 2,570 aa (5.8 h training on 7g.80gb). Not run by decision; popEVE
   is the reference.
