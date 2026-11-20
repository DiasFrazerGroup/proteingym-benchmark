# Shared by all alignment models: the EVE-protocol MSA filter (`*_filtered` folders)

Applies to `eve_filtered`, `evmutation_filtered`, `evmutation2_filtered`, `gemme_filtered`, `rsalor_filtered` and
`structuredca_filtered`. Full rules and sources: `models_long_proteins.md` §3.5–3.6.

## Why

The archives carry the ColabFold MSA, up to 1.16M rows. In the zoo, no alignment model handles that on long proteins
(memory, O(N²) weighting), and some score unscored positions as the wild type. EVE and popEVE ran the human proteome
at full length with a filtered alignment. User rule: start from the EVE/popEVE implementation, give every alignment
model the same alignment, and do not improve the models.

## Changes (identical in every `*_filtered` folder)

- **`msa_filter.py`** (one byte-identical copy per folder, git blob `d9bd2d4b`, also vendored on the data branch).
  The filter has five steps:
  1. **Match states.** Keep match columns only (insertions removed).
  2. **Fragments.** Drop rows that are more than 50% gaps.
  3. **Non-standard residues.** Drop rows containing any.
  4. **Bit-score sweep.** Sweep thresholds in MMseqs2 bits per residue, read from the a3m headers.
  5. **Lowest threshold that qualifies.** The first rule that some threshold satisfies decides:
     - strict: L_cov ≥ 0.8 L and 10 L ≤ N ≤ 100k (EVE);
     - relaxed: L_cov ≥ 0.7 L, N ≤ 200k (EVE);
     - low_depth (ours, the paper's "best available").

  All columns are kept, as in popEVE.
- **`filtered_model.py`:** `class XFiltered(MSAFilterMixin, X)`. Filter parameters are popped from the card; the
  zoo class gets everything else unchanged.
- **`__init__.py`:** a header hook that attaches the a3m bit scores to the sequences. evedesign drops headers, so the
  hook is needed; the model sees the same sequences as before.
- **Unscored variants become NaN** instead of the wrappers' wild type/0. This covers EVmutation gap-excluded columns,
  RSALOR positions without RSA, GEMME NA cells and EVE non-focus columns.
- Each model keeps its own weighting and column handling (zoo cards unchanged).
- `tests/test_msa_filter.py` tests the filter on synthetic alignments.
- `msa_filter: false` gives the zoo model exactly. That is the T1 test.

## Data side (user decision 2026-11-19)

The final archives (`$WORK/genebass/final_le2700/`) store only the filtered rows. That is 868 MB instead of 7.1 GB,
built with the same `msa_filter.py`. On them the model-side filter is a no-op: it is idempotent, with a max
difference of 4e-16 on RSALOR. Data branch `clinical/genebass-data` 22c1d1b, `docs/clinical/genebass.md`
"Final dataset".

## Deviations from EVE/popEVE (forced by the data, none is an improvement)

- ColabFold/MMseqs2 instead of jackhmmer.
- Bit scores on the MMseqs2 scale.
- Thresholds are cuts of one search instead of re-searches.
- The low_depth step: the paper says "best available" without defining it.

15 of 63 proteins end below N = 10 L. They keep their scores and carry the weak-alignment flag `msa_quality`.

## Runner additions (validation scripts, not model code)

- `scripts/clinical/models/msa_filter/validate_msa_filter.{py,sbatch}`: one build plus score per job in the zoo
  image, with the package mounted. It writes time, peak RSS (self and children) and GPU peak, and keeps a record of
  failures.
- `run_genebass.sbatch`: scores the union of a protein's Genebass mutants over its archives; card overrides via
  `HP`.
- `export_genebass.py`: writes the `v2__*` exports and `status.tsv`.
