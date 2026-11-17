---
# Model identifier used for referencing this model in the benchmark system
name: eve_filtered

# Import string for the evedesign model class to run. Accepts entry-point style
# ("package.module:ClassName") or fully dotted style ("package.module.ClassName").
# The class is instantiated with hyper_parameters, then built and scored.
model_class: "proteingym.models.eve_filtered.filtered_model:EVEFiltered"

tags: ["zero-shot", "embeddable"]

# Multi-target support flag
multi_y: false

# Hyperparameters are forwarded verbatim to the evedesign EVE wrapper
# constructor (evedesign.models.eve.EVE), so every key here must be a valid
# EVE.__init__ argument.
hyper_parameters:
    # --- clinical: EVE-protocol MSA filter (msa_filter.py, docs/clinical/models_long_proteins.md §3.6) ---
    # false = the zoo model unchanged (except NaN for unscored variants)
    msa_filter: true
    # EVE fragment filter: keep rows with gap fraction <= 1 - min_seq_coverage (rows with X/B/Z/J... also
    # dropped, EVE's indeterminate-residue rule); query always kept
    min_seq_coverage: 0.5
    # a column counts towards L_cov if >= this fraction of the kept rows have a residue there (L_cov only selects
    # the threshold; all columns are kept, popEVE-like)
    col_occupancy: 0.7
    # bit-score threshold: lowest one with L_cov >= lcov_target * L and min_depth_per_length * L <= N <= max_depth,
    # relaxed as in the EVE paper to lcov_relaxed and max_depth_relaxed, then "best available" (msa_filter.py)
    lcov_target: 0.8
    lcov_relaxed: 0.7
    min_depth_per_length: 10.0
    max_depth: 100000
    max_depth_relaxed: 200000
    # random subsample to max_depth_relaxed only if no threshold gives N <= max_depth_relaxed
    seed: 42
    # bits/residue = MMseqs2 bit score (a3m header) / L
    bits_thresholds: [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5]
    # bit scores come from the dataset's a3m headers; proteingym-base archives rewrite headers to seq_<i>, so for
    # archives the selection is skipped ("no_bitscores") unless this points at the original ColabFold a3m files
    # (matched by query sequence, verified row by row; validation only, see §3.6)
    a3m_headers_dir: null
    # per-protein filter stats appended to <stats_dir>/msa_filter_stats.jsonl (next to the predictions)
    stats_dir: "/opt/program/output"
    # --- zoo model hyper-parameters (unchanged) ---
    # VAE encoder/decoder/training hyperparameters (path to a JSON file or a
    # dict). null uses EVE's default_model_params.json (400,000 training steps)
    model_hyperparameters: null
    # Path to an EVE VAE checkpoint. If it exists it is loaded instead of
    # training; otherwise the trained VAE is saved there. null keeps the
    # trained VAE in memory only (use the --model-file option to persist it)
    model_checkpoint_path: null
    # MSA sequence re-weighting identity threshold (EVE default 0.2)
    theta: 0.2
    # Use weights attached to the dataset MSA where available, otherwise
    # compute EVE weights
    use_weights: true
    # Location to load/save EVE-computed MSA weights (.npy). null keeps them
    # in a temporary working directory
    msa_weights_path: null
    # EVE MSA pre-processing (drop fragment sequences and low-occupancy
    # columns)
    preprocess_msa: true
    # Number of samples used to approximate the ELBO when scoring (EVE uses
    # 20000)
    num_samples: 20000
    # Maximum number of sequences scored concurrently
    batch_size: 2048
    # Random seed for the VAE
    random_seed: 42
    # Device to run the model on ("cpu", "cuda", or "mps"). Training EVE for
    # the default 400,000 steps is only practical on a GPU
    device: "cuda"
    # Maximum number of non-focus MSA sequences to train on; deeper MSAs are
    # randomly subsampled (re-weighting is O(depth^2)). null disables
    max_msa_depth: 200000
    # Fraction of gaps above which a wild-type-covered column is excluded from
    # training. 1.0 keeps every covered column; 0.3 is the original EVE value
    threshold_focus_cols_frac_gaps: 1.0
---

# Model Card for eve_filtered (clinical: EVE-protocol alignment)

Same model as the zoo's `eve` (wrapper unchanged), trained on the dataset's MSA after the EVE/popEVE alignment
filter: insertions removed, fragments (> 50% gaps) and rows with indeterminate residues dropped (EVE code), a
bit-score threshold (MMseqs2 bits/residue from the ColabFold a3m headers) chosen for column coverage L_cov >= 0.8 L
and 10 L <= N <= 100,000, relaxed as in the EVE paper to L_cov >= 0.7 L and N <= 200,000, then "best available"
(documented order); all columns are kept, as in the popEVE release, which scores every position. Every `*_filtered` model uses the same filter and parameters, so all alignment models see the same alignment.
Variants at positions the model does not score are NaN (the zoo wrapper scores them as wild type). Per-protein filter
statistics (rows in / after coverage, threshold, N, N/L, L_cov, relaxation step, quality flag) are logged and
appended to `msa_filter_stats.jsonl` in the output folder. See `docs/clinical/models_long_proteins.md` §3.6.

---


# Model Card for EVE (evedesign-wrapped)

This model scores protein sequences with the [evedesign](https://github.com/evedesignbio/evedesign)
`EVE` wrapper. EVE (Evolutionary model of Variant Effect) is a family-specific
Bayesian variational autoencoder trained on a multiple sequence alignment of
homologs of the target protein. Scoring is zero-shot: the VAE is trained on
the dataset's MSA, and each sequence is assigned its evidence lower bound
(ELBO), averaged over `num_samples` samples, as an approximation to its
log-likelihood under the model.

The wrapper uses EVE code vendored from
[OATML-Markslab/EVE](https://github.com/OATML-Markslab/EVE) (MIT license).
It also exposes `transform()` (latent embeddings) and `generate()` (sampling
from the VAE prior), so this card is tagged `embeddable`.

EVE is a fixed-length model: variants with insertions, deletions or
non-standard residues cannot be scored and are left without a prediction.

For detailed information on the method, please refer to the
[EVE paper](https://doi.org/10.1038/s41586-021-04043-8) (Frazer et al.,
*Nature* 2021).

The benchmark entrypoint converts the ProteinGym subset into evedesign
datatypes with `evedesign.proteingym.dataset_to_evedesign` (which uses the
first MSA of the dataset), builds the wrapper on the resulting `System` (MSA
processing, sequence re-weighting and VAE training), and scores all variants.

## Network access

The Docker/Singularity build fetches Python packages (PyPI, and
`https://github.com/ProteinGym/proteingym-base.git` via GET) and the `uv`
binary. Nothing is fetched at `train` time.
