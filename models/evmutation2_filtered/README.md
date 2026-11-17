---
# Model identifier used for referencing this model in the benchmark system
name: evmutation2_filtered

# Import string for the evedesign model class to run. Accepts entry-point style
# ("package.module:ClassName") or fully dotted style ("package.module.ClassName").
# The class is instantiated with hyper_parameters, then built and scored
model_class: "proteingym.models.evmutation2_filtered.filtered_model:EVmutation2Filtered"

tags: ["zero-shot", "embeddable"]

# Multi-target support flag
multi_y: false

# Hyperparameters are forwarded verbatim to the evedesign EVmutation2 wrapper
# constructor (evedesign.models.evmutation2.EVmutation2)
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
    # (fetching https://huggingface.co/thomashopf/evmutation2/resolve/main/
    # msa-only-small.ckpt at build() time)
    model_file_path: "/opt/program/data/msa-only-small.ckpt"
    # Number of encoder samples (single/pair representations) drawn from the
    # MSA encoder and averaged over
    encoder_num_samples: 16
    # Recycling steps run when computing each MSA encoding.
    encoder_num_recycling_steps: 4
    # Number of sequences sampled from the MSA when computing the encoding.
    encoder_max_num_msa: 2048
    # Maximum number of sequences decoded concurrently.
    decoder_batch_size: 64
    # Sampled decoding orders averaged over when computing full-sequence
    # scores
    decoder_num_full_samples: 16
    # Sampled decoding orders averaged over when computing single-mutant/
    # conditional scores
    decoder_num_mutant_samples: 16
    # Device to run the model on ("cpu", "cuda", or "mps")
    device: "cpu"
---

# Model Card for evmutation2_filtered (clinical: EVE-protocol alignment)

Same model as the zoo's `evmutation2` (wrapper unchanged), trained on the dataset's MSA after the EVE/popEVE alignment
filter: insertions removed, fragments (> 50% gaps) and rows with indeterminate residues dropped (EVE code), a
bit-score threshold (MMseqs2 bits/residue from the ColabFold a3m headers) chosen for column coverage L_cov >= 0.8 L
and 10 L <= N <= 100,000, relaxed as in the EVE paper to L_cov >= 0.7 L and N <= 200,000, then "best available"
(documented order); all columns are kept, as in the popEVE release, which scores every position. Every `*_filtered` model uses the same filter and parameters, so all alignment models see the same alignment.
Variants at positions the model does not score are NaN (the zoo wrapper scores them as wild type). Per-protein filter
statistics (rows in / after coverage, threshold, N, N/L, L_cov, relaxation step, quality flag) are logged and
appended to `msa_filter_stats.jsonl` in the output folder. See `docs/clinical/models_long_proteins.md` §3.6.

---


# Model Card for EVmutation2 (evedesign-wrapped)

This model scores protein sequences with the [evedesign](https://github.com/evedesignbio/evedesign)
`EVmutation2` wrapper. EVmutation2 encodes the evolutionary context of a target protein
family from a multiple sequence alignment (an axial-attention-style encoder
producing single/pair representations), then autoregressively decodes 
sequence and mutant log likelihoods from those representations.
Scoring is zero-shot - each sequence is assigned a log-likelihood-derived
score.

For detailed information on the model architecture and training data, please
refer to the [accompanying preprint](https://doi.org/10.64898/2026.03.17.712115).

## Citation

```
@article{Hopf2026evedesign,
  doi = {10.64898/2026.03.17.712115},
  url = {https://doi.org/10.64898/2026.03.17.712115},
  year = {2026},
  publisher = {bioRxiv},
  title = {evedesign: accessible biosequence design with a unified framework},
  author = {Hopf, Thomas A. and Belahsen, Khaoula and others}
}
```

## License

Both the `evmutation2` PyPI package's source
([evedesignbio/evmutation2](https://github.com/evedesignbio/evmutation2)) and
the [`thomashopf/evmutation2`](https://huggingface.co/thomashopf/evmutation2)
checkpoint on HuggingFace Hub are released under the **MIT license** - see
`LICENSE` in this directory.
