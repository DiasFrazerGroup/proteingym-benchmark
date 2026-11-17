---
# Model identifier used for referencing this model in the benchmark system
name: structuredca_filtered

# Import string for the evedesign model class to run. Accepts entry-point style
# ("package.module:ClassName") or fully dotted style ("package.module.ClassName").
# The class is instantiated with hyper_parameters, then built and scored.
model_class: "proteingym.models.structuredca_filtered.filtered_model:StructureDCAFiltered"

tags: ["zero-shot"]

# Multi-target support flag
multi_y: false

# Hyperparameters are forwarded verbatim to the structuredca StructureDCA
# wrapper constructor, so every key here must be a valid StructureDCA.__init__
# argument. This is a CPU-only model
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
    # Use StructureDCA's RSA-complement reweighting for all scores
    reweight_by_rsa: false
    # Distance threshold (Angstrom) for a residue-residue contact to be kept in
    # the sparse DCA model
    distance_cutoff: 8.0
    # Remove contacts at positions with low structure confidence
    use_contacts_plddt_filter: false
    contacts_plddt_cutoff: 70.0
    contacts_plddt_keep_window: 1
    # Drop contacts at MSA positions with gap-ratio above this threshold
    contacts_gap_cutoff: null
    # L2 regularisation strength on fields h_i and couplings J_ij
    lambda_h: 1.0
    lambda_J: 1.0
    # Asymptotic correction to the L2 regularisation (as Neff -> +inf)
    lambda_asymptotic: 0.001
    # Frequency-level regularisation, only used to initialise fields h
    theta_regularization: 0.10
    # If true, gap is excluded as a modelled state (cannot be scored then)
    exclude_gaps: true
    # Discard MSA sequences below this identity to the target (null = no filter)
    min_seqid: 0.25
    # Sequence-identity clustering threshold for reweighting (null = no reweighting)
    weights_seqid: 0.90
    count_target_sequence: true
    # Solver settings
    num_threads: 4
    max_iterations: 2000
    use_sparse_J: true
    # Caching paths (null = no caching, always refit)
    distance_cache_path: null
    rsa_cache_path: null
    weights_cache_path: null
    dca_cache_path: null
    # Logging
    verbose: false
    log_gd_steps: false
    disable_warnings: false
    disable_log_colors: false
    disable_solver_logs: false
---

# Model Card for structuredca_filtered (clinical: EVE-protocol alignment)

Same model as the zoo's `structuredca` (wrapper unchanged), trained on the dataset's MSA after the EVE/popEVE alignment
filter: insertions removed, fragments (> 50% gaps) and rows with indeterminate residues dropped (EVE code), a
bit-score threshold (MMseqs2 bits/residue from the ColabFold a3m headers) chosen for column coverage L_cov >= 0.8 L
and 10 L <= N <= 100,000, relaxed as in the EVE paper to L_cov >= 0.7 L and N <= 200,000, then "best available"
(documented order); all columns are kept, as in the popEVE release, which scores every position. Every `*_filtered` model uses the same filter and parameters, so all alignment models see the same alignment.
Variants at positions the model does not score are NaN (the zoo wrapper scores them as wild type). Per-protein filter
statistics (rows in / after coverage, threshold, N, N/L, L_cov, relaxation step, quality flag) are logged and
appended to `msa_filter_stats.jsonl` in the output folder. See `docs/clinical/models_long_proteins.md` §3.6.

---


# Model Card for StructureDCA (evedesign-wrapped)

This model scores protein sequences with the [structuredca](.)
`StructureDCA` wrapper around the
[StructureDCA](https://github.com/3BioCompBio/StructureDCA) package, which fits a
sparse, structure-informed pairwise Potts model of a protein family: couplings 
between residue pairs that are not in structural contact (per a residue distance 
cutoff on a provided 3D structure) are constrained to zero. Scoring is zero-shot: 
each sequence is assigned its (negative) statistical energy under the fitted model, 
no labelled training fold is required.

## License

StructureDCA's code (the `structuredca` PyPI package, depended on directly, not
vendored) is released under the MIT license, Copyright (c), 2026, Matsvei
Tsishyn, Hugo Talibart, Marianne Rooman and Fabrizio Pucci (see `LICENSE`,
copied verbatim from the upstream repository). There are no model weights to
separately license
