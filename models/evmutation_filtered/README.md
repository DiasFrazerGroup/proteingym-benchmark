---
# Model identifier used for referencing this model in the benchmark system
name: evmutation_filtered

# Import string for the evedesign model class to run. Accepts entry-point style
# ("package.module:ClassName") or fully dotted style ("package.module.ClassName").
# The class is instantiated with hyper_parameters, then built and scored.
model_class: "proteingym.models.evmutation_filtered.filtered_model:EVmutationFiltered"

tags: ["zero-shot"]

# Multi-target support flag
multi_y: false

# Hyperparameters are forwarded verbatim to the evedesign EVcouplingsPLM wrapper
# constructor (evedesign.models.evcouplings.EVcouplingsPLM), so every key here
# must be a valid EVcouplingsPLM.__init__ argument. This is a CPU-only model-
# there is no "device" hyperparameter.
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
    # bit scores are read from the dataset MSA's own a3m headers (raw bytes in the archive); archives built with
    # seq_<i> headers have none and the selection is skipped ("no_bitscores"). Test fallback only (never in
    # benchmark runs): a directory of the original ColabFold a3m files, matched by query, verified row by row
    a3m_headers_dir: null
    # per-protein filter stats appended to <stats_dir>/msa_filter_stats.jsonl (next to the predictions)
    stats_dir: "/opt/program/output"
    # --- zoo model hyper-parameters (unchanged) ---
    # Alignment columns whose (unweighted) gap frequency strictly exceeds this
    # threshold are excluded from the fitted model
    max_gap_fraction: 0.5
    # Sequence reweighting identity threshold, used as a fallback when the
    # provided Sequences carry no precomputed weights. Leave null to require
    # precomputed weights on the input. (will fail if no weights + null)
    theta: 0.8
    # L2 regularisation strength on fields h_i
    lambda_h: 0.01
    # L2 regularisation strength on couplings J_ij, before scaling
    lambda_J: 0.01
    # Scale lambda_J by the number of states and modelled positions, as in
    # standard EVcouplings/plmc practice
    lambda_J_times_Lq: true
    # Group L1 regularisation strength on couplings (null = plmc default)
    lambda_group: null
    # Scale weights of sequence clusters by this value (null = plmc default)
    scale_clusters: null
    # Maximum L-BFGS iterations
    iterations: 100
    # If true, exclude gaps from parameter inference (gaps then cannot be scored)
    ignore_gaps: false
    # If true, fit an independent-site (no couplings) model instead
    independent_model: false
    # Path to / name of the plmc binary; "plmc" resolves via PATH
    plmc_binary: "plmc"
    # Number of CPU cores for plmc to use ("max" or an integer). Requires an
    # OpenMP-compiled plmc, which is how this model's container builds it
    cpu: "max"
---

# Model Card for evmutation_filtered (clinical: EVE-protocol alignment)

Same model as the zoo's `evmutation` (wrapper unchanged), trained on the dataset's MSA after the EVE/popEVE alignment
filter: insertions removed, fragments (> 50% gaps) and rows with indeterminate residues dropped (EVE code), a
bit-score threshold (MMseqs2 bits/residue from the ColabFold a3m headers) chosen for column coverage L_cov >= 0.8 L
and 10 L <= N <= 100,000, relaxed as in the EVE paper to L_cov >= 0.7 L and N <= 200,000, then "best available"
(documented order); all columns are kept, as in the popEVE release, which scores every position. Every `*_filtered` model uses the same filter and parameters, so all alignment models see the same alignment.
Variants at positions the model does not score are NaN (the zoo wrapper scores them as wild type). Per-protein filter
statistics (rows in / after coverage, threshold, N, N/L, L_cov, relaxation step, quality flag) are logged and
appended to `msa_filter_stats.jsonl` in the output folder. See `docs/clinical/models_long_proteins.md` §3.6.

---


# Model Card for EVmutation (evedesign-wrapped)

This model scores protein sequences with the [evedesign](https://github.com/evedesignbio/evedesign)
`EVcouplingsPLM` wrapper, which fits a pairwise Potts model of a protein family
from a multiple sequence alignment using the plmc pseudo-likelihood solver. 
Scoring is zero-shot: each sequence is assigned its statistical energy 
(Hamiltonian) under the fitted model.

The wrapper also exposes `single_mutation_scan` and `score_conditional` for 
mutation-effect and conditional scoring, but does not implement the `Transformer` 
interface, so this card is not tagged `embeddable`.

For detailed information on the underlying method, please refer to
the [EVmutation paper](https://doi.org/10.1038/nbt.3769) (Hopf et al., *Nature
Biotechnology* 2017) and the [EVcouplings framework paper](https://doi.org/10.1093/bioinformatics/bty862)
(Hopf et al., *Bioinformatics* 2019).

The wrapper implements two fitting engines: `EVcouplingsMeanField` (mean-field
DCA, pure Python/numpy, no external dependency) and `EVcouplingsPLM` (used by
this card), which fits with the pseudo-likelihood solver via the external
`plmc` binary ([debbiemarkslab/plmc](https://github.com/debbiemarkslab/plmc)).

Regularisation strengths (`lambda_h`, `lambda_J`, `lambda_group`) and sequence
reweighting (`theta`, or precomputed weights on the input) govern the
statistical quality of the fit; `max_gap_fraction` controls which alignment
columns are modelled at all. `iterations` and `independent_model` trade fit
quality/runtime against a simpler, couplings-free baseline.
