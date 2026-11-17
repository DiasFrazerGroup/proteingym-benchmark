---
# Model identifier used for referencing this model in the benchmark system
name: rsalor_filtered

# Import string for the evedesign model class to run. Accepts entry-point style
# ("package.module:ClassName") or fully dotted style ("package.module.ClassName").
# The class is instantiated with hyper_parameters, then built and scored.
model_class: "proteingym.models.rsalor_filtered.filtered_model:RSALORFiltered"

tags: ["zero-shot"]

# Multi-target support flag
multi_y: false

# This is a CPU-only model; there is no "device"/GPU hyperparameter
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
    # Regularization term for LOR/LR at the amino-acid frequencies level
    theta_regularization: 0.01
    # Regularization term for LOR/LR at the amino-acid counts level
    n_regularization: 0.0
    # Count the target/WT sequence itself in the frequency computation
    count_target_sequence: true
    # Pre-process the MSA to remove exactly-duplicate sequences before
    # weighting
    remove_redundant_sequences: true
    # Sequence-identity threshold for clustering sequences into weight groups
    # (EVE/EVcouplings-style reweighting, computed by RSALOR's own C++
    # backend)
    seqid_weights: 0.80
    # Discard MSA sequences whose identity to the target sequence is below
    # this threshold before weighting
    min_seqid: 0.35
    # Number of CPU threads for RSALOR's C++ sequence-weighting backend
    num_threads: 4
    # RSA solver run against the entity's structure. "biopython" needs no 
    # external binary; "DSSP"/"MuSiC" are supported but require executable 
    # to be installed, which this wrapper/container does not provide
    rsa_solver: "biopython"
    # "LOR" (log odd ratio) or "LR" (log ratio)
    metric: "LOR"
    # Multiply the per-position metric by the RSA-derived weight before
    # taking differences ("RSA*LOR", the method this wrapper is named/
    # submitted for). Structure is required regardless of this flag
    use_rsa_factor: true
    # Key into the entity's structures mapping to select a structure, if more
    # than one is present leave null 
    structure_name: null
    verbose: false
---

# Model Card for rsalor_filtered (clinical: EVE-protocol alignment)

Same model as the zoo's `rsalor` (wrapper unchanged), trained on the dataset's MSA after the EVE/popEVE alignment
filter: insertions removed, fragments (> 50% gaps) and rows with indeterminate residues dropped (EVE code), a
bit-score threshold (MMseqs2 bits/residue from the ColabFold a3m headers) chosen for column coverage L_cov >= 0.8 L
and 10 L <= N <= 100,000, relaxed as in the EVE paper to L_cov >= 0.7 L and N <= 200,000, then "best available"
(documented order); all columns are kept, as in the popEVE release, which scores every position. Every `*_filtered` model uses the same filter and parameters, so all alignment models see the same alignment.
Variants at positions the model does not score are NaN (the zoo wrapper scores them as wild type). Per-protein filter
statistics (rows in / after coverage, threshold, N, N/L, L_cov, relaxation step, quality flag) are logged and
appended to `msa_filter_stats.jsonl` in the output folder. See `docs/clinical/models_long_proteins.md` §3.6.

---


# Model Card for rsalor (evedesign-wrapped)

This model scores protein sequences with the
`proteingym.models.rsalor.model.rsalor` wrapper around
[RSALOR](https://github.com/3BioCompBio/RSALOR), which combines two
CPU-side, structure/alignment-derived quantities into a  per-position +
per-substitution score: (regularized, sequence-identity-reweighted) MSA
column's Log Odd Ratio (**LOR**) between the wt and mutant amino acid,
scaled by Relative Solvent Accessibility (**RSA**) weight computed from a 3D
structure (`RSA_factor = 1 - min(RSA, 100) / 100` Scoring is zero-shot,
 wrapper's `build()` computes this "RSA*LOR" matrix once from the alignment and
structure, and `score()`/`single_mutation_scan()`/`score_conditional()` just
read off

**Sign convention note:** upstream's own `MSA.eval_mutations()`/`get_scores()`
deliberately produce a score that is *positive for destabilizing/disruptive*
mutations ("large positive values predict highly destabilizing / disruptive
mutations" - see upstream README). This wrapper flips that sign so
*higher = more fit/tolerated* and the WT always scores 0 (see `model.py`'s
module and `build()` docstrings for exact formula)

## Architecture notes

- **No GPU, no pretrained checkpoint.** RSALOR has no learned weights -
  everything (sequence reweighting, frequency estimation, RSA) is computed
  from the MSA + structure at `build()` time.
- **Real CPU parallelism.** RSALOR's C++ sequence-weighting backend
  takes `num_threads` argument and parallelizes the pairwise
  sequence-identity computation - `supports_cpu_parallel=True`.
- **No indel support.** RSALOR's `Mutation` class only accepts standard
  amino-acid one-letter codes for the wild-type and mutant symbol (no
  gap); it is a fixed per-position profile (one score per aligned MSA
  column), with no notion of a position that doesn't exist in the target at
  all. `requires_fixed_length=True`, `handles_deletions=False`,
  `handles_insertions=False`
- **Cost driver.** The dominant cost of `build()` is an O(depth^2 x length)
  all-pairs sequence-identity computation used for sequence reweighting
  RSA computation is linear in structure size and comparatively cheap

## Citation

```
@article{Tsishyn2025RSALOR,
  author = {Tsishyn, Matsvei and Hermans, Pauline and Rooman, Marianne and Pucci, Fabrizio},
  title = {Residue conservation and solvent accessibility are (almost) all you need for predicting mutational effects in proteins},
  journal = {Bioinformatics},
  volume = {41},
  number = {6},
  pages = {btaf322},
  year = {2025},
  doi = {10.1093/bioinformatics/btaf322}
}

@article{Hermans2025Exploring,
  author = {Hermans, Pauline and Tsishyn, Matsvei and Schwersensky, Martin and Rooman, Marianne and Pucci, Fabrizio},
  title = {Exploring evolution to uncover insights into protein mutational stability},
  journal = {Molecular Biology and Evolution},
  volume = {42},
  number = {1},
  pages = {msae267},
  year = {2025},
  doi = {10.1093/molbev/msae267}
}
```

## License

RSALOR's code (the `rsalor` PyPI package, depended on directly, not vendored)
is released under the MIT license, Copyright (c) 2025 Matsvei Tsishyn (see
`LICENSE`, copied verbatim from the upstream repository). There are no model
weights to separately license
