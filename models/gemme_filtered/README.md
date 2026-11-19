---
# Model identifier used for referencing this model in the benchmark system
name: gemme_filtered

# Import string for the evedesign model class to run. Accepts entry-point style
# ("package.module:ClassName") or fully dotted style ("package.module.ClassName").
# The class is instantiated with hyper_parameters, then built and scored.
model_class: "proteingym.models.gemme_filtered.filtered_model:GEMMEFiltered"

tags: ["zero-shot"]

# Multi-target support flag
multi_y: false

# Hyperparameters are forwarded verbatim to the GEMME wrapper constructor
# (proteingym.models.gemme.model.GEMME). CPU-only; there is no "device" option.
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
    # Number of JET2 (iJET) resampling iterations used to compute conservation
    # (10 matches GEMME's published usage)
    n_iter: 10
    # Maximum number of alignment sequences (including the query) passed to JET2
    max_sequences: 20000
    # Path to the GEMME directory (computePred.R, pred.R, default.conf, ...)
    gemme_path: "/opt/program/GEMME"
    # Path to the JET2 directory (jet/JET.class, jet/extLibs/vecmath.jar)
    jet_path: "/opt/program/JET2"
    # Path to / name of the java binary used to run JET2
    java_binary: "java"
    # Path to / name of the Rscript binary (needs seqinr and RColorBrewer)
    rscript_binary: "Rscript"
    # Path to / name of muscle 3.8.x (muscle 5.x does not work with JET2)
    muscle_binary: "muscle"
    # -Xmx value passed to the JET2 JVM
    java_max_heap: "4096m"
    # clinical: score = computePred.R output "evolEpi" (zoo wrapper) or "evolCombi" (GEMME combined model, as
    # ProteinGym reports it). Both are benchmarked as separate prediction sets (user 2026-11-19)
    output: "evolEpi"
---

# Model Card for gemme_filtered (clinical: EVE-protocol alignment)

Same model as the zoo's `gemme` (wrapper unchanged), trained on the dataset's MSA after the EVE/popEVE alignment
filter: insertions removed, fragments (> 50% gaps) and rows with indeterminate residues dropped (EVE code), a
bit-score threshold (MMseqs2 bits/residue from the ColabFold a3m headers) chosen for column coverage L_cov >= 0.8 L
and 10 L <= N <= 100,000, relaxed as in the EVE paper to L_cov >= 0.7 L and N <= 200,000, then "best available"
(documented order); all columns are kept, as in the popEVE release, which scores every position. Every `*_filtered` model uses the same filter and parameters, so all alignment models see the same alignment.
Variants at positions the model does not score are NaN (the zoo wrapper scores them as wild type). Per-protein filter
statistics (rows in / after coverage, threshold, N, N/L, L_cov, relaxation step, quality flag) are logged and
appended to `msa_filter_stats.jsonl` in the output folder. See `docs/clinical/models_long_proteins.md` §3.6.

---


# Model Card for GEMME (evedesign-wrapped)

GEMME (Global Epistatic Model for predicting Mutational Effects) is an MSA-based
model that combines per-position conservation from JET2 with a "global epistasis"
model of evolutionary distance to homologs. It needs a multiple sequence alignment
of the target family; no structure is used. Each substitution gets GEMME's evolEpi
score, and a sequence is scored as the sum over its residues. Substitutions only;
insertions and deletions are not scored. CPU-only.

For details on the method, see the [paper](https://doi.org/10.1093/molbev/msz179).

## Citation

```
@article{Laine2019GEMME,
  author = {Laine, Elodie and Karami, Yasaman and Carbone, Alessandra},
  title = {GEMME: A Simple and Fast Global Epistatic Model Predicting Mutational Effects},
  journal = {Molecular Biology and Evolution},
  volume = {36},
  number = {11},
  pages = {2604--2619},
  year = {2019},
  doi = {10.1093/molbev/msz179}
}
```

## License

GEMME is released under the MIT license (see `LICENSE`). JET2 and muscle are
downloaded at build time and are not redistributed in this folder.
