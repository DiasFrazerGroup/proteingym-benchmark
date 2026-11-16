---
# Model identifier used for referencing this model in the benchmark system
name: esm2_window

# Import string for the evedesign model class to run. Accepts entry-point style
# ("package.module:ClassName") or fully dotted style ("package.module.ClassName").
# The class is instantiated with hyper_parameters, then built and scored.
model_class: "proteingym.models.esm2_window.model:ESM2Window"

tags: ["zero-shot", "embeddable"]

# Multi-target support flag
multi_y: false

# Hyperparameters are forwarded verbatim to the evedesign ESM2 wrapper
# constructor (evedesign.models.esm2.ESM2), so every key here must be a valid
# ESM2.__init__ argument.
hyper_parameters:
    # clinical (long proteins): residues per window, one window per mutated position, masked marginal
    window: 1022
    # "termini": <cls>/<eos> only at the true protein ends; "every_window": around every window
    special_tokens: "termini"
    # ESM-2 checkpoint identifier (loaded from facebook/<model_name> on the Hub)
    model_name: "esm2_t33_650M_UR50D"
    # Number of sequences scored per forward pass
    batch_size: 64
    # Device to run the model on ("cpu", "cuda", or "mps")
    device: "cpu"
---

# Model Card for esm2_window (clinical: long proteins)

Same model and checkpoint as the zoo's `esm2`, with windowed masked-marginal scoring so that proteins of any length can be scored: for each mutated position, a `window`-residue window centred on it (shifted at the protein ends) with the position masked; score = log p(mut) - log p(wt), summed over substitutions. For proteins no longer than the window this equals the zoo's masked marginal exactly. See `docs/clinical/models_long_proteins.md` (section 3.1) for the method, the `special_tokens` option and the validation.

---

# Model Card for ESM-2 (evedesign-wrapped)

This model scores protein sequences with the [evedesign](https://github.com/evedesignbio/evedesign)
`ESM2` wrapper, which loads an ESM-2 masked language model via HuggingFace
`transformers`. Scoring is zero-shot: each sequence is assigned the sum of the
per-position (unmasked) log-likelihoods of its residues under the model.

The benchmark entrypoint converts the ProteinGym subset into evedesign
datatypes with `evedesign.proteingym.dataset_to_evedesign`, builds the wrapper
on the resulting `System`, and scores the requested test fold. Because the
wrapper exposes the full evedesign model API, additional capabilities such as
`transform()` (embeddings), `generate()` (Gibbs sampling), and mutation scoring
are also available.

For detailed information on the model architecture and training data,
please refer to the [accompanying paper](https://www.biorxiv.org/content/10.1101/2022.07.20.500902v2).

| Checkpoint name | Num layers | Num parameters |
|------------------------------|----|----------|
| [esm2_t48_15B_UR50D](https://huggingface.co/facebook/esm2_t48_15B_UR50D) | 48 | 15B     |
| [esm2_t36_3B_UR50D](https://huggingface.co/facebook/esm2_t36_3B_UR50D) | 36 | 3B      |
| [esm2_t33_650M_UR50D](https://huggingface.co/facebook/esm2_t33_650M_UR50D) | 33 | 650M    |
| [esm2_t30_150M_UR50D](https://huggingface.co/facebook/esm2_t30_150M_UR50D) | 30 | 150M    |
| [esm2_t12_35M_UR50D](https://huggingface.co/facebook/esm2_t12_35M_UR50D) | 12 | 35M     |
| [esm2_t6_8M_UR50D](https://huggingface.co/facebook/esm2_t6_8M_UR50D)  | 6  | 8M      |
