---
# Model identifier used for referencing this model in the benchmark system
name: proteinmpnn_window

# Import string for the evedesign model class to run. Accepts entry-point style
# ("package.module:ClassName") or fully dotted style ("package.module.ClassName").
# The class is instantiated with hyper_parameters, then built and scored.
model_class: "proteingym.models.proteinmpnn_window.window_model:ProteinMPNNWindow"

tags: ["zero-shot"]

# Multi-target support flag
multi_y: false

# Forwarded to evedesign.models.mpnn.LigandMPNN (plus `seed`, see model.py)
hyper_parameters:
    # ProteinMPNN vanilla weights, 48 neighbours, 0.20 A training noise (upstream default)
    model_name: "proteinmpnn_v_48_020"
    model_file_path: "/opt/program/data/proteinmpnn_v_48_020.pt"
    # one decoding order drawn at build (torch seeded with `seed`) and reused for every sequence
    fix_full_decoding_order: true
    vary_decoding_order_per_instance: false
    seed: 0
    # protein-only system: no ligand atoms
    use_ligand_context: false
    # keep the weights in memory between build and score
    keep_model_after_build: true
    batch_size: 1
    device: "cpu"
---

# Model Card for proteinmpnn_window (clinical: long proteins)

Same model as the zoo's `proteinmpnn`, for proteins of any length. With one structure (every protein up to 2,700 aa
in AlphaFold DB) it is the zoo model unchanged. With several structures (AlphaFold DB fragments of proteins longer
than 2,700 aa) each fragment is placed on the wild type by sequence, the zoo model is built on each fragment, and
each mutant is scored in the fragment where it is most central, as the summed log-likelihood difference to the wild
type over that fragment. See `docs/clinical/models_long_proteins.md` (section 3.7).

---

# Model Card for ProteinMPNN (evedesign-wrapped)

ProteinMPNN is an inverse-folding message-passing network: given a protein backbone it predicts the amino acid at
every position, decoding the sequence autoregressively in a random order over a k-nearest-neighbour graph
(48 neighbours). It needs a structure (`entity.structures`) and no MSA.

Scoring is zero-shot: each sequence is scored as the mean log-likelihood of its residues given the wild-type
backbone (one decoding order for all sequences, fixed by `seed`), as in ProteinGym's ProteinMPNN baseline (which
draws a new random order per mutant). Higher scores mean fitter variants. Substitutions only.

The model runs through evedesign's `LigandMPNN` wrapper (`model_type` protein_mpnn), which vendors the LigandMPNN
code (MIT). Weights: `v_48_020` from [dauparas/ProteinMPNN](https://github.com/dauparas/ProteinMPNN) (commit
8907e66, sha256 pinned), downloaded at image build time.

```
@article{dauparas2022proteinmpnn,
  title={Robust deep learning-based protein sequence design using ProteinMPNN},
  author={Dauparas, J. and Anishchenko, I. and Bennett, N. and others},
  journal={Science},
  volume={378},
  pages={49--56},
  year={2022},
  doi={10.1126/science.add2187}
}
```
