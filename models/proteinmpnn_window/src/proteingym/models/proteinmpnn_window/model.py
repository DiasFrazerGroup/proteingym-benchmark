"""ProteinMPNN through evedesign's LigandMPNN wrapper (`evedesign.models.mpnn:LigandMPNN`, model_type protein_mpnn).

The wrapper scores a sequence as the mean log-likelihood of its residues given the backbone, decoded in one random
order (ProteinMPNN `score` with `use_sequence=True`), the same quantity as ProteinGym's ProteinMPNN baseline
(`proteingym/baselines/protein_mpnn/compute_fitness.py`: mutant-sequence global score, one random order per mutant).
With `fix_full_decoding_order=True` the wrapper draws one order at build time and reuses it for every sequence;
`torch.randn` is unseeded there, so this subclass only seeds torch before build (`seed`), making scores reproducible.
Nothing else is changed.
"""
from evedesign.models.mpnn import LigandMPNN

try:
    import torch
except ImportError:  # pragma: no cover
    torch = None


class ProteinMPNN(LigandMPNN):
    name: str = "ProteinMPNN"
    citations: list[str] = ["doi: 10.1126/science.add2187"]

    def __init__(self, seed: int | None = 0, **kwargs):
        super().__init__(**kwargs)
        self.seed = seed

    def build(self, system, data=None, status_callback=None):
        if self.seed is not None:
            torch.manual_seed(int(self.seed))
        return super().build(system, data, status_callback)
