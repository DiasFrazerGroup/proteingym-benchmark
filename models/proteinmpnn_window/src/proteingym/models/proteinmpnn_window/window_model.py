"""ProteinMPNN for long proteins: AlphaFold DB fragments (docs/clinical/models_long_proteins.md §3.7).

- One structure (every protein <= 2,700 aa): the zoo model unchanged (score = mean log-likelihood of the sequence
  given the full-length backbone).
- Several structures (AFDB fragments of proteins > 2,700 aa): every fragment is placed on the wild type by sequence
  (fragments.place) and renumbered in wild-type coordinates; the zoo model is built on each fragment separately (its
  own backbone, its own decoding order drawn from `seed`). Each mutant goes to the fragment where its (mean) mutated
  position is most central and is scored there as
      score = n_frag * (mean_ll(mutant | fragment) - mean_ll(wild type | fragment)),
  i.e. the summed log-likelihood difference over the fragment's residues (same rule as esm-if1_window).
  A mutant not covered by any fragment is NaN.
The wild type over each fragment is the entity's wild type (residues where the structure differs are scored as the
wild type, as the zoo wrapper does for any structure).
"""
import numpy as np
from evedesign.model import assign_scores_to_instances
from evedesign.proteingym import seqs_to_instances
from evedesign.utils import ensure_sequence

from . import fragments
from .model import ProteinMPNN


def _chain(structure):
    return list(ensure_sequence(structure))[0]


class ProteinMPNNWindow(ProteinMPNN):
    name: str = "ProteinMPNNWindow"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._kwargs = dict(kwargs)
        self._frames = None   # list of (Placement, fragment model, n residues)

    def build(self, system, data=None, status_callback=None):
        target = system[0]
        self._frames = None
        if len(target.structures) == 1:
            return super().build(system, data, status_callback)
        self.can_model_or_raise(system, data)
        self._system = system
        structures = target.structures
        frames = []
        try:
            for pl in fragments.placements_for(target, _chain):
                chain = _chain(structures[pl.name]).copy()
                chain.atom_array.res_id = chain.atom_array.res_id + pl.shift + target.first_index
                target.structures = {pl.name: chain}
                m = ProteinMPNN(**self._kwargs).build(system, data)
                frames.append((pl, m, len(m._native_seq)))
        finally:
            target.structures = structures
        self._frames = frames
        return self

    @property
    def ready(self) -> bool:
        return self._system is not None and (self._frames is not None or self._feature_dict is not None)

    def score(self, instances, status_callback=None):
        if self._frames is None:
            return super().score(instances, status_callback)
        self.ready_or_raise()
        self._validate_instances(instances)
        ref = "".join(self.system[0].rep)
        groups = {}
        scores = np.full(len(instances), np.nan)
        for k, x in enumerate(instances):
            s = "".join(x[0].rep)
            diff = [i for i, (a, b) in enumerate(zip(ref, s)) if a != b]
            c = float(np.mean(diff)) if diff else len(ref) / 2
            cand = [(min(c - pl.start, pl.end - c), f) for f, (pl, _, _) in enumerate(self._frames)
                    if all(pl.start <= i <= pl.end for i in diff)]
            if cand:
                groups.setdefault(max(cand)[1], []).append(k)
        wt = seqs_to_instances([ref])[0]
        for f, idx in groups.items():
            _, m, n = self._frames[f]
            v = np.array([float(x.score) for x in m.score([wt] + [instances[i] for i in idx])])
            scores[idx] = n * (v[1:] - v[0])
        return assign_scores_to_instances(instances, scores)

