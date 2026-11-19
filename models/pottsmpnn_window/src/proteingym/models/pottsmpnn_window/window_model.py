"""PottsMPNN for long proteins: AlphaFold DB fragments (docs/clinical/models_long_proteins.md §3.7).

- One structure (every protein <= 2,700 aa): the zoo model unchanged (score = E(wild type) - E(variant) under the
  Potts model of the full-length structure).
- Several structures (AFDB fragments of proteins > 2,700 aa): every fragment is placed on the wild type by sequence
  (fragments.place); the zoo's energy table is computed for each fragment (same model, same code path as the zoo
  build), and each mutant is scored with the energy table of the fragment where its (mean) mutated position is most
  central: score = E(wild type) - E(variant) over that fragment (energies are additive, so no rescaling).
  A mutant not covered by any fragment is NaN.
The reference sequence of each fragment is the entity's wild type over its span (residues where the AFDB fragment
differs from the wild type, e.g. APOB Y1422C, are scored as the wild type).
"""
import numpy as np
import torch
from evedesign.model import assign_scores_to_instances
from evedesign.utils import ensure_sequence

from . import fragments
from .model import PottsMPNN, _structure_to_pdb_dict, get_etab


def _chain(structure):
    return list(ensure_sequence(structure))[0]


class PottsMPNNWindow(PottsMPNN):
    name: str = "PottsMPNNWindow"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._frames = None   # list of (Placement, etab, E_idx, wt sequence over the fragment's residues)

    def build(self, system, data=None, status_callback=None):
        target = system[0]
        self._frames = None
        if self.structure_name is not None or len(target.structures) == 1:
            return super().build(system, data, status_callback)
        self.can_model_or_raise(system, data)
        self._system = system
        wt = "".join(target.rep)
        model = self._load_model()
        frames = []
        for pl in fragments.placements_for(target, _chain):
            chain = _chain(target.structures[pl.name])
            pdb_dict = _structure_to_pdb_dict(chain)
            lo = int(chain.atom_array.res_id.min()) + pl.shift          # wild-type index of the dict's first residue
            n = len(pdb_dict["seq"])
            with torch.no_grad():
                etab, E_idx, _ = get_etab(model, [pdb_dict], self._cfg(), None)
            frames.append((pl, lo, etab.double().cpu(), E_idx.cpu(), wt[lo:lo + n]))
        del model
        self._frames = frames
        self._etab = frames[0][2]          # marks the model as ready
        return self

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
            cand = [(min(c - fr[0].start, fr[0].end - c), f) for f, fr in enumerate(self._frames)
                    if all(fr[0].start <= i <= fr[0].end for i in diff)]
            if cand:
                groups.setdefault(max(cand)[1], []).append(k)
        for f, idx in groups.items():
            _, lo, etab, E_idx, wt_frag = self._frames[f]
            seqs = []
            for i in idx:
                s = list(wt_frag)
                for p, (a, b) in enumerate(zip(ref, "".join(instances[i][0].rep))):
                    if a != b:
                        s[p - lo] = b.upper()
                seqs.append("".join(s))
            self._etab, self._E_idx, self._wt_seq = etab, E_idx, wt_frag
            scores[idx] = -self._ddg(seqs)
        return assign_scores_to_instances(instances, scores)
