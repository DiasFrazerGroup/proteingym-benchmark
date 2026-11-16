"""ESM-IF1 for long proteins: structure fragments + tiles (docs/clinical/models_long_proteins.md §3.2).

Zoo ESM-IF1 score: mean autoregressive log-likelihood of the full mutant sequence given the backbone, one unbatched
pass per mutant, on one structure numbered like the wild type.

Here:
- Frames: every structure of the entity is placed on the wild type by sequence (fragments.place). A frame is one
  structure's backbone over the wild-type span it covers. With one structure numbered like the wild type there is
  one frame over the whole protein. Fragments are separate frames (their coordinates are not in a common frame).
- Each mutant goes to the frame where its (mean) mutated position is most central, then to the tile of `window`
  residues (tiles every window//2 inside the frame's span) where it is most central, and is scored as
      score = sum_tile log p(mutant | tile backbone) - sum_tile log p(wt | tile backbone),
  batched (`batch_size` sequences per forward pass).
- For one full-length structure and L <= window this is L * (zoo(mut) - zoo(wt)): a fixed linear transform of the
  zoo score within a protein (same ranking, same Pearson/Spearman per protein).
Residues that differ between a structure and the wild type are allowed (the zoo raises); the wild type is scored.
"""
import numpy as np
import torch
import torch.nn.functional as F
import biotite.structure as struc
import esm.inverse_folding.util as if_util
from evedesign.model import assign_scores_to_instances
from evedesign.utils import ensure_sequence

from . import fragments
from .model import ESMIF1


def tiles(start: int, end: int, window: int) -> list[tuple[int, int]]:
    """[a, b) tiles of `window` residues every window//2 covering [start, end] (inclusive end)."""
    n = end - start + 1
    if n <= window:
        return [(start, end + 1)]
    step = window // 2
    s = list(range(start, end - window + 2, step))
    if s[-1] != end - window + 1:
        s.append(end - window + 1)
    return [(a, a + window) for a in s]


def _chain(structure):
    return list(ensure_sequence(structure))[0]


class ESMIF1Window(ESMIF1):
    name: str = "ESMIF1Window"

    def __init__(self, window: int = 1022, batch_size: int = 16, **kwargs):
        super().__init__(**kwargs)
        self.window, self.batch_size = int(window), int(batch_size)
        self._frames = None   # list of (name, start, end, coords (L, 3, 3) indexed by wild-type position)

    def build(self, system, data=None, status_callback=None):
        target = system[0]
        if self.structure_name is not None or len(target.structures) == 1:
            super().build(system, data, status_callback)
            if self._multichain:
                raise ValueError("ESMIF1Window supports single-chain structures only")
            L = len(target.rep)
            self._frames = [("single", 0, L - 1, self._coords)]
            return self
        # several structures (fragments): place each, extract its backbone in wild-type coordinates
        self.can_model_or_raise(system, data)
        self._system = system
        self._model, self._alphabet = self._load_model()
        self._multichain = False
        L = len(target.rep)
        self._frames = []
        for pl in fragments.placements_for(target, _chain):
            atoms = _chain(target.structures[pl.name]).atom_array
            bb = atoms[struc.filter_backbone(atoms)]
            res_ids, _ = struc.get_residues(bb)
            coords, _ = if_util.extract_coords_from_structure(bb)
            full = np.full((L, 3, 3), np.nan, dtype=np.float32)
            for i, rid in enumerate(res_ids):
                p = int(rid) + pl.shift
                if 0 <= p < L:
                    full[p] = coords[i]
            self._frames.append((pl.name, pl.start, pl.end, full))
        self._coords = None
        self._scoring_coords = None
        return self

    @property
    def ready(self) -> bool:
        return self._system is not None and self._frames is not None

    def _ll_sums(self, coords, seqs):
        conv = if_util.CoordBatchConverter(self._alphabet)
        out = []
        for b0 in range(0, len(seqs), self.batch_size):
            batch = [(coords, None, s) for s in seqs[b0:b0 + self.batch_size]]
            coords_t, conf, _, tokens, pad = conv(batch, device=self.device)
            prev, target = tokens[:, :-1], tokens[:, 1:]
            logits, _ = self._model.forward(coords_t, pad, conf, prev)
            loss = F.cross_entropy(logits, target, reduction="none")
            keep = target != self._alphabet.padding_idx
            out.extend((-(loss * keep).sum(dim=1)).detach().cpu().numpy().tolist())
        return np.array(out)

    def score(self, instances, status_callback=None):
        self.ready_or_raise()
        self._validate_instances(instances)
        ref = "".join(self.system[0].rep)
        seqs = ["".join(x[0].rep) for x in instances]
        units = {}   # (frame index, tile) -> instance indices
        for k, s in enumerate(seqs):
            diff = [i for i, (a, b) in enumerate(zip(ref, s)) if a != b]
            c = float(np.mean(diff)) if diff else len(ref) / 2
            fr = [(min(c - st, en - c), f) for f, (_, st, en, _) in enumerate(self._frames)
                  if all(st <= i <= en for i in diff)]
            if not fr:
                raise ValueError(f"no structure covers all mutated positions of instance {k}")
            f = max(fr)[1]
            _, st, en, _ = self._frames[f]
            T = tiles(st, en, self.window)
            cand = [t for t in T if all(t[0] <= i < t[1] for i in diff)] or T
            t = min(cand, key=lambda x: abs((x[0] + x[1] - 1) / 2 - c))
            units.setdefault((f, t), []).append(k)
        scores = np.zeros(len(seqs))
        with torch.no_grad():
            for (f, (a, b)), idx in units.items():
                coords = self._frames[f][3][a:b]
                ll = self._ll_sums(coords, [ref[a:b]] + [seqs[i][a:b] for i in idx])
                scores[idx] = ll[1:] - ll[0]
        return assign_scores_to_instances(instances, scores)
