"""SaProt for long proteins: structure fragments + windowed masked marginal (docs/clinical/models_long_proteins.md §3.2).

Differences from the zoo's SaProt (everything else is inherited):
- Structures: every structure of the entity is placed on the wild type by sequence (fragments.place) and each
  position takes its 3Di token from the structure in which it is most central (fragments.assign). With one
  structure numbered like the wild type this is exactly the zoo's token string. Residues that differ between a
  structure and the wild type (e.g. UniProt version changes) keep the structure's 3Di token; the zoo raises instead.
- Scoring: each masked position is scored inside a `window`-residue window of (residue, 3Di) tokens centred on it,
  with <cls>/<eos> per `special_tokens` (as esm2_window). For L <= window this is the zoo's computation.
"""
import numpy as np
import torch
from evedesign.utils import ensure_sequence

from . import fragments
from .model import FOLDSEEK_SEQ_VOCAB, FOLDSEEK_STRUC_VOCAB, MASK, SaProt, foldseek_3di

SPECIAL_TOKENS = ("termini", "every_window")


def window_bounds(pos: int, length: int, window: int) -> tuple[int, int]:
    """[start, end) of the window for 0-based `pos`: centred on it, shifted to fit inside the sequence."""
    if length <= window:
        return 0, length
    start = max(0, pos - window // 2)
    end = start + window
    if end > length:
        end, start = length, length - window
    return start, end


def _chain(structure):
    chain = list(ensure_sequence(structure))[0]
    if len(chain.chains()) != 1:
        raise ValueError(f"Expected a single-chain structure, got chains={chain.chains()}")
    return chain


class SaProtWindow(SaProt):
    name: str = "SaProtWindow"

    def __init__(self, window: int = 1022, special_tokens: str = "termini", **kwargs):
        super().__init__(**kwargs)
        if special_tokens not in SPECIAL_TOKENS:
            raise ValueError(f"special_tokens must be one of {SPECIAL_TOKENS}")
        self.window, self.special_tokens = int(window), special_tokens
        self.structure_assignment = None   # per position: structure name (for reporting)

    def _structure_tokens(self, target, rep: str) -> str:
        if self.structure_name is not None:          # explicit choice: zoo behaviour
            return super()._structure_tokens(target, rep)
        placements = fragments.placements_for(target, _chain)
        assigned = fragments.assign(placements, len(rep))
        tokens = [MASK] * len(rep)
        for pl in placements:
            res_ids, struc_seq = foldseek_3di(_chain(target.structures[pl.name]), self.foldseek_path)
            for rid, s in zip(res_ids, struc_seq):
                p = int(rid) + pl.shift
                if 0 <= p < len(rep) and assigned[p] == pl.name:
                    tokens[p] = s if s in FOLDSEEK_STRUC_VOCAB else MASK
        n_missing = tokens.count(MASK)
        if n_missing and not self.allow_missing_structure:
            raise ValueError(f"Structures do not cover {n_missing}/{len(rep)} target positions")
        self.structure_assignment = assigned
        return "".join(tokens)

    def _masked_group_logp(self, seqs, masked):
        if all(len(s) <= self.window for s in seqs):
            return super()._masked_group_logp(seqs, masked)
        self._load_model()
        vocab = self._tokenizer.get_vocab()
        cls_id, eos_id = self._tokenizer.cls_token_id, self._tokenizer.eos_token_id
        struc = self.structure_tokens
        items = []
        for seq, mpos in zip(seqs, masked):
            if len(mpos) != 1:
                raise ValueError("SaProtWindow scores one masked position per sequence on long proteins")
            p, L = mpos[0], len(seq)
            s, e = window_bounds(p, L, self.window)
            aas = [aa if aa in FOLDSEEK_SEQ_VOCAB else MASK for aa in seq.upper()[s:e]]
            aas[p - s] = MASK
            start_tok = self.special_tokens == "every_window" or s == 0
            end_tok = self.special_tokens == "every_window" or e == L
            ids = [vocab[a + t] for a, t in zip(aas, struc[s:e])]
            items.append(([cls_id] * start_tok + ids + [eos_id] * end_tok, p - s + start_tok))
        out = np.empty((len(items), self._group_index.shape[0]), dtype=np.float32)
        for b0 in range(0, len(items), self.batch_size):
            chunk = items[b0:b0 + self.batch_size]
            width = max(len(c[0]) for c in chunk)
            pad = self._tokenizer.pad_token_id
            ids = torch.tensor([c[0] + [pad] * (width - len(c[0])) for c in chunk], device=self.device)
            att = torch.tensor([[1] * len(c[0]) + [0] * (width - len(c[0])) for c in chunk], device=self.device)
            with torch.no_grad():
                logits = self._model(input_ids=ids, attention_mask=att).logits.float()
                for i, (_, idx) in enumerate(chunk):
                    row = torch.log_softmax(logits[i, idx], dim=-1)
                    out[b0 + i] = torch.logsumexp(row[self._group_index], dim=-1).cpu().numpy()
        return [out[i:i + 1] for i in range(len(items))]
