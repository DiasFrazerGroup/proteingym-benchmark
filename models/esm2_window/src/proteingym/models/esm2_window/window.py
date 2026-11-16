"""Windowed masked-marginal scoring for fixed-context protein language models on long proteins.

Shared by the `esm2_window` and `esmc_window` model folders (copied verbatim into each, so that every model folder
stays self-contained, as in the rest of the zoo). Documented in docs/clinical/models_long_proteins.md (§3.1).

Scoring (masked marginal, as evedesign's `ESM2.score_mutants`): for every unique mutated position p, the model sees a
window of `window` residues containing p, with p masked; a mutant's score is the sum over its substitutions of
log p(mut | window) - log p(wt | window). The window is centred on p and shifted to fit inside the protein, so for
L <= window it is the whole protein and the result equals the zoo's `score_mutants` (test T1).

Special tokens (`special_tokens`):
  "termini"       <cls> only if the window starts at residue 1, <eos> only if it ends at the last residue: the model
                  is told where the protein really starts and ends (dark_energy convention, esm1v_core.py).
  "every_window"  <cls> ... <eos> around every window, as the tokenizer adds by default (ProteinGym v1 ESM scripts,
                  2026-11-06 prototype).
Both are identical for L <= window.
"""
from typing import Sequence

import torch
from evedesign.model import assign_scores_to_instances
from evedesign.system import Mutation, System, SystemInstance

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


class WindowedMaskedMarginal:
    """Mixin for an evedesign ESM2-style class (attributes tokenizer, model, device, batch_size, _load_model).

    Removes the length cap and scores with per-position windows. Put it first in the bases.
    """

    window: int = 1022
    special_tokens: str = "termini"

    def _init_window(self, window: int, special_tokens: str) -> None:
        if special_tokens not in SPECIAL_TOKENS:
            raise ValueError(f"special_tokens must be one of {SPECIAL_TOKENS}, got {special_tokens!r}")
        if window < 1:
            raise ValueError("window must be positive")
        self.window = int(window)
        self.special_tokens = special_tokens
        self.max_seq_length = None

    # --- no length cap: same checks as the zoo's ESMC, minus the length ---------------------------------------------
    @classmethod
    def can_model(cls, system: System, data: None = None) -> tuple[bool, str]:
        if data is not None:
            return False, "Model does not support data parameter (must be None)"
        if len(system) != 1 or system[0].type != "protein":
            return False, "Can only handle single-component protein system"
        if not system[0].defined_sequence():
            return False, "Entity must have defined rep sequence"
        return True, ""

    def build(self, system, data=None, status_callback=None):
        self.can_model_or_raise(system, data)
        self._system = system
        self.encoding = None
        self.token_ids = None
        return self

    def _validate_instances_and_max_length(self, instances: Sequence[SystemInstance]) -> None:
        for instance in instances:
            self.system.valid_instance(instance, validate_reps=True, require_reps=True, fixed_length=True,
                                       allow_deletions=False, raise_invalid=True)

    # --- scoring -----------------------------------------------------------------------------------------------------
    def score(self, instances, status_callback=None):
        """Full-length mutant sequences -> masked-marginal scores relative to the wild type."""
        self.ready_or_raise()
        self._validate_instances_and_max_length(instances)
        target = self.system[0]
        ref_seq = "".join(target.rep)
        mutants = [[Mutation(entity=0, pos=target.first_index + i, ref=r, to=a)
                    for i, (r, a) in enumerate(zip(ref_seq, "".join(inst[0].rep))) if a != r]
                   for inst in instances]
        scored = self.score_mutants(self.system.rep_to_instance(), mutants, status_callback)
        return assign_scores_to_instances(instances, [inst.score for inst in scored])

    def _window_inputs(self, seq: str, positions: list[int]):
        """Token ids (padded) and the index of the masked token for each position."""
        tok = self.tokenizer
        L = len(seq)
        rows, mask_idx = [], []
        for p in positions:
            s, e = window_bounds(p, L, self.window)
            ids = tok(seq[s:e], add_special_tokens=False)["input_ids"]
            if len(ids) != e - s:
                raise ValueError(f"tokenizer is not one token per residue ({len(ids)} tokens for {e - s} residues)")
            ids[p - s] = tok.mask_token_id
            start_tok = self.special_tokens == "every_window" or s == 0
            end_tok = self.special_tokens == "every_window" or e == L
            row = ([tok.cls_token_id] if start_tok else []) + ids + ([tok.eos_token_id] if end_tok else [])
            rows.append(row)
            mask_idx.append(p - s + (1 if start_tok else 0))
        width = max(len(r) for r in rows)
        input_ids = torch.full((len(rows), width), tok.pad_token_id, dtype=torch.long)
        attention = torch.zeros((len(rows), width), dtype=torch.long)
        for i, r in enumerate(rows):
            input_ids[i, :len(r)] = torch.tensor(r)
            attention[i, :len(r)] = 1
        return input_ids.to(self.device), attention.to(self.device), mask_idx

    def masked_logprobs(self, seq: str, positions: list[int]) -> dict[int, torch.Tensor]:
        """0-based position -> log-softmax over the vocabulary at that position, masked, in its window."""
        self._load_model()
        out = {}
        with torch.no_grad():
            for b in range(0, len(positions), self.batch_size):
                batch = positions[b:b + self.batch_size]
                input_ids, attention, mask_idx = self._window_inputs(seq, batch)
                logits = self.model(input_ids=input_ids, attention_mask=attention).logits
                for i, p in enumerate(batch):
                    out[p] = torch.log_softmax(logits[i, mask_idx[i]].float(), dim=-1).cpu()
        return out

    def score_mutants(self, instance, mutants, status_callback=None):
        self.ready_or_raise()
        self.system.valid_mutants(instance, mutants, deletions=False, insertions=False, raise_invalid=True)
        target = self.system[0]
        seq = "".join(instance[0].rep)
        positions = sorted({s.pos - target.first_index for m in mutants for s in m})
        logp = self.masked_logprobs(seq, positions)
        tid = self.tokenizer.convert_tokens_to_ids
        scores = []
        for m in mutants:
            t = 0.0
            for sub in m:
                p = sub.pos - target.first_index
                if seq[p] == sub.to:
                    continue
                t += (logp[p][tid(sub.to)] - logp[p][tid(seq[p])]).item()
            scores.append(t)
        out = self.system.mutate(instance, mutants)
        for inst, sc in zip(out, scores):
            inst.score, inst.confidence = sc, None
        return out
