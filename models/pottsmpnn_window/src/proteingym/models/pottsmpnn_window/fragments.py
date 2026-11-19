"""Structure fragments for long proteins: place each structure on the wild type, pick one structure per position.

AlphaFold DB splits proteins longer than 2,700 aa into fragments F1..Fn (1,400 aa, starting every 200 aa), each
numbered from 1. An archive can therefore carry several structures of the same protein, each covering a span of the
wild type with local residue numbers. Here every structure is placed on the wild type by its residue sequence (the
shift with the most identical residues), and every wild-type position gets the structure in which it is most
central (largest distance to the structure's ends). Single full-length structures numbered like the wild type are
the trivial case: shift 0, one structure for every position.

Copied verbatim into each structure-model folder (self-contained zoo folders). Documented in
docs/clinical/models_long_proteins.md (§3.2).
"""
from dataclasses import dataclass

import numpy as np
from biotite.sequence import ProteinSequence

MIN_IDENTITY = 0.9


@dataclass
class Placement:
    name: str            # key in entity.structures
    shift: int           # wild-type 0-based index = residue id + shift
    start: int           # first covered wild-type index (0-based)
    end: int             # last covered wild-type index (0-based, inclusive)
    identity: float      # fraction of the structure's residues identical to the wild type


def chain_residues(chain) -> tuple[np.ndarray, str]:
    """Residue ids and one-letter residues of the CA atoms of a single chain (evedesign Structure)."""
    aa = chain.atom_array
    ca = aa[(aa.atom_name == "CA") & ~aa.hetero]
    letters = "".join(ProteinSequence.convert_letter_3to1(n) if n in ProteinSequence._dict_3to1 else "X"
                      for n in ca.res_name)
    return np.asarray(ca.res_id, dtype=int), letters


def place(name: str, res_ids: np.ndarray, letters: str, wt: str) -> Placement:
    """Shift that maps residue ids onto the wild type with the most identical residues."""
    w = np.frombuffer(wt.encode(), dtype=np.uint8)
    s = np.frombuffer(letters.encode(), dtype=np.uint8)
    best = (-1, 0)
    lo, hi = -int(res_ids.min()), len(wt) - int(res_ids.max()) - 1   # shifts keeping every residue inside wt
    for k in range(lo, hi + 1):
        n = int((w[res_ids + k] == s).sum())
        if n > best[0]:
            best = (n, k)
    n, k = best
    ident = n / len(res_ids)
    if ident < MIN_IDENTITY:
        raise ValueError(f"structure {name}: best placement on the wild type has identity {ident:.2f}")
    return Placement(name, k, int(res_ids.min()) + k, int(res_ids.max()) + k, ident)


def assign(placements: list[Placement], length: int) -> list[str | None]:
    """Per wild-type position (0-based): the structure in which it is most central, None if uncovered."""
    out = []
    for p in range(length):
        cand = [(min(p - x.start, x.end - p), x.name) for x in placements if x.start <= p <= x.end]
        out.append(max(cand)[1] if cand else None)
    return out


def placements_for(entity, chain_of) -> list[Placement]:
    """Place every structure of an evedesign entity; `chain_of(structure)` returns its single chain."""
    wt = "".join(entity.rep)
    out = []
    for name, st in entity.structures.items():
        res_ids, letters = chain_residues(chain_of(st))
        out.append(place(name, res_ids, letters, wt))
    return out
