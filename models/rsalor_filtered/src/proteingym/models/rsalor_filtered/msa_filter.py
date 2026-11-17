"""EVE-protocol MSA filter for the alignment models (docs/clinical/models_long_proteins.md §3.6).

Shared by the `*_filtered` model folders; this file is copied verbatim into each folder (folders stay
self-contained). The base alignment is the dataset's MSA (ColabFold/MMseqs2 a3m); every MSA model gets the same
filtered alignment from the same card parameters. Benchmark rule: follow EVE (Frazer et al. 2021 preprint Methods;
EVE code `utils/data_utils.py`, commit 460d70e as vendored in evedesign) / popEVE as published; where they differ,
take the choice with the higher coverage (popEVE); where they are silent, take the simplest option (marked "ours").

Rows (EVE code, MSA_processing):
1. insertions (lowercase, '.') removed (EVE: a2m match states); rows whose match states do not have the query length
   are dropped (malformed);
2. fragments dropped: gap fraction > 1 - `min_seq_coverage` (EVE: threshold_sequence_frac_gaps = 0.5; paper: "align
   to at least 50% of the target");
3. rows with an indeterminate residue (anything but the 20 standard amino acids and '-') dropped (EVE:
   remove_sequences_with_indeterminate_AA_in_focus_cols; every column is a focus column here, see Columns).
   The query (row 0) is always kept.
Bit-score selection (paper: thresholds explored around 0.3 bits/residue; prefer L_cov >= 0.8 L and 10 L <= N <=
100,000, relaxed "sequentially down to L_cov >= 0.7 L and N <= 200,000", then "best available"): bits/residue =
MMseqs2 bit score (a3m header field 2) / L (deviation: MMseqs2 scale, not jackhmmer/HMMER, so `bits_thresholds` is
a sweep 0-1.5, the 0.3 reference does not transfer). N(t) = rows with bits/residue >= t (query included); L_cov(t) =
columns with >= `col_occupancy` residue occupancy among them. The first rule some threshold satisfies picks the
lowest such t (deepest alignment; ours: the paper does not say which qualifying threshold):
     step 0 "strict":    L_cov >= lcov_target  * L and min_depth_per_length * L <= N <= max_depth           (paper)
     step 1 "relaxed":   L_cov >= lcov_relaxed * L and min_depth_per_length * L <= N <= max_depth_relaxed   (paper)
     step 2 "low_depth": L_cov >= lcov_relaxed * L and 2 <= N <= max_depth_relaxed      (ours: depth floor dropped)
     step 3 "low_lcov":  2 <= N <= max_depth_relaxed                                    (ours: deepest alignment)
     step 4 "subsample": no threshold gives N <= max_depth_relaxed: highest t, then a uniform random subsample of
                         max_depth_relaxed rows (query first, `seed`, original order) (ours)
   Rows without a bit score are kept only at the lowest threshold. Without any bit scores (headers not available)
   the selection is skipped ("no_bitscores": all rows passing 1-3, subsampled as in step 4 if needed).
Columns: ALL columns are kept and modelled (popEVE choice: the popEVE release scores 100% of positions, including
positions with gap frequency up to 1.0; EVE code's default threshold_focus_cols_frac_gaps = 0.3 would leave ~24% of
positions unscored). L_cov is a selection criterion only. Column handling inside a model stays the zoo model's own
(EVE: zoo card threshold_focus_cols_frac_gaps = 1.0 = popEVE-like; EVmutation: its own > 50% gap exclusion).

Bit scores: proteingym-base keeps the raw a3m bytes (`MSA.file`) but evedesign's `msa_to_sequences` drops the
headers (`Sequence(id=None)`). `install_header_hook()` (called from the model package's __init__, i.e. before the
runner converts the dataset) wraps `evedesign.proteingym.msa_to_sequences` so the returned `Sequences` carries the
bit scores as an attribute (`pg_bitscores`, float array aligned with `.seqs`) and the query header
(`pg_query_name`). Sequence ids stay None, so nothing changes for the wrapped models.

Unscored variants: `nan_unscored()` sets the score to NaN for an instance that mutates a position outside the model's
`positions()` (EVmutation: gap-excluded columns; EVE: non-focus columns) or hits a NaN cell of the model's
per-position matrix `_pssm` (RSALOR: no RSA; GEMME: NA from pred.R). The upstream wrappers score these as WT (0).
"""
from __future__ import annotations

import copy
import json
import os
import re
import resource
import sys
import time
from typing import Any

import numpy as np

STANDARD_AA = "ACDEFGHIKLMNPQRSTVWY"
DEFAULT_BITS_THRESHOLDS = [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                           1.25, 1.5]
FILTER_DEFAULTS: dict[str, Any] = dict(
    msa_filter=True,
    min_seq_coverage=0.5,
    col_occupancy=0.7,
    lcov_target=0.8,
    lcov_relaxed=0.7,
    min_depth_per_length=10.0,
    max_depth=100_000,
    max_depth_relaxed=200_000,
    seed=42,
    bits_thresholds=DEFAULT_BITS_THRESHOLDS,
    a3m_headers_dir=None,
    stats_dir="/opt/program/output",
)

_DEL_INSERTS = str.maketrans("", "", "abcdefghijklmnopqrstuvwxyz.")
_DEL_STANDARD_GAP = str.maketrans("", "", STANDARD_AA + "-")
_HEADER_RE = re.compile(rb"^>([^\n]*)", re.MULTILINE)


# ---------------------------------------------------------------------------------------------------------------
# bit scores from the a3m headers
# ---------------------------------------------------------------------------------------------------------------
def bitscores_from_a3m(raw: bytes) -> tuple[np.ndarray, str]:
    """Bit score of every record in file order (NaN when the header has no tab-separated score) + query name."""
    out, query = [], None
    for m in _HEADER_RE.finditer(raw):
        h = m.group(1)
        if query is None:
            query = h.decode(errors="replace").split()[0] if h.strip() else ""
        f = h.split(b"\t", 2)
        try:
            out.append(float(f[1]) if len(f) > 1 else np.nan)
        except ValueError:
            out.append(np.nan)
    return np.asarray(out, dtype=float), (query or "")


def install_header_hook() -> None:
    """Wrap evedesign.proteingym.msa_to_sequences so the Sequences keep the a3m bit scores (idempotent)."""
    import evedesign.proteingym as pg

    if getattr(pg.msa_to_sequences, "_pg_msa_filter_hook", False):
        return
    original = pg.msa_to_sequences

    def msa_to_sequences(dataset):
        seqs = original(dataset)
        if seqs is None or not dataset.msas or getattr(dataset.msas[0], "file", None) is None:
            return seqs
        try:
            bits, query = bitscores_from_a3m(dataset.msas[0].file.getvalue())
        except Exception as e:  # never break the conversion because of the hook
            print(f"msa_filter: could not read a3m headers ({type(e).__name__}: {e})", flush=True)
            return seqs
        if len(bits) == len(seqs.seqs):
            seqs.pg_bitscores = bits
            seqs.pg_query_name = query
        else:
            print(f"msa_filter: {len(bits)} a3m headers for {len(seqs.seqs)} sequences; bit scores not attached",
                  flush=True)
        return seqs

    msa_to_sequences._pg_msa_filter_hook = True
    msa_to_sequences.__wrapped__ = original
    pg.msa_to_sequences = msa_to_sequences


# ---------------------------------------------------------------------------------------------------------------
# selection
# ---------------------------------------------------------------------------------------------------------------
def choose_threshold(L: int, N, lcov, *, lcov_target=0.8, lcov_relaxed=0.7, min_depth_per_length=10.0,
                     max_depth=100_000, max_depth_relaxed=200_000):
    """(index of the chosen threshold, step, rule name); see module docstring."""
    N, lcov = np.asarray(N), np.asarray(lcov)
    deep = N >= min_depth_per_length * L
    rules = [
        ("strict", (lcov >= lcov_target * L) & deep & (N <= max_depth)),
        ("relaxed", (lcov >= lcov_relaxed * L) & deep & (N <= max_depth_relaxed)),
        ("low_depth", (lcov >= lcov_relaxed * L) & (N >= 2) & (N <= max_depth_relaxed)),
        ("low_lcov", (N >= 2) & (N <= max_depth_relaxed)),
    ]
    for step, (name, mask) in enumerate(rules):
        idx = np.nonzero(mask)[0]
        if len(idx):
            return int(idx[0]), step, name
    return len(N) - 1, 4, "subsample"


def _occupied(m: str) -> np.ndarray:
    return np.frombuffer(m.encode("ascii", "replace"), dtype=np.uint8) != 45   # '-'


def filter_sequences(sequences, L: int, *, min_seq_coverage=0.5, col_occupancy=0.7, lcov_target=0.8,
                     lcov_relaxed=0.7, min_depth_per_length=10.0, max_depth=100_000, max_depth_relaxed=200_000,
                     seed=42, bits_thresholds=None, **_):
    """Apply the filter to an evedesign Sequences (a3m, query first). Returns (new Sequences, stats)."""
    from evedesign.sequence import Sequence, Sequences

    t0 = time.time()
    thresholds = np.asarray(DEFAULT_BITS_THRESHOLDS if bits_thresholds is None else bits_thresholds, dtype=float)
    if len(thresholds) == 0 or np.any(np.diff(thresholds) <= 0):
        raise ValueError("bits_thresholds must be a non-empty, strictly increasing list")
    seqs = sequences.seqs
    bits = getattr(sequences, "pg_bitscores", None)
    have_bits = bits is not None and len(bits) == len(seqs) and bool(np.isfinite(bits[1:]).any())
    weights = sequences.weights

    query = seqs[0].seq.translate(_DEL_INSERTS)
    if len(query) != L:
        raise ValueError(f"MSA query has {len(query)} match states, target length is {L}")
    min_res = min_seq_coverage * L
    T = len(thresholds) if have_bits else 1
    n_bucket = np.zeros(T, dtype=np.int64)
    occ = np.zeros((T, L), dtype=np.int64)
    kept_idx, kept_seq, kept_bucket = [0], [query], [T - 1]   # query: in every bucket's cumulative set
    n_bad = n_frag = n_indet = 0
    for i in range(1, len(seqs)):
        m = seqs[i].seq.translate(_DEL_INSERTS)
        if len(m) != L:
            n_bad += 1
            continue
        if L - m.count("-") < min_res:           # EVE: gap fraction > 0.5 = fragment
            n_frag += 1
            continue
        if m.translate(_DEL_STANDARD_GAP):        # EVE: indeterminate residue (X, B, Z, J, U, O, ...)
            n_indet += 1
            continue
        if have_bits:
            b = bits[i] / L
            k = 0 if np.isnan(b) else int(np.searchsorted(thresholds, b, side="right")) - 1
            if k < 0:   # below the lowest threshold
                continue
        else:
            k = 0
        kept_idx.append(i); kept_seq.append(m); kept_bucket.append(k)
        n_bucket[k] += 1
        occ[k] += _occupied(m)
    n_rows = len(kept_idx)   # incl. query
    N = np.cumsum(n_bucket[::-1])[::-1] + 1                     # rows with bucket >= k, + query
    occ = np.cumsum(occ[::-1], axis=0)[::-1] + _occupied(query)
    lcov = (occ >= col_occupancy * N[:, None]).sum(axis=1)

    if have_bits:
        k, step, step_name = choose_threshold(L, N, lcov, lcov_target=lcov_target, lcov_relaxed=lcov_relaxed,
                                              min_depth_per_length=min_depth_per_length, max_depth=max_depth,
                                              max_depth_relaxed=max_depth_relaxed)
    else:
        k, step, step_name = 0, -1, "no_bitscores"
    sel = np.nonzero(np.asarray(kept_bucket) >= k)[0]           # positions in kept_*; sel[0] == 0 (query)
    subsampled = False
    if max_depth_relaxed is not None and len(sel) > max_depth_relaxed:
        rng = np.random.default_rng(seed)
        rest = np.sort(rng.choice(sel[1:], size=max_depth_relaxed - 1, replace=False))
        sel = np.concatenate([[0], rest])
        subsampled = True

    new = Sequences(
        seqs=[Sequence(seq=kept_seq[j], id=None, type=seqs[kept_idx[j]].type_) for j in sel],
        aligned=True,
        type=sequences.type_,
        weights=None if weights is None else [weights[kept_idx[j]] for j in sel],
        format="a3m",
    )
    n_final = len(sel)
    lcov_final = int(lcov[k]) if not subsampled else int(
        (np.sum([_occupied(kept_seq[j]) for j in sel], axis=0) >= col_occupancy * n_final).sum())
    flags = [f for f, bad in (("low_depth", n_final < min_depth_per_length * L),
                              ("low_lcov", lcov_final < lcov_relaxed * L)) if bad]
    stats = dict(
        query=getattr(sequences, "pg_query_name", None), L=L, rows_in=len(seqs), rows_bad_length=n_bad,
        rows_fragments=n_frag, rows_indeterminate=n_indet, rows_after_row_filters=n_rows, bitscores=bool(have_bits),
        bitscore_source=getattr(sequences, "pg_bitscore_source", None) if have_bits else None,
        threshold=float(thresholds[k]) if have_bits else None,
        N_at_threshold=int(N[k]), N=n_final, N_over_L=round(n_final / L, 3), Lcov=lcov_final,
        Lcov_frac=round(lcov_final / L, 4), relax_step=step, relax_rule=step_name, subsampled=subsampled,
        quality="+".join(flags) or "ok",
        sweep=[dict(thr=float(t), N=int(n), Lcov=int(c)) for t, n, c in zip(thresholds[:T], N, lcov)],
        t_filter_s=round(time.time() - t0, 2),
    )
    return new, stats


def bitscores_from_dir(sequences, wt: str, a3m_dir: str):
    """Bit scores from the original ColabFold a3m in `a3m_dir` whose query equals the wild type, verified row by row.

    Needed because proteingym-base's MSA.dump writes `seq_<i>` headers, so dataset archives carry no bit scores.
    Returns (bits, file name) or (None, reason).
    """
    import glob

    path = None
    for f in sorted(glob.glob(os.path.join(a3m_dir, "*.a3m"))):
        with open(f) as fh:
            fh.readline()
            if fh.readline().strip() == wt:
                path = f
                break
    if path is None:
        return None, "no a3m with this query in a3m_headers_dir"
    seqs = sequences.seqs
    bits, n, ok = [], 0, True
    with open(path) as fh:
        hdr = None
        for line in fh:
            if line.startswith(">"):
                hdr = line[1:].rstrip("\n")
                continue
            if n >= len(seqs) or seqs[n].seq.translate(_DEL_INSERTS) != line.strip().translate(_DEL_INSERTS):
                ok = False
                break
            f = hdr.split("\t", 2)
            try:
                bits.append(float(f[1]) if len(f) > 1 else np.nan)
            except ValueError:
                bits.append(np.nan)
            n += 1
    if not ok or n != len(seqs):
        return None, f"{os.path.basename(path)}: rows differ from the dataset MSA (row {n})"
    return np.asarray(bits, dtype=float), os.path.basename(path)


def filter_system(system, params: dict):
    """Filtered copy of a single-entity system (entity 0's MSA replaced) + stats. Input system is unchanged."""
    from evedesign.system import System

    target = system[0]
    L = len(target.rep)
    params = dict(params)
    a3m_dir = params.pop("a3m_headers_dir", None)
    sq = target.sequences
    bits = getattr(sq, "pg_bitscores", None)
    if bits is not None and len(bits) > 1 and np.isfinite(bits[1:]).any():
        sq.pg_bitscore_source = "dataset a3m headers"
    elif a3m_dir:
        bits, src = bitscores_from_dir(sq, "".join(target.rep), a3m_dir)
        print(f"msa_filter: bit scores from a3m_headers_dir: {src}", flush=True)
        if bits is not None:
            sq.pg_bitscores, sq.pg_bitscore_source = bits, f"a3m_headers_dir/{src}"
    new_seqs, stats = filter_sequences(sq, L, **params)
    entity = copy.copy(target)
    entity.sequences = new_seqs
    return System([entity] + [system[i] for i in range(1, len(system))]), stats


def write_stats(stats: dict, stats_dir: str | None, model_name: str) -> None:
    short = {k: v for k, v in stats.items() if k != "sweep"}
    print(f"msa_filter [{model_name}]: " + " ".join(f"{k}={v}" for k, v in short.items()), flush=True)
    if not stats_dir or not os.path.isdir(stats_dir) or not os.access(stats_dir, os.W_OK):
        return
    with open(os.path.join(stats_dir, "msa_filter_stats.jsonl"), "a") as fh:
        fh.write(json.dumps(dict(model=model_name, **stats)) + "\n")


# ---------------------------------------------------------------------------------------------------------------
# unscored variants -> NaN
# ---------------------------------------------------------------------------------------------------------------
def _as_bytes(rep) -> np.ndarray:
    return np.frombuffer("".join(str(c) for c in rep).encode("ascii", "replace"), dtype=np.uint8)


def nan_unscored(model, instances, scored):
    """Set NaN scores for instances that mutate an unmodelled position or a NaN `_pssm` cell (in place)."""
    target = model.system[0]
    wt = _as_bytes(target.rep)
    first = target.first_index
    try:
        modelled = np.zeros(len(wt), dtype=bool)
        for _, p in model.positions():
            modelled[int(p) - first] = True
    except Exception:
        modelled = None
    pssm = getattr(model, "_pssm", None)
    aa_idx = {a: i for i, a in enumerate(STANDARD_AA)}
    n = 0
    for inst, out in zip(instances, scored):
        rep = _as_bytes(inst[0].rep)
        if len(rep) != len(wt):
            continue
        for p in np.nonzero(rep != wt)[0]:
            bad = modelled is not None and not modelled[p]
            aa = chr(rep[p])
            if not bad and pssm is not None and aa in aa_idx:
                bad = bool(np.isnan(pssm[p, aa_idx[aa]]))
            if bad:
                out.score = float("nan")
                n += 1
                break
    if n:
        print(f"msa_filter: {n} of {len(scored)} instances unscored (NaN)", flush=True)
    return scored


# ---------------------------------------------------------------------------------------------------------------
# mixin for the filtered model classes
# ---------------------------------------------------------------------------------------------------------------
class MSAFilterMixin:
    """Put before the zoo class: class XFiltered(MSAFilterMixin, X). Filter keys are popped from the kwargs."""

    def _pop_filter_params(self, kwargs: dict) -> dict:
        p = {k: kwargs.pop(k, v) for k, v in FILTER_DEFAULTS.items()}
        self.msa_filter_params = p
        self.msa_filter_stats = None
        return kwargs

    def build(self, system, data=None, status_callback=None):
        p = dict(self.msa_filter_params)
        if p.pop("msa_filter"):
            system, stats = filter_system(system, p)
            stats["maxrss_gb"] = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20, 2)
            self.msa_filter_stats = stats
            write_stats(stats, p.get("stats_dir"), type(self).__name__)
        return super().build(system, data=data, status_callback=status_callback)

    def score(self, instances, status_callback=None):
        scored = super().score(instances, status_callback=status_callback)
        return nan_unscored(self, instances, scored)


if __name__ == "__main__":   # quick check on a local a3m: python msa_filter.py <file.a3m>
    from evedesign.sequence import Sequences

    s = Sequences.from_file(sys.argv[1], format="a3m")
    with open(sys.argv[1], "rb") as fh:
        s.pg_bitscores, s.pg_query_name = bitscores_from_a3m(fh.read())
    L = len(s.seqs[0].seq.translate(_DEL_INSERTS))
    _, st = filter_sequences(s, L)
    print(json.dumps({k: v for k, v in st.items() if k != "sweep"}, indent=1))
