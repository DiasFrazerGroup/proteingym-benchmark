"""Unit tests of the shared EVE-protocol MSA filter (synthetic alignments, no model run)."""
import numpy as np
from evedesign.sequence import Sequence, Sequences

from proteingym.models.structuredca_filtered.msa_filter import choose_threshold, filter_sequences

Q = "ACDEFGHIKL"


def _msa(rows, bits):
    s = Sequences(seqs=[Sequence(seq=r) for r in [Q] + rows], aligned=True, format="a3m")
    s.pg_bitscores = np.array([np.nan] + bits, dtype=float)
    return s


def test_row_filters():
    rows = ["ACDEFGHIKL", "ACDEF-----", "ACDE------", "XXXXXXXXXX", "ACDEFGHIKX", "ACdeDEFGHIKL", "ACD"]
    new, st = filter_sequences(_msa(rows, [10.0] * 7), len(Q), bits_thresholds=[0.0], min_depth_per_length=0)
    seqs = [x.seq for x in new.seqs]
    assert seqs[0] == Q
    # 50% residues kept, 40% dropped (fragment), X rows dropped (indeterminate), inserts removed, bad length dropped
    assert seqs[1:] == ["ACDEFGHIKL", "ACDEF-----", "ACDEFGHIKL"]
    assert (st["rows_bad_length"], st["rows_fragments"], st["rows_indeterminate"]) == (1, 1, 2)


def test_threshold_and_subsample():
    rows = ["ACDEFGHIKL"] * 40 + ["ACDEF-----"] * 40
    bits = [9.0] * 40 + [1.0] * 40       # bits/residue 0.9 and 0.1
    new, st = filter_sequences(_msa(rows, bits), 10, bits_thresholds=[0.0, 0.5])
    # at 0.0 half the rows cover only 5 columns -> L_cov 5/10 < 0.7; at 0.5 all 10 columns covered, N 41 < 10 L
    assert st["threshold"] == 0.5 and st["N"] == 41 and st["relax_rule"] == "low_depth"
    new, st = filter_sequences(_msa(rows, bits), 10, bits_thresholds=[0.0, 0.5], max_depth_relaxed=20, seed=1)
    assert st["relax_rule"] == "subsample" and st["subsampled"] and len(new.seqs) == 20 and new.seqs[0].seq == Q


def test_relaxation_order():
    L = 100
    assert choose_threshold(L, [5000, 2000, 1500], [60, 85, 100]) == (1, 0, "strict")
    assert choose_threshold(L, [150_000, 2000, 500], [75, 75, 100]) == (0, 1, "relaxed")
    assert choose_threshold(L, [5000, 500, 50], [60, 75, 100]) == (1, 2, "low_depth")
    assert choose_threshold(L, [5000, 500, 1], [60, 65, 100]) == (0, 3, "low_lcov")
    assert choose_threshold(L, [300_000, 250_000, 210_000], [85, 90, 95]) == (2, 4, "subsample")


def test_no_bitscores():
    s = Sequences(seqs=[Sequence(seq=r) for r in [Q, Q, "AC--------"]], aligned=True, format="a3m")
    new, st = filter_sequences(s, 10)
    assert st["relax_rule"] == "no_bitscores" and len(new.seqs) == 2


def test_nan_bitscores_are_no_bitscores():
    s = Sequences(seqs=[Sequence(seq=r) for r in [Q, Q, "AC--------"]], aligned=True, format="a3m")
    s.pg_bitscores = np.array([np.nan, np.nan, np.nan])   # archive headers seq_<i>
    new, st = filter_sequences(s, 10)
    assert st["relax_rule"] == "no_bitscores" and not st["bitscores"]
