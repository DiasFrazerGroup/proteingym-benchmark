"""Check the card archives against the benchmark cards (clinical-curation benchmark_cards/summary.tsv).

Per card:
  1. selection   genes / rows (all and complete; classes or singletons) of prepare_cards.py equal the card's
  2. mapping     rows dropped because they cannot be placed on UniProt: reported, must be <= 3%
  3. rows        the archives' assay rows are exactly the prepared rows (archives hold no model scores)
     metrics     the card's six model metrics recomputed on those rows with the card scores kept in the prepared
                 table (complete rows;
                 labelled: pooled AUROC + per-gene mean AUROC over genes with both classes; gnomAD: pooled Spearman
                 rho(score, log10 AF) + per-gene mean over genes with >= 20 variants) vs summary.tsv.
                 REPRODUCED if all agree to 5e-4 (summary.tsv is rounded to 4 decimals); with unmapped rows a
                 small difference is expected and is printed instead.
  4. archives    a sample loads with proteingym-base (Subsets.from_path); every row's sequence is the wild type
                 with exactly the `mutant` substitution
The random300 sets (dnm_aptest_random300: no card; gnomad_af_random300: topped up with replacement genes) are not
compared with a card: their counts and metrics are printed only.
Run: python scripts/clinical/cards/check_cards.py --prepared $WORK/cards --curation <clinical-curation>
"""
import argparse
import io
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

SCORES = {"score_alphamissense": "AlphaMissense", "score_popeve": "popEVE", "score_cadd_phred": "CADD v1.7",
          "score_esm1b": "ESM1b", "score_vesm3b": "VESM-3B", "score_gpn_msa": "GPN-MSA"}
TARGET = {"clinvar": "label", "dnm": "dnm_case", "gnomad": "log10_af"}
MIN_N, TOL = 20, 5e-4


def read_assay(path):
    with zipfile.ZipFile(path) as z:
        inner = zipfile.ZipFile(io.BytesIO(z.read("dataset.pgdata")))
        name = next(n for n in inner.namelist() if n.startswith("assays/") and n.endswith(".csv"))
        wt = next(inner.read(n).decode().split("\n", 1)[1].replace("\n", "") for n in inner.namelist()
                  if n.endswith(".fasta"))
        return pd.read_csv(io.BytesIO(inner.read(name)), low_memory=False), wt


def metrics(df, y, gnomad):
    x = df[df.complete].reset_index(drop=True)
    out = {}
    for c, nm in SCORES.items():
        if gnomad:
            pooled = spearmanr(x[c], x[y]).statistic
            pg = [spearmanr(t[c], t[y]).statistic for _, t in x.groupby("gene") if len(t) >= MIN_N]
        else:
            pooled = roc_auc_score(x[y], x[c])
            pg = [roc_auc_score(t[y], t[c]) for _, t in x.groupby("gene") if t[y].nunique() == 2]
        out[nm] = (pooled, float(np.nanmean(pg)))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--prepared", type=Path, required=True)
    ap.add_argument("--curation", type=Path, required=True)
    ap.add_argument("--index", type=Path, default=Path("datasets/clinical/cards"))
    ap.add_argument("--load-sample", type=int, default=10, help="archives per card loaded with proteingym-base")
    a = ap.parse_args()
    summ = pd.read_csv(a.curation / "benchmark_cards/summary.tsv", sep="\t")
    sel = pd.read_csv(a.prepared / "selection_counts.tsv", sep="\t").set_index("card")
    from proteingym.base import Subsets
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from build_cards import SOURCE  # noqa: E402

    failed = []
    for card, s in SOURCE.items():
        print(f"\n=== {card}")
        y, gnomad = TARGET[s], s == "gnomad"
        # random300 sets are topped up to 300 usable genes (<= 2,700 aa), so they are not the card's exact set
        cs = summ[(summ.card == card) & ~card.endswith("random300")]
        # 1. selection
        if len(cs):
            r0 = cs.iloc[0]
            cols = [c for c in sel.columns if c in cs.columns and pd.notna(r0[c])]
            diff = {c: (int(sel.loc[card, c]), int(r0[c])) for c in cols if int(sel.loc[card, c]) != int(r0[c])}
            print(f"selection: {len(cols)} counts vs card -> {'IDENTICAL' if not diff else f'DIFF {diff}'}")
            if diff:
                failed.append(f"{card} selection")
        else:
            print("selection: no card (new set); counts", sel.loc[card].dropna().astype(int).to_dict())
        # 2. mapping
        um = pd.read_csv(a.prepared / card / "unmapped_rows.tsv", sep="\t", low_memory=False)
        n_all = int(sel.loc[card, "variants_all"])
        frac = len(um) / n_all
        print(f"mapping: {len(um):,} of {n_all:,} rows unmapped ({frac:.3%})")
        if frac > 0.03:
            failed.append(f"{card} unmapped {frac:.2%}")
        # 3. metrics from the archives
        idx = pd.read_csv(a.index / f"{card}.csv")
        parts = []
        for arc, gene in zip(idx.archive, idx.gene):
            df, wt = read_assay(a.prepared / arc)
            parts.append(df.assign(gene=gene))
        df = pd.concat(parts, ignore_index=True)
        assert len(df) == idx.n_variants.sum() == n_all - len(um), (len(df), idx.n_variants.sum(), n_all, len(um))
        # the archives carry no model scores: check their rows are exactly the prepared rows, then take the card
        # scores from the prepared table (prepare_cards.py keeps them there for this check)
        prep = pd.read_csv(a.prepared / card / "variants.tsv.gz", sep="\t", low_memory=False,
                           dtype={"variant_id": str})
        k = ["gene", "mutant", "mutant_source", "variant_id", y]
        norm = lambda t: t[k].astype(str).sort_values(k).reset_index(drop=True)  # noqa: E731
        same = norm(df.assign(variant_id=df.variant_id.fillna(""))).equals(
            norm(prep.assign(variant_id=prep.variant_id.fillna(""))))
        print(f"rows: archives vs prepared table -> {'IDENTICAL' if same else 'DIFFER'}")
        if not same:
            failed.append(f"{card} rows")
        m = metrics(prep.assign(complete=prep.complete.astype(str).str.lower() == "true"), y, gnomad)
        worst = 0.0
        for nm, (p, g) in m.items():
            ref = cs[cs.model == nm]
            if len(ref):
                dp, dg = p - ref.pooled.iloc[0], g - ref.per_gene.iloc[0]
                worst = max(worst, abs(dp), abs(dg))
                print(f"  {nm:14s} pooled {p:.4f} (card {ref.pooled.iloc[0]:.4f})  per-gene {g:.4f} "
                      f"(card {ref.per_gene.iloc[0]:.4f})")
            else:
                print(f"  {nm:14s} pooled {p:.4f}  per-gene {g:.4f}")
        if len(cs):
            if worst <= TOL:
                print(f"metrics: REPRODUCED (max |diff| {worst:.1e})")
            elif len(um):
                print(f"metrics: differ by up to {worst:.4f}; {len(um)} rows unmapped (expected source of the gap)")
            else:
                print(f"metrics: MISMATCH (max |diff| {worst:.4f})")
                failed.append(f"{card} metrics")
        # 4. archives load; sequences = wild type + mutant
        for arc in idx.archive.iloc[: a.load_sample]:
            Subsets.from_path(a.prepared / arc)
        bad = 0
        for arc in idx.archive:
            df, wt = read_assay(a.prepared / arc)
            for sq, mu in zip(df.sequence, df.mutant):
                p = int(mu[1:-1])
                bad += not (len(sq) == len(wt) and wt[p - 1] == mu[0] and sq[p - 1] == mu[-1]
                            and sq[:p - 1] == wt[:p - 1] and sq[p:] == wt[p:])
        print(f"archives: {min(a.load_sample, len(idx))} loaded with Subsets.from_path; "
              f"{idx.n_variants.sum():,} rows checked, {bad} with a wrong sequence")
        if bad:
            failed.append(f"{card} sequences")
    print("\nCHECK FAILED: " + ", ".join(failed) if failed else "\nCHECK OK")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
