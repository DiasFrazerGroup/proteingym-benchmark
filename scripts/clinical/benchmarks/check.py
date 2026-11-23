"""Check the ClinVar, de novo and gnomAD benchmark archives, and report the reference models' results on them.

Per benchmark:
  1. mapping    variants dropped because they cannot be placed on UniProt: reported, must be <= 3%
  2. rows       the archives' assay rows are exactly the prepared variants
  3. archives   a sample loads with proteingym-base (Subsets.from_path); every variant's sequence is the wild type
                with exactly the `mutant` substitution
  4. results    six reference models (clinical-curation scores, kept in the prepared table, never archived) on the
                variants every model scores. Headline = pooled metric over all variants: AUROC (clinvar, denovo),
                Spearman rho(score, log10 AF) (gnomad, lower = better). Also the per-gene mean (genes with both
                classes; gnomAD: genes with >= 20 variants). De novo also reports a gene-only baseline: whether the
                variant's gene has another case variant in the ap_test_set (gene identity alone separates many cases
                from controls there).
ClinVar is also compared with clinical-curation's reference numbers (benchmark_cards/summary.tsv, clinvar_v0):
selection counts identical, and metrics REPRODUCED to 5e-4 or a gap explained by unmapped variants.
Run: bash scripts/clinical/benchmarks/submit.sh check
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
TARGET = {"clinvar": "label", "denovo": "dnm_case", "gnomad": "log10_af"}
MIN_N, TOL, MAX_UNMAPPED = 20, 5e-4, 0.03


def read_assay(path):
    with zipfile.ZipFile(path) as z:
        inner = zipfile.ZipFile(io.BytesIO(z.read("dataset.pgdata")))
        name = next(n for n in inner.namelist() if n.startswith("assays/") and n.endswith(".csv"))
        wt = next(inner.read(n).decode().split("\n", 1)[1].replace("\n", "") for n in inner.namelist()
                  if n.endswith(".fasta"))
        return pd.read_csv(io.BytesIO(inner.read(name)), low_memory=False), wt


def metric(x, y, c, gnomad):
    if gnomad:
        pooled = spearmanr(x[c], x[y]).statistic
        pg = [spearmanr(t[c], t[y]).statistic for _, t in x.groupby("gene") if len(t) >= MIN_N]
    else:
        pooled = roc_auc_score(x[y], x[c])
        pg = [roc_auc_score(t[y], t[c]) for _, t in x.groupby("gene") if t[y].nunique() == 2]
    return pooled, float(np.nanmean(pg)) if pg else np.nan, len(pg)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--prepared", type=Path, required=True)
    ap.add_argument("--curation", type=Path, required=True)
    ap.add_argument("--index", type=Path, default=Path("datasets/clinical"))
    ap.add_argument("--load-sample", type=int, default=10, help="archives per benchmark loaded with proteingym-base")
    ap.add_argument("--results", type=Path, default=None, help="write the results table here (TSV)")
    a = ap.parse_args()
    from proteingym.base import Subsets
    sel = pd.read_csv(a.prepared / "selection_counts.tsv", sep="\t").set_index("benchmark")
    ref = pd.read_csv(a.curation / "benchmark_cards/summary.tsv", sep="\t")
    ref = ref[ref.card == "clinvar_v0"]
    full_dnm = pd.read_csv(a.curation / "data/processed/benchmarks/dnm_aptest_scores.tsv", sep="\t",
                           low_memory=False)
    dnm_cases = full_dnm.groupby("gene").dnm_case.sum()

    failed, results = [], []
    for bench, y in TARGET.items():
        print(f"\n=== {bench}")
        gnomad = bench == "gnomad"
        # 1. mapping
        um = pd.read_csv(a.prepared / bench / "unmapped_rows.tsv", sep="\t", low_memory=False)
        n_all = int(sel.loc[bench, "variants_all"])
        print(f"mapping: {len(um):,} of {n_all:,} variants unmapped ({len(um) / n_all:.2%})")
        if len(um) / n_all > MAX_UNMAPPED:
            failed.append(f"{bench} unmapped")
        # 2. rows: archives == prepared variants
        idx = pd.read_csv(a.index / bench / f"{bench}.csv")
        parts, bad = [], 0
        for arc, gene in zip(idx.archive, idx.gene):
            df, wt = read_assay(a.prepared / arc)
            parts.append(df.assign(gene=gene))
            for sq, mu in zip(df.sequence, df.mutant):  # 3. sequences
                p = int(mu[1:-1])
                bad += not (len(sq) == len(wt) and wt[p - 1] == mu[0] and sq[p - 1] == mu[-1]
                            and sq[:p - 1] == wt[:p - 1] and sq[p:] == wt[p:])
        df = pd.concat(parts, ignore_index=True)
        prep = pd.read_csv(a.prepared / bench / "variants.tsv.gz", sep="\t", low_memory=False,
                           dtype={"variant_id": str})
        k = ["gene", "mutant", "mutant_source", "variant_id", y]
        norm = lambda t: t[k].fillna("").astype(str).sort_values(k).reset_index(drop=True)  # noqa: E731
        same = len(df) == n_all - len(um) and norm(df).equals(norm(prep))
        print(f"rows: {len(idx)} archives, {len(df):,} variants; archives vs prepared -> "
              f"{'IDENTICAL' if same else 'DIFFER'}")
        for arc in idx.archive.iloc[: a.load_sample]:
            Subsets.from_path(a.prepared / arc)
        print(f"archives: {min(a.load_sample, len(idx))} loaded with Subsets.from_path; {bad} wrong sequences; "
              f"without structure {(idx.n_structures == 0).sum()}, without MSA {idx.msa_n_sequences.isna().sum()}")
        if not same:
            failed.append(f"{bench} rows")
        if bad:
            failed.append(f"{bench} sequences")
        # 4. results on the variants every reference model scores
        x = prep[prep.complete.astype(str).str.lower() == "true"].reset_index(drop=True)
        print(f"results on {len(x):,} variants, {x.gene.nunique()} genes "
              f"({'pooled Spearman, lower = better' if gnomad else 'pooled AUROC'} | per-gene mean):")
        for c, nm in SCORES.items():
            p, g, ng = metric(x, y, c, gnomad)
            results.append(dict(benchmark=bench, model=nm, pooled=p, per_gene=g, n_genes_per_gene=ng,
                                n_variants=len(x), n_genes=x.gene.nunique()))
            print(f"  {nm:14s} {p:7.4f} | {g:7.4f} ({ng} genes)")
        if bench == "denovo":
            other = x.gene.map(dnm_cases).fillna(0) - x[y]
            x["gene_only"] = (other > 0).astype(float)
            p = roc_auc_score(x[y], x.gene_only)
            results.append(dict(benchmark=bench, model="gene-only (other case in gene)", pooled=p, per_gene=np.nan,
                                n_genes_per_gene=0, n_variants=len(x), n_genes=x.gene.nunique()))
            print(f"  {'gene-only':14s} {p:7.4f}   (does the gene have another case variant in the ap_test_set)")
        if bench == "clinvar":
            r0 = ref.iloc[0]
            cols = [c for c in sel.columns if c in ref.columns and pd.notna(r0[c]) and pd.notna(sel.loc[bench, c])]
            diff = {c: (int(sel.loc[bench, c]), int(r0[c])) for c in cols if int(sel.loc[bench, c]) != int(r0[c])}
            worst = max(max(abs(r["pooled"] - ref[ref.model == r["model"]].pooled.iloc[0]),
                            abs(r["per_gene"] - ref[ref.model == r["model"]].per_gene.iloc[0]))
                        for r in results if r["benchmark"] == "clinvar")
            print(f"vs clinical-curation clinvar_v0: selection {'IDENTICAL' if not diff else f'DIFF {diff}'}; "
                  f"metrics {'REPRODUCED' if worst <= TOL else 'differ'} (max |diff| {worst:.1e})")
            if diff or (worst > TOL and not len(um)):
                failed.append("clinvar vs reference")
    if a.results:
        pd.DataFrame(results).to_csv(a.results, sep="\t", index=False, float_format="%.4f")
    print("\nCHECK FAILED: " + ", ".join(failed) if failed else "\nCHECK OK")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
