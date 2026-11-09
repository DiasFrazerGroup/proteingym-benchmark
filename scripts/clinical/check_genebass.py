"""Check the Genebass Fig. 6b archives: read them back with proteingym-base, join the VESM paper's own scores and
reproduce the paper's printed Fig. 6b numbers (DJZN25-Fig6.ipynb: mean -log10 P overall, cell 26; per biomarker
group, cell 31).

Per archive: Pearson r(beta, score) over the assay rows and its two-sided P, as the paper does (score = paper model
score, higher = more damaging; the sign of r does not matter for -log10 P). Paper scores are taken per row from the betas table (the paper's input), joined on
(variant_id, hgvsp). Rows are read as stored (no averaging of
the 7 duplicate protein changes, which the paper counts twice); the effect of Dataset.to_df() averaging is reported.

Run: python scripts/clinical/check_genebass.py --fig6 <Fig6_data> [--out datasets/clinical/genebass]
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from proteingym.base import Subsets
from scipy.stats import pearsonr

PAPER_ALL = {"VESM++": 2.47, "VESM (3B)": 2.42, "VESM (650M)": 2.40, "VESM (150M)": 2.39, "AlphaMissense": 2.37,
             "VESM3": 2.35, "VESM (35M)": 2.18, "ESM1b": 2.17, "ESM2 (650M)": 2.10, "ESM2 (3B)": 2.05,
             "ESM2 (15B)": 1.99, "ESM3_struct": 1.78}
GROUPS = ["Lipid Metabolism", "Liver Function", "Renal Function", "Other"]
PAPER_GROUP = {  # cell 31 output
    "VESM++": [3.118714, 2.148785, 1.638412, 2.278652], "VESM3": [2.834545, 2.102417, 1.519999, 2.405236],
    "VESM (3B)": [3.150827, 2.057727, 1.699880, 2.011209], "VESM (150M)": [3.068647, 2.047556, 1.698348, 2.071367],
    "VESM (650M)": [3.086480, 2.047345, 1.697052, 2.090657], "AlphaMissense": [3.143235, 1.861097, 1.531538, 2.285358],
    "VESM (35M)": [2.839852, 1.824261, 1.662825, 1.791446], "ESM2 (650M)": [2.645675, 1.869973, 1.376365, 1.851277],
    "ESM1b": [2.872824, 1.835658, 1.336772, 1.838134], "ESM2 (3B)": [2.809681, 1.602843, 1.463730, 1.634046],
    "ESM3_struct": [1.971691, 1.596456, 1.302508, 2.083540], "ESM2 (15B)": [2.785609, 1.516043, 1.441094, 1.554455]}
MODELS = list(PAPER_ALL)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fig6", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=Path("datasets/clinical/genebass"))
    a = ap.parse_args()
    idx = pd.read_csv(a.out / "genebass_fig6b.csv")
    bet = pd.read_csv(a.fig6 / "GeneBass_BiochemAssays_top_hits_all_variants_betas.csv.gz")
    hits = pd.read_csv(a.fig6 / "GeneBass_BiochemAssays_top_hits.csv")

    # independent path: the notebook's own code on the raw table (cells 5, 23-25)
    b = bet[bet.gene_pheno.isin(hits[hits.pval_mis <= 1e-10].gene_pheno) & (bet.pval < 0.05)]
    k = b.groupby("gene_pheno").hgvsp.count()
    b = b[b.gene_pheno.isin(k[k > 25].index)]
    nb = pd.DataFrame([dict(group=x.biomarker_group.iloc[0], model=m, logp=-np.log10(pearsonr(x.beta, x[m])[1]))
                       for _, x in b.groupby(["gene", "phenotype"]) for m in MODELS])
    nb = nb.groupby("model").logp.mean().to_frame("all").join(nb.pivot_table(index="model", columns="group",
                                                                             values="logp")[GROUPS])

    rows, n_avg, missing = [], [], []
    for r in idx.itertuples():
        ds = Subsets.from_path(a.out / r.archive).dataset
        assert len(ds.assays) == 1 and ds.name == r.dataset
        assay = ds.assays[0]
        cols = [f.name for f in assay.fields]
        df = pd.DataFrame([[x.value if hasattr(x, "value") else x for x in rec] for rec in assay.records])
        df.columns = [assay.sequence_feature_name] + cols if df.shape[1] == len(cols) + 1 else cols
        ref = bet.loc[bet.gene_pheno == r.gene_pheno, ["variant_id", "hgvsp"] + MODELS]
        df = df.merge(ref, on=["variant_id", "hgvsp"], how="left", validate="one_to_one")
        assert len(df) == r.n_variants
        if df[MODELS].isna().any().any():                     # rows the paper does not have (e.g. API-built set)
            missing.append((r.dataset, int(df[MODELS[0]].isna().sum())))
            df = df.dropna(subset=MODELS)
        n_avg.append(len(ds.to_df()))
        for m in MODELS:
            rr, p = pearsonr(df.beta.astype(float), df[m].astype(float))
            rows.append(dict(dataset=r.dataset, group=r.biomarker_group, model=m, r=rr, logp=-np.log10(p)))
    c = pd.DataFrame(rows)
    res = c.groupby("model").logp.mean().to_frame("all").join(c.pivot_table(index="model", columns="group",
                                                                            values="logp")[GROUPS])
    ref = pd.DataFrame(PAPER_GROUP, index=GROUPS).T.join(pd.Series(PAPER_ALL, name="all"))
    diff = (res - ref[res.columns]).abs()
    dnb = (res - nb).abs().max().max()
    pd.set_option("display.width", 200)
    print(res.loc[MODELS].round(3).to_string())
    print(f"\nmax |ours - paper|: overall {diff['all'].max():.4f} (paper printed 2 decimals), "
          f"per group {diff[GROUPS].max().max():.2e}")
    print(f"max |archives - notebook code on raw table|: {dnb:.2e}")
    if missing:
        print(f"rows without paper scores (dropped for this check): {sum(n for _, n in missing)} in "
              f"{len(missing)} archives: {missing[:10]}")
    print(f"archives {len(idx)}, rows {idx.n_variants.sum()}, rows after Dataset.to_df() averaging {sum(n_avg)}")
    ok = diff["all"].max() <= 0.005 and diff[GROUPS].max().max() < 1e-5 and dnb < 1e-9
    print("REPRODUCED" if ok else "MISMATCH")


if __name__ == "__main__":
    main()
