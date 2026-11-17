"""Summarise validate_msa_filter.py runs (local env: pandas, scipy). Tags are <model>__<UNIPROT>__<variant>.

  python compare_msa_filter.py <val_dir>            -> runs.tsv (one row per run), pairs.tsv (comparisons)
Comparisons per (model, protein): every variant vs `filtered` (Pearson, Spearman on all scored mutants and on the
archive's Genebass mutants, max abs diff), and for eve also vs the popEVE release `EVE` column (protein found by exact
sequence in the popEVE sequence index; read-only).
"""
import io
import json
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

POPEVE_ZIP = "/no_backup/diasfrazer/db/popeve/raw/popEVE_ukbb_20250312.zip"
POPEVE_SEQ = "/no_backup/diasfrazer/egalpern_pgb_clinical/layers/popeve_sequences.parquet"
MSAS = "/no_backup/diasfrazer/egalpern_pgb_clinical/genebass_msas"


def corr(x, y):
    ok = ~(np.isnan(x) | np.isnan(y))
    if ok.sum() < 3:
        return np.nan, np.nan, int(ok.sum()), np.nan
    return (pearsonr(x[ok], y[ok])[0], spearmanr(x[ok], y[ok])[0], int(ok.sum()),
            float(np.max(np.abs(x[ok] - y[ok]))))


def popeve_eve(uniprot):
    with open(f"{MSAS}/{uniprot}.a3m") as f:
        f.readline(); wt = f.readline().strip()
    idx = pd.read_parquet(POPEVE_SEQ)
    hit = idx[idx.sequence == wt]
    if hit.empty:
        return None, None
    npid = hit.np.iloc[0]
    z = zipfile.ZipFile(POPEVE_ZIP)
    name = next(n for n in z.namelist() if n.endswith(f"/{npid}.csv"))
    d = pd.read_csv(io.BytesIO(z.read(name)), usecols=["mutant", "EVE", "gap frequency"])
    return npid, d


def main():
    val = Path(sys.argv[1])
    runs = []
    for j in sorted(val.glob("*.json")):
        if j.name.endswith(".timeout.json"):
            r = json.loads(j.read_text()); r["status"] = "timeout"
        else:
            r = json.loads(j.read_text()); r["status"] = "ok"
            for k, v in (r.pop("msa_filter_stats", None) or {}).items():
                r[f"f_{k}"] = v
        r["model"], r["uniprot"], r["variant"] = r["tag"].split("__")
        runs.append(r)
    runs = pd.DataFrame(runs)
    runs.to_csv(val / "runs.tsv", sep="\t", index=False)
    pairs = []
    for (m, u), g in runs[runs.status == "ok"].groupby(["model", "uniprot"]):
        ref = val / f"{m}__{u}__filtered.csv"
        if not ref.exists():
            ref = val / f"{m}__{u}__off.csv"   # T1-only proteins
        if not ref.exists():
            continue
        f = pd.read_csv(ref).set_index("mutant")
        for v in g.variant:
            if f"{m}__{u}__{v}.csv" == ref.name:
                continue
            o = pd.read_csv(val / f"{m}__{u}__{v}.csv").set_index("mutant").reindex(f.index)
            for subset, mask in (("all", np.ones(len(f), bool)), ("archive", f.in_archive.values)):
                r, rho, n, mx = corr(f.pred.values[mask], o.pred.values[mask])
                pairs.append(dict(model=m, uniprot=u, a=ref.stem.split("__")[2], b=v, subset=subset, pearson=r, spearman=rho,
                                  n=n, max_abs_diff=mx, nan_a=int(np.isnan(f.pred.values[mask]).sum()),
                                  nan_b=int(np.isnan(o.pred.values[mask]).sum())))
        t1z, t1o = val / f"{m}__{u}__zoo.csv", val / f"{m}__{u}__off.csv"
        if t1z.exists() and t1o.exists():   # T1: filter disabled vs the zoo model
            z, o = pd.read_csv(t1z).set_index("mutant"), pd.read_csv(t1o).set_index("mutant")
            o = o.reindex(z.index)
            r, rho, n, mx = corr(z.pred.values, o.pred.values)
            pairs.append(dict(model=m, uniprot=u, a="off", b="zoo", subset="all", pearson=r, spearman=rho, n=n,
                              max_abs_diff=mx, nan_a=int(o.pred.isna().sum()), nan_b=int(z.pred.isna().sum())))
        if m == "eve":
            npid, d = popeve_eve(u)
            if d is None:
                continue
            d = d.set_index("mutant").reindex(f.index)
            for v in g.variant:
                o = pd.read_csv(val / f"{m}__{u}__{v}.csv").set_index("mutant").reindex(f.index)
                for subset, mask in (("all", np.ones(len(f), bool)), ("archive", f.in_archive.values)):
                    r, rho, n, _ = corr(o.pred.values[mask], d.EVE.values[mask])
                    pairs.append(dict(model=m, uniprot=u, a=v, b=f"popEVE_release_EVE({npid})", subset=subset,
                                      pearson=r, spearman=rho, n=n,
                                      nan_a=int(np.isnan(o.pred.values[mask]).sum()),
                                      nan_b=int(np.isnan(d.EVE.values[mask]).sum())))
    pairs = pd.DataFrame(pairs)
    pairs.to_csv(val / "pairs.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
    cols = [c for c in ["tag", "status", "L", "msa_rows", "f_rows_after_row_filters", "f_threshold", "f_N",
                        "f_N_over_L", "f_Lcov_frac", "f_relax_rule", "f_quality", "t_load_s", "t_build_s",
                        "t_score_s", "maxrss_self_gb", "maxrss_children_gb", "n_scored", "n_nan"] if c in runs]
    print(runs[cols].to_string(index=False))
    if len(pairs):
        print(pairs.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
