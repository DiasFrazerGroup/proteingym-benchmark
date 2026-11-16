"""Export model-API scores ($WORK/models/<family>/<model>/scores.parquet) as per-protein prediction CSVs for the
Genebass evaluation (proteingym_clinical/genebass/scripts/16_eval_genebass_source.py reads predictions/<name>/<UNIPROT>.csv,
columns mutant,pred; pred = higher is fitter) and compare per variant with reference predictions.
Usage: python export_scores.py   (paths below)
"""
from pathlib import Path

import numpy as np
import pandas as pd

W = Path("/no_backup/diasfrazer/egalpern_pgb_clinical")
PRED = Path("/users/diasfrazer/egalpern/project/proteingym_clinical/genebass/predictions")
idx = pd.read_csv(W / "genebass_api/build/genebass_fig6b.csv")
entry2u = dict(zip(idx.entry_name, idx.uniprot_id))

runs = {"plm_windows/esm2": "esm2_window", "plm_windows/esmc": "esmc_window",
        "struct_windows/saprot": "saprot_window", "struct_windows/esmif1": "esm-if1_window"}
REFS = {"esm2_window": ["dark_energy__ESM2_650M_masked", "esm2_window"], "esmc_window": ["esmc_window"],
        "saprot_window": ["saprot_clinical"], "esm-if1_window": ["esm_if1_clinical"]}
for sub, model in runs.items():
    f = W / "models" / sub / "scores.parquet"
    if not f.exists():
        continue
    s = pd.read_parquet(f)
    s["uniprot_id"] = s.entry.map(entry2u)
    cols = [c for c in s.columns if c in ("score",) or c.startswith(("termini_", "every_window_"))]
    for c in cols:
        name = f"v2__{model}" + ("" if c == "score" else f"__{c}")
        d = PRED / name
        d.mkdir(exist_ok=True)
        for u, g in s.groupby("uniprot_id"):
            g[["mutant", c]].rename(columns={c: "pred"}).to_csv(d / f"{u}.csv", index=False)
        for ref in REFS[model]:
            r = pd.concat([pd.read_csv(x).assign(uniprot_id=x.stem) for x in (PRED / ref).glob("*.csv")])
            j = s[["uniprot_id", "mutant", "L", c]].merge(r, on=["uniprot_id", "mutant"]).dropna()
            for lab, m in [("L<=1022", j.L <= 1022), ("L>1022", j.L > 1022)]:
                if m.sum() > 2:
                    x, y = j.loc[m, c], j.loc[m, "pred"]
                    print(f"{name:42s} vs {ref:30s} {lab:8s} n={m.sum():5d} r={np.corrcoef(x, y)[0, 1]:.5f} "
                          f"max|diff|={np.abs(x - y).max():.3g}")
    print(f"exported {model}: {len(s)} variants, columns {cols}")
