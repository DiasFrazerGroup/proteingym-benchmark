"""Export *_filtered Genebass runs (run_genebass.sbatch) and write a per-protein status table (local env).

  python export_genebass.py <run_dir> <predictions_dir> <model folder> [<model folder> ...]
For each model folder (e.g. rsalor_filtered) and each protein of <run_dir>/proteins.tsv:
- predictions: <predictions_dir>/v2__<model folder>/<UNIPROT>.csv with mutant,pred (pred = score(variant) -
  score(WT), higher = fitter), NaN rows dropped (= not scored by the model);
- status: <predictions_dir>/v2__<model folder>/status.tsv: n_mutants, n_scored, n_unscored, status (ok / failed /
  missing), failure reason, time, peak RSS, filter stats.
Reuse: a run in <run_dir>/reuse/<model folder>__<UNIPROT>__genebass.{csv,json} (e.g. an all-singles validation run
with the same card) is used when there is no own run; it is restricted to the protein's Genebass mutants.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    run, preds = Path(sys.argv[1]), Path(sys.argv[2])
    proteins = pd.read_csv(run / "proteins.tsv", sep="\t")
    for folder in sys.argv[3:]:
        out = preds / f"v2__{folder}"
        out.mkdir(parents=True, exist_ok=True)
        rows = []
        for p in proteins.itertuples():
            u = p.uniprot_id
            tag = f"{folder}__{u}__genebass"
            rec = dict(uniprot=u, gene=p.gene, L=p.L, msa_rows=p.msa, n_mutants=np.nan)
            src = None
            for d in (run, run / "reuse"):
                if (d / f"{tag}.json").exists():
                    src = d
                    break
            if src is None:
                failed = run / f"{tag}.failed.json"
                if failed.exists():
                    j = json.loads(failed.read_text())
                    rec.update(status="failed", reason=j.get("error"), t_total_min=round(j["t_total_s"] / 60, 1),
                               maxrss_gb=max(j.get("maxrss_self_gb", 0), j.get("maxrss_children_gb", 0)))
                else:
                    reason = (run / f"{tag}.reason.txt").read_text().strip() if (run / f"{tag}.reason.txt").exists() \
                        else "no run output"
                    rec.update(status="missing", reason=reason)
                rows.append(rec)
                continue
            j = json.loads((src / f"{tag}.json").read_text())
            df = pd.read_csv(src / f"{tag}.csv")
            if src.name == "reuse":
                df = df[df.in_archive]
            f = j.get("msa_filter_stats") or {}
            rec.update(status="ok", reason="reused: " + j["tag"] if src.name == "reuse" else "",
                       n_mutants=len(df), n_scored=int(df.pred.notna().sum()), n_unscored=int(df.pred.isna().sum()),
                       t_build_min=round(j["t_build_s"] / 60, 1), t_total_min=round(j["t_total_s"] / 60, 1),
                       maxrss_gb=max(j["maxrss_self_gb"], j["maxrss_children_gb"]),
                       **{f"filter_{k}": f.get(k) for k in ("bitscore_source", "rows_after_row_filters", "threshold",
                                                            "N", "N_over_L", "Lcov_frac", "relax_rule", "quality")})
            df[df.pred.notna()][["mutant", "pred"]].to_csv(out / f"{u}.csv", index=False)
            rows.append(rec)
        st = pd.DataFrame(rows)
        st.to_csv(out / "status.tsv", sep="\t", index=False)
        ok = st[st.status == "ok"]
        print(f"{folder}: {len(ok)}/{len(st)} proteins ok, {int(ok.n_scored.sum())} scored, "
              f"{int(ok.n_unscored.sum())} unscored mutants; failed/missing: "
              f"{', '.join(st[st.status != 'ok'].gene)}")


if __name__ == "__main__":
    main()
