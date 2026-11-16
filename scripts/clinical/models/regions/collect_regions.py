"""Map region predictions back to UniProt mutants; T5 summary; export for the Genebass benchmark.

Reads the region layout (make_regions.py: <ARCH>/<UNIPROT>/W<w>/{regions.json, variants.csv}) and the region
scores (run_regions.py: <PREDS>/<model>/<tag>/<UNIPROT>/W<w>/{<region>.csv, status.jsonl}). Every UniProt mutant
gets the prediction of its region; NaN if the region failed or the model did not score it (the `unscored` flag of
run_regions.py: evmutation gap-excluded columns, eve non-focus columns, NaN cells of rsalor/gemme). For
evmutation the flag is cross-checked against the gap rule recomputed from the region MSA (regions.json).

  collect_regions.py summary  --models rsalor evmutation --tag t5          -> T5 table (tsv + stdout)
  collect_regions.py compare  --model rsalor --tag t5 --uniprot Q8NBP7 --windows 250 500 10000
  collect_regions.py export   --model rsalor --tag full --window 500 [--dest .../predictions/v2__rsalor]
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

WORK = Path("/no_backup/diasfrazer/egalpern_pgb_clinical")
ARCH = WORK / "models/regions/archives"
PREDS = WORK / "models/regions/preds"
PRED_DEST = Path("/users/diasfrazer/egalpern/project/proteingym_clinical/genebass/predictions")


def mapped(model: str, tag: str, u: str, w: int) -> tuple[pd.DataFrame, list[dict]]:
    """One row per UniProt mutant of the protein: region, pred (NaN if unscored/failed), reason."""
    lay = ARCH / u / f"W{w}"
    reg = json.loads((lay / "regions.json").read_text())
    var = pd.read_csv(lay / "variants.csv")
    pd_dir = PREDS / model / tag / u / f"W{w}"
    status = {}
    if (pd_dir / "status.jsonl").exists():
        for line in (pd_dir / "status.jsonl").read_text().splitlines():
            r = json.loads(line)
            status[r["region"]] = r            # last attempt wins
    parts = []
    for r in reg["regions"]:
        rid = f"{u}_{r['region']}"
        f = pd_dir / f"{rid}.csv"
        sub = var[var.region == r["region"]][["mutant", "region", "region_mutant"]].copy()
        if f.exists():
            p = pd.read_csv(f, keep_default_na=False)   # empty `unscored` = scored
            p["pred"] = pd.to_numeric(p["pred"], errors="coerce")
            p["unscored"] = p["unscored"].astype(str)
            sub = sub.merge(p[["region_mutant", "pred", "unscored"]], on="region_mutant", how="left")
            sub["reason"] = np.where(sub.unscored != "", "unscored:" + sub.unscored, "")
            sub.loc[sub.unscored != "", "pred"] = np.nan
            if model == "evmutation":   # cross-check the flag with the gap rule recomputed from the region MSA
                excl = set(r["evmutation_gap_excluded"])
                rule = sub.region_mutant.map(lambda m: any(int(x[1:-1]) in excl for x in m.split(":")))
                flag = sub.unscored != ""
                if (rule != flag).any():
                    raise AssertionError(f"{rid}: evmutation unscored flag != gap rule for {(rule != flag).sum()}")
            sub = sub.drop(columns=["unscored"])
        else:
            st = status.get(rid, {})
            sub["pred"] = np.nan
            sub["reason"] = "region_failed:" + st.get("error", "not run")[:80] if st else "region_not_run"
        parts.append(sub)
    return pd.concat(parts, ignore_index=True), list(status.values())


def summary(models, tag, out):
    rows = []
    for m in models:
        for d in sorted((PREDS / m / tag).glob("*/W*")):
            u, w = d.parent.name, int(d.name[1:])
            if not (ARCH / u / d.name / "regions.json").exists():
                continue
            mp, st = mapped(m, tag, u, w)
            reg = json.loads((ARCH / u / d.name / "regions.json").read_text())
            ok = [s for s in st if s.get("ok")]
            fails = [s for s in st if not s.get("ok")]
            rows.append(dict(
                model=m, uniprot=u, L=reg["length"], W=w, regions=len(reg["regions"]), regions_ok=len(ok),
                regions_failed=len(fails),
                depth_kept_median=int(np.median([r["depth_kept"] for r in reg["regions"]])),
                n_mutants=len(mp), coverage=round(mp.pred.notna().mean(), 3),
                n_unscored=int(mp.reason.str.startswith("unscored").sum()),
                n_region_failed=int(mp.reason.str.startswith("region_").sum()),
                time_total_min=round(sum(s.get("t_total_s", 0) for s in st) / 60, 1),
                time_max_region_min=round(max([s.get("t_total_s", 0) for s in ok] or [0]) / 60, 1),
                peak_rss_gb=max([max(s.get("maxrss_self_gb", 0), s.get("maxrss_children_gb", 0)) for s in st] or [0]),
                errors="; ".join(sorted({s.get("error", "")[:60] for s in fails}))))
    t = pd.DataFrame(rows)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 30)
    print(t.to_string(index=False))
    if out:
        t.to_csv(out, sep="\t", index=False)


def compare(model, tag, u, windows):
    tabs = {w: mapped(model, tag, u, w)[0].set_index("mutant").pred for w in windows}
    df = pd.DataFrame(tabs)
    print(f"{model} {u}: coverage " + ", ".join(f"W{w} {df[w].notna().mean():.3f}" for w in windows))
    for i, a in enumerate(windows):
        for b in windows[i + 1:]:
            x = df[[a, b]].dropna()
            r = x[a].corr(x[b]) if len(x) > 2 else np.nan
            rs = x[a].corr(x[b], method="spearman") if len(x) > 2 else np.nan
            print(f"  W{a} vs W{b}: n={len(x)} pearson {r:.4f} spearman {rs:.4f}")
    return df


def export(model, tag, w, dest):
    dest.mkdir(parents=True, exist_ok=True)
    n = 0
    for d in sorted((PREDS / model / tag).glob(f"*/W{w}")):
        u = d.parent.name
        mp, _ = mapped(model, tag, u, w)
        mp[["mutant", "pred"]].dropna().to_csv(dest / f"{u}.csv", index=False)
        n += 1
    print(f"exported {n} proteins to {dest}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    s = sp.add_parser("summary")
    s.add_argument("--models", nargs="+", required=True)
    s.add_argument("--tag", required=True)
    s.add_argument("--out", type=Path)
    c = sp.add_parser("compare")
    c.add_argument("--model", required=True)
    c.add_argument("--tag", required=True)
    c.add_argument("--uniprot", required=True)
    c.add_argument("--windows", type=int, nargs="+", required=True)
    e = sp.add_parser("export")
    e.add_argument("--model", required=True)
    e.add_argument("--tag", required=True)
    e.add_argument("--window", type=int, required=True)
    e.add_argument("--dest", type=Path)
    a = ap.parse_args()
    if a.cmd == "summary":
        summary(a.models, a.tag, a.out)
    elif a.cmd == "compare":
        compare(a.model, a.tag, a.uniprot, a.windows)
    else:
        export(a.model, a.tag, a.window, a.dest or PRED_DEST / f"v2__{a.model}")


if __name__ == "__main__":
    main()
