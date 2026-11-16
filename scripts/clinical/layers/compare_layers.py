"""Compare two layered runs of the same archives (e.g. dbNSFP-based OLD vs own-source NEW): per-layer coverage
(rows, archives, proteins), what is still missing, agreement where both have a value, and optionally the
AlphaMissense layer against per-protein reference predictions (CSV mutant,pred; pred = -pathogenicity, ProteinGym
orientation, e.g. the VESM paper's AlphaMissense files).

Run: python compare_layers.py OLD_DIR NEW_DIR [--am-ref DIR]
"""
import argparse
import io
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import toml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from add_layers import read_archive  # noqa: E402

SCORES = ["score_alphamissense", "score_popeve", "score_cadd_phred"]


def load(d):
    out = []
    for p in sorted(Path(d).glob("*.pgdata")) + sorted(Path(d).glob("*.pgarchive")):
        _, inner = read_archive(p)
        man = toml.loads(inner["manifest.lock"][0].decode())
        wt = next(s for s in man["sequences"] if s.get("type", "wild_type") == "wild_type")
        for a in man["assays"]:
            df = pd.read_csv(io.BytesIO(inner[a["path"]][0]), dtype=str, keep_default_na=False)
            out.append(df.assign(archive=p.name, uniprot=(wt.get("uniprot_id") or "").split("-")[0]))
    return pd.concat(out, ignore_index=True)


def coverage(d, c):
    has = d[c] != ""
    full = d.groupby("archive")[c].apply(lambda s: (s != "").all()).sum()
    none = sorted(set(d.uniprot) - set(d[has].uniprot))
    return (f"rows {has.sum()}/{len(d)} ({has.mean():.1%}, missing {(~has).sum()}); archives with any "
            f"{d[has].archive.nunique()}/{d.archive.nunique()}, complete {full}; proteins with any "
            f"{d[has].uniprot.nunique()}/{d.uniprot.nunique()}; proteins with none: {' '.join(none) or '-'}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("old", type=Path)
    ap.add_argument("new", type=Path)
    ap.add_argument("--am-ref", type=Path)
    a = ap.parse_args()
    old, new = load(a.old), load(a.new)
    assert len(old) == len(new) and (old.archive.values == new.archive.values).all(), "different archives/rows"
    for c in SCORES:
        print(f"\n## {c}\n  old: {coverage(old, c)}\n  new: {coverage(new, c)}")
        miss = new[new[c] == ""]
        if len(miss):
            print("  still missing per protein:", miss.groupby("uniprot").size().to_dict())
            cols = [x for x in ["uniprot", "mutant", "variant_id", "uniprot_pos", "score_snv_source"] if x in miss]
            print(miss[cols].head(15).to_string())
        both = (old[c] != "") & (new[c] != "")
        o, n = old.loc[both, c].astype(float), new.loc[both, c].astype(float)
        dif = (o - n).abs()
        print(f"  where both have a value: {both.sum()} rows, identical {(dif < 1e-9).sum()}, |diff| > 1e-3 "
              f"{(dif > 1e-3).sum()}, max |diff| {dif.max():.4g}, Pearson {np.corrcoef(o, n)[0, 1]:.6f}")
        if (dif > 1e-3).any():
            x = new.loc[both].assign(old=o, new=n)[dif > 1e-3]
            print("  differing rows per protein:", x.groupby("uniprot").size().to_dict())
            print(x[["uniprot", "mutant", "variant_id", "old", "new"]].head(15).to_string())
    for c in ["score_alphamissense_source", "score_snv_source"]:
        if c in new:
            print(f"\n{c}: {new[c].value_counts().to_dict()}")
    if a.am_ref:
        rows = []
        for u, g in new.groupby("uniprot"):
            f = a.am_ref / f"{u}.csv"
            if not f.exists():
                continue
            ref = pd.read_csv(f).drop_duplicates("mutant").set_index("mutant").pred
            key = [f"{m[0]}{p}{m[-1]}" if p and ":" not in p else "" for m, p in zip(g.mutant, g.uniprot_pos)]
            rows.append(pd.DataFrame({"uniprot": u, "mutant": g.mutant.values, "am": g.score_alphamissense.values,
                                      "ref": ref.reindex(key).values}))
        r = pd.concat(rows)
        print(f"\n## AlphaMissense vs {a.am_ref}: reference files for {r.uniprot.nunique()} proteins, "
              f"{len(r)} rows, reference has {r.ref.notna().sum()}, new layer has {(r.am != '').sum()}")
        r = r[(r.am != "") & r.ref.notna()]
        am, ref = r.am.astype(float), r.ref.astype(float)
        dif = (am + ref).abs()
        print(f"  both: {len(r)} rows; Pearson(new, ref) {np.corrcoef(am, ref)[0, 1]:.6f}; new == -ref "
              f"{(dif < 1e-9).sum()}, |new + ref| > 1e-3 {(dif > 1e-3).sum()}, max {dif.max():.4g}")
        if (dif > 1e-3).any():
            print(r[dif > 1e-3].head(15).to_string())


if __name__ == "__main__":
    main()
