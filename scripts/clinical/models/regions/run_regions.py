"""Score region archives with an unchanged zoo model, inside the model's container (one container start).

Per region archive (make_regions.py): the zoo runner's own steps (Subsets.from_path -> dataset_to_evedesign ->
model class from the card -> build -> score), with two additions that do not change the model:
- the region wild type is scored together with the variants and pred = score(variant) - score(wild type).
  eve, evmutation, evmutation2 and structuredca return absolute sequence scores (ELBO, Hamiltonian,
  log-likelihood, -energy) whose offset differs between regions; gemme and rsalor are already WT-relative
  (score(WT) = 0, so the subtraction is a no-op for them).
- unscored variants are flagged (NaN after mapping back), never left as WT:
  * any mutated position outside model.positions() (evmutation: columns excluded for > 50% gaps; eve: non-focus
    columns), and
  * a NaN cell in the model's per-position matrix `_pssm` at (position, mutant aa) (rsalor: no RSA at that
    position; gemme: NA cell of pred.R). The wrappers add 0 for these.
One failing region is logged and skipped; the others are still scored.

Usage (inside the container): run_regions.py <card.md> <out_dir> <region archives...> [--hp '{"key": value}']
Writes <out_dir>/<UNIPROT>/W<window>/<region>.csv (region_mutant, uniprot_mutant, score, score_wt, pred, unscored) and appends one
JSON line per region to <out_dir>/status.jsonl (time, peak RSS, error).
"""
import argparse
import importlib
import io
import json
import os
import resource
import sys
import time
import traceback
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

from evedesign.proteingym import dataset_to_evedesign, seqs_to_instances
from proteingym.base import Subsets
from proteingym.base.model import ModelCard

VALID_AA_SORTED = "ACDEFGHIKLMNPQRSTVWY"


def import_class(s: str):
    mod, _, attr = s.partition(":")
    if not attr:
        mod, _, attr = s.rpartition(".")
    return getattr(importlib.import_module(mod), attr)


def maxrss_gb():
    return (round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20, 2),
            round(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 2**20, 2))


def assay_table(path: Path) -> pd.DataFrame:
    """The region assay CSV with all columns (Dataset.to_df keeps only sequence and targets)."""
    with zipfile.ZipFile(path) as z, zipfile.ZipFile(io.BytesIO(z.read("dataset.pgdata"))) as inner:
        name = next(n for n in inner.namelist() if n.startswith("assays/"))
        return pd.read_csv(io.BytesIO(inner.read(name)))


def unscored_flags(model, mutants: list[str]) -> list[str]:
    """Reason a region mutant is not scored by the model ('' = scored). Positions 1-based in the region."""
    try:
        modelled = {int(p) for _, p in model.positions()}
    except Exception:
        modelled = None
    pssm = getattr(model, "_pssm", None)
    out = []
    for m in mutants:
        why = []
        for x in m.split(":"):
            p, aa = int(x[1:-1]), x[-1]
            if modelled is not None and p not in modelled:
                why.append(f"pos{p}:not_modelled")
            elif pssm is not None and aa in VALID_AA_SORTED and np.isnan(pssm[p - 1, VALID_AA_SORTED.index(aa)]):
                why.append(f"pos{p}:nan_cell")
        out.append(";".join(why))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("card", type=Path)
    ap.add_argument("out", type=Path)
    ap.add_argument("archives", type=Path, nargs="+")
    ap.add_argument("--hp", default=os.environ.get("REGIONS_HP") or "{}",
                    help="JSON dict overriding card hyper_parameters (default: $REGIONS_HP)")
    ap.add_argument("--redo", action="store_true")
    a = ap.parse_args()

    card = ModelCard.from_path(a.card)
    hp = {**card.hyper_parameters, **json.loads(a.hp)}
    Model = import_class(card.model_class)
    a.out.mkdir(parents=True, exist_ok=True)
    print(f"model {card.model_class} hp overrides {a.hp}; {len(a.archives)} regions", flush=True)
    for f in a.archives:
        rid = f.name.removesuffix(".splits.pgdata")
        # region archives live in <UNIPROT>/W<window>/: keep that layout (same span can occur in several windows)
        od = a.out / f.parent.parent.name / f.parent.name if f.parent.name.startswith("W") else a.out
        od.mkdir(parents=True, exist_ok=True)
        dst = od / f"{rid}.csv"
        if dst.exists() and not a.redo:
            continue
        t = time.time()
        rec = dict(region=rid, archive=str(f))
        try:
            subsets = Subsets.from_path(f)
            target = subsets.dataset.assay_targets[0].name
            system, _, test = dataset_to_evedesign(subsets, split=None, target=target, test_fold=None)
            rec["n_msa"] = len(system[0].sequences.seqs) if system[0].sequences is not None else 0
            rec["has_structure"] = bool(system[0].structures)
            t_load = time.time()
            model = Model(**hp).build(system, data=None)
            t_build = time.time()
            instances, _, _, _ = test.select(name=target, drop_missing=False)
            wt = "".join(system[0].rep)
            scored = model.score(list(instances) + seqs_to_instances([wt]))
            t_score = time.time()
            s = np.array([float(x.score) for x in scored])
            # instances come in Dataset.to_df() order: match them to assay rows by sequence (unique per mutant)
            df = assay_table(f).drop_duplicates("sequence").set_index("sequence")
            if "uniprot_mutant" not in df:   # an original (non-region) archive, for comparisons
                df["uniprot_mutant"] = df["mutant"]
            df = df.loc[["".join(x[0].rep) for x in instances]]
            muts = df["mutant"].tolist()
            res = pd.DataFrame(dict(region_mutant=muts, uniprot_mutant=df["uniprot_mutant"].tolist(),
                                    score=s[:-1], score_wt=s[-1], pred=s[:-1] - s[-1],
                                    unscored=unscored_flags(model, muts)))
            res.to_csv(dst, index=False)
            rec.update(ok=True, n=len(res), n_unscored=int((res.unscored != "").sum()),
                       t_load_s=round(t_load - t, 1), t_build_s=round(t_build - t_load, 1),
                       t_score_s=round(t_score - t_build, 1))
            print(f"[{time.strftime('%H:%M:%S')}] OK {rid} n={len(res)} unscored={rec['n_unscored']} "
                  f"msa={rec['n_msa']} build {t_build - t_load:.0f}s score {t_score - t_build:.0f}s "
                  f"rss {maxrss_gb()}", flush=True)
        except Exception as e:
            rec.update(ok=False, error=f"{type(e).__name__}: {str(e)[:300]}")
            (od / f"{rid}.error.log").write_text(traceback.format_exc())
            print(f"[{time.strftime('%H:%M:%S')}] FAILED {rid} {rec['error']}", flush=True)
        rec["t_total_s"] = round(time.time() - t, 1)
        rec["maxrss_self_gb"], rec["maxrss_children_gb"] = maxrss_gb()
        with open(od / "status.jsonl", "a") as fh:
            fh.write(json.dumps(rec) + "\n")
    print("RUN REGIONS DONE", flush=True)


if __name__ == "__main__":
    sys.exit(main())
