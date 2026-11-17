"""Validate the *_filtered alignment models (docs/clinical/models_long_proteins.md §3.6, §9): one build + score per call.

Runs inside a zoo model's image with the filtered package mounted on PYTHONPATH (validate_msa_filter.sbatch). The
zoo runner's own steps (Subsets.from_path -> dataset_to_evedesign -> model_class(**hyper_parameters) -> build ->
score), on one archive, scoring: the archive's variants (Genebass mutants), all single substitutions (or a seeded
sample of `--singles` of them) and the wild type. pred = score - score(WT) (absolute scores of EVE / EVmutation /
StructureDCA are not WT-relative). Writes <out>/<tag>.csv (mutant, in_archive, score, pred) and <out>/<tag>.json
(times, peak RSS self/children, MSA depth, filter stats, NaN count).

  python validate_msa_filter.py run <card.md> <archive> <out_dir> <tag> [--hp JSON] [--class X] [--singles all|N]
"""
import argparse
import importlib
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

AA = "ACDEFGHIKLMNPQRSTVWY"


def import_class(s: str):
    mod, _, attr = s.partition(":")
    if not attr:
        mod, _, attr = s.rpartition(".")
    return getattr(importlib.import_module(mod), attr)


def maxrss_gb():
    return (round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20, 2),
            round(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 2**20, 2))


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg} rss={maxrss_gb()}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["run"])
    ap.add_argument("card", type=Path)
    ap.add_argument("archive", type=Path)
    ap.add_argument("out", type=Path)
    ap.add_argument("tag")
    ap.add_argument("--hp", default="{}", help="JSON dict overriding card hyper_parameters")
    ap.add_argument("--class", dest="cls", default=None, help="model class import string (default: card)")
    ap.add_argument("--singles", default="all", help="'all' or a number of single substitutions to sample")
    a = ap.parse_args()

    from proteingym.base.model import ModelCard
    card = ModelCard.from_path(a.card)
    Model = import_class(a.cls or card.model_class)   # imports the package first: installs the header hook
    from evedesign.proteingym import dataset_to_evedesign, seqs_to_instances
    from proteingym.base import Subsets

    hp = {**card.hyper_parameters, **json.loads(a.hp)}
    a.out.mkdir(parents=True, exist_ok=True)
    rec = dict(tag=a.tag, archive=a.archive.name, model_class=a.cls or card.model_class,
               hp_overrides=json.loads(a.hp))
    t0 = time.time()
    subsets = Subsets.from_path(a.archive)
    target = subsets.dataset.assay_targets[0].name
    system, _, test = dataset_to_evedesign(subsets, split=None, target=target, test_fold=None)
    wt = "".join(system[0].rep)
    L = len(wt)
    rec.update(L=L, msa_rows=len(system[0].sequences.seqs), t_load_s=round(time.time() - t0, 1),
               bitscores_attached=hasattr(system[0].sequences, "pg_bitscores"))
    log(f"loaded {a.archive.name}: L={L} rows={rec['msa_rows']} in {rec['t_load_s']}s")

    # mutants: archive variants + singles + WT
    inst, _, _, _ = test.select(name=target, drop_missing=False)
    arch = set()
    for x in inst:
        s = "".join(x[0].rep)
        d = [f"{w}{i + 1}{m}" for i, (w, m) in enumerate(zip(wt, s)) if w != m]
        if d and len(s) == L:
            arch.add(":".join(d))
    singles = [f"{w}{i + 1}{m}" for i, w in enumerate(wt) for m in AA if m != w]
    if a.singles != "all":
        rng = np.random.default_rng(0)
        singles = list(rng.choice(singles, size=min(int(a.singles), len(singles)), replace=False))
    muts = sorted(set(singles) | arch, key=lambda m: (int(m.split(":")[0][1:-1]), m))

    def apply(m):
        s = list(wt)
        for x in m.split(":"):
            s[int(x[1:-1]) - 1] = x[-1]
        return "".join(s)

    seqs = [apply(m) for m in muts] + [wt]
    rec["n_scored"] = len(muts)

    t1 = time.time()
    model = Model(**hp).build(system, data=None)
    rec["t_build_s"] = round(time.time() - t1, 1)
    rec["msa_filter_stats"] = {k: v for k, v in (getattr(model, "msa_filter_stats", None) or {}).items()
                               if k != "sweep"} or None
    log(f"built in {rec['t_build_s']}s")

    last = [0.0]

    def progress(status, pct=None, msg=None, *args):
        if time.time() - last[0] > 300:
            last[0] = time.time()
            log(f"score progress {pct} {msg}")

    t2 = time.time()
    try:
        scored = model.score(seqs_to_instances(seqs), status_callback=progress)
    except TypeError:
        scored = model.score(seqs_to_instances(seqs))
    rec["t_score_s"] = round(time.time() - t2, 1)
    s = np.array([np.nan if x.score is None else float(x.score) for x in scored])
    df = pd.DataFrame(dict(mutant=muts, in_archive=[m in arch for m in muts], score=s[:-1], pred=s[:-1] - s[-1]))
    df.to_csv(a.out / f"{a.tag}.csv", index=False)
    rec.update(score_wt=float(s[-1]), n_nan=int(np.isnan(df.pred).sum()), n_archive=len(arch),
               t_total_s=round(time.time() - t0, 1))
    rec["maxrss_self_gb"], rec["maxrss_children_gb"] = maxrss_gb()
    (a.out / f"{a.tag}.json").write_text(json.dumps(rec, indent=1))
    log(f"DONE {a.tag}: {json.dumps(rec)}")


if __name__ == "__main__":
    sys.exit(main())
