"""Validate and run the structure models for long proteins (saprot_window / esm-if1_window) on Genebass archives.

Runs inside the model's zoo image with the *_window package on PYTHONPATH.
  T1  proteins with one structure and L <= 1022: window model == zoo model from the same package
      saprot: identical scores; esm-if1: window == L * (zoo(mut) - zoo(wt))   (first 3 proteins, <= 40 mutants each)
  T4  proteins with several structures (AFDB fragments): placement and per-position assignment are logged
Every archive is scored; writes per-variant scores (deduplicated per protein) to --out/scores.parquet (entry, mutant, L, n_structures, score).
Usage: python validate_struct_windows.py saprot|esmif1 ARCHIVE_DIR CARD --out DIR [--limit N]
"""
import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd
from evedesign.proteingym import dataset_to_evedesign
from proteingym.base import Subsets
from proteingym.base.model import ModelCard

ap = argparse.ArgumentParser()
ap.add_argument("which", choices=["saprot", "esmif1"])
ap.add_argument("archives", type=Path)
ap.add_argument("card", type=Path)
ap.add_argument("--out", type=Path, required=True)
ap.add_argument("--limit", type=int, default=0)
a = ap.parse_args()
a.out.mkdir(parents=True, exist_ok=True)

if a.which == "saprot":
    from proteingym.models.saprot_window.window_model import SaProtWindow as Win
    from proteingym.models.saprot_window.model import SaProt as Zoo
else:
    from proteingym.models.esm_if1_window.window_model import ESMIF1Window as Win
    from proteingym.models.esm_if1_window.model import ESMIF1 as Zoo
hp = dict(ModelCard.from_path(a.card).hyper_parameters)
hp["device"] = "cuda"
zoo_hp = {k: v for k, v in hp.items() if k not in ("window", "special_tokens", "batch_size", "context")}
if a.which == "saprot" and "batch_size" in hp:
    zoo_hp["batch_size"] = hp["batch_size"]

seen, rows, t1 = set(), [], []
files = sorted(a.archives.glob("*.splits.pgdata"), key=lambda f: f.stat().st_size)
for f in files:
    sub = Subsets.from_path(f)
    system, _, test = dataset_to_evedesign(sub, split=None, target="beta", test_fold=None)
    ref = "".join(system[0].rep)
    first = ref not in seen                      # every archive is scored (pairs of one protein differ in variants)
    seen.add(ref)
    if a.limit and len(seen) > a.limit:
        break
    inst, _, _, _ = test.select(name="beta", drop_missing=False)
    L, ns = len(ref), len(system[0].structures)
    entry = f.name.split("_GENEBASS")[0]
    t = time.time()
    try:
        m = Win(**hp).build(system)
        sc = np.array(m.score(inst).scores if hasattr(m.score(inst), "scores") else [x.score for x in m.score(inst)])
    except Exception as e:
        print(f"{entry} L={L} structures={ns}: FAILED {type(e).__name__}: {str(e)[:200]}", flush=True)
        continue
    muts = ["".join(f"{r}{i + 1}{b}" for i, (r, b) in enumerate(zip(ref, "".join(x[0].rep))) if b != r) for x in inst]
    rows.append(pd.DataFrame(dict(entry=entry, mutant=muts, L=L, n_structures=ns, score=sc)))
    msg = f"{entry} L={L} structures={ns} n={len(inst)} {time.time() - t:.1f}s"
    if ns > 1 and a.which == "saprot" and first:
        asg = m.structure_assignment
        msg += f" | {len(set(asg))} structures used, uncovered {sum(x is None for x in asg)}"
    print(msg, flush=True)
    if first and ns == 1 and L <= hp.get("window", 1022) and len(t1) < 3:
        sub_inst = inst[:40]
        z = Zoo(**zoo_hp).build(system)
        zs = np.array([x.score for x in z.score(sub_inst)])
        ws = sc[:len(sub_inst)]
        if a.which == "saprot":
            d = float(np.abs(ws - zs).max())
        else:
            wt = np.array([x.score for x in z.score([system.rep_to_instance()])])[0]
            d = float(np.abs(ws - L * (zs - wt)).max())
        t1.append(d)
        print(f"  T1 {entry}: max diff vs zoo = {d:.2e}", flush=True)

out = pd.concat(rows, ignore_index=True)
spread = out.groupby(["entry", "mutant"]).score.agg(lambda x: x.max() - x.min()).max()   # same variant, other archive
print(f"same variant scored in several archives: max spread {spread:.2e}")
out.drop_duplicates(["entry", "mutant"]).to_parquet(a.out / "scores.parquet")
print(f"\nT1 ({len(t1)} proteins): max diff = {max(t1) if t1 else float('nan'):.2e} "
      f"{'PASS' if t1 and max(t1) < 1e-3 else 'CHECK'}")
