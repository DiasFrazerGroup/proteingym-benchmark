"""Validate and run the windowed PLMs (esm2_window / esmc_window) on the Genebass archives.

Runs inside the zoo container of the base model (esm2 or esmc image) with the window package on PYTHONPATH.
Tests (docs/clinical/models_long_proteins.md §4):
  T1  proteins with L <= 1022: window model == the zoo's own evedesign ESM2.score_mutants (both special_tokens modes)
  T2  every protein: special_tokens termini vs every_window (comparison with reference scores, e.g. dark energy,
      in compare_plm_windows.py, outside the container)
  T3  (--windows) other window sizes, e.g. 2048 and full length (0 = whole protein)
Writes per-variant scores for every mode to --out/scores.parquet (uniprot_id, mutant, L, column per mode).
Usage: python validate_plm_windows.py esm2|esmc ARCHIVE_DIR CARD --out DIR [--windows 2048,0] [--limit N]
"""
import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd
from evedesign.models.esm2 import ESM2
from evedesign.proteingym import dataset_to_evedesign
from evedesign.system import Mutation
from proteingym.base import Subsets
from proteingym.base.model import ModelCard

ap = argparse.ArgumentParser()
ap.add_argument("which", choices=["esm2", "esmc"])
ap.add_argument("archives", type=Path)
ap.add_argument("card", type=Path)
ap.add_argument("--out", type=Path, required=True)
ap.add_argument("--windows", default="")
ap.add_argument("--limit", type=int, default=0, help="first N proteins only (smoke test)")
ap.add_argument("--batch-size", type=int, default=16, help="windows per forward pass (L <= 2048)")
a = ap.parse_args()
a.out.mkdir(parents=True, exist_ok=True)

if a.which == "esm2":
    from proteingym.models.esm2_window.model import ESM2Window as Cls
else:
    from proteingym.models.esmc_window.model import ESMCWindow as Cls
hp = dict(ModelCard.from_path(a.card).hyper_parameters)
hp.pop("window", None), hp.pop("special_tokens", None)
hp["device"] = "cuda"

# one archive per protein; union of its pairs' mutants
by_prot = {}
for f in sorted(a.archives.glob("*.splits.pgdata")):
    sub = Subsets.from_path(f)
    system, _, test = dataset_to_evedesign(sub, split=None, target="beta", test_fold=None)
    inst, _, _, _ = test.select(name="beta", drop_missing=False)
    ref = "".join(system[0].rep)
    d = by_prot.setdefault(ref, dict(system=system, muts={}, name=f.name))
    for x in inst:
        m = [(i, r, b) for i, (r, b) in enumerate(zip(ref, "".join(x[0].rep))) if b != r]
        d["muts"]["".join(f"{r}{i + 1}{b}" for i, r, b in m)] = m
print(f"{len(by_prot)} proteins, {sum(len(d['muts']) for d in by_prot.values())} unique mutants", flush=True)

modes = [("termini", 1022), ("every_window", 1022)] + [("termini", int(w)) for w in a.windows.split(",") if w]
rows, t1 = [], []
for k, (ref, d) in enumerate(sorted(by_prot.items(), key=lambda x: len(x[0]))):
    if a.limit and k >= a.limit:
        break
    system, L = d["system"], len(ref)
    fi = system[0].first_index
    names = list(d["muts"])
    muts = [[Mutation(entity=0, pos=fi + i, ref=r, to=b) for i, r, b in d["muts"][n]] for n in names]
    t = time.time()
    res = dict(mutant=names)
    for st, w in modes:
        w_eff = w if w > 0 else 10 ** 7
        col = f"{st}_{w if w > 0 else 'full'}"
        try:
            bs = a.batch_size if (L <= 2048 or 0 < w <= 2048) else 2
            m = Cls(window=w_eff, special_tokens=st, **dict(hp, batch_size=bs))
            m.build(system)
            res[col] = [x.score for x in m.score_mutants(system.rep_to_instance(), muts)]
        except Exception as e:  # e.g. out of memory at full length
            print(f"  {d['name']} L={L} {col}: FAILED {str(e)[:150]}", flush=True)
            res[col] = [np.nan] * len(names)
    if L <= 1022:   # T1: the zoo's own code path on the same model object
        zoo = np.array([x.score for x in ESM2.score_mutants(m, system.rep_to_instance(), muts)])
        diff = max(np.nanmax(np.abs(np.array(res[c]) - zoo)) for c in ("termini_1022", "every_window_1022"))
        t1.append(diff)
    uid = d["name"].split("_GENEBASS")[0]
    df = pd.DataFrame(res).assign(entry=uid, L=L)
    rows.append(df)
    print(f"[{k + 1}/{len(by_prot)}] {uid} L={L} n={len(names)} {time.time() - t:.1f}s", flush=True)

s = pd.concat(rows, ignore_index=True)
s.to_parquet(a.out / "scores.parquet")
if t1:
    print(f"\nT1 (L<=1022, {len(t1)} proteins): max |window - zoo score_mutants| = {max(t1):.2e} "
          f"{'PASS' if max(t1) < 1e-4 else 'FAIL'}")
long = s.L > 1022
print(f"T2 termini vs every_window, L>1022: r = {s[long][['termini_1022', 'every_window_1022']].corr().iloc[0, 1]:.5f}, "
      f"max|diff| = {np.abs(s[long].termini_1022 - s[long].every_window_1022).max():.3f}")
print(f"T2 L<=1022: max|termini - every_window| = {np.abs(s[~long].termini_1022 - s[~long].every_window_1022).max():.2e}")
