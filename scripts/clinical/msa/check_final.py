"""Checks of the final archives (filter_msas.py) against their sources.

For every archive in the capped index:
1. outer entries other than dataset.pgdata, and inner entries other than the MSA, are byte-identical to the source
   archive (assays, sequences, structures, layers, slices); exactly one structure; seq_length <= the cap;
2. the MSA is an ordered subset of the original ColabFold a3m: every record (header and sequence, byte for byte)
   occurs in the original, in the same order, query first;
3. the model-side filter is idempotent on it: filtering the stored MSA again keeps every row (same bytes);
4. the per-protein table matches (rows, sha256).

  check_final.py --index datasets/clinical/genebass/genebass_le2700.csv --final $WORK/genebass/final_le2700 \
      --source $WORK/genebass/splits_layers --msa-dir $WORK/genebass_msas \
      --filter-table datasets/clinical/genebass/genebass_le2700_msa_filter.tsv [--workers N]
"""
import argparse
import hashlib
import io
import os
import sys
import tomllib
import zipfile
from multiprocessing import Pool
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import filter_msas as fm   # noqa: E402


def inner(z):
    return zipfile.ZipFile(io.BytesIO(z.read("dataset.pgdata")))


def check_archive(args):
    name, final, source, max_len = args
    with zipfile.ZipFile(final / name) as zf, zipfile.ZipFile(source / name) as zs:
        assert [i.filename for i in zf.infolist()] == [i.filename for i in zs.infolist()], "outer entries"
        for i in zs.infolist():
            if i.filename != "dataset.pgdata":
                assert zf.read(i.filename) == zs.read(i.filename), f"outer {i.filename}"
        fi, si = inner(zf), inner(zs)
        assert [i.filename for i in fi.infolist()] == [i.filename for i in si.infolist()], "inner entries"
        man = tomllib.loads(fi.read("manifest.lock").decode())
        entry = man["msas"][0]["path"]
        for i in si.infolist():
            if i.filename != entry:
                assert fi.read(i.filename) == si.read(i.filename), f"inner {i.filename}"
        assert len(man.get("structures", [])) == 1, "structures"
        wt = b"".join(fi.read(man["sequences"][0]["path"]).split(b"\n")[1:]).strip()
        assert len(wt) <= max_len, "length"
        return name, man["msas"][0]["reference_sequence_name"], hashlib.sha256(fi.read(entry)).hexdigest(), len(wt)


def check_protein(args):
    uniprot, L, msa_bytes, msa_dir = args
    orig = fm.records((msa_dir / f"{uniprot}.a3m").read_bytes())
    new = fm.records(msa_bytes)
    assert new[0][0] == orig[0][0], "query first"
    j = 0
    for rec in new:                                  # ordered subset, exact records
        while j < len(orig) and orig[j][0] != rec[0]:
            j += 1
        assert j < len(orig), "record not found in order in the original a3m"
        j += 1
    again, stats, idx = fm.filter_a3m(msa_bytes, L)
    assert again == msa_bytes and len(idx) == len(new), f"not idempotent: {len(idx)} of {len(new)} rows kept"
    return uniprot, len(new), hashlib.sha256(msa_bytes).hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--index", type=Path, required=True)
    ap.add_argument("--final", type=Path, required=True)
    ap.add_argument("--source", type=Path, required=True)
    ap.add_argument("--msa-dir", type=Path, required=True)
    ap.add_argument("--filter-table", type=Path, required=True)
    ap.add_argument("--max-length", type=int, default=2700)
    ap.add_argument("--workers", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 4)))
    a = ap.parse_args()
    idx = pd.read_csv(a.index)
    files = sorted(p.name for p in a.final.glob("*.splits.pgdata"))
    assert files == sorted(Path(p).name for p in idx.archive), "archive set differs from the index"
    with Pool(a.workers) as pool:
        arch = pool.map(check_archive, [(n, a.final, a.source, a.max_length) for n in files])
    print(f"1. {len(arch)} archives: non-MSA entries byte-identical to the source, one structure, L <= "
          f"{a.max_length}: OK", flush=True)
    by_prot = {}
    for name, u, sha, L in arch:
        by_prot.setdefault(u, set()).add((sha, L))
    assert all(len(v) == 1 for v in by_prot.values()), "archives of one protein carry different MSAs"
    jobs = []
    for u, v in by_prot.items():
        sha, L = next(iter(v))
        name = next(n for n, uu, s, _ in arch if uu == u)
        with zipfile.ZipFile(a.final / name) as z:
            fi = inner(z)
            man = tomllib.loads(fi.read("manifest.lock").decode())
            jobs.append((u, L, fi.read(man["msas"][0]["path"]), a.msa_dir))
    jobs.sort(key=lambda j: -len(j[2]))
    with Pool(a.workers) as pool:
        prot = pool.map(check_protein, jobs)
    print(f"2-3. {len(prot)} proteins: MSA = ordered subset of the ColabFold a3m (exact records), filter idempotent: OK")
    t = pd.read_csv(a.filter_table, sep="\t").set_index("uniprot")
    for u, n, sha in prot:
        assert t.loc[u, "N"] == n and t.loc[u, "a3m_sha256"] == sha, f"{u}: filter table mismatch"
    assert set(t.index) == set(by_prot), "filter table proteins"
    print("4. filter table rows and sha256 match: OK")
    print("CHECK FINAL OK")


if __name__ == "__main__":
    main()
