"""EVE-protocol alignment selection on ColabFold a3m files: sweep bits/residue thresholds (docs/clinical/models_long_proteins.md §3.5).

EVE (Frazer et al. 2021; reused by popEVE): alignments at several bit-score thresholds (reference 0.3 bits/residue),
keep sequences aligned to >= 50% of the target, columns with >= 70% occupancy count as covered (L_cov); choose an
alignment with L_cov >= 0.8 L and 10 L <= N <= 100,000, relaxing when none qualifies.
ColabFold/MMseqs2 a3m headers carry: name, bitscore, seqid, evalue, qstart, qend, qlen, tstart, tend, tlen. Here
bits/residue = MMseqs2 bitscore / query length (MMseqs2 bit scores are not on HMMER's scale: thresholds are swept).
Reads each a3m once (streamed); writes one row per (protein, threshold).
"""
import sys, glob, os
from multiprocessing import Pool
import numpy as np
import pandas as pd

THR = [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5]
DEL = bytes(range(97, 123))   # lowercase = insertions


def sweep(path):
    u = os.path.basename(path).split(".")[0]
    T = np.array(THR)
    with open(path, "rb") as f:
        f.readline(); q = f.readline().strip()
        L = len(q)
        n = np.zeros(len(T), dtype=np.int64)            # sequences per bucket (highest threshold passed)
        occ = np.zeros((len(T), L), dtype=np.int64)
        n_all = n_cov = 0
        hdr = None
        for line in f:
            if line.startswith(b">"):
                hdr = line; continue
            n_all += 1
            s = line.strip().translate(None, DEL)
            if len(s) != L:
                continue
            ng = np.frombuffer(s, dtype=np.uint8) != 45   # '-'
            if ng.sum() < 0.5 * L:
                continue
            n_cov += 1
            bits = float(hdr.split(b"\t")[1]) / L
            b = np.searchsorted(T, bits, side="right") - 1
            n[b] += 1; occ[b] += ng
    n = np.cumsum(n[::-1])[::-1]; occ = np.cumsum(occ[::-1], axis=0)[::-1]
    rows = []
    for k, t in enumerate(T):
        lcov = int((occ[k] >= 0.7 * max(n[k], 1)).sum()) if n[k] else 0
        rows.append(dict(uniprot=u, L=L, thr=t, N=int(n[k]), N_over_L=n[k] / L, Lcov=lcov, Lcov_frac=lcov / L,
                         n_rows=n_all, n_cov50=n_cov))
    return rows


if __name__ == "__main__":
    files = sorted(glob.glob(sys.argv[1] + "/*.a3m"), key=os.path.getsize, reverse=True)
    with Pool(int(os.environ.get("SLURM_CPUS_PER_TASK", 4))) as p:
        out = [r for rs in p.imap_unordered(sweep, files) for r in rs]
    df = pd.DataFrame(out).sort_values(["L", "uniprot", "thr"])
    df.to_csv(sys.argv[2], sep="\t", index=False)
    ok = df[(df.Lcov_frac >= 0.8) & (df.N <= 100_000) & (df.N >= 10 * df.L)]
    pick = ok.sort_values("thr").groupby("uniprot").head(1)
    print(f"{df.uniprot.nunique()} proteins; strict EVE criteria met by {pick.uniprot.nunique()}")
    print(pick[["uniprot", "L", "thr", "N", "Lcov_frac", "n_rows"]].to_string(index=False))
