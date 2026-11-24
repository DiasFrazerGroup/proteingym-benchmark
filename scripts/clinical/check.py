"""Check the clinical benchmark archives against prepare.py's variants, and report the reference models.

Per benchmark:
  rows       the index lists every dataset; the archives' assay rows are exactly the prepared variants
  sequences  every assay sequence is the wild type with exactly its `mutant` substitution
  content    no model score columns; the MSA's first record is the wild type; one structure per protein <= 2,700 aa
  loading    the smallest archives load with proteingym-base (Subsets.from_path)
  reference  clinvar / denovo / gnomad: six published predictors (clinical-curation's scores, kept beside the
             variants, never archived) on the variants all six score. Pooled metric over all variants and the mean
             per gene: AUROC (clinvar, denovo; per gene needs both classes), Spearman rho with log10 AF (gnomad,
             lower = better; per gene needs >= 20 variants). De novo also gets the gene-only baseline (does the gene
             have another case variant in the ap_test_set), which beats every predictor there.
Writes <work>/reference_results.tsv. Exits 1 if a check fails.
"""
import argparse
import io
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from proteingym.base import Subsets
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

from common import MAX_LEN

TARGET = {"clinvar": "label", "denovo": "dnm_case", "gnomad": "log10_af", "genebass": "beta"}
SCORES = {"score_alphamissense": "AlphaMissense", "score_popeve": "popEVE", "score_cadd_phred": "CADD v1.7",
          "score_esm1b": "ESM1b", "score_vesm3b": "VESM-3B", "score_gpn_msa": "GPN-MSA"}


def read_archive(path):
    """(assay table, wild type, first MSA record, number of structures) without parsing the MSA."""
    with zipfile.ZipFile(path) as z, zipfile.ZipFile(io.BytesIO(z.read("dataset.pgdata"))) as d:
        names = d.namelist()
        assay = pd.read_csv(d.open(next(n for n in names if n.startswith("assays/"))), low_memory=False)
        wt = "".join(d.read(next(n for n in names if n.startswith("sequences/"))).decode().split("\n")[1:])
        with d.open(next(n for n in names if n.startswith("msas/"))) as f:
            f.readline()
            first = f.readline().decode().strip()
        return assay, wt, first, sum(n.startswith("structures/") for n in names)


def reference(b, v, dnm_cases):
    y, rows = TARGET[b], []
    x = v.dropna(subset=list(SCORES))
    for c, model in SCORES.items():
        if b == "gnomad":
            pooled = spearmanr(x[c], x[y]).statistic
            per = [spearmanr(t[c], t[y]).statistic for _, t in x.groupby("gene") if len(t) >= 20]
        else:
            pooled = roc_auc_score(x[y], x[c])
            per = [roc_auc_score(t[y], t[c]) for _, t in x.groupby("gene") if t[y].nunique() == 2]
        rows.append(dict(benchmark=b, model=model, pooled=pooled, per_gene=np.nanmean(per), n_genes_per_gene=len(per),
                         n_variants=len(x), n_genes=x.gene.nunique()))
    if b == "denovo":
        other = (x.gene.map(dnm_cases).fillna(0) - x[y] > 0).astype(float)
        rows.append(dict(benchmark=b, model="gene-only baseline", pooled=roc_auc_score(x[y], other),
                         n_variants=len(x), n_genes=x.gene.nunique()))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("benchmarks", nargs="*", default=list(TARGET))
    ap.add_argument("--work", type=Path, required=True, help="prepare.py / build.py output folder")
    ap.add_argument("--index", type=Path, default=Path("datasets/clinical"))
    ap.add_argument("--curation", type=Path, default=None, help="clinical-curation repo (de novo gene-only baseline)")
    ap.add_argument("--load", type=int, default=5, help="archives per benchmark loaded with proteingym-base")
    a = ap.parse_args()
    dnm_cases = pd.Series(dtype=float)
    if a.curation:
        dnm_cases = pd.read_csv(a.curation / "data/processed/benchmarks/dnm_aptest_scores.tsv", sep="\t",
                                low_memory=False).groupby("gene").dnm_case.sum()
    failed, results = [], []
    for b in a.benchmarks:
        y = TARGET[b]
        v = pd.read_csv(a.work / b / "variants.tsv.gz", sep="\t", low_memory=False)
        idx = pd.read_csv(a.index / b / f"{b}.csv")
        bad = dict(rows=set(idx.dataset) != set(v.dataset), sequences=0, content=0)
        for r in idx.itertuples():
            assay, wt, first, n_struct = read_archive(a.work / b / "splits" / r.archive)
            x = v[v.dataset == r.dataset]
            bad["rows"] |= not (assay.mutant.tolist() == x.mutant.tolist() and np.allclose(assay[y], x[y]))
            for s, m in zip(assay.sequence, assay.mutant):
                p = int(m[1:-1])
                bad["sequences"] += not (wt[p - 1] == m[0] != m[-1] and s == wt[:p - 1] + m[-1] + wt[p:])
            bad["content"] += (any(c.startswith("score_") for c in assay) or first != wt
                               or (len(wt) <= MAX_LEN and n_struct != 1) or n_struct != r.n_structures)
        for arc in idx.sort_values("msa_n_sequences").archive[:a.load]:
            Subsets.from_path(a.work / b / "splits" / arc)
        print(f"{b}: {len(idx)} archives, {idx.uniprot_id.nunique()} proteins, {idx.n_variants.sum():,} variants; "
              f"longest protein {idx.seq_length.max():,} aa; MSA depth median {int(idx.msa_n_sequences.median()):,} "
              f"(min {idx.msa_n_sequences.min():,}); " + ", ".join(f"{k} {'FAILED' if n else 'ok'}"
                                                                    for k, n in bad.items()), flush=True)
        failed += [f"{b} {k}" for k, n in bad.items() if n]
        if b != "genebass":
            results += reference(b, v, dnm_cases)
    if results:
        res = pd.DataFrame(results)
        res.to_csv(a.work / "reference_results.tsv", sep="\t", index=False, float_format="%.4f")
        print(res.to_string(index=False, float_format="%.3f"))
    print("CHECK FAILED: " + ", ".join(failed) if failed else "CHECK OK")
    sys.exit(bool(failed))


if __name__ == "__main__":
    main()
