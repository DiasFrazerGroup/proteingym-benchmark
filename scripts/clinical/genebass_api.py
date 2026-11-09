"""Genebass benchmark tables straight from Genebass (public API behind app.genebass.org), same rules as VESM Fig. 6.

Replaces the VESM authors' extraction (figshare Fig6_data) as the source of association data. Writes tables in the
same schema as the paper's, so build_genebass.py can use them (--hits / --betas):

  pairs   Genebass portal top pLoF gene-phenotype associations (SKAT-O P < 1e-6) for the "Blood biochemistry"
          phenotypes (the paper's 332 = this list; checked against the paper and reported)
  gene    pval_/beta_ {plof, mis, syn} = SKAT-O Pvalue / BETA_Burden of burden sets pLoF, missense|LC, synonymous
  variant single-variant BETA / Pvalue (+ test AC/AF) from variant-associations, joined to variant-qc annotation
          (consequence, hgvsp, cohort allele_count/number/frequency); missense_variant rows only
  group   biomarker_group per phenotype: the paper's grouping (a paper definition, not Genebass data)

Inputs: --raw   local API download (db/genebass/raw, read-only); pairs missing there are fetched into --fetch
        --fig6  paper Fig6_data (biomarker groups + comparison only)
Outputs (--out): genebass_api_hits.csv, genebass_api_betas.csv.gz, compare_vs_paper.md
"""
import argparse
import gzip
import json
import time
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

API = "https://main.genebass.org/api"
SETS = {"plof": "pLoF", "mis": "missense_LC", "syn": "synonymous"}


def jload(p):
    with gzip.open(p, "rt") as f:
        return json.load(f)


def api_get(path):
    url = f"{API}/{path}"
    for attempt in range(6):
        try:
            req = urllib.request.Request(url, headers={"Accept-Encoding": "gzip"})
            with urllib.request.urlopen(req, timeout=300) as r:
                b = r.read()
            time.sleep(2)                                           # keep requests sequential and paced
            return json.loads(gzip.decompress(b) if b[:2] == b"\x1f\x8b" else b)
        except Exception as e:  # noqa: BLE001
            print(f"  retry {attempt + 1} {url}: {e}", flush=True)
            time.sleep(15 * (attempt + 1))
    raise RuntimeError(f"failed: {url}")


def cached(raw, fetch, rel, path):
    for base in (raw, fetch):
        if (base / rel).exists():
            return jload(base / rel)
    d = api_get(path)
    (fetch / rel).parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(fetch / rel, "wt") as f:
        json.dump(d, f)
    print(f"fetched {rel}", flush=True)
    return d


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--raw", type=Path, default=Path("/no_backup/diasfrazer/db/genebass/raw"))
    ap.add_argument("--fetch", type=Path, required=True)
    ap.add_argument("--fig6", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    a.fetch.mkdir(parents=True, exist_ok=True)

    phen = pd.DataFrame(jload(a.raw / "meta/phenotypes.json.gz"))
    bio = phen[phen.astype(str).apply(lambda r: "Blood biochemistry" in " ".join(r), axis=1)]
    desc = bio.set_index("analysis_id").description.to_dict()
    top = pd.DataFrame(jload(a.raw / "meta/top_associations__pLoF.json.gz"))
    pairs = top[top.analysis_id.isin(desc) & (top.Pvalue < 1e-6)].copy()
    pairs["code"] = pairs.analysis_id.str.split("-").str[1]
    pairs["gene_pheno"] = pairs.gene_symbol + "_" + pairs.code
    pairs = pairs.drop_duplicates("gene_pheno")

    # gene-level results per burden set
    for k, s in SETS.items():
        vals = {}
        for aid in pairs.analysis_id.unique():
            g = {x["gene_id"]: x for x in jload(a.raw / f"gene_burden/{aid}__{s}.json.gz")}
            for r in pairs[pairs.analysis_id == aid].itertuples():
                x = g.get(r.gene_id)
                vals[r.gene_pheno] = (x["Pvalue"], x["BETA_Burden"]) if x else (np.nan, np.nan)
        pairs[f"pval_{k}"] = pairs.gene_pheno.map(lambda gp: vals[gp][0])
        pairs[f"beta_{k}"] = pairs.gene_pheno.map(lambda gp: vals[gp][1])
    hits = pd.DataFrame({"gene_pheno": pairs.gene_pheno, "Description": pairs.analysis_id.map(desc),
                         "Phenotype": pairs.code.astype(int), "Gene": pairs.gene_symbol,
                         **{c: pairs[c] for c in ["pval_plof", "pval_mis", "pval_syn", "beta_plof", "beta_mis",
                                                  "beta_syn"]},
                         "analysis_id": pairs.analysis_id, "ensg": pairs.gene_id})
    hits.to_csv(a.out / "genebass_api_hits.csv", index=False)

    # variants
    paper_b = pd.read_csv(a.fig6 / "GeneBass_BiochemAssays_top_hits_all_variants_betas.csv.gz")
    group = paper_b.drop_duplicates("phenotype").set_index("phenotype").biomarker_group.to_dict()
    rows = []
    for r in hits.itertuples():
        va = pd.DataFrame(cached(a.raw, a.fetch, f"variant_assoc/{r.analysis_id}/{r.ensg}.json.gz",
                                 f"gene/{r.ensg}/variant-associations/{r.analysis_id}"))
        an = pd.DataFrame(cached(a.raw, a.fetch, f"variant_annotation/{r.ensg}.json.gz",
                                 f"gene/{r.ensg}/variant-qc"))
        if va.empty or an.empty:
            print(f"WARNING no variant data for {r.gene_pheno}", flush=True)
            continue
        an = an[an.consequence == "missense_variant"]
        v = va.rename(columns={"BETA": "beta", "Pvalue": "pval", "AC": "assoc_ac", "AF": "assoc_af"}).merge(
            an[["variant_id", "consequence", "hgvsp", "allele_count", "allele_number", "allele_frequency",
                "homozygote_count"]], on="variant_id")
        v = v.assign(gene_pheno=r.gene_pheno, gene=r.Gene, phenotype=r.Description, code=r.Phenotype, ensg=r.ensg,
                     biomarker_group=group.get(r.Description, "unassigned"))
        rows.append(v)
    bet = pd.concat(rows, ignore_index=True)
    bet.to_csv(a.out / "genebass_api_betas.csv.gz", index=False)

    # comparison with the paper's tables
    ph = pd.read_csv(a.fig6 / "GeneBass_BiochemAssays_top_hits.csv")
    P, A = set(ph.gene_pheno), set(hits.gene_pheno)
    m = paper_b.merge(bet, on=["gene_pheno", "variant_id"], how="outer", suffixes=("_paper", "_api"),
                      indicator=True)
    both = m[m._merge == "both"]
    gm = ph.merge(hits, on="gene_pheno", suffixes=("_paper", "_api"))
    gm = gm.dropna(subset=["pval_mis_paper"])

    def sel(b, h):
        b = b[b.gene_pheno.isin(h[h.pval_mis <= 1e-10].gene_pheno) & (b.pval < 0.05)]
        n = b.groupby("gene_pheno").size()
        return b[b.gene_pheno.isin(n[n > 25].index)]

    f_p, f_a = sel(paper_b, ph), sel(bet, hits)
    rep = [
        "# Genebass API vs VESM paper tables", "",
        f"- pairs: paper {len(P)}, API {len(A)}, common {len(P & A)}; only paper {sorted(P - A)}; "
        f"only API {sorted(A - P)}",
        f"- gene-level missense SKAT-O P (pairs where the paper has it, n={len(gm)}): max |rel diff| "
        f"{(np.abs(gm.pval_mis_api - gm.pval_mis_paper) / gm.pval_mis_paper).max():.2e}",
        f"- variant rows: paper {len(paper_b)}, API missense {len(bet)}, common {len(both)}, "
        f"only paper {(m._merge == 'left_only').sum()}, only API {(m._merge == 'right_only').sum()}",
        f"- common rows: beta identical {np.isclose(both.beta_paper, both.beta_api, rtol=1e-9, atol=0).sum()}, "
        f"P identical {np.isclose(both.pval_paper, both.pval_api, rtol=1e-9, atol=0).sum()}, "
        f"paper allele_count == API annotation AC {(both.allele_count_paper == both.allele_count_api).sum()}, "
        f"== API test AC {(both.allele_count_paper == both.assoc_ac).sum()}",
        f"- Fig. 6b selection: paper {f_p.gene_pheno.nunique()} pairs / {len(f_p)} rows; "
        f"API {f_a.gene_pheno.nunique()} pairs / {len(f_a)} rows; pairs only paper "
        f"{sorted(set(f_p.gene_pheno) - set(f_a.gene_pheno))}, only API {sorted(set(f_a.gene_pheno) - set(f_p.gene_pheno))}",
        f"- unassigned biomarker groups (phenotypes not in the paper): "
        f"{sorted(bet.loc[bet.biomarker_group == 'unassigned', 'phenotype'].unique())}",
    ]
    (a.out / "compare_vs_paper.md").write_text("\n".join(rep) + "\n")
    print("\n".join(rep), flush=True)
    print("GENEBASS_API OK")


if __name__ == "__main__":
    main()
