"""Select the variants of the four clinical benchmarks and place them on UniProt proteins.

  clinvar   ClinVar missense (Nov 2026): >= 2 review stars, germline, pathogenic / likely pathogenic vs benign /
            likely benign; genes with >= 10 of each (clinical-curation clinvar_v0)              target label
  denovo    de novo missense variants of the popEVE paper (Orenbuch et al. 2025, Suppl. Table 5, ap_test_set):
            severe developmental disorder cases vs controls; 300 random genes                   target dnm_case
  gnomad    gnomAD v4.1 exome rare missense (PASS, MANE Select, non-UKB AF < 0.1%, AN >= 50% of max);
            300 random genes                                                                    target log10_af
  genebass  Genebass single-variant associations with blood biomarkers, VESM Fig. 6b selection: gene-phenotype
            pairs with missense SKAT-O P <= 1e-10, variants with P < 0.05, pairs with > 25 variants;
            proteins <= 2,700 aa                                                                target beta

Placing on UniProt, the same for every benchmark: each variant's SNV is looked up in the AlphaMissense SNV table
(build_snv_index.py), which gives the protein change on a UniProt canonical isoform; it is kept when the wild-type
residue is the Swiss-Prot one. A gene's protein is the Swiss-Prot entry named after the gene if its SNVs map to it,
else the entry most of its SNVs map to. Variants that do not land on that protein are dropped (counted in the log).
Genebass first uses the VESM paper's own mapping, so the paper's variant set is kept.

Random genes: one draw of 300 (numpy default_rng(0), without replacement) from the sorted eligible genes, i.e.
genes whose protein is <= 2,700 aa and has at least 1 (denovo) or 20 (gnomad) placed variants.

Output: <out>/<benchmark>/variants.tsv.gz, one row per variant: dataset (archive name), gene, uniprot_id, mutant,
the target and the other assay columns. clinvar / denovo / gnomad also keep clinical-curation's reference model
scores (score_*), which check.py reports and build.py never archives.
"""
import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl

from common import MAX_LEN, read_sprot, snv_to_protein

N_GENES, MIN_GNOMAD, RARE = 300, 20, 1e-3
SUFFIX = {"clinvar": "CLINVAR_2026", "denovo": "DENOVO_2025", "gnomad": "GNOMAD_2024", "genebass": "GENEBASS_2022"}
SCORES = {"am": "score_alphamissense", "popeve": "score_popeve", "cadd": "score_cadd_phred", "esm1b": "score_esm1b",
          "vesm": "score_vesm3b", "gpn_msa": "score_gpn_msa"}


def vid(d):
    return d.chrom.astype(str).str.removeprefix("chr") + "-" + d.pos.astype(str) + "-" + d.ref + "-" + d.alt


def place(d, index, sp):
    """d: gene, variant_id (SNVs, ';'-joined if several). Adds uniprot_id and mutant; drops unplaced variants."""
    e = d[["gene"]].assign(row=d.index, variant_id=d.variant_id.str.split(";")).explode("variant_id")
    h = e.merge(snv_to_protein(index, e.variant_id, sp), on="variant_id").drop_duplicates(["row", "uniprot_id"])
    n = h.groupby(["gene", "uniprot_id"]).size().rename("n").reset_index()
    n["named"] = [sp[u]["gene"] == g for g, u in zip(n.gene, n.uniprot_id)]
    acc = n.sort_values(["named", "n", "uniprot_id"], ascending=[False, False, True]).drop_duplicates("gene")
    acc = acc.drop_duplicates("uniprot_id")           # a protein belongs to one gene (drops e.g. readthrough genes)
    h = h[h.uniprot_id == h.gene.map(acc.set_index("gene").uniprot_id)].set_index("row")
    out = d.drop(columns=["mutant"], errors="ignore").join(h[["uniprot_id", "mutant"]], how="inner")
    print(f"  placed {len(out):,} of {len(d):,} variants ({1 - len(out) / len(d):.2%} dropped), "
          f"{out.gene.nunique():,} of {d.gene.nunique():,} genes", flush=True)
    return out


def draw(d, sp, min_variants):
    """One draw of N_GENES genes from the eligible ones (protein <= MAX_LEN aa, >= min_variants placed variants)."""
    n = d.groupby("gene").size()
    length = d.drop_duplicates("gene").set_index("gene").uniprot_id.map(lambda u: len(sp[u]["sequence"]))
    eligible = sorted(n.index[(n >= min_variants) & (length[n.index] <= MAX_LEN)])
    genes = np.random.default_rng(0).choice(eligible, N_GENES, replace=False)
    print(f"  drew {N_GENES} of {len(eligible):,} eligible genes", flush=True)
    return d[d.gene.isin(set(genes))]


def clinvar(a, sp):
    b = a.curation / "data/processed/benchmarks"
    d = pd.read_csv(b / "clinvar_v0_scores.tsv", sep="\t").merge(
        pd.read_csv(b / "clinvar_v0_extra_scores.tsv", sep="\t"), on=["gene", "hgvs_p"], how="left", validate="1:1")
    snv = pd.read_csv(b / "clinvar_v0_snvs.tsv", sep="\t")                 # the SNVs causing each protein change
    snv = snv.assign(v=vid(snv)).groupby(["gene", "hgvs_p"]).v.agg(lambda s: ";".join(sorted(set(s))))
    d = d.merge(snv.rename("variant_id").reset_index(), on=["gene", "hgvs_p"], how="left", validate="1:1")
    d = place(d.rename(columns={"hgvs_p": "hgvsp"}), a.snv_index / "am_hg38.parquet", sp)
    return d, ["label", "variant_id", "hgvsp", "clnsig", "max_stars", "n_records"]


def denovo(a, sp):
    b = a.curation / "data/processed/benchmarks"
    k = ["chrom", "pos", "ref", "alt"]
    s = {"chrom": str, "pos": str}
    d = pd.read_csv(b / "dnm_aptest_scores.tsv", sep="\t", dtype=s).merge(
        pd.read_csv(b / "dnm_aptest_extra_scores.tsv", sep="\t", dtype=s), on=k, how="left", validate="m:1")
    d["variant_id"] = vid(d)                                               # the source coordinates are GRCh37
    d = draw(place(d, a.snv_index / "am_hg19.parquet", sp), sp, 1)
    # GRCh38 SNV of the same protein change, where unique
    g38 = pl.scan_parquet(a.snv_index / "am_hg38.parquet").filter(
        pl.col("uniprot_id").is_in(sorted(set(d.uniprot_id)))).collect().to_pandas()
    g38 = g38.assign(variant_id=vid(g38), mutant=g38.protein_variant).drop_duplicates(
        ["uniprot_id", "mutant", "ref", "alt"], keep=False)
    d = d.rename(columns={"variant_id": "variant_id_grch37"}).merge(
        g38[["uniprot_id", "mutant", "ref", "alt", "variant_id"]], on=["uniprot_id", "mutant", "ref", "alt"],
        how="left")
    return d, ["dnm_case", "variant_id", "variant_id_grch37", "n_sdd", "n_control"]


def gnomad(a, sp):
    b = a.curation / "data/interim"
    k = ["chrom", "pos", "ref", "alt"]
    d = pd.read_csv(b / "gnomad_missense_scores_all.tsv.gz", sep="\t",
                    usecols=["gene", *k, "AC_non_ukb", "AN_non_ukb", "am", "popeve", "cadd"])
    d = d[(d.AC_non_ukb > 0) & (d.AN_non_ukb >= 0.5 * d.AN_non_ukb.max())]
    d = d.assign(af=d.AC_non_ukb / d.AN_non_ukb)
    d = d[d.af < RARE].merge(pd.read_csv(b / "gnomad_missense_extra_scores_all.tsv.gz", sep="\t"), on=k, how="left")
    d = d.assign(log10_af=np.log10(d.af), variant_id=vid(d))
    d = draw(place(d, a.snv_index / "am_hg38.parquet", sp), sp, MIN_GNOMAD)
    return d, ["log10_af", "variant_id", "af", "AC_non_ukb", "AN_non_ukb"]


def genebass(a, sp):
    def fig6b(hits, betas):  # the VESM Fig. 6b selection (DJZN25-Fig6.ipynb, cells 23-24)
        v = betas[betas.gene_pheno.isin(hits.gene_pheno[hits.pval_mis <= 1e-10]) & (betas.pval < 0.05)]
        n = v.groupby("gene_pheno").size()
        return v[v.gene_pheno.isin(n.index[n > 25])].copy()

    hits = pd.read_csv(a.genebass / "genebass_api_hits.csv")               # Genebass's own tables (genebass_api.py)
    d = fig6b(hits, pd.read_csv(a.genebass / "genebass_api_betas.csv.gz"))
    # the paper's own tables: which pairs and rows are in its Fig. 6b, and its variant -> UniProt mapping
    p = "GeneBass_BiochemAssays_top_hits"
    paper = fig6b(pd.read_csv(a.fig6 / f"{p}.csv"), pd.read_csv(a.fig6 / f"{p}_all_variants_betas.csv.gz"))
    d["in_vesm_paper"] = pd.MultiIndex.from_frame(d[["gene_pheno", "variant_id"]]).isin(
        pd.MultiIndex.from_frame(paper[["gene_pheno", "variant_id"]]))
    m = pd.read_csv(a.fig6 / f"{p}_all_variants_scores.csv.gz",
                    usecols=["variant_id", "hgvsp", "Gene", "uniprot_id", "uniprot_mut"])
    m = m.rename(columns={"Gene": "gene", "uniprot_mut": "mutant"}).dropna().drop_duplicates(
        ["variant_id", "hgvsp", "gene"])
    ok = [u in sp and (x := re.fullmatch(r"([A-Z])(\d+)[A-Z]", w)) is not None
          and sp[u]["sequence"][int(x[2]) - 1:int(x[2])] == x[1] for u, w in zip(m.uniprot_id, m.mutant)]
    by_paper = d.merge(m[ok], on=["variant_id", "hgvsp", "gene"])
    rest = d.merge(by_paper[["gene_pheno", "variant_id", "hgvsp"]], how="left", indicator=True)
    rest = place(rest[rest._merge == "left_only"].drop(columns="_merge"), a.snv_index / "am_hg38.parquet", sp)
    rest = rest[rest.gene == rest.uniprot_id.map(lambda u: sp[u]["gene"])]
    protein = by_paper.drop_duplicates("gene").set_index("gene").uniprot_id    # one protein per gene: the paper's
    d = pd.concat([by_paper, rest[rest.uniprot_id == rest.gene.map(protein).fillna(rest.uniprot_id)]])
    d = d[d.uniprot_id.map(lambda u: len(sp[u]["sequence"])) <= MAX_LEN]
    d = d.merge(hits.rename(columns={"analysis_id": "genebass_analysis_id", "pval_mis": "skato_pval_missense"})[
        ["gene_pheno", "genebass_analysis_id", "skato_pval_missense"]], on="gene_pheno", validate="m:1")
    d = d.assign(phenotype_code=d.code, pair_in_vesm_paper=d.gene_pheno.isin(set(paper.gene_pheno)))
    d["pheno"] = "_" + d.phenotype.str.upper().str.replace(r"[^A-Z0-9]+", "_", regex=True).str.strip("_")
    return d, ["beta", "variant_id", "hgvsp", "pval", "allele_count", "allele_number", "allele_frequency",
               "in_vesm_paper", "phenotype", "phenotype_code", "biomarker_group", "genebass_analysis_id",
               "skato_pval_missense", "pair_in_vesm_paper"]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("benchmarks", nargs="*", default=list(SUFFIX))
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--curation", type=Path, required=True, help="clinical-curation repo")
    ap.add_argument("--genebass", type=Path, required=True, help="genebass_api.py output")
    ap.add_argument("--fig6", type=Path, required=True, help="VESM paper Fig6_data")
    ap.add_argument("--sprot", type=Path, required=True, help="uniprot_sprot_human.dat.gz")
    ap.add_argument("--snv-index", type=Path, required=True, help="folder with am_hg38.parquet, am_hg19.parquet")
    a = ap.parse_args()
    sp = read_sprot(a.sprot)
    for b in a.benchmarks:
        print(b, flush=True)
        d, cols = globals()[b](a, sp)
        d["dataset"] = d.uniprot_id.map(lambda u: sp[u]["entry_name"]) + "_" + SUFFIX[b] + d.get("pheno", "")
        d["pos"] = d.mutant.str[1:-1].astype(int)
        d = d.sort_values(["dataset", "pos", "mutant", "variant_id"]).rename(columns=SCORES)
        d = d[["dataset", "gene", "uniprot_id", "mutant", *cols, *[c for c in SCORES.values() if c in d]]]
        (a.out / b).mkdir(parents=True, exist_ok=True)
        d.to_csv(a.out / b / "variants.tsv.gz", sep="\t", index=False)
        print(f"  {d.dataset.nunique():,} datasets, {d.uniprot_id.nunique():,} proteins, {len(d):,} variants",
              flush=True)


if __name__ == "__main__":
    main()
