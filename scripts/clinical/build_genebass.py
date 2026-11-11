"""Build the Genebass clinical datasets (VESM paper Fig. 6b setting) as ProteinGym archives.

One archive per gene–phenotype pair, as PG2 DMS do one archive per (protein, condition), e.g.
APOB_HUMAN_GENEBASS_2022_LDL_DIRECT.splits.pgdata. Each archive holds:
  assay      the Fig. 6b variants of the pair: one row per (variant, protein change), target `beta`
             (Genebass single-variant effect), non-targets mutant / pval / allele counts / variant ids
  variables  phenotype_code, biomarker_group (Lipid Metabolism | Liver Function | Renal Function | Other)
  sequence   wild type = the UniProt sequence the paper mapped variants onto (compared with Swiss-Prot 2026_03)
  structures AlphaFold DB v6 model(s); proteins > 2,700 aa come as AFDB fragments F1..Fn, each listed with its
             residue span in `metadata` (span_start/span_end, 1-based on the wild type)
  (no MSA yet: alignment-based models are out of scope for now)

Fig. 6b selection, exactly as DJZN25-Fig6.ipynb cells 23-24 (on all variant rows of the betas table):
  1. pairs with missense SKAT-O P <= 1e-10 (Genebass gene-level result, top-hits table)
  2. variants with single-variant P < 0.05
  3. pairs with > 25 variants left
-> 103 pairs, 61 proteins, 6,977 variant rows; every row is then mapped onto a UniProt wild type:
  - the paper's own mapping (scores table: uniprot_id, uniprot_mut, Sequence), wild-type residue checked;
  - otherwise (UGT1A5 only: the paper gives no UniProt mapping, its hgvsp uses ENSP00000303174, 2 residues
    shorter at the N-terminus) the gene's Swiss-Prot entry, at the single position offset under which every
    variant's wild-type residue matches. Recorded per row in the `mapping` non-target.
Evaluation (not here): per pair, Pearson r(beta, score) and its two-sided P; Fig. 6b top = mean -log10 P over all
103 pairs, bottom = mean per biomarker group (the paper's heatmap shows Lipid / Liver / Renal; the 16 "Other"
pairs count only in the overall mean).

Inputs (all public):
  --fig6     VESM paper figure data, Fig6_data/ (figshare 10.6084/m9.figshare.31162471): top-hits table,
             single-variant betas, variant -> UniProt mapping
  --sprot    uniprot_sprot_human.dat.gz (entry name, Pfam, sequence check)
  --afdb     folder with AF-<UNIPROT>-F<k>-model_v6.pdb (AlphaFold DB v6 human proteome, UP000005640)
Outputs (--out, default datasets/clinical/genebass):
  genebass_fig6b.csv              index: one row per archive (pair metadata, gene-level stats, counts, structures)
  splits/<name>.splits.pgdata     the archives (one trivial 2-fold split; zero-shot models ignore it)
  build/                          staged manifests, CSVs, per-protein inputs (rebuildable, not committed)

Run (Python env with proteingym-base):
  python scripts/clinical/build_genebass.py --fig6 <Fig6_data> --sprot <dat.gz> --afdb <dir>
"""
import argparse
import gzip
import re
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
from proteingym.base import Dataset, Manifest
from proteingym.base.splits import KFoldSplitter

PVAL_MIS, PVAL_VAR, NVAR = 1e-10, 0.05, 25
AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H",
       "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V"}
NON_TARGETS = ["mutant", "pval", "allele_count", "allele_number", "allele_frequency", "variant_id", "hgvsp",
               "mapping"]
GENEBASS_DOI = "10.1016/j.xgen.2022.100168"
VESM_DOI = "10.1101/2025.03.09.642148"
AA3P = {k.capitalize(): v for k, v in AA3.items()}


def read_sprot(path):
    """Swiss-Prot human: accession -> dict(entry_name, gene, pfam_ids, sequence), and gene symbol -> accession."""
    out, rec = {}, None
    with gzip.open(path, "rt") as f:
        for line in f:
            tag = line[:2]
            if tag == "ID":
                rec = dict(entry_name=line.split()[1], acc=None, gene=None, pfam_ids=[], seq=[], in_seq=False)
            elif tag == "AC" and rec["acc"] is None:
                rec["acc"] = line[5:].split(";")[0].strip()
            elif tag == "GN" and rec["gene"] is None and (g := re.search(r"Name=([^;{ ]+)", line)):
                rec["gene"] = g[1]
            elif tag == "DR" and line[5:].startswith("Pfam;"):
                rec["pfam_ids"].append(line[5:].split(";")[1].strip())
            elif tag == "SQ":
                rec["in_seq"] = True
            elif tag == "//":
                out[rec["acc"]] = dict(entry_name=rec["entry_name"], gene=rec["gene"], pfam_ids=rec["pfam_ids"],
                                       sequence="".join(rec["seq"]))
            elif rec is not None and rec["in_seq"]:
                rec["seq"].append(line.replace(" ", "").strip())
    return out, {v["gene"]: k for k, v in out.items() if v["gene"]}


def map_by_offset(hgvsp, wt):
    """Map ENSP p.Xaa123Yaa changes onto `wt` with one constant offset under which every wild-type residue matches.
    Returns (offset, [X(123+offset)Y, ...]) or None."""
    ch = [(AA3P[g[1]], int(g[2]), AA3P[g[3]]) for h in hgvsp
          if (g := re.search(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$", h))]
    if len(ch) < len(hgvsp):
        return None
    for o in sorted(range(-100, 101), key=abs):
        if all(0 < p + o <= len(wt) and wt[p + o - 1] == a for a, p, _ in ch):
            return o, [f"{a}{p + o}{b}" for a, p, b in ch]
    return None


def pdb_info(path):
    """(sequence from CA atoms, mean pLDDT from the CA B-factors)."""
    seq, b = [], []
    for line in open(path):
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            seq.append(AA3.get(line[17:20], "X"))
            b.append(float(line[60:66]))
    return "".join(seq), float(np.mean(b))


def structures_for(u, wt, afdb):
    """AFDB v6 fragments of `u`, each placed on the wild type at its AFDB offset (fragment k starts at 200*(k-1)+1)."""
    files = sorted(afdb.glob(f"AF-{u}-F*-model_v6.pdb"), key=lambda p: int(re.search(r"-F(\d+)-", p.name)[1]))
    out = []
    for p in files:
        k = int(re.search(r"-F(\d+)-", p.name)[1])
        s, plddt = pdb_info(p)
        start = 1 if len(files) == 1 else 200 * (k - 1) + 1
        ref = wt[start - 1:start - 1 + len(s)]
        mism = [f"{a}{start + i}{b}" for i, (a, b) in enumerate(zip(ref, s)) if a != b]
        assert len(ref) == len(s) and len(mism) <= 0.01 * len(s), (u, k, len(ref), len(s), len(mism))
        out.append(dict(path=p, fragment=k, n_fragments=len(files), span_start=start, span_end=start + len(s) - 1,
                        n_mismatch=len(mism), mismatches=",".join(mism) or "none", mean_plddt=round(plddt, 2)))
    return out


MSA_DESC = ("MSA generated with ColabFold (MMseqs2) on UniRef100 (UniProt 2026_01) with all post-alignment "
            "filtering disabled.")


def a3m_check(path, wt):
    """First a3m record must be the wild type (match states = upper case + '-'); returns the depth."""
    n, first, cur = 0, None, []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                n += 1
                if n == 2:
                    first = "".join(cur)
                cur = []
            elif n == 1:
                cur.append(line.strip())
    first = first if first is not None else "".join(cur)
    assert first.replace("-", "") == wt and first == first.upper(), f"{path}: first record is not the wild type"
    return dict(n_sequences=n)


def slug(x):
    return re.sub(r"[^A-Z0-9]+", "_", x.upper()).strip("_")


def mutate(s, m):
    p = int(m[1:-1])
    assert s[p - 1] == m[0]
    return s[:p - 1] + m[-1] + s[p:]


def toml_str(x):
    return '"' + str(x).replace("\\", "\\\\").replace('"', '\\"') + '"'


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fig6", type=Path, required=True)
    ap.add_argument("--hits", type=Path, default=None,
                    help="gene-phenotype table (default: paper top hits; genebass_api.py writes one from Genebass)")
    ap.add_argument("--betas", type=Path, default=None,
                    help="single-variant table (default: paper betas; genebass_api.py writes one from Genebass)")
    ap.add_argument("--sprot", type=Path, required=True)
    ap.add_argument("--afdb", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=Path("datasets/clinical/genebass"))
    ap.add_argument("--msa-dir", type=Path, default=None,
                    help="<UNIPROT>.a3m per protein (scripts/clinical/msa: ColabFold on UniRef100 2026_01, no "
                         "post-alignment filtering); without it the archives carry no MSA")
    a = ap.parse_args()
    build, splits = a.out / "build", a.out / "splits"
    for d in (build, splits):
        shutil.rmtree(d, ignore_errors=True)
        d.mkdir(parents=True)

    # --- raw tables ---------------------------------------------------------------------------------------------
    hits = pd.read_csv(a.hits or a.fig6 / "GeneBass_BiochemAssays_top_hits.csv")
    bet = pd.read_csv(a.betas or a.fig6 / "GeneBass_BiochemAssays_top_hits_all_variants_betas.csv.gz")
    sc = pd.read_csv(a.fig6 / "GeneBass_BiochemAssays_top_hits_all_variants_scores.csv.gz",
                     usecols=["variant_id", "hgvsp", "Gene", "uniprot_id", "uniprot_mut", "Sequence"])
    sc = sc.rename(columns={"Gene": "gene"}).drop_duplicates(["variant_id", "hgvsp", "gene"])
    seqs = sc.dropna(subset=["uniprot_id"]).groupby("uniprot_id").Sequence.agg(set)
    assert seqs.map(len).max() == 1
    seq = seqs.map(lambda s: next(iter(s))).to_dict()
    sp, gene2acc = read_sprot(a.sprot)

    # --- 1-3. Fig. 6b selection (paper cells 23-24, on all rows) -------------------------------------------------
    pairs = hits.rename(columns={"Gene": "gene", "Description": "phenotype", "Phenotype": "phenotype_code"})
    var = bet[bet.gene_pheno.isin(pairs.loc[pairs.pval_mis <= PVAL_MIS, "gene_pheno"])]
    var = var[var.pval < PVAL_VAR]
    n = var.groupby("gene_pheno").hgvsp.count()
    var = var[var.gene_pheno.isin(n[n > NVAR].index)].copy()

    # --- map every selected row onto a UniProt wild type ---------------------------------------------------------
    var = var.merge(sc[["variant_id", "hgvsp", "gene", "uniprot_id", "uniprot_mut"]],
                    on=["variant_id", "hgvsp", "gene"], how="left", validate="many_to_one")

    def on_wt(u, m):
        x = re.fullmatch(r"([A-Z])(\d+)([A-Z])", m) if isinstance(m, str) else None
        return bool(x and u in seq and int(x[2]) <= len(seq[u]) and seq[u][int(x[2]) - 1] == x[1] and x[1] != x[3])

    ok = pd.Series([on_wt(u, m) for u, m in zip(var.uniprot_id, var.uniprot_mut)], index=var.index)
    var["mapping"] = np.where(ok, "paper_uniprot", "")
    unmapped = []
    for g, rows in var[~ok].groupby("gene"):
        acc = gene2acc.get(g)
        ref = seq.get(acc, sp[acc]["sequence"]) if acc else None      # keep the paper wild type if it has one
        hit = map_by_offset(rows.hgvsp.tolist(), ref) if acc else None
        if not hit:                                                     # cannot be placed on a UniProt wild type
            print(f"WARNING {g}: {len(rows)} rows cannot be mapped onto Swiss-Prot ({acc}); dropped "
                  f"(listed in build/unmapped_rows.tsv)", flush=True)
            unmapped.append(rows)
            continue
        o, muts = hit
        seq.setdefault(acc, ref)
        var.loc[rows.index, ["uniprot_id", "uniprot_mut", "mapping"]] = [[acc, m, f"swissprot_offset{o:+d}"]
                                                                         for m in muts]
        print(f"{g}: {len(rows)} rows without a paper UniProt mapping -> {acc} ({sp[acc]['entry_name']}), "
              f"ENSP position {o:+d}, all wild-type residues match")
    if unmapped:
        um = pd.concat(unmapped)
        um.to_csv(build / "unmapped_rows.tsv", sep="\t", index=False)
        var = var.drop(um.index)
    var = var.rename(columns={"uniprot_mut": "mutant"})
    pairs = pairs[pairs.gene_pheno.isin(var.gene_pheno)].merge(
        var[["gene_pheno", "uniprot_id", "biomarker_group"]].drop_duplicates(), on="gene_pheno", validate="1:1")
    assert var.groupby("gene_pheno").uniprot_id.nunique().max() == 1

    # --- per-protein input store: wild type + structures -----------------------------------------------------------
    prot = {}
    for u in sorted(pairs.uniprot_id.unique()):
        d = build / "proteins" / u
        d.mkdir(parents=True)
        wt, info = seq[u], sp[u]
        (d / f"{u}.fasta").write_text(f">{u} {info['entry_name']}\n{wt}\n")
        st = structures_for(u, wt, a.afdb)
        assert st, f"no AFDB structure for {u}"
        for s in st:
            shutil.copy(s["path"], d / s["path"].name)
        msa = None
        if a.msa_dir and (a.msa_dir / f"{u}.a3m").exists():
            src = a.msa_dir / f"{u}.a3m"
            msa = dict(path=d / f"{u}_colabfold_uniref100.a3m", **a3m_check(src, wt))
            shutil.copy(src, msa["path"])
        elif a.msa_dir:
            print(f"WARNING no MSA for {u} in {a.msa_dir}: archive built without MSA", flush=True)
        prot[u] = dict(info, structures=st, msa=msa, sprot_identical=(info["sequence"] == wt), length=len(wt))

    # --- one archive per pair ------------------------------------------------------------------------------------
    rows = []
    for p in pairs.sort_values(["gene", "phenotype"]).itertuples():
        u, pr = p.uniprot_id, prot[p.uniprot_id]
        name = f"{pr['entry_name']}_GENEBASS_2022_{slug(p.phenotype)}"
        v = var[var.gene_pheno == p.gene_pheno].copy()
        v["mutant_pos"] = v.mutant.str[1:-1].astype(int)
        v = v.sort_values(["mutant_pos", "mutant", "variant_id"])
        d = build / "datasets" / name
        (d / "assays").mkdir(parents=True)
        assay = pd.DataFrame({"sequence": [mutate(seq[u], m) for m in v.mutant], "beta": v.beta.to_numpy()})
        for c in NON_TARGETS:
            assay[c] = v[c].to_numpy()
        assay.to_csv(d / "assays" / f"{name}.csv", index=False)
        group = p.biomarker_group
        m = [f"name = {toml_str(name)}",
             f"description = {toml_str(f'Genebass (UK Biobank exomes) single-variant associations of {p.gene} missense variants with {p.phenotype} (field {p.phenotype_code}); VESM Fig. 6b setting')}",
             'version = "0.1.0"', "",
             "[[assay_variables]]", 'name = "phenotype_code"', 'description = "UK Biobank field code"', "",
             "[[assay_variables]]", 'name = "biomarker_group"',
             'description = "VESM Fig. 6b biomarker group: Lipid Metabolism, Liver Function, Renal Function, Other"',
             "",
             "[[assay_targets]]", 'name = "beta"',
             'description = "Genebass single-variant effect size (SAIGE, inverse-rank-normalised trait); sign = '
             'direction on the trait, not on fitness"', 'unit = "SD"', "",
             "[[assays]]", f"name = {toml_str(name)}", f'path = "assays/{name}.csv"',
             f"description = {toml_str(f'Genebass {p.analysis_id}; gene-level missense SKAT-O P = {p.pval_mis:.3g}; variants with single-variant P < 0.05')}",
             'sequence_alphabet = "AA"', 'target_phenotype = "other"', f"number_of_variants = {len(assay)}", "",
             "[[assays.targets]]", 'name = "beta"', "",
             *[x for c in NON_TARGETS for x in ("[[assays.non_targets]]", f'name = "{c}"', "")],
             "[assays.variables]", f"phenotype_code = {int(p.phenotype_code)}",
             f"biomarker_group = {toml_str(group)}", "",
             "[[sequences]]", 'type = "wild_type"', 'alphabet = "AA"', f'path = "../../proteins/{u}/{u}.fasta"',
             f'uniprot_id = "{u}"', f"pfam_ids = [{', '.join(toml_str(x) for x in pr['pfam_ids'])}]",
             'taxon_root = "Eukaryota"', f'molecule_name = "{pr["entry_name"]}"', 'organism = "Homo sapiens"', ""]
        for s in pr["structures"]:
            sdesc = (f"AlphaFold DB v6 model, fragment {s['fragment']} of {s['n_fragments']}, "
                     f"wild-type residues {s['span_start']}-{s['span_end']}")
            m += ["[[structures]]", f'path = "../../proteins/{u}/{s["path"].name}"',
                  f'name = "{s["path"].stem}"', f"description = {toml_str(sdesc)}", "",
                  "[structures.metadata]", f'afdb_entry_id = "{s["path"].stem.removesuffix("-model_v6")}"',
                  'afdb_version = "v6"',
                  *[f'{k} = "{s[k]}"' for k in ["fragment", "n_fragments", "span_start", "span_end", "n_mismatch",
                                               "mismatches", "mean_plddt"]], ""]
        if pr["msa"]:
            mp = pr["msa"]
            m += ["[[msas]]", f'path = "../../proteins/{u}/{mp["path"].name}"', f'name = "{name}_colabfold"',
                  'format = "a3m"', f'reference_sequence_name = "{u}"',
                  f"description = {toml_str(MSA_DESC)}", "",
                  "[msas.metadata]", f'n_sequences = "{mp["n_sequences"]}"', 'database = "UniRef100 2026_01"',
                  'method = "colabfold_search monomer, filter 0, no cluster expansion, no envdb"', ""]
        m += ["[publication]", 'title = "Systematic single-variant and gene-based association testing of thousands of phenotypes in 394,841 UK Biobank exomes"',
              'author = "Karczewski, Konrad J. and others"', 'journal = "Cell Genomics"', 'volume = "2"', 'number = "9"', 'year = "2022"',
              f'doi = "{GENEBASS_DOI}"', ""]
        (d / "manifest.toml").write_text("\n".join(m))

        ds = Dataset.from_manifest(Manifest.from_path(d / "manifest.toml"))
        out = KFoldSplitter(2).split(ds, targets=["beta"]).dump(path=splits)
        rows.append(dict(dataset=name, gene=p.gene, uniprot_id=u, entry_name=pr["entry_name"],
                         phenotype=p.phenotype, phenotype_code=int(p.phenotype_code), biomarker_group=group,
                         genebass_analysis_id=p.analysis_id, gene_pheno=p.gene_pheno,
                         skato_pval_missense=p.pval_mis, burden_beta_missense=p.beta_mis,
                         skato_pval_plof=p.pval_plof, burden_beta_plof=p.beta_plof, seq_length=pr["length"],
                         n_variants=len(assay), n_unique_mutants=assay["mutant"].nunique(),
                         n_structures=len(pr["structures"]),
                         msa_n_sequences=pr["msa"]["n_sequences"] if pr["msa"] else None,
                         mean_plddt=round(np.mean([s["mean_plddt"] for s in pr["structures"]]), 2),
                         wt_equals_sprot_2026_03=pr["sprot_identical"], archive=f"splits/{out.name}"))
        print(f"{name:55s} {len(assay):4d} variants  {group}", flush=True)

    idx = pd.DataFrame(rows)
    idx.to_csv(a.out / "genebass_fig6b.csv", index=False)
    print(f"\n{len(idx)} archives, {idx.uniprot_id.nunique()} proteins, {idx.n_variants.sum():,} variant rows; "
          f"groups {idx.biomarker_group.value_counts().to_dict()}; "
          f"wild type != Swiss-Prot 2026_03: {(~idx.drop_duplicates('uniprot_id').wt_equals_sprot_2026_03).sum()}")


if __name__ == "__main__":
    main()
