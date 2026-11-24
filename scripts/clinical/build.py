"""Build one clinical benchmark as PG2 archives, one per dataset, from prepare.py's variants table.

Every archive <DATASET>.splits.pgdata is laid out like the DMS archives:
  assays/<DATASET>.csv             sequence (full-length mutated sequence), the target, then mutant (UniProt
                                   numbering) and the source's own columns. No model scores
  sequences/<DATASET>.fasta        wild type: Swiss-Prot 2026_03 canonical sequence
  structures/<DATASET>.pdb         AlphaFold DB v6 model (proteins > 2,700 aa: fragments <DATASET>_F<k>.pdb, with
                                   their residue span in the manifest)
  msas/<DATASET>_colabfold.a3m     ColabFold search on UniRef100 2026_04, no post-alignment filtering (msa.sh);
                                   the search output byte for byte, so the hit headers (bit scores) are kept
  split                            trivial 2-fold (zero-shot models ignore it)
Also writes the index CSV, one row per archive.

Run: python build.py <benchmark> --variants <variants.tsv.gz> --out <dir> --index <csv> --sprot ... --afdb ... --msas ...
"""
import argparse
import shutil
import tempfile
import zipfile
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from proteingym.base import Dataset, Manifest
from proteingym.base.splits import KFoldSplitter

from proteingym.base.publication import DoiField

from common import a3m_depth, afdb_models, mutate, read_sprot

DoiField.resolve = lambda self, doi: {}  # publications are given in full below: no DOI lookup, the build is offline

MSA_DESC = ("MSA generated with ColabFold (MMseqs2) on UniRef100 (UniProt 2026_04) with all post-alignment "
            "filtering disabled.")
SPEC = {  # target, its description [, unit]; dataset description; per-dataset variables; publication
    "clinvar": dict(
        target="label", target_desc="1 = pathogenic / likely pathogenic, 0 = benign / likely benign",
        desc="ClinVar germline missense variants of {gene} (November 2026 release, >= 2 review stars): "
             "pathogenic / likely pathogenic vs benign / likely benign",
        pub=dict(title="ClinVar: improving access to variant interpretations and supporting evidence",
                 author="Landrum, Melissa J. and others", journal="Nucleic Acids Research", volume="46",
                 number="D1", year="2018", pages="D1062-D1067", doi="10.1093/nar/gkx1153")),
    "denovo": dict(
        target="dnm_case", target_desc="1 = de novo variant seen in severe developmental disorder cases, "
                                       "0 = seen only in controls",
        desc="De novo missense variants of {gene} in severe developmental disorder cases and controls "
             "(Orenbuch et al. 2025, Suppl. Table 5, ap_test_set)",
        pub=dict(title="Proteome-wide model for human disease genetics", author="Orenbuch, Rose and others",
                 journal="Nature Genetics", volume="57", number="12", year="2025", pages="3165-3174",
                 doi="10.1038/s41588-025-02400-1")),
    "gnomad": dict(
        target="log10_af", unit="log10",
        target_desc="log10 allele frequency in gnomAD v4.1 exomes excluding UK Biobank (rare: AF < 1e-3)",
        desc="gnomAD v4.1 rare missense variants of {gene} (exomes, PASS, MANE Select, non-UKB AF < 0.1%)",
        pub=dict(title="A genomic mutational constraint map using variation in 76,156 human genomes",
                 author="Chen, Siwei and others", journal="Nature", volume="625", number="7993", year="2024",
                 pages="92-100", doi="10.1038/s41586-023-06045-0")),
    "genebass": dict(
        target="beta", unit="SD",
        target_desc="Genebass single-variant effect size (SAIGE, inverse-rank-normalised trait); sign = direction "
                    "on the trait, not on fitness",
        desc="Genebass (UK Biobank exomes) single-variant associations of {gene} missense variants with "
             "{phenotype} (field {phenotype_code}), variants with P < 0.05; VESM Fig. 6b setting",
        variables={"phenotype": "Genebass phenotype", "phenotype_code": "UK Biobank field code",
                   "biomarker_group": "VESM Fig. 6b biomarker group", "genebass_analysis_id": "Genebass analysis",
                   "skato_pval_missense": "gene-level missense SKAT-O P value",
                   "pair_in_vesm_paper": "gene-phenotype pair is in the VESM paper's Fig. 6b"},
        pub=dict(title="Systematic single-variant and gene-based association testing of thousands of phenotypes "
                       "in 394,841 UK Biobank exomes", author="Karczewski, Konrad J. and others",
                 journal="Cell Genomics", volume="2", number="9", year="2022", pages="100168",
                 doi="10.1016/j.xgen.2022.100168")),
}


def q(x):
    if isinstance(x, (bool, np.bool_)):
        return str(bool(x)).lower()
    if isinstance(x, (int, float, np.number)):
        return str(x)
    return '"' + str(x).replace("\\", "\\\\").replace('"', '\\"') + '"'


def toml(table, items):
    return [table, *[f"{k} = {q(v)}" for k, v in items.items()], ""]


def build(job):
    name, v, a, spec, sp = job
    u, gene = v.uniprot_id.iloc[0], v.gene.iloc[0]
    assert v.uniprot_id.nunique() == v.gene.nunique() == 1, name
    wt, entry = sp[u]["sequence"], sp[u]["entry_name"]
    variables = {k: v[k].iloc[0] for k in spec.get("variables", {})}
    target, non_targets = spec["target"], [c for c in v.columns if c not in
                                           {"dataset", "gene", "uniprot_id", spec["target"], *variables}
                                           and not c.startswith("score_")]
    models = afdb_models(u, wt, a.afdb)
    msa = a.msas / f"{u}.a3m"
    depth = a3m_depth(msa, wt) if msa.exists() else None
    if not models or not depth:
        return f"{name} ({u}): " + ("no AlphaFold DB model matching the wild type" if not models else
                                    f"no MSA with the wild type as first record ({msa})")

    with tempfile.TemporaryDirectory(dir=a.out) as tmp:
        d = Path(tmp)
        (d / "assays").mkdir()
        assay = pd.DataFrame({"sequence": [mutate(wt, m) for m in v.mutant], target: v[target].to_numpy()})
        for c in non_targets:
            assay[c] = v[c].to_numpy()
        assay.to_csv(d / "assays" / f"{name}.csv", index=False)
        (d / f"{name}.fasta").write_text(f">{name}\n{wt}\n")
        (d / "stub.a3m").write_text(f">{u}\n{wt}\n")      # the real MSA is streamed into the archive below
        desc = spec["desc"].format(gene=gene, **variables)
        tdesc = {"name": target, "description": spec["target_desc"], **({"unit": spec["unit"]} if "unit" in spec
                                                                         else {})}
        m = toml("", dict(name=name, description=desc, version="1.0.0"))[1:]
        for k, kd in {"gene": "HGNC gene symbol", **spec.get("variables", {})}.items():
            m += toml("[[assay_variables]]", dict(name=k, description=kd))
        m += toml("[[assay_targets]]", tdesc)
        m += toml("[[assays]]", dict(name=name, path=f"assays/{name}.csv", description=desc, sequence_alphabet="AA",
                                     target_phenotype="other", number_of_variants=len(assay)))
        m += toml("[[assays.targets]]", dict(name=target))
        for c in non_targets:
            m += toml("[[assays.non_targets]]", dict(name=c))
        m += toml("[assays.variables]", {"gene": gene, **variables})
        m += ["[[sequences]]", 'type = "wild_type"', 'alphabet = "AA"', f'path = "{name}.fasta"',
              f'uniprot_id = "{u}"', f"pfam_ids = [{', '.join(q(x) for x in sp[u]['pfam_ids'])}]",
              'taxon_root = "Eukaryota"', f'molecule_name = "{entry}"', 'organism = "Homo sapiens"', ""]
        for i, (pdb, start, end, plddt) in enumerate(models, 1):
            sname = name if len(models) == 1 else f"{name}_F{i}"
            shutil.copy(pdb, d / f"{sname}.pdb")
            m += toml("[[structures]]", dict(path=f"{sname}.pdb", name=sname,
                                             description=f"AlphaFold structure prediction for {entry}"))
            m += toml("[structures.metadata]", dict(afdb_entry_id=pdb.stem.removesuffix("-model_v6"),
                                                    afdb_version="v6", mean_plddt=str(plddt),
                                                    span_start=str(start), span_end=str(end)))
        m += toml("[[msas]]", dict(path="stub.a3m", name=f"{name}_colabfold", format="a3m",
                                   reference_sequence_name=name, description=MSA_DESC))
        m += toml("[publication]", spec["pub"])
        (d / "manifest.toml").write_text("\n".join(m))
        ds = Dataset.from_manifest(Manifest.from_path(d / "manifest.toml"))
        raw = KFoldSplitter(2).split(ds, targets=[target]).dump(path=d)
        # rewrite compressed (proteingym-base stores), with the search's a3m in place of the stub
        with zipfile.ZipFile(raw) as z:
            (d / "in.pgdata").write_bytes(z.read("dataset.pgdata"))
            slices = z.read("slices.json")
        with zipfile.ZipFile(d / "in.pgdata") as zi, zipfile.ZipFile(d / "dataset.pgdata", "w",
                                                                     zipfile.ZIP_DEFLATED) as zo:
            for i in zi.infolist():
                if i.filename.endswith(".a3m"):
                    zo.write(msa, i.filename)
                else:
                    zo.writestr(i.filename, zi.read(i))
        with zipfile.ZipFile(d / "out", "w", zipfile.ZIP_DEFLATED) as zo:
            zo.write(d / "dataset.pgdata", "dataset.pgdata", zipfile.ZIP_STORED)
            zo.writestr("slices.json", slices)
        (d / "out").replace(a.out / raw.name)
    print(f"{name}\t{len(assay):,} variants", flush=True)
    return dict(dataset=name, gene=gene, uniprot_id=u, entry_name=entry, **variables, seq_length=len(wt),
                n_variants=len(assay), n_unique_mutants=assay.mutant.nunique(), n_structures=len(models),
                mean_plddt=round(float(np.mean([x[3] for x in models])), 2), msa_n_sequences=depth,
                archive=raw.name)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("benchmark", choices=list(SPEC))
    ap.add_argument("--variants", type=Path, required=True, help="prepare.py output")
    ap.add_argument("--out", type=Path, required=True, help="folder for the archives")
    ap.add_argument("--index", type=Path, required=True, help="index CSV to write")
    ap.add_argument("--sprot", type=Path, required=True)
    ap.add_argument("--afdb", type=Path, required=True, help="folder with AF-<UNIPROT>-F<k>-model_v6.pdb")
    ap.add_argument("--msas", type=Path, required=True, help="folder with <UNIPROT>.a3m")
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--first", type=int, default=None, help="only the first N datasets (smoke test)")
    a = ap.parse_args()
    shutil.rmtree(a.out, ignore_errors=True)
    a.out.mkdir(parents=True)
    sp, spec = read_sprot(a.sprot), SPEC[a.benchmark]
    v = pd.read_csv(a.variants, sep="\t", low_memory=False)
    jobs = [(n, x, a, spec, {u: sp[u] for u in x.uniprot_id.unique()}) for n, x in v.groupby("dataset")][:a.first]
    with Pool(a.jobs) as pool:
        res = list(pool.imap_unordered(build, jobs))
    failed = [r for r in res if isinstance(r, str)]
    pd.DataFrame([r for r in res if isinstance(r, dict)]).sort_values("dataset").to_csv(a.index, index=False)
    print(f"{a.benchmark}: {len(res) - len(failed)} archives in {a.out}, index {a.index}")
    if failed:
        raise SystemExit("NOT BUILT:\n  " + "\n  ".join(sorted(failed)))


if __name__ == "__main__":
    main()
