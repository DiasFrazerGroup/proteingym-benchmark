"""Build the ClinVar, de novo and gnomAD benchmarks as PG2 archives, one per gene, from prepare.py output.

Each archive <ENTRY>_<SOURCE>_<YEAR>.splits.pgdata (e.g. BRCA1_HUMAN_CLINVAR_2026) holds:
  assay      the gene's variants: full-length mutated `sequence`, the target (label | dnm_case | log10_af) and
             non-targets: mutant (UniProt numbering), mutant_source, variant_id (GRCh38), mapping, source fields.
             No model scores: the archives are benchmark data only
  sequence   Swiss-Prot 2026_03 canonical wild type
  structures AlphaFold DB v6 model(s), fragments F1..Fn above 2,700 aa (as build_genebass.py)
  msas       ColabFold on UniRef100 2026_01, no post-alignment filtering (scripts/clinical/msa)
  split      trivial 2-fold (zero-shot models ignore it)
A protein with no usable AFDB model or MSA is built without it (WARNING, and the index records it).

Outputs: <out>/<benchmark>/splits/*.splits.pgdata and the index datasets/clinical/<benchmark>/<benchmark>.csv
(one row per archive).
Run: bash scripts/clinical/benchmarks/submit.sh build
"""
import argparse
import io
import shutil
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from proteingym.base import Dataset, Manifest
from proteingym.base.splits import KFoldSplitter

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from build_genebass import MSA_DESC, a3m_check, mutate, read_sprot, structures_for, toml_str  # noqa: E402

COMMON = ["mutant", "mutant_source", "variant_id", "mapping"]
SPEC = {
    "clinvar": dict(
        suffix="CLINVAR_2026", target="label", unit=None,
        extra=["hgvsp", "source_transcript", "clnsig", "max_stars", "n_records"],
        target_desc="1 = pathogenic / likely pathogenic, 0 = benign / likely benign (ClinVar, >=2 review stars)",
        desc="ClinVar germline missense variants of {gene} (Nov 2026 release; >=2 review stars, pathogenic / "
             "likely pathogenic vs benign / likely benign)"),
    "denovo": dict(
        suffix="DENOVO_2025", target="dnm_case", unit=None,
        extra=["variant_id_grch37", "source_protein", "n_sdd", "n_control"],
        target_desc="1 = de novo variant seen in severe developmental disorder cases, 0 = only in controls",
        desc="De novo missense variants of {gene} in severe developmental disorder cases and controls "
             "(Orenbuch et al. 2025, Suppl. Table 5, ap_test_set)",
        doi="10.1038/s41588-025-02400-1"),
    "gnomad": dict(
        suffix="GNOMAD_2024", target="log10_af", unit="log10",
        extra=["af", "AC_non_ukb", "AN_non_ukb", "AC", "AN", "enst", "mane_select"],
        target_desc="log10 allele frequency in gnomAD v4.1 exomes excluding UK Biobank (rare: AF < 1e-3)",
        desc="gnomAD v4.1 rare missense variants of {gene} (exomes, PASS, MANE Select, non-UKB AF < 0.1%)"),
}


def protein_inputs(u, sp, store, afdb, msa_dirs, cache):
    """Stage wild type, AFDB structures and MSA of `u` once; tolerate missing or mismatching inputs."""
    if u in cache:
        return cache[u]
    d = store / u
    d.mkdir(parents=True, exist_ok=True)
    wt, info = sp[u]["sequence"], sp[u]
    (d / f"{u}.fasta").write_text(f">{u} {info['entry_name']}\n{wt}\n")
    note = []
    try:
        st = structures_for(u, wt, afdb)
    except AssertionError as e:
        st = []
        note.append(f"AFDB model does not match the wild type ({e})")
    if not st and not note:
        note.append("no AFDB model")
    for s in st:
        if not (d / s["path"].name).exists():
            shutil.copy(s["path"], d / s["path"].name)
    msa = None
    src = next((m / f"{u}.a3m" for m in msa_dirs if (m / f"{u}.a3m").exists()), None)
    if src:
        try:
            msa = dict(path=d / f"{u}_colabfold_uniref100.a3m", **a3m_check(src, wt))
            if not msa["path"].exists():
                shutil.copy(src, msa["path"])
        except AssertionError:
            msa = None
            note.append(f"MSA first record is not the wild type ({src})")
    elif msa_dirs:
        note.append("no MSA")
    for n in note:
        print(f"WARNING {u} ({info['entry_name']}): {n}; archive built without it", flush=True)
    cache[u] = dict(info, wt=wt, structures=st, msa=msa, note="; ".join(note))
    return cache[u]


def compress(path):
    """proteingym-base writes archives uncompressed (MSAs make them GBs): deflate the inner dataset.pgdata members
    and store it in the outer zip. Same content, readable by Subsets.from_path."""
    tmp = path.with_name(path.name + ".tmp")
    with zipfile.ZipFile(path) as z, zipfile.ZipFile(tmp, "w") as o:
        for i in z.infolist():
            data = z.read(i)
            if i.filename.endswith(".pgdata"):
                b = io.BytesIO()
                with zipfile.ZipFile(io.BytesIO(data)) as zi, zipfile.ZipFile(b, "w", zipfile.ZIP_DEFLATED,
                                                                              compresslevel=6) as zo:
                    for j in zi.infolist():
                        zo.writestr(j.filename, zi.read(j))
                data = b.getvalue()
                o.writestr(i.filename, data, compress_type=zipfile.ZIP_STORED)
            else:
                o.writestr(i.filename, data, compress_type=zipfile.ZIP_DEFLATED)
    tmp.replace(path)


def manifest(name, bench, gene, u, pr, n_rows, non_targets):
    sp = SPEC[bench]
    m = [f"name = {toml_str(name)}", f"description = {toml_str(sp['desc'].format(gene=gene))}",
         'version = "0.1.0"', "",
         "[[assay_variables]]", 'name = "gene"', 'description = "HGNC gene symbol in the source table"', "",
         "[[assay_targets]]", f'name = "{sp["target"]}"', f"description = {toml_str(sp['target_desc'])}",
         *([f'unit = "{sp["unit"]}"'] if sp["unit"] else []), "",
         "[[assays]]", f"name = {toml_str(name)}", f'path = "assays/{name}.csv"',
         f"description = {toml_str(sp['desc'].format(gene=gene))}",
         'sequence_alphabet = "AA"', 'target_phenotype = "other"', f"number_of_variants = {n_rows}", "",
         "[[assays.targets]]", f'name = "{sp["target"]}"', "",
         *[x for c in non_targets for x in ("[[assays.non_targets]]", f'name = "{c}"', "")],
         "[assays.variables]", f"gene = {toml_str(gene)}", "",
         "[[sequences]]", 'type = "wild_type"', 'alphabet = "AA"', f'path = "../../proteins/{u}/{u}.fasta"',
         f'uniprot_id = "{u}"', f"pfam_ids = [{', '.join(toml_str(x) for x in pr['pfam_ids'])}]",
         'taxon_root = "Eukaryota"', f'molecule_name = "{pr["entry_name"]}"', 'organism = "Homo sapiens"', ""]
    for st in pr["structures"]:
        sdesc = (f"AlphaFold DB v6 model, fragment {st['fragment']} of {st['n_fragments']}, "
                 f"wild-type residues {st['span_start']}-{st['span_end']}")
        m += ["[[structures]]", f'path = "../../proteins/{u}/{st["path"].name}"', f'name = "{st["path"].stem}"',
              f"description = {toml_str(sdesc)}", "", "[structures.metadata]",
              f'afdb_entry_id = "{st["path"].stem.removesuffix("-model_v6")}"', 'afdb_version = "v6"',
              *[f'{k} = "{st[k]}"' for k in ["fragment", "n_fragments", "span_start", "span_end", "n_mismatch",
                                            "mismatches", "mean_plddt"]], ""]
    if pr["msa"]:
        mp = pr["msa"]
        m += ["[[msas]]", f'path = "../../proteins/{u}/{mp["path"].name}"', f'name = "{name}_colabfold"',
              'format = "a3m"', f'reference_sequence_name = "{u}"', f"description = {toml_str(MSA_DESC)}", "",
              "[msas.metadata]", f'n_sequences = "{mp["n_sequences"]}"', 'database = "UniRef100 2026_01"',
              'method = "colabfold_search monomer, filter 0, no cluster expansion, no envdb"', ""]
    if sp.get("doi"):
        m += ["[publication]", f'doi = "{sp["doi"]}"', ""]
    return "\n".join(m)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--prepared", type=Path, required=True, help="prepare.py --out")
    ap.add_argument("--out", type=Path, default=None, help="default: --prepared")
    ap.add_argument("--benchmarks", nargs="*", default=list(SPEC))
    ap.add_argument("--genes", type=int, default=None, help="smoke test: first N genes per benchmark")
    ap.add_argument("--sprot", type=Path,
                    default=Path("/no_backup/diasfrazer/db/uniprot_sprot_human/uniprot_sprot_human.dat.gz"))
    ap.add_argument("--afdb", type=Path,
                    default=Path("/users/diasfrazer/egalpern/project/human_proteome_de/data/structures"))
    ap.add_argument("--msa-dir", type=Path, nargs="*", default=[], help="folders of <UNIPROT>.a3m, in order")
    ap.add_argument("--index-out", type=Path, default=Path("datasets/clinical"))
    a = ap.parse_args()
    out = a.out or a.prepared
    sp, _ = read_sprot(a.sprot)
    for bench in a.benchmarks:
        spec = SPEC[bench]
        y = spec["target"]
        v = pd.read_csv(a.prepared / bench / "variants.tsv.gz", sep="\t", low_memory=False,
                        dtype={"variant_id": str})
        stage, splits = out / "build" / bench, out / bench / "splits"
        for d in (stage, splits):
            shutil.rmtree(d, ignore_errors=True)
            d.mkdir(parents=True)
        store, cache = stage / "proteins", {}
        genes = sorted(v.gene.unique())[: a.genes]
        nt = COMMON + spec["extra"]
        rows = []
        for gene in genes:
            x = v[v.gene == gene].copy()
            u = x.uniprot_id.iloc[0]
            assert x.uniprot_id.nunique() == 1, gene
            pr = protein_inputs(u, sp, store, a.afdb, a.msa_dir, cache)
            name = f"{pr['entry_name']}_{spec['suffix']}"
            x["pos_"] = x.mutant.str[1:-1].astype(int)
            x = x.sort_values(["pos_", "mutant", "variant_id"])
            d = stage / "datasets" / name
            (d / "assays").mkdir(parents=True)
            assay = pd.DataFrame({"sequence": [mutate(pr["wt"], m) for m in x.mutant], y: x[y].to_numpy()})
            for c in nt:
                assay[c] = x[c].to_numpy()
            assay.to_csv(d / "assays" / f"{name}.csv", index=False)
            (d / "manifest.toml").write_text(manifest(name, bench, gene, u, pr, len(assay), nt))
            ds = Dataset.from_manifest(Manifest.from_path(d / "manifest.toml"))
            arc = KFoldSplitter(2).split(ds, targets=[y]).dump(path=splits)
            compress(arc)
            shutil.rmtree(d)  # staged assay CSV (full-length sequences) is in the archive now
            r = dict(dataset=name, gene=gene, uniprot_id=u, entry_name=pr["entry_name"],
                     seq_length=len(pr["wt"]), n_variants=len(assay), n_unique_mutants=assay.mutant.nunique())
            if bench == "gnomad":
                r.update(n_singletons=int((assay.AC_non_ukb == 1).sum()))
            else:
                r.update(n_pos=int((assay[y] == 1).sum()), n_neg=int((assay[y] == 0).sum()))
            r.update(mapping=";".join(f"{k}:{n}" for k, n in assay.mapping.value_counts().items()),
                     n_structures=len(pr["structures"]),
                     mean_plddt=round(np.mean([t["mean_plddt"] for t in pr["structures"]]), 2) if pr["structures"]
                     else None,
                     msa_n_sequences=pr["msa"]["n_sequences"] if pr["msa"] else None, inputs_note=pr["note"],
                     archive=f"{bench}/splits/{arc.name}")
            rows.append(r)
            print(f"{bench:8s} {name:40s} {len(assay):7,} variants", flush=True)
        shutil.rmtree(stage)  # staged inputs (MSA copies) are in the archives now
        idx = pd.DataFrame(rows)
        (a.index_out / bench).mkdir(parents=True, exist_ok=True)
        idx.to_csv(a.index_out / bench / f"{bench}.csv", index=False)
        print(f"{bench}: {len(idx)} archives, {idx.n_variants.sum():,} variants; without structure "
              f"{(idx.n_structures == 0).sum()}, without MSA {idx.msa_n_sequences.isna().sum()}", flush=True)
    print("BUILD OK")


if __name__ == "__main__":
    main()
