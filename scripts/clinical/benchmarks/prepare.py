"""Prepare the ClinVar, de novo and gnomAD benchmarks: variant selection, protein choice, UniProt mapping.

The variant tables come from the clinical-curation repo. This selects each benchmark's variants, then places every
variant on a UniProt (Swiss-Prot canonical) wild type, so build.py can write one PG2 archive per gene.

Benchmarks:
  clinvar  ClinVar missense (Nov 2026 release): >= 2 review stars, germline, pathogenic / likely pathogenic vs
           benign / likely benign, genes with >= 10 of each                                         target label
  denovo   de novo missense variants of the popEVE paper (Orenbuch et al. 2025, Suppl. Table 5), ap_test_set:
           severe developmental disorder cases vs controls; 300 random genes*                       target dnm_case
  gnomad   gnomAD v4.1 exome rare missense (PASS, MANE Select, non-UKB AF < 0.1%, AN >= 50% of max);
           300 random genes*                                                                         target log10_af
  * Exactly 300 genes, each with a protein <= 2,700 aa and mapped variants (gnomAD: >= 20). De novo: drawn from all
    ap_test_set genes (seed 0). gnomAD: the 300 genes drawn in clinical-curation (seed 0, from genes with >= 20
    scored rare variants). Genes that fail the rule are replaced by genes from the same pool in a fixed random order
    (seed 1). Gene lists with each gene's origin: datasets/clinical/<benchmark>/genes.tsv.

Protein per gene: the Swiss-Prot accession that most of the gene's SNVs map to in the AlphaMissense tables (UniProt
canonical isoforms); genes without SNV hits fall back to Swiss-Prot by gene symbol.
Mapping per variant (column `mapping`):
  am_snv_index        the variant's SNV(s) -> protein change on that accession, wild-type residue checked. ClinVar and
                      gnomAD are GRCh38 (AlphaMissense hg38 index); the de novo variants are GRCh37 (AlphaMissense hg19
                      table, same isoforms), and their GRCh38 SNV is recovered from the hg38 index (variant_id; the
                      source coordinate is kept as variant_id_grch37)
  source_offset+k     source protein numbering (RefSeq / MANE) shifted by a constant offset k, fitted per (gene,
                      source isoform) on that isoform's SNV-mapped variants, or else on its remaining variants:
                      accepted when >= 90% of >= 10 variants match the wild-type residue, applied only to variants
                      whose residue matches. Isoforms with < 10 variants reuse an offset accepted for the gene if all
                      their residues match. (Needed for genes absent from the AlphaMissense canonical set.)
  otherwise dropped, listed in <out>/<benchmark>/unmapped_rows.tsv

Outputs (--out):
  <benchmark>/variants.tsv.gz   mapped variants: gene, uniprot_id, mutant (UniProt numbering), target, provenance;
                                plus clinical-curation's reference scores, used only by check.py (never archived)
  <benchmark>/unmapped_rows.tsv variants that could not be placed
  selection_counts.tsv          per benchmark, genes / variants before mapping
  proteins.tsv                  one row per protein: gene, accession, entry, length, how chosen, reusable MSA
  queries.fasta                 Swiss-Prot wild types of proteins with no reusable MSA (input of the MSA search)

Run (proteingym_clinical env; ~30 GB RAM for the gnomAD table): bash scripts/clinical/benchmarks/submit.sh prepare
"""
import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.dataset as pads

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from build_genebass import read_sprot  # noqa: E402

# clinical-curation reference scores (check.py only)
MODELS = ["am", "popeve", "cadd", "esm1b", "vesm", "gpn_msa"]
SCORE_COL = {"am": "score_alphamissense", "popeve": "score_popeve", "cadd": "score_cadd_phred",
             "esm1b": "score_esm1b", "vesm": "score_vesm3b", "gpn_msa": "score_gpn_msa"}
BASE_MODELS = ["am", "popeve", "cadd"]          # clinical-curation's gnomAD random-300 eligibility models
MIN_N, RARE, OFFSET_MIN, OFFSET_FRAC = 20, 1e-3, 10, 0.9
N_GENES, MAX_LEN, REPLACE_BUFFER = 300, 2700, 100
MUT = re.compile(r"^([A-Z])(\d+)([A-Z])$")

BENCHMARKS = {  # name: (source table, target)
    "clinvar": ("clinvar", "label"),
    "denovo": ("dnm", "dnm_case"),
    "gnomad": ("gnomad", "log10_af"),
}
EXTRA = {  # source-specific non-target columns kept in the assay
    "clinvar": ["hgvsp", "source_transcript", "clnsig", "max_stars", "n_records"],
    "dnm": ["variant_id_grch37", "source_protein", "n_sdd", "n_control"],
    "gnomad": ["af", "AC_non_ukb", "AN_non_ukb", "AC", "AN", "enst", "mane_select"],
}
KEYS = {"clinvar": ["gene", "hgvsp"], "dnm": ["variant_id", "gene", "mutant_source"],
        "gnomad": ["variant_id", "gene", "enst"]}


def split_mut(m):
    x = MUT.match(m) if isinstance(m, str) else None
    return (x[1], int(x[2]), x[3]) if x else (None, np.nan, None)


def load_sources(cur):
    bm = cur / "data/processed/benchmarks"
    # ClinVar: one row per (gene, protein change); SNVs from clinvar_v0_snvs.tsv (several per change possible)
    cv = pd.read_csv(bm / "clinvar_v0_scores.tsv", sep="\t", low_memory=False)
    cv = cv.merge(pd.read_csv(bm / "clinvar_v0_extra_scores.tsv", sep="\t"), on=["gene", "hgvs_p"], how="left",
                  validate="1:1")
    snv = pd.read_csv(bm / "clinvar_v0_snvs.tsv", sep="\t", dtype={"chrom": str, "pos": str})
    snv["vid"] = snv.chrom.str.removeprefix("chr") + "-" + snv.pos + "-" + snv.ref + "-" + snv.alt
    cv = cv.merge(snv.groupby(["gene", "hgvs_p"]).vid.agg(lambda s: ";".join(sorted(set(s)))).rename("variant_id")
                  .reset_index(), on=["gene", "hgvs_p"], how="left", validate="1:1")
    cv = cv.rename(columns={"hgvs_p": "hgvsp", "tx": "source_transcript", "mutant": "mutant_source"})

    # de novo variants: ap_test_set table (GRCh37)
    dn = pd.read_csv(bm / "dnm_aptest_scores.tsv", sep="\t", low_memory=False, dtype={"chrom": str, "pos": str})
    dn = dn.merge(pd.read_csv(bm / "dnm_aptest_extra_scores.tsv", sep="\t", dtype={"chrom": str, "pos": str}),
                  on=["chrom", "pos", "ref", "alt"], how="left", validate="m:1")
    dn["variant_id"] = dn.chrom.str.removeprefix("chr") + "-" + dn.pos + "-" + dn.ref + "-" + dn.alt
    dn = dn.rename(columns={"protein": "source_protein", "mutant": "mutant_source"})

    # gnomAD rare missense
    usecols = ["gene", "chrom", "pos", "ref", "alt", "AC", "AN", "AC_non_ukb", "AN_non_ukb", "enst", "mane_select",
               "prot_pos", "wt", "mt", "am", "popeve", "cadd"]
    g = pd.read_csv(cur / "data/interim/gnomad_missense_scores_all.tsv.gz", sep="\t", usecols=usecols,
                    dtype={"chrom": str, "pos": str, "prot_pos": str})
    g = g.merge(pd.read_csv(cur / "data/interim/gnomad_missense_extra_scores_all.tsv.gz", sep="\t",
                            dtype={"chrom": str, "pos": str}), on=["chrom", "pos", "ref", "alt"], how="left")
    g = g[(g.AC_non_ukb > 0) & (g.AN_non_ukb >= 0.5 * g.AN_non_ukb.max())].copy()
    g["af"] = g.AC_non_ukb / g.AN_non_ukb
    g = g[g.af < RARE].copy()
    g["log10_af"] = np.log10(g.af)
    g["variant_id"] = g.chrom.str.removeprefix("chr") + "-" + g.pos + "-" + g.ref + "-" + g.alt
    g["mutant_source"] = g.wt + g.prot_pos + g.mt
    return cv, dn, g


def replacement_order(pool, orig, n=REPLACE_BUFFER):
    """Replacement candidates: the pool minus the original draw, in a fixed random order (seed 1)."""
    rest = np.array(sorted(set(pool) - set(orig)))
    return list(np.random.default_rng(1).permutation(rest)[:n])


def select(cur, cv, dn, g):
    """Benchmark -> variants (before mapping); for the random-gene benchmarks also (original draw, replacement
    candidates). Their variants cover both; finalise_random() picks the 300 genes after mapping."""
    pool = np.array(sorted(dn.gene.unique()))
    orig = sorted(np.random.default_rng(0).choice(pool, N_GENES, replace=False))
    rand = {"denovo": (orig, replacement_order(pool, orig))}
    orig = sorted(pd.read_csv(cur / "benchmark_cards/gnomad_af_random300_genes.tsv", sep="\t").gene)
    nper = g.dropna(subset=BASE_MODELS).groupby("gene").size()       # the pool it was drawn from
    rand["gnomad"] = (orig, replacement_order(np.array(sorted(nper[nper >= MIN_N].index)), orig))
    b = {"clinvar": cv.copy(),
         "denovo": dn[dn.gene.isin(set().union(*rand["denovo"]))].copy(),
         "gnomad": g[g.gene.isin(set().union(*rand["gnomad"]))].copy()}
    return b, rand


def finalise_random(bench, rand, mp, prots, genes_out):
    """Keep the original draw's usable genes, top up to N_GENES with replacement candidates in order.
    Usable = protein <= MAX_LEN aa and enough mapped variants (gnomAD: >= MIN_N; de novo: >= 1).
    Writes <genes_out>/<benchmark>/genes.tsv (gene, origin, protein length, mapped variants)."""
    length = prots.set_index("gene").length
    for n, (orig, repl) in rand.items():
        s, _ = BENCHMARKS[n]
        d = bench[n]
        k = s + ":" + d[KEYS[s]].astype(str).agg("|".join, axis=1)
        nmap = d[k.isin(mp.index)].groupby("gene").size()
        need = MIN_N if s == "gnomad" else 1

        def usable(gn):
            return gn in length.index and length[gn] <= MAX_LEN and nmap.get(gn, 0) >= need
        keep = [x for x in orig if usable(x)]
        add = [x for x in repl if usable(x)][: N_GENES - len(keep)]
        assert len(keep) + len(add) == N_GENES, (n, len(keep), len(add))
        why = {str(x): (f"protein > {MAX_LEN} aa" if x in length.index and length[x] > MAX_LEN else
                        f"< {need} mapped variants") for x in orig if x not in keep}
        print(f"{n}: kept {len(keep)} of the original {N_GENES}, replaced {len(why)} {why} by "
              f"{[str(x) for x in add]}", flush=True)
        bench[n] = d[d.gene.isin(set(keep) | set(add))].copy()
        (genes_out / n).mkdir(parents=True, exist_ok=True)
        pd.DataFrame({"gene": keep + add, "origin": ["original_draw"] * len(keep) + ["replacement"] * len(add),
                      "protein_length": [int(length[x]) for x in keep + add],
                      "n_mapped_variants": [int(nmap[x]) for x in keep + add]}).sort_values("gene").to_csv(
            genes_out / n / "genes.tsv", sep="\t", index=False)


def counts(name, d):
    src, y = BENCHMARKS[name]
    out = {"benchmark": name}
    for s, x in [("all", d), ("complete", d.dropna(subset=MODELS))]:
        if src == "gnomad":
            n = x.groupby("gene").size()
            out.update({f"genes_{s}": len(n), f"genes_min_{s}": int((n >= MIN_N).sum()), f"variants_{s}": len(x),
                        f"singletons_{s}": int((x.AC_non_ukb == 1).sum())})
        else:
            gg = x.groupby("gene")[y].agg(["sum", "size"])
            out.update({f"genes_{s}": len(gg),
                        f"genes_both_{s}": int(((gg["sum"] > 0) & (gg["size"] > gg["sum"])).sum()),
                        f"variants_{s}": len(x), f"pos_{s}": int(x[y].sum()), f"neg_{s}": int((x[y] == 0).sum())})
    return out


def _positions(vids):
    v = pd.Series(sorted(vids)).str.split("-", expand=True)
    return set(v[(v[2].str.len() == 1) & (v[3].str.len() == 1)][1].astype(int))


def snv_lookup(index_path, vids):
    """GRCh38 SNV ids (chrom-pos-ref-alt, no 'chr') -> DataFrame(variant_id, uniprot_id, protein_variant)."""
    t = pads.dataset(index_path).to_table(filter=pads.field("pos").isin(sorted(_positions(vids)))).to_pandas()
    t["variant_id"] = t.chrom.astype(str) + "-" + t.pos.astype(str) + "-" + t.ref + "-" + t.alt
    return t[t.variant_id.isin(set(vids))][["variant_id", "uniprot_id", "protein_variant"]]


def snv_lookup_hg19(am_hg19, vids):
    """GRCh37 SNV ids -> DataFrame(variant_id, uniprot_id, protein_variant), from AlphaMissense_hg19.tsv.gz
    (same UniProt canonical isoforms as the hg38 table)."""
    pos, out = _positions(vids), []
    for ch in pd.read_csv(am_hg19, sep="\t", skiprows=3, chunksize=5_000_000, dtype={"POS": int},
                          usecols=["#CHROM", "POS", "REF", "ALT", "uniprot_id", "protein_variant"]):
        ch = ch[ch.POS.isin(pos)]
        out.append(ch.assign(variant_id=ch["#CHROM"].str.removeprefix("chr") + "-" + ch.POS.astype(str) + "-"
                             + ch.REF + "-" + ch.ALT))
    t = pd.concat(out)
    return t[t.variant_id.isin(set(vids))][["variant_id", "uniprot_id", "protein_variant"]]


def grch38_ids(index_path, mapped):
    """GRCh38 SNV of each (uniprot_id, mutant, ref, alt) from the hg38 index; '' when absent or ambiguous."""
    t = pads.dataset(index_path).to_table(
        filter=pads.field("uniprot_id").isin(sorted(set(mapped.uniprot_id)))).to_pandas()
    t["vid38"] = t.chrom.astype(str) + "-" + t.pos.astype(str) + "-" + t.ref + "-" + t.alt
    t = t.groupby(["uniprot_id", "protein_variant", "ref", "alt"]).vid38.agg(
        lambda s: s.iloc[0] if len(s) == 1 else "").reset_index()
    m = mapped.merge(t.rename(columns={"protein_variant": "mutant"}), on=["uniprot_id", "mutant", "ref", "alt"],
                     how="left")
    return m.vid38.fillna("").to_numpy()


def best_offset(aa, pos, wt, max_shift=500):
    """Offset in -max_shift..max_shift under which most source residues (aa at pos) match wt.
    Returns (offset, fraction matched)."""
    w = np.frombuffer(wt.encode(), dtype="S1")
    a = np.array(aa, dtype="S1")
    p = np.asarray(pos, dtype=int)
    best = (0, -1.0)
    for o in sorted(range(-max_shift, max_shift + 1), key=abs):
        q = p + o
        ok = (q >= 1) & (q <= len(w))
        f = (w[q[ok] - 1] == a[ok]).sum() / len(p)
        if f > best[1]:
            best = (o, f)
    return best


def shift(li, o, wt, acc):
    """Variants of li whose source residue matches wt at source position + o, placed on UniProt numbering."""
    ok = [0 < p + o <= len(wt) and wt[int(p + o) - 1] == a and a != b
          for a, p, b in zip(li.src_wt, li.src_pos, li.src_mt)]
    x = li[ok]
    return x.assign(uniprot_id=acc, mapping=f"source_offset{o:+d}",
                    mutant=x.src_wt + (x.src_pos + o).astype(int).astype(str) + x.src_mt)


def map_rows(rows, sp, gene2acc, ex):
    """rows: DataFrame(key, gene, mutant_source, src_iso); ex: SNV hits DataFrame(key, uniprot_id, protein_variant).
    Returns (proteins, mapped, unmapped)."""
    r = rows.copy()
    r[["src_wt", "src_pos", "src_mt"]] = pd.DataFrame([split_mut(m) for m in r.mutant_source], index=r.index)
    ex = ex[ex.uniprot_id.isin(sp)]
    ex = ex.merge(r[["key", "gene"]], on="key")
    prots, mapped, unmapped = [], [], []
    for gene, gr in r.groupby("gene", sort=True):
        e = ex[ex.gene == gene]
        if len(e):
            acc, how = e.drop_duplicates(["key", "uniprot_id"]).uniprot_id.value_counts().index[0], "am_snv_majority"
        elif gene in gene2acc:
            acc, how = gene2acc[gene], "sprot_gene_symbol"
        else:
            unmapped.append(gr.assign(reason="no UniProt entry"))
            continue
        wt = sp[acc]["sequence"]
        prots.append(dict(gene=gene, uniprot_id=acc, entry_name=sp[acc]["entry_name"], length=len(wt),
                          protein_choice=how, sprot_gene=sp[acc]["gene"]))
        e = e[e.uniprot_id == acc]
        hit = {}
        for k, pv in zip(e.key, e.protein_variant):
            a, p, _ = split_mut(pv)
            if k not in hit and a and p <= len(wt) and wt[int(p) - 1] == a:
                hit[k] = pv
        g1 = gr[gr.key.isin(hit)].assign(uniprot_id=acc, mutant=lambda x: x.key.map(hit), mapping="am_snv_index")
        left = gr[~gr.key.isin(hit)]
        # variants without an SNV hit: constant offset between the source isoform's numbering and UniProt, one per
        # (gene, source isoform), from that isoform's SNV-mapped variants when there are enough, else from its
        # remaining variants; accepted when >= OFFSET_FRAC of >= OFFSET_MIN match the wild-type residue, applied
        # only to variants whose residue matches
        new, accepted, small = [], [], []
        for iso, li in left.groupby("src_iso", dropna=False, sort=False):
            li = li[li.src_pos.notna()]
            gi = g1[(g1.src_iso == iso) & g1.src_pos.notna()]
            o = None
            if len(gi) >= OFFSET_MIN:
                d = gi.mutant.map(lambda m: split_mut(m)[1]) - gi.src_pos
                if d.value_counts(normalize=True).iloc[0] >= OFFSET_FRAC:
                    o = int(d.mode().iloc[0])
            elif len(li) >= OFFSET_MIN:
                o, f = best_offset(li.src_wt, li.src_pos, wt)
                o = o if f >= OFFSET_FRAC else None
            elif len(li):
                small.append(li)
                continue
            if o is not None:
                accepted.append(o)
                new.append(shift(li, o, wt, acc))
        # isoforms with too few variants to fit an offset: reuse an offset accepted for this gene on another
        # isoform, only when every variant's residue matches under it
        for li in small:
            for o in dict.fromkeys(accepted):
                s = shift(li, o, wt, acc)
                if len(s) == len(li):
                    new.append(s)
                    break
        if new:
            g1 = pd.concat([g1, *new])
            left = left[~left.key.isin(g1.key)]
        if len(left):
            unmapped.append(left.assign(uniprot_id=acc, reason="no SNV hit on the protein and no source offset"))
        mapped.append(g1)
    um = pd.concat(unmapped) if unmapped else rows.iloc[:0]
    return pd.DataFrame(prots), pd.concat(mapped), um


def msa_ok(path, wt):
    """Reusable a3m: first record (match states) equals the wild type."""
    seq, n = [], 0
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                n += 1
                if n == 2:
                    break
            else:
                seq.append(line.strip())
    s = "".join(seq)
    return s.replace("-", "") == wt and s == s.upper()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--curation", type=Path, required=True, help="clinical-curation repo")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--sprot", type=Path,
                    default=Path("/no_backup/diasfrazer/db/uniprot_sprot_human/uniprot_sprot_human.dat.gz"))
    ap.add_argument("--snv-index", type=Path,
                    default=Path("/no_backup/diasfrazer/egalpern_pgb_clinical/layers/am_hg38_snv_index.parquet"))
    ap.add_argument("--am-hg19", type=Path,
                    default=Path("/no_backup/diasfrazer/db/alphamissense/raw/AlphaMissense_hg19.tsv.gz"),
                    help="AlphaMissense hg19 table, to map the GRCh37 de novo coordinates")
    ap.add_argument("--msa-dirs", type=Path, nargs="*", default=[],
                    help="folders of <UNIPROT>.a3m to reuse (searched in order)")
    ap.add_argument("--genes-out", type=Path, default=Path("datasets/clinical"))
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    cv, dn, g = load_sources(a.curation)
    bench, rand = select(a.curation, cv, dn, g)

    # map all variants together (one protein per gene), then split per benchmark
    sp, gene2acc = read_sprot(a.sprot)
    hits = []
    for n, d in bench.items():
        s = BENCHMARKS[n][0]
        d["key"] = s + ":" + d[KEYS[s]].astype(str).agg("|".join, axis=1)
        e = d[["key", "variant_id"]].dropna().assign(variant_id=lambda x: x.variant_id.str.split(";")).explode(
            "variant_id")
        look = (snv_lookup_hg19(a.am_hg19, set(e.variant_id)) if s == "dnm" else
                snv_lookup(a.snv_index, set(e.variant_id)))
        hits.append(e.merge(look, on="variant_id")[["key", "uniprot_id", "protein_variant"]])
        print(f"{n}: {e.variant_id.nunique():,} SNVs ({'GRCh37' if s == 'dnm' else 'GRCh38'}) -> "
              f"{len(look):,} (SNV, protein) hits", flush=True)
    iso_col = {"clinvar": "source_transcript", "dnm": "source_protein", "gnomad": "enst"}
    pooled = pd.concat([d[["key", "gene", "mutant_source"]].assign(
        src_iso=BENCHMARKS[n][0] + ":" + d[iso_col[BENCHMARKS[n][0]]].astype(str)) for n, d in bench.items()],
        ignore_index=True)
    prots, mapped, um = map_rows(pooled, sp, gene2acc, pd.concat(hits, ignore_index=True))
    mp = mapped.set_index("key")[["uniprot_id", "mutant", "mapping"]]
    print(f"mapped {len(mapped):,} of {len(pooled):,} variants; {len(um):,} unmapped; "
          f"{mapped.mapping.value_counts().to_dict()}", flush=True)

    finalise_random(bench, rand, mp, prots, a.genes_out)
    sel = pd.DataFrame([counts(n, d) for n, d in bench.items()])
    sel.to_csv(a.out / "selection_counts.tsv", sep="\t", index=False)
    print(sel.to_string(index=False), flush=True)

    for n, d in bench.items():
        s, y = BENCHMARKS[n]
        d = d.copy()
        d["complete"] = d[MODELS].notna().all(axis=1)
        m = d.join(mp, on="key", how="inner")
        if s == "dnm":  # GRCh37 source coordinates: keep them, and recover the GRCh38 SNV of the mapped change
            m["variant_id_grch37"] = m.variant_id
            m["variant_id"] = grch38_ids(a.snv_index, m[["uniprot_id", "mutant", "ref", "alt"]])
        out = m[["gene", "uniprot_id", "mutant", y, "mutant_source", "variant_id", "mapping", "complete",
                 *EXTRA[s], *MODELS]].rename(columns=SCORE_COL)
        (a.out / n).mkdir(exist_ok=True)
        out.to_csv(a.out / n / "variants.tsv.gz", sep="\t", index=False)
        u = d[~d.key.isin(mp.index)]
        u.drop(columns=["key"]).to_csv(a.out / n / "unmapped_rows.tsv", sep="\t", index=False)
        print(f"{n:8s} {out.gene.nunique():4d} genes  {len(out):8,} variants mapped  {len(u):5,} unmapped "
              f"({len(u) / len(d):.2%})", flush=True)

    # proteins: reusable MSAs and the query FASTA for the rest
    used = set().union(*[set(pd.read_csv(a.out / n / "variants.tsv.gz", sep="\t", usecols=["uniprot_id"]).uniprot_id)
                         for n in bench])
    prots = prots[prots.uniprot_id.isin(used)].copy()
    reuse = []
    for u in prots.uniprot_id:
        hit = next((d / f"{u}.a3m" for d in a.msa_dirs
                    if (d / f"{u}.a3m").exists() and msa_ok(d / f"{u}.a3m", sp[u]["sequence"])), None)
        reuse.append(str(hit) if hit else "")
    prots["msa_reused"] = reuse
    prots.to_csv(a.out / "proteins.tsv", sep="\t", index=False)
    q = prots[prots.msa_reused == ""].drop_duplicates("uniprot_id").sort_values("length")
    with open(a.out / "queries.fasta", "w") as f:
        for u in q.uniprot_id:
            f.write(f">{u}\n{sp[u]['sequence']}\n")
    print(f"{prots.uniprot_id.nunique()} proteins ({prots.gene.nunique()} genes); MSA reused for "
          f"{(prots.msa_reused != '').sum()}; queries.fasta: {len(q)} proteins")


if __name__ == "__main__":
    main()
