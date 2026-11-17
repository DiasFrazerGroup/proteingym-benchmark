"""Prepare the clinical benchmark-card variant sets for archiving: selection, protein choice, UniProt mapping.

The cards are defined in clinical-curation (scripts/16_benchmark_cards.py -> benchmark_cards/). This re-derives
each card's variant rows from the same curation tables with the same rules, then places every row on a UniProt
(Swiss-Prot canonical) wild type, so build_cards.py can write one PG2 archive per gene and card.

Cards (one archive set each; `complete` = scored by all six card models, the card's performance set):
  clinvar_v0                 ClinVar v0 (>=2 stars, germline, P/LP vs B/LB, genes >=10 of each)   target label
  dnm_aptest_ge1             popEVE ST5 ap_test_set DNMs, genes >=1 case & >=1 control             target dnm_case
  dnm_aptest_ge2             same, >=2 & >=2 (= DNM v1)                                           target dnm_case
  dnm_aptest_random300       same, 300 random ap_test_set genes (seed 0; NEW, no card yet)         target dnm_case
  gnomad_af_clinvar_dnm_ge1  gnomAD v4.1 rare missense (non-UKB AF < 0.1%, AN >= 50% of max),
                             genes of clinvar_v0 U dnm_aptest_ge1                                  target log10_af
  gnomad_af_clinvar_dnm_ge2  same, genes of clinvar_v0 U dnm_aptest_ge2                            target log10_af
  gnomad_af_random300        same, benchmark_cards/gnomad_af_random300_genes.tsv                   target log10_af
The reference cards dnm_aptest_full and gnomad_af_full are not archived.

Protein per gene: the Swiss-Prot accession that most of the gene's SNVs map to in the AlphaMissense hg38 SNV index
(UniProt canonical isoforms); genes without SNV hits fall back to Swiss-Prot by gene symbol.
Mapping per row (column `mapping`):
  am_snv_index        the row's SNV(s) -> protein change on that accession, wild-type residue checked. ClinVar and
                      gnomAD are GRCh38 (AlphaMissense hg38 index); the popEVE DNMs are GRCh37 (AlphaMissense hg19
                      table, same isoforms), and their GRCh38 SNV is then recovered from the hg38 index
                      (variant_id; the source coordinate is kept as variant_id_grch37)
  source_offset+k     source numbering (RefSeq / MANE protein) shifted by a constant offset k, fitted per (gene,
                      source isoform) on that isoform's SNV-mapped rows, or else on its remaining rows: accepted when
                      >= 90% of >= 10 rows match the wild-type residue, applied only to rows whose residue matches.
                      Isoforms with < 10 rows reuse an offset accepted for the gene if all their residues match.
                      (Needed for genes absent from the AlphaMissense canonical set, e.g. CACNA1C, COG4, TTC8.)
  otherwise dropped, listed in <out>/<card>/unmapped_rows.tsv

Outputs (--out):
  <card>/variants.tsv.gz   mapped rows: gene, uniprot_id, mutant (UniProt numbering), target, provenance, scores
  <card>/unmapped_rows.tsv rows that could not be placed
  selection_counts.tsv     per card, genes / rows before mapping (compare with benchmark_cards/summary.tsv)
  proteins.tsv             one row per protein: gene, accession, entry, length, how chosen, reusable MSA
  queries.fasta            Swiss-Prot wild types of proteins with no reusable MSA (input of the MSA search)
and --genes-out/dnm_aptest_random300_genes.tsv (the new random DNM gene list, committed).

Run (proteingym_clinical env; ~30 GB RAM for the gnomAD table): see cards.sbatch.
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

MODELS = ["am", "popeve", "cadd", "esm1b", "vesm", "gpn_msa"]
SCORE_COL = {"am": "score_alphamissense", "popeve": "score_popeve", "cadd": "score_cadd_phred",
             "esm1b": "score_esm1b", "vesm": "score_vesm3b", "gpn_msa": "score_gpn_msa"}
MIN_N, RARE, OFFSET_MIN, OFFSET_FRAC = 20, 1e-3, 10, 0.9
MUT = re.compile(r"^([A-Z])(\d+)([A-Z])$")

CARDS = {  # name: (source, target)
    "clinvar_v0": ("clinvar", "label"),
    "dnm_aptest_ge1": ("dnm", "dnm_case"),
    "dnm_aptest_ge2": ("dnm", "dnm_case"),
    "dnm_aptest_random300": ("dnm", "dnm_case"),
    "gnomad_af_clinvar_dnm_ge1": ("gnomad", "log10_af"),
    "gnomad_af_clinvar_dnm_ge2": ("gnomad", "log10_af"),
    "gnomad_af_random300": ("gnomad", "log10_af"),
}
EXTRA = {  # card-specific non-target columns kept in the assay
    "clinvar": ["hgvsp", "source_transcript", "clnsig", "max_stars", "n_records"],
    "dnm": ["variant_id_grch37", "source_protein", "n_sdd", "n_control"],
    "gnomad": ["af", "AC_non_ukb", "AN_non_ukb", "AC", "AN", "enst", "mane_select"],
}


def gene_filter(d, y, k):
    g = d.groupby("gene")[y].agg(["sum", "size"])
    return d[d.gene.isin(g[(g["sum"] >= k) & (g["size"] - g["sum"] >= k)].index)].copy()


def split_mut(m):
    x = MUT.match(m) if isinstance(m, str) else None
    return (x[1], int(x[2]), x[3]) if x else (None, np.nan, None)


def load_sources(cur):
    bm = cur / "data/processed/benchmarks"
    # ClinVar v0: one row per (gene, protein change); SNVs from clinvar_v0_snvs.tsv (several per change possible)
    cv = pd.read_csv(bm / "clinvar_v0_scores.tsv", sep="\t", low_memory=False)
    cv = cv.merge(pd.read_csv(bm / "clinvar_v0_extra_scores.tsv", sep="\t"), on=["gene", "hgvs_p"], how="left",
                  validate="1:1")
    snv = pd.read_csv(bm / "clinvar_v0_snvs.tsv", sep="\t", dtype={"chrom": str, "pos": str})
    snv["vid"] = snv.chrom.str.removeprefix("chr") + "-" + snv.pos + "-" + snv.ref + "-" + snv.alt
    cv = cv.merge(snv.groupby(["gene", "hgvs_p"]).vid.agg(lambda s: ";".join(sorted(set(s)))).rename("variant_id")
                  .reset_index(), on=["gene", "hgvs_p"], how="left", validate="1:1")
    cv = cv.rename(columns={"hgvs_p": "hgvsp", "tx": "source_transcript", "mutant": "mutant_source"})

    # DNMs: ap_test_set table (popEVE = gnomAD version, as the cards)
    dn = pd.read_csv(bm / "dnm_aptest_scores.tsv", sep="\t", low_memory=False, dtype={"chrom": str, "pos": str})
    dn = dn.merge(pd.read_csv(bm / "dnm_aptest_extra_scores.tsv", sep="\t", dtype={"chrom": str, "pos": str}),
                  on=["chrom", "pos", "ref", "alt"], how="left", validate="m:1")
    dn["variant_id"] = dn.chrom.str.removeprefix("chr") + "-" + dn.pos + "-" + dn.ref + "-" + dn.alt
    dn = dn.rename(columns={"protein": "source_protein", "mutant": "mutant_source"})

    # gnomAD rare missense (popEVE = UKBB), filters exactly as scripts/16
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


def select(cur, cv, dn, g, genes_out):
    """Card name -> rows (before mapping), as scripts/16_benchmark_cards.py."""
    c = {"clinvar_v0": cv.copy(),
         "dnm_aptest_ge1": gene_filter(dn, "dnm_case", 1),
         "dnm_aptest_ge2": gene_filter(dn, "dnm_case", 2)}
    rnd = np.random.default_rng(0).choice(np.array(sorted(dn.gene.unique())), 300, replace=False)
    pd.Series(sorted(rnd), name="gene").to_csv(genes_out / "dnm_aptest_random300_genes.tsv", sep="\t", index=False)
    c["dnm_aptest_random300"] = dn[dn.gene.isin(rnd)].copy()
    cvg = set(cv.gene)
    r300 = set(pd.read_csv(cur / "benchmark_cards/gnomad_af_random300_genes.tsv", sep="\t").gene)
    for name, genes in [("gnomad_af_clinvar_dnm_ge1", cvg | set(c["dnm_aptest_ge1"].gene)),
                        ("gnomad_af_clinvar_dnm_ge2", cvg | set(c["dnm_aptest_ge2"].gene)),
                        ("gnomad_af_random300", r300)]:
        c[name] = g[g.gene.isin(genes)].copy()
    return c


def counts(name, d):
    src, y = CARDS[name]
    out = {"card": name}
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
    (same UniProt canonical isoforms as the hg38 table). The popEVE DNM coordinates are GRCh37."""
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
    """Rows of li whose source residue matches wt at source position + o, placed on UniProt numbering."""
    ok = [0 < p + o <= len(wt) and wt[int(p + o) - 1] == a and a != b
          for a, p, b in zip(li.src_wt, li.src_pos, li.src_mt)]
    x = li[ok]
    return x.assign(uniprot_id=acc, mapping=f"source_offset{o:+d}",
                    mutant=x.src_wt + (x.src_pos + o).astype(int).astype(str) + x.src_mt)


def map_rows(rows, sp, gene2acc, ex):
    """rows: DataFrame(key, gene, mutant_source); ex: SNV hits DataFrame(key, uniprot_id, protein_variant).
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
        # rows without an SNV hit: constant offset between the source isoform's numbering and UniProt, one per
        # (gene, source isoform), from that isoform's SNV-mapped rows when there are enough, else from its
        # remaining rows; accepted when >= OFFSET_FRAC of >= OFFSET_MIN rows match the wild-type residue, applied
        # only to rows whose residue matches
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
        # isoforms with too few rows to fit an offset: reuse an offset accepted for this gene on another isoform,
        # only when every row's residue matches under it
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
                    help="AlphaMissense hg19 table, to map the GRCh37 DNM coordinates")
    ap.add_argument("--msa-dirs", type=Path, nargs="*", default=[],
                    help="folders of <UNIPROT>.a3m to reuse (searched in order)")
    ap.add_argument("--genes-out", type=Path, default=Path("datasets/clinical/cards"))
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    a.genes_out.mkdir(parents=True, exist_ok=True)

    cv, dn, g = load_sources(a.curation)
    cards = select(a.curation, cv, dn, g, a.genes_out)
    sel = pd.DataFrame([counts(n, d) for n, d in cards.items()])
    sel.to_csv(a.out / "selection_counts.tsv", sep="\t", index=False)
    print(sel.to_string(index=False), flush=True)

    # map each source once (on the union of the rows its cards use), then subset per card
    sp, gene2acc = read_sprot(a.sprot)
    keys = {"clinvar": ["gene", "hgvsp"], "dnm": ["variant_id", "gene", "mutant_source"],
            "gnomad": ["variant_id", "gene", "enst"]}
    src_rows = {}
    for n, d in cards.items():
        s = CARDS[n][0]
        src_rows[s] = pd.concat([src_rows.get(s, d.iloc[:0]), d]).drop_duplicates(keys[s])
    hits = []
    for s, d in src_rows.items():
        d["key"] = s + ":" + d[keys[s]].astype(str).agg("|".join, axis=1)
        e = d[["key", "variant_id"]].dropna().assign(variant_id=lambda x: x.variant_id.str.split(";")).explode(
            "variant_id")
        look = snv_lookup_hg19(a.am_hg19, set(e.variant_id)) if s == "dnm" else snv_lookup(a.snv_index,
                                                                                             set(e.variant_id))
        hits.append(e.merge(look, on="variant_id")[["key", "uniprot_id", "protein_variant"]])
        print(f"{s}: {e.variant_id.nunique():,} SNVs ({'GRCh37' if s == 'dnm' else 'GRCh38'}) -> "
              f"{len(look):,} (SNV, protein) hits", flush=True)
    # one protein per gene across all sources: map the pooled rows of every source together
    iso_col = {"clinvar": "source_transcript", "dnm": "source_protein", "gnomad": "enst"}
    pooled = pd.concat([d[["key", "gene", "mutant_source"]].assign(src_iso=s + ":" + d[iso_col[s]].astype(str))
                        for s, d in src_rows.items()], ignore_index=True)
    prots, mapped, um = map_rows(pooled, sp, gene2acc, pd.concat(hits, ignore_index=True))
    mp = mapped.set_index("key")[["uniprot_id", "mutant", "mapping"]]
    print(f"mapped {len(mapped):,} of {len(pooled):,} source rows; {len(um):,} unmapped; "
          f"{mapped.mapping.value_counts().to_dict()}", flush=True)

    for n, d in cards.items():
        s, y = CARDS[n]
        d = d.copy()
        d["key"] = s + ":" + d[keys[s]].astype(str).agg("|".join, axis=1)
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
        print(f"{n:28s} {out.gene.nunique():4d} genes  {len(out):8,} rows mapped  {len(u):5,} unmapped "
              f"({len(u) / len(d):.2%})", flush=True)

    # proteins: reusable MSAs and the query FASTA for the rest
    used = set().union(*[set(pd.read_csv(a.out / n / "variants.tsv.gz", sep="\t", usecols=["uniprot_id"]).uniprot_id)
                         for n in cards])
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
          f"{(prots.msa_reused != '').sum()}; queries.fasta: {len(q)} proteins, {q.length.sum():,} residues, "
          f"longest {q.nlargest(5, 'length')[['gene', 'length']].values.tolist()}")


if __name__ == "__main__":
    main()
