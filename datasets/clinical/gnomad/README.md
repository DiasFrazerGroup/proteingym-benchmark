# gnomAD rare variants

Do model scores track how rare a missense variant is in the population? Damaging variants are kept rare by
selection, so more damaging scores should go with lower allele frequencies.

| | |
|---|---|
| Variants | gnomAD v4.1 exome missense SNVs: PASS, MANE Select, rare (allele frequency < 0.1% outside UK Biobank), sites with allele number ≥ 50% of the maximum |
| Genes | 300 random genes (see below) |
| Size | 127,261 variants |
| Target | `log10_af`: log10 allele frequency in gnomAD v4.1 exomes excluding UK Biobank |
| Metric | **Spearman ρ(score, log10 AF) over all variants (headline; lower = better)**; also the mean ρ per gene (genes with ≥ 20 variants) |

**Genes.** The 300 genes drawn at random (seed 0) in clinical-curation, from genes with ≥ 20 rare variants. Genes
whose protein is longer than 2,700 aa or with fewer than 20 variants that map onto the UniProt protein are replaced
by genes from the same pool (seed 1). `genes.tsv` lists each gene and its origin (2 replacements).

**Caution: mutation rate.** Allele frequency also depends on how mutable a site is (e.g. CpG sites), so this metric
partly rewards predicting mutability; a DNA model (GPN-MSA) does best here for that reason.

One archive per gene: `<ENTRY>_GNOMAD_2024.splits.pgdata`. Layout and columns:
[docs/clinical/clinvar_denovo_gnomad.md](../../../docs/clinical/clinvar_denovo_gnomad.md). `gnomad.csv` is the index.
