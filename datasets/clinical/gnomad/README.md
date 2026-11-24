# gnomAD rare variants

Do model scores track how rare a missense variant is in the population? Damaging variants are kept rare by
selection, so more damaging scores should go with lower allele frequencies.

| | |
|---|---|
| Variants | gnomAD v4.1 exome missense SNVs: PASS, MANE Select, rare (allele frequency < 0.1% outside UK Biobank), sites with allele number ≥ 50% of the maximum |
| Genes | 300 random genes (see below) |
| Size | 130,198 variants |
| Target | `log10_af`: log10 allele frequency in gnomAD v4.1 exomes excluding UK Biobank |
| Metric | **Spearman ρ(score, log10 AF) over all variants (headline; lower = better)**; also the mean ρ per gene (genes with ≥ 20 variants) |

**Genes.** One draw of 300 genes (`numpy` `default_rng(0)`, without replacement) from the 17,564 eligible genes:
those whose protein is at most 2,700 aa and has at least 20 rare variants placed on it.

**Caution: mutation rate.** Allele frequency also depends on how mutable a site is (e.g. CpG sites), so this metric
partly rewards predicting mutability; a DNA model (GPN-MSA) does best here for that reason.

One archive per gene: `<ENTRY>_GNOMAD_2024.splits.pgdata`. Layout and columns:
[docs/clinical/benchmarks.md](../../../docs/clinical/benchmarks.md). `gnomad.csv` is the index.
