# De novo variants

Can a model tell de novo missense variants of children with severe developmental disorders from those of controls?

| | |
|---|---|
| Variants | de novo missense variants from the popEVE paper (Orenbuch et al. 2025, Nat Genet, Suppl. Table 5), `ap_test_set` |
| Genes | 300 random genes (see below) |
| Size | 722 variants: 479 cases, 249 controls |
| Target | `dnm_case`: 1 = seen in a case, 0 = seen only in controls |
| Metric | **AUROC over all variants (headline)**; also the mean AUROC per gene (only 35 genes have both classes) |

**Genes.** 300 genes drawn at random (seed 0) from all `ap_test_set` genes. Genes whose protein is longer than
2,700 aa or with no variant that maps onto the UniProt protein are replaced by genes drawn from the same pool
(seed 1). `genes.tsv` lists each gene and its origin (9 replacements).

**Caution: gene identity.** In the `ap_test_set` every case is in a known developmental-disorder gene, while
controls are in any gene. So most genes have only cases or only controls, and knowing the gene alone separates
much of the data: a gene-only baseline (does the gene have another case variant) scores a higher AUROC than any
variant effect predictor. Report it next to every model.

One archive per gene: `<ENTRY>_DENOVO_2025.splits.pgdata`. The coordinates in the source are GRCh37
(`variant_id_grch37`); `variant_id` is the GRCh38 SNV. Layout and columns:
[docs/clinical/clinvar_denovo_gnomad.md](../../../docs/clinical/clinvar_denovo_gnomad.md). `denovo.csv` is the index.
