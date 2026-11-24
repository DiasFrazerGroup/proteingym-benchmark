# De novo variants

Can a model tell de novo missense variants of children with severe developmental disorders from those of controls?

| | |
|---|---|
| Variants | de novo missense variants from the popEVE paper (Orenbuch et al. 2025, Nat Genet, Suppl. Table 5), `ap_test_set` |
| Genes | 300 random genes (see below) |
| Size | 802 variants: 543 cases, 259 controls |
| Target | `dnm_case`: 1 = seen in a case, 0 = seen only in controls |
| Metric | **AUROC over all variants (headline)**; also the mean AUROC per gene (only 42 genes have both classes) |

**Genes.** One draw of 300 genes (`numpy` `default_rng(0)`, without replacement) from the 3,977 eligible genes:
`ap_test_set` genes whose protein is at most 2,700 aa and has at least one variant placed on it.

**Caution: gene identity.** In the `ap_test_set` every case is in a known developmental-disorder gene, while
controls are in any gene. So most genes have only cases or only controls, and knowing the gene alone separates
much of the data: a gene-only baseline (does the gene have another case variant) scores a higher AUROC than any
variant effect predictor. Report it next to every model.

One archive per gene: `<ENTRY>_DENOVO_2025.splits.pgdata`. The coordinates in the source are GRCh37
(`variant_id_grch37`); `variant_id` is the GRCh38 SNV. Layout and columns:
[docs/clinical/benchmarks.md](../../../docs/clinical/benchmarks.md). `denovo.csv` is the index.
