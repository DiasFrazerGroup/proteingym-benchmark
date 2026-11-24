# ClinVar

Can a model tell pathogenic missense variants from benign ones?

| | |
|---|---|
| Variants | ClinVar germline missense variants (November 2026 release), ≥ 2 review stars, pathogenic / likely pathogenic vs benign / likely benign |
| Genes | 158 genes with ≥ 10 variants of each class |
| Size | 14,237 variants |
| Target | `label`: 1 = pathogenic / likely pathogenic, 0 = benign / likely benign |
| Metric | **AUROC over all variants (headline)**; also the mean AUROC per gene |

One archive per gene: `<ENTRY>_CLINVAR_2026.splits.pgdata`, e.g. `BRCA1_HUMAN_CLINVAR_2026`. The archive layout and
columns are described in [docs/clinical/benchmarks.md](../../../docs/clinical/benchmarks.md).
`clinvar.csv` is the index (one row per archive).
