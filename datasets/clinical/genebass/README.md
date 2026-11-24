# Genebass

How strongly does a model's score track a missense variant's effect on a blood biomarker?

| | |
|---|---|
| Variants | [Genebass](https://doi.org/10.1016/j.xgen.2022.100168) single-variant associations in UK Biobank exomes, selected as in the VESM paper's Fig. 6b: gene–phenotype pairs with missense SKAT-O P ≤ 1e-10, variants with P < 0.05, pairs with more than 25 variants |
| Proteins | up to 2,700 aa |
| Size | 100 gene–phenotype pairs, 63 proteins, 5,995 variants |
| Target | `beta`: single-variant effect on the trait (SD units); its sign is the direction on the trait, not on fitness |
| Metric | per pair, Pearson r(beta, score) and its P value; **headline = mean −log10 P over the pairs**, also per biomarker group |

One archive per gene–phenotype pair: `<ENTRY>_GENEBASS_2022_<PHENOTYPE>.splits.pgdata`, e.g.
`APOA1_HUMAN_GENEBASS_2022_HDL_CHOLESTEROL`. `pair_in_vesm_paper` (pair) and `in_vesm_paper` (variant) flag the VESM
paper's own set: 94 pairs, 5,252 variants.

Layout and columns: [docs/clinical/benchmarks.md](../../../docs/clinical/benchmarks.md). Source and selection:
[docs/clinical/genebass.md](../../../docs/clinical/genebass.md). `genebass.csv` is the index.
