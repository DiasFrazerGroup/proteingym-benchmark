# Clinical benchmarks

Status: 2026-11-25. All 858 archives are built and pass `check.py`.

Four clinical variant-effect benchmarks in the PG2 format. Each is a set of `.splits.pgdata` archives laid out like
the DMS archives, with benchmark data only: no model scores.

| Benchmark | Question | Archives | Variants | Target | Headline metric |
|---|---|---|---|---|---|
| [`clinvar`](../../datasets/clinical/clinvar/) | pathogenic vs benign (ClinVar, ≥ 2 stars) | 158 genes | 14,237 | `label` | AUROC, all variants |
| [`denovo`](../../datasets/clinical/denovo/) | de novo variants in severe developmental disorder cases vs controls (popEVE paper, `ap_test_set`) | 300 random genes | 802 | `dnm_case` | AUROC, all variants |
| [`gnomad`](../../datasets/clinical/gnomad/) | rarity in the population (gnomAD v4.1 rare missense) | 300 random genes | 130,198 | `log10_af` | Spearman ρ, all variants (lower = better) |
| [`genebass`](../../datasets/clinical/genebass/) | effect on blood biomarkers (Genebass, VESM Fig. 6b setting) | 100 gene–phenotype pairs, 63 proteins | 5,995 | `beta` | mean −log10 P of Pearson r per pair |

The selection rules are in each benchmark's README. `datasets/clinical/<benchmark>/<benchmark>.csv` is the index,
one row per archive.

## Random genes

`denovo` and `gnomad` use 300 random genes, drawn once:

1. **Eligible genes:** the gene's protein is at most 2,700 aa and has at least 1 (`denovo`) or 20 (`gnomad`)
   variants placed on it. This gives 3,977 genes for `denovo` and 17,564 for `gnomad`.
2. **One draw:** `numpy.random.default_rng(0).choice(sorted(eligible), 300, replace=False)`.

There is no second draw and no replacement of genes.

## Proteins and variants

Every variant is placed on a Swiss-Prot (2026_03) canonical protein the same way:

- The variant's SNV is looked up in the AlphaMissense SNV table, which gives the protein change on a UniProt
  canonical isoform. It is kept when the wild-type residue is the Swiss-Prot one.
- A gene's protein is the Swiss-Prot entry named after the gene, if its SNVs map to it; otherwise the entry most of
  its SNVs map to.
- Variants that don't land on that protein are dropped: 0.9% (ClinVar), 2.7% (de novo) and 3.8% (gnomAD) of the
  source variants, before the draw.
- The de novo source coordinates are GRCh37 and go through the AlphaMissense hg19 table.
- Genebass first uses the VESM paper's own mapping, so the paper's variants are all kept (94 pairs, 5,252 rows,
  flagged `in_vesm_paper`); its other variants are placed by SNV.

Protein length: `denovo`, `gnomad` and `genebass` only have proteins up to 2,700 aa, which have one full-length
AlphaFold DB model. ClinVar keeps all its genes; 28 of its proteins are longer and come with the AlphaFold DB
fragments.

## Archives

`<ENTRY>_<SOURCE>_<YEAR>.splits.pgdata`: `BRCA1_HUMAN_CLINVAR_2026`, `ABCA4_HUMAN_DENOVO_2025`,
`MITOS_HUMAN_GNOMAD_2024`, `APOA1_HUMAN_GENEBASS_2022_HDL_CHOLESTEROL`.

| In the archive | Content |
|---|---|
| `assays/<DATASET>.csv` | `sequence` (full-length mutated sequence), the target, `mutant` (UniProt numbering, e.g. `R44W`), then the source's columns |
| `sequences/<DATASET>.fasta` | wild type: Swiss-Prot 2026_03 canonical sequence |
| `structures/<DATASET>.pdb` | AlphaFold DB v6 model (ClinVar proteins over 2,700 aa: `<DATASET>_F<k>.pdb` fragments, with their residue span in the manifest) |
| `msas/<DATASET>_colabfold.a3m` | ColabFold search on UniRef100 2026_04, no post-alignment filtering. The file is the search output byte for byte, so the hit headers (bit scores) are there |
| split | trivial 2-fold (zero-shot models ignore it) |

Source columns in the assay:

| Benchmark | Columns after `mutant` |
|---|---|
| `clinvar` | `variant_id` (GRCh38 `chrom-pos-ref-alt`; all SNVs causing the change, `;`-joined), `hgvsp`, `clnsig`, `max_stars`, `n_records` |
| `denovo` | `variant_id` (GRCh38), `variant_id_grch37` (the source coordinate), `n_sdd`, `n_control` |
| `gnomad` | `variant_id`, `af`, `AC_non_ukb`, `AN_non_ukb` |
| `genebass` | `variant_id`, `hgvsp`, `pval`, `allele_count`, `allele_number`, `allele_frequency`, `in_vesm_paper` |

Genebass archives also carry the pair's `phenotype`, `phenotype_code`, `biomarker_group`,
`genebass_analysis_id`, `skato_pval_missense` and `pair_in_vesm_paper` as assay variables.

Alignments are not filtered in the archives: models that need a filtered alignment filter it themselves.

## MSAs

One MSA per protein (809 proteins), searched on UniRef100 2026_04 (227 million sequences) with
[scripts/clinical/msa.sh](../../scripts/clinical/msa.sh): the ColabFold monomer search with no post-alignment
filtering. The two prefilters of MMseqs2 were used, by protein length:

| Proteins | Prefilter | Where |
|---|---|---|
| up to 1,250 aa (653) | GPU, ungapped (`GPU=1`); at most 524,288 candidates per protein | MareNostrum 5, 4 × H100 per node |
| over 1,250 aa (156) | CPU, k-mer (`-s 8`); at most 1,000,000 candidates | CRG cluster |

The steps after the prefilter are the same. The split is practical: after the GPU prefilter, the alignment of one
protein runs on one CPU thread and its time grows with the square of the protein length (about 8 minutes per pass
at 1,000 aa, four passes), which didn't fit the 2-hour jobs available on MareNostrum 5.

MSA depth, median (minimum) sequences: ClinVar 179,394 (1,764); de novo 99,294 (286); gnomAD 71,501 (114);
Genebass 187,656 (3,042). 17 proteins have fewer than 1,000 sequences (2 in de novo, 15 in gnomAD). Each archive's
depth is in the index column `msa_n_sequences`.

## Reference results

Six published predictors (clinical-curation's scores, not stored in the archives), on the variants all six score.
Pooled over all variants = headline; the mean per gene is in brackets.

| Model | ClinVar AUROC | De novo AUROC | gnomAD ρ (lower = better) |
|---|---|---|---|
| AlphaMissense | 0.963 (0.972) | 0.755 (0.671) | −0.128 (−0.134) |
| VESM-3B | **0.964** (0.963) | 0.719 (0.697) | −0.075 (−0.089) |
| ESM1b | 0.945 (0.944) | 0.718 (0.661) | −0.064 (−0.074) |
| CADD v1.7 | 0.909 (0.936) | 0.691 (0.582) | −0.051 (−0.049) |
| GPN-MSA | 0.889 (0.926) | 0.742 (0.646) | **−0.166** (−0.176) |
| popEVE | 0.886 (0.957) | **0.760** (0.681) | −0.099 (−0.083) |
| *gene-only baseline* | | ***0.854*** | |

- **ClinVar:** 9,697 variants in 122 genes scored by all six.
- **De novo:** 768 variants in 289 genes; the per-gene mean covers the 41 genes with both classes. **Gene identity
  alone beats every model:** every case in the `ap_test_set` is in a known developmental-disorder gene, while
  controls are in any gene. The baseline (does the gene have another case variant) must be reported next to every
  model.
- **gnomAD:** 120,225 variants in 287 genes. Allele frequency also tracks how mutable a site is, which is why the
  DNA model GPN-MSA leads.

Written by `check.py` to `$WORK/reference_results.tsv`. Genebass has no reference scores here.

## Where the archives are

The archives are generated, not committed (31.6 GB):

- CRG: `/no_backup/diasfrazer/cpugh/pgb_clinical/<benchmark>/splits/`
- MareNostrum 5: `/gpfs/projects/ehpc706/cpugh/pgb_clinical/<benchmark>/splits/`

They load with proteingym-base 0.2.0 (`Subsets.from_path`).

## Rebuild

From the repo root. Paths and tools are set in [scripts/clinical/env.sh](../../scripts/clinical/env.sh).

```shell
bash scripts/clinical/run.sh prepare   # select variants, place them on UniProt, list proteins without an MSA
bash scripts/clinical/run.sh msa       # MSAs for those proteins (GPU=1 for the GPU search)
bash scripts/clinical/run.sh build     # archives and index CSVs
bash scripts/clinical/run.sh check     # verify the archives, report the reference models
```

See [scripts/clinical/README.md](../../scripts/clinical/README.md) for what each step reads and writes.
