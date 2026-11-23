# ClinVar, de novo and gnomAD benchmarks

Status: 2026-11-23.

Three clinical variant-effect benchmarks in the PG2 format, one `.splits.pgdata` archive per gene, built like the
Genebass archives ([genebass.md](genebass.md)). Each archive holds benchmark data only:
- the wild type, AFDB structure and MSA;
- the mutated sequences;
- the target;
- provenance.

There are no model scores.

| Benchmark | Question | Genes | Variants | Target | Headline metric |
|---|---|---|---|---|---|
| [`clinvar`](../../datasets/clinical/clinvar/) | pathogenic vs benign (ClinVar, ≥ 2 stars) | 158 | 14,267 | `label` | AUROC, all variants |
| [`denovo`](../../datasets/clinical/denovo/) | de novo variants in severe developmental disorder cases vs controls (popEVE paper, `ap_test_set`) | 300 random | 722 | `dnm_case` | AUROC, all variants |
| [`gnomad`](../../datasets/clinical/gnomad/) | rarity in the population (gnomAD v4.1 rare missense) | 300 random | 127,261 | `log10_af` | Spearman ρ, all variants (lower = better) |

Each benchmark also reports the mean metric per gene. The variant tables come from the clinical-curation repo; the
selection rules are in each benchmark's README.

## Archives

`<ENTRY>_<SOURCE>_<YEAR>.splits.pgdata`, e.g. `BRCA1_HUMAN_CLINVAR_2026`, `ABCA4_HUMAN_DENOVO_2025`,
`MITOS_HUMAN_GNOMAD_2024`.

- **Wild type:** Swiss-Prot 2026_03 canonical sequence.
- **Structure:** AlphaFold DB v6. ClinVar has 28 proteins over 2,700 aa, given as AFDB fragments with their residue
  spans. The two random-gene benchmarks have none: the longest is 2,426 aa.
- **MSA:** ColabFold on UniRef100 2026_01, no post-alignment filtering (the Genebass recipe,
  `scripts/clinical/msa`). Every protein has one.
- **Split:** trivial 2-fold (zero-shot models ignore it).

| Assay column | Content |
|---|---|
| `sequence` | full-length mutated sequence |
| `label` / `dnm_case` / `log10_af` | **target** |
| `mutant` | UniProt numbering, e.g. `R44W` |
| `mutant_source` | the change in the source's own protein numbering (RefSeq / MANE) |
| `variant_id` | GRCh38 `chrom-pos-ref-alt` (ClinVar: all SNVs causing the change, `;`-joined) |
| `mapping` | how the variant was placed on UniProt: `am_snv_index` or `source_offset±k` |
| source fields | ClinVar `hgvsp`, `source_transcript`, `clnsig`, `max_stars`, `n_records`; de novo `variant_id_grch37`, `source_protein`, `n_sdd`, `n_control`; gnomAD `af`, `AC_non_ukb`, `AN_non_ukb`, `AC`, `AN`, `enst`, `mane_select` |

**Placing variants on UniProt:**
- **By SNV (99.6%):** through the AlphaMissense tables. The de novo source coordinates are GRCh37 and go through the
  hg19 table; their GRCh38 position is recovered.
- **Otherwise:** a constant offset on the source protein numbering, fitted per source isoform.
- **Unplaced:** 0.7% (ClinVar), 0.8% (de novo) and 1.4% (gnomAD) of variants can't be placed, because they lie on
  exons the UniProt canonical protein doesn't have. They are listed in the build output.

The archives are generated, not committed: `/no_backup/diasfrazer/cpugh/pgb_clinical_benchmarks/benchmarks/<benchmark>/splits/`
(34 GB). Each folder's `<benchmark>.csv` is the index, one row per archive.

## Reference results

Six published predictors (clinical-curation's scores, not stored in the archives), on the variants all six score.
Pooled = headline.

| Model | ClinVar AUROC (per gene) | De novo AUROC (per gene, 35 genes) | gnomAD ρ (per gene) |
|---|---|---|---|
| AlphaMissense | 0.963 (0.972) | 0.706 (0.618) | −0.118 (−0.129) |
| VESM-3B | **0.964** (0.963) | 0.685 (0.633) | −0.070 (−0.086) |
| ESM1b | 0.945 (0.945) | 0.698 (0.578) | −0.056 (−0.070) |
| CADD v1.7 | 0.909 (0.936) | 0.692 (0.659) | −0.046 (−0.043) |
| GPN-MSA | 0.889 (0.926) | **0.709** (0.482) | **−0.169** (−0.171) |
| popEVE | 0.886 (0.957) | 0.704 (0.569) | −0.077 (−0.081) |
| *gene-only baseline* | *0.780* | ***0.876*** | *−0.089* |

- **ClinVar:** 9,697 variants scored by all six. Identical to clinical-curation's numbers.
- **De novo:** 700 variants. **Gene identity alone beats every model.** Every case in the `ap_test_set` is in a known
  developmental-disorder gene, and controls are in any gene; only 35 of the 300 genes have both cases and controls.
  The baseline ("does the gene have another case variant") must be reported next to every model.
- **gnomAD:** 124,197 variants. Allele frequency also tracks how mutable a site is, which is why the DNA model
  GPN-MSA leads. The gene-only baselines for ClinVar and gnomAD are clinical-curation's (leave-one-out case rate;
  gnomAD missense o/e).

## Rebuild

From the repo root, on the CRG cluster (settings in `scripts/clinical/benchmarks/env.sh`):

```shell
bash scripts/clinical/benchmarks/submit.sh prepare   # select variants, place them on UniProt, list MSA queries
bash scripts/clinical/benchmarks/submit.sh msa       # MSAs for proteins without one, then build + check
bash scripts/clinical/benchmarks/submit.sh build     # build (one job per benchmark) + check
```

`check.py` verifies the archives:
- the rows match the prepared variants;
- every sequence carries exactly its mutation;
- the archives load with proteingym-base;
- ClinVar matches clinical-curation's reference.

It then writes the results table above.
