# Genebass clinical datasets (VESM Fig. 6b setting)

Single-variant associations between missense variants and blood biomarkers in UK Biobank exomes
([Genebass](https://doi.org/10.1016/j.xgen.2022.100168)). The selection follows the VESM paper's Fig. 6b
exactly ([bioRxiv 10.1101/2025.03.09.642148](https://doi.org/10.1101/2025.03.09.642148), notebook
`DJZN25-Fig6.ipynb`). The task: for each gene–phenotype pair, how strongly does a model's variant score
correlate with the variants' effect on the trait?

Status: data files only. Long-protein and sparse scoring are handled separately in the model API, and the
pooled evaluation stage isn't written yet (see [docs/clinical_benchmark.md](../../../docs/clinical_benchmark.md)).

## Contents

| | |
|---|---|
| Pairs (one archive each) | 103 |
| Proteins | 61 (18 longer than 1022 aa; APOB, PKD1 and LRP2 longer than 2700 aa) |
| Variant rows | 6,977 (missense, single-variant P < 0.05) |
| Biomarker groups | Lipid Metabolism 39, Liver Function 35, Renal Function 13, Other 16 |

- [genebass_fig6b.csv](genebass_fig6b.csv) is the index: one row per archive, with the pair, UniProt ID,
  biomarker group, Genebass gene-level results (missense SKAT-O P, burden betas), counts, structures and
  archive path.
- The archives themselves (`splits/<name>.splits.pgdata`, 258 MB) and `build/` are generated, not committed
  (see [Hosting](#hosting)). `build/` holds the staged manifests, CSVs and a per-protein input store.

## Selection

From the paper's own tables, in the paper's order (cells 23–24):

1. pairs with gene-level missense SKAT-O P ≤ 1e-10;
2. variants with single-variant P < 0.05;
3. pairs with more than 25 variants left.

## Archive layout

There is one archive per gene–phenotype pair, as PG2 does one archive per (protein, condition), e.g.
`P53_HUMAN_EBERT_2019_ACTIVITY_DMSO`. Names are `<ENTRY>_GENEBASS_2022_<PHENOTYPE>`, e.g.
`APOB_HUMAN_GENEBASS_2022_LDL_DIRECT`.

- **Assay:** one per archive, one row per (variant, protein change), `sequence` = full-length mutated sequence.
  - Target `beta`: the Genebass single-variant effect (SAIGE, inverse-rank-normalised trait, SD units). Its sign
    is the direction on the trait, not on fitness.
  - Non-targets: `mutant` (UniProt numbering, e.g. `R3527Q`), `pval`, `allele_count`, `allele_number`,
    `allele_frequency`, `variant_id` (GRCh38 chr-pos-ref-alt), `hgvsp` (Ensembl protein) and `mapping`.
- **Assay variables:** `phenotype_code` (UK Biobank field) and `biomarker_group`.
- **Wild type:** the UniProt sequence the paper mapped variants onto, with Swiss-Prot metadata (entry name, Pfam).
- **Structures:** AlphaFold DB v6 models.
  - Proteins longer than 2,700 aa come as AFDB fragments F1..Fn (1,400 aa each, one starting every 200 aa).
    All fragments are included, each with its residue span on the wild type in `metadata` (`span_start`,
    `span_end`), plus `n_mismatch` and mean pLDDT.
  - Choosing a fragment per position is left to the model API.
- **MSA:** ColabFold MSA on UniRef100 2026_01 with no post-alignment filtering, the PG2 protocol (see
  [scripts/clinical/msa](../../../scripts/clinical/msa/README.md)). One per protein, the same file in all its
  archives, named `<archive>_colabfold`. Added when the builder gets `--msa-dir`.
- **Split:** one trivial 2-fold split, because the zero-shot runner needs a `.splits.pgdata` and ignores the folds.

## Evaluation (the paper's metric, to be implemented as the clinical stage)

- **Per archive:** Pearson r(beta, score) over the assay rows, and its two-sided P.
- **Fig. 6b top:** mean −log10 P over all 103 pairs.
- **Fig. 6b bottom:** mean −log10 P per biomarker group. The paper's heatmap shows Lipid, Liver and Renal; the
  16 `Other` pairs count only in the overall mean.

## Caveats

- **UGT1A5 mapping.** The paper's tables give no UniProt mapping for UGT1A5's 72 rows: their `hgvsp` uses
  ENSP00000303174, which is 2 residues shorter at the N-terminus. These rows are mapped onto Swiss-Prot P35504
  at position +2, where all 140 UGT1A5 variants' wild-type residues match (`mapping = swissprot_offset+2`). Every
  other row uses the paper's mapping (`paper_uniprot`).
- **APOB wild type.** It is the paper's sequence, which differs from Swiss-Prot 2026_03 at one residue (Y1422C).
- **Duplicate protein changes.** 7 protein changes occur twice within a pair (different codons, different betas).
  - Both rows are kept, as the paper counts both.
  - `Dataset.to_df()` would average them (6,970 rows), so the evaluation should read the assay rows.

## Build and check

```shell
# Python env with proteingym-base (requirements.txt)
python scripts/clinical/build_genebass.py \
  --fig6  <VESM figshare 10.6084/m9.figshare.31162471>/Fig6_data \
  --sprot uniprot_sprot_human.dat.gz \
  --afdb  <AlphaFold DB v6 human proteome PDBs>
python scripts/clinical/check_genebass.py --fig6 <...>/Fig6_data
```

The check reads the archives back with `proteingym.base.Subsets`, joins the paper's own model scores per row,
and recomputes Fig. 6b. Results:

- matches the notebook's code on the raw table to 1e-15;
- matches the paper's printed per-group values (cell 31) to 5e-7;
- matches the printed overall values (2 decimals) within rounding.

Build time is about 2.5 min.

## Hosting

`*.pgdata` files are Git LFS (`.gitattributes`). Following the DMS convention (full collection on Hugging Face,
a few examples in the repo), the archives are not committed here yet. Decide between Hugging Face and LFS
examples before pushing them.
