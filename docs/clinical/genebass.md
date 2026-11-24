# Genebass clinical benchmark (VESM Fig. 6b setting)

Status: 2026-11-11. Branch `clinical/genebass-data` of
[DiasFrazerGroup/proteingym-benchmark](https://github.com/DiasFrazerGroup/proteingym-benchmark).

## Task

**Question:** do model scores for missense variants track the variants' measured effects on blood biomarkers in
UK Biobank?

**Evaluation unit:** a gene–phenotype pair. For each pair, compute the Pearson r between each variant's
single-variant effect size (`beta`) and the model score, and its two-sided P.

**Headline metric:** mean −log10 P over all pairs (VESM paper Fig. 6b, top). It is also reported per biomarker
group: Lipid Metabolism, Liver Function and Renal Function. That is the paper's Fig. 6b bottom panel; pairs in the
fourth group, "Other", count only in the overall mean.

Reference: Dinh, Jang, Zaitlen, Ntranos (VESM), bioRxiv 10.1101/2025.03.09.642148, notebook `DJZN25-Fig6.ipynb`.

## Data source and selection

**Source:** Genebass (Karczewski et al., Cell Genomics 2022), 394,841 UK Biobank exomes, through the public API
behind app.genebass.org.

- Local copy: `/no_backup/diasfrazer/db/genebass/raw`, downloaded 2026-11-02.
- Pairs missing from that copy are fetched by `scripts/clinical/genebass_api.py`.

The VESM figshare tables are used only for the biomarker groups (the paper's own grouping) and for comparison.

**Selection rules** (they reproduce the paper's pair set):

1. **Phenotypes:** the paper's 27 blood-biochemistry phenotypes. That is the 31 Genebass "Blood biochemistry"
   phenotypes minus C-reactive protein, oestradiol, rheumatoid factor and the custom TG/HDL ratio. All are
   inverse-rank-normalised continuous traits.
2. **Candidate pairs:** the Genebass portal's top pLoF associations for those phenotypes, i.e. pairs where any pLoF
   gene-level test (SKAT-O, SKAT or burden) has P < 1e-6. This gives 344 pairs: the paper's 332 plus 12 it does not
   have (G6PC ×4, SLCO1B3-SLCO1B7 ×3, APOC3, CST3, GPLD1, PACC1, SHBG), flagged by `in_vesm_paper`.
3. **Fig. 6b filter:**
   - missense (`missense|LC`) SKAT-O P ≤ 1e-10;
   - then missense variants with single-variant P < 0.05;
   - then pairs with more than 25 such variants.

### Genebass vs the paper's tables

From `compare_vs_paper.md`:

- **All 174,073 of the paper's variant rows are in Genebass,** with identical beta and P. The paper's allele counts
  are Genebass's annotation (cohort) counts.
- **Gene-level missense SKAT-O P agrees to within 2e-11.**
- **Genebass has 5,729 more missense rows,** mostly in the 12 extra pairs.
- **Fig. 6b selection: Genebass 111 pairs / 7,924 rows; paper 103 / 6,977.**
  - All 103 paper pairs are included.
  - New pairs: CST3, GPLD1, SHBG, SLCO1B3-SLCO1B7 ×2 (from the 12 extra), and ALPL, GPT, SLC22A12 (paper pairs
    that pass with the extra rows).
- **Rows added to existing paper pairs:** LPA × 30790 (131 vs 60), UGT1A5 × 30660 (50 vs 31), UGT1A5 × 30840
  (59 vs 41), STAB1 × 30860 (108 vs 107).
- **Reproduction:** on the paper's own rows, the archives give the paper's Fig. 6b values exactly, overall and per
  group, e.g. VESM++ 2.472; per group 3.1187 / 2.1488 / 1.6384. Test: `check_genebass.py --paper-rows-only`
  (`REPRODUCED`, max difference 8.9e-16 against the notebook's code). The plain check reports MISMATCH on the
  Genebass set, as expected, because of the extra pairs and rows.

## Mapping variants to proteins

Each variant must sit on a Swiss-Prot wild type:

1. **The paper's mapping table**, where it covers the variant, with the wild-type residue checked.
2. **Otherwise the variant's own GRCh38 SNV**, looked up in the AlphaMissense hg38 table, on the gene's protein.
3. **Otherwise the row is dropped.**

The SLCO1B3-SLCO1B7 readthrough gene has no Swiss-Prot entry of its own, so its 2 pairs are dropped. LPA isn't in
the AlphaMissense hg38 table, so it keeps only the paper's 60 rows.

## Data files

**100 archives (pairs), 63 proteins, 5,995 variant rows.** By group: Lipid 34, Liver 38, Renal 11, Other 17.

- **Proteins up to 2,700 aa only.** This drops APOB, LRP2 and PKD1 (9 pairs). Every remaining protein has one
  full-length AlphaFold DB model.
- **Paper set:** 94 pairs and 5,252 rows are the VESM paper's own. They are flagged by the assay variable
  `pair_in_vesm_paper` and the assay column `in_vesm_paper`.
- **Layout:** the same as the other clinical benchmarks, see [benchmarks.md](benchmarks.md). The archives have no
  score layers, and the MSA is the unfiltered ColabFold search output with its headers.
- **Index:** [`genebass.csv`](../../datasets/clinical/genebass/genebass.csv), one row per archive.
- **Duplicate protein changes:** a few protein changes occur twice within a pair (different SNVs, different betas).
  Both rows are kept. `Dataset.to_df()` would average them, so the evaluation should read the assay rows.

Build: `BENCHMARKS=genebass bash scripts/clinical/run.sh prepare`, then `build` and `check`
([scripts/clinical/README.md](../../scripts/clinical/README.md)).

The earlier build (109 pairs, score layers, filtered alignments, the Fig. 6b reproduction check) is on branch
`clinical/genebass-data`.
