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

Each variant must sit on a UniProt wild type. Rules, in order:

1. **The paper's mapping table**, where it covers the variant: UniProt ID, position and sequence. The wild-type
   residue is checked.
2. **The variant's own GRCh38 SNV**, looked up in the AlphaMissense hg38 table (canonical isoforms), with the
   wild-type residue checked against the sequence.
3. **A constant Ensembl-to-UniProt offset**, only when at least 10 rows of the gene all match under it.
4. **Otherwise the row is dropped** and listed in `build/unmapped_rows.tsv`.

In the current build:

- 109 of the 111 pairs survive. The SLCO1B3-SLCO1B7 readthrough gene has no canonical Swiss-Prot entry, so its 2
  pairs (110 rows) are dropped.
- 23 more rows across 12 genes can't be placed and are dropped.
- The column `mapping` (`paper_uniprot`, `am_snv_index`, `swissprot_offset±k`) records the route per row.

Result: **109 archives, 66 proteins, 7,791 variant rows.** By group: Lipid 39, Liver 38, Renal 15, Other 17.

## Data files

### Archive layout

One archive per gene–phenotype pair, as PG2 does one per (protein, condition):
`<ENTRY>_GENEBASS_2022_<PHENOTYPE>.splits.pgdata`, e.g. `APOB_HUMAN_GENEBASS_2022_LDL_DIRECT`.

```
<name>.splits.pgdata = zip(dataset.pgdata, slices.json)       # trivial 2-fold split; zero-shot ignores it
dataset.pgdata       = zip(manifest.lock, assays/<name>.csv, sequences/<UNIPROT>.fasta,
                           structures/AF-<UNIPROT>-F<k>-model_v6.pdb [all AFDB fragments if > 2,700 aa],
                           msas/<UNIPROT>_colabfold_uniref100.a3m)
```

### Assay CSV

One row per (variant, protein change). Rows are not merged: 7 protein changes appear twice in a pair (different
codons), as in the paper.

| Column | Kind | Content |
|---|---|---|
| `sequence` | input | full-length mutated sequence |
| `beta` | **target** | Genebass single-variant effect (SAIGE; SD of the inverse-rank-normalised trait); its sign is the direction on the trait |
| `mutant` | id | e.g. `P4T`, UniProt numbering |
| `pval`, `allele_count`, `allele_number`, `allele_frequency` | association | Genebass single-variant P; cohort allele counts |
| `variant_id`, `hgvsp` | id | GRCh38 `chr-pos-ref-alt`; Ensembl protein HGVS |
| `mapping` | provenance | how the variant was placed on UniProt |
| `in_vesm_paper` | subset flag | row is in the paper's exact Fig. 6b set (same pair and SNV); `False` for rows only Genebass has |
| `uniprot_pos` | layer | UniProt position(s) |
| `mask_plddt`, `mask_plddt_source`, `mask_protein_length` | mask layer | pLDDT of the site (AFDB v6, see below); full canonical length |
| `score_alphamissense`, `score_alphamissense_source` | score layer | AlphaMissense hg38 table (Zenodo 10.5281/zenodo.8208688) by UniProt + the row's SNV, else by protein change |
| `score_popeve` | score layer | popEVE v1.1 (UKBB 2025-03-12) for the UniProt entry's RefSeq protein (from Swiss-Prot DR lines; NP used is in the column description) |
| `score_cadd_phred`, `score_snv_source` | score layer | CADD from dbNSFP 5.3.1a for the row's SNV; PrimateAI-3D once the licensed file is available |

### Manifest and inputs

- **Assay variables:** `phenotype_code` (UK Biobank field) and `biomarker_group`.
- **Wild type:** the paper's UniProt sequence, or Swiss-Prot 2026_03 for proteins the paper does not map. The APOB
  wild type differs from Swiss-Prot at Y1422C.
- **Structures:** AlphaFold DB v6, each with its residue span; proteins over 2,700 aa (APOB, LRP2, PKD1) come as
  their AFDB fragments F1..Fn. Decided 2026-11-16: keep AFDB v6 with fragments, no AF2 re-folding.
- **MSA:** the PG2 protocol, ColabFold on UniRef100 2026_01 with no post-alignment filtering. This is the
  `colabfold_search` monomer search with filtering off. Because UniRef100 is not clustered there is no cluster
  expansion, no environmental database, and the hit limit is raised to 1M. Depths reach 1.16M sequences.

### Index file

`genebass_fig6b.csv` has one row per archive: pair, gene, UniProt ID, phenotype, group, Genebass gene-level results
(missense and pLoF SKAT-O P, burden betas), counts, structures and MSA depth. `in_vesm_paper` marks the paper's
103 pairs and `n_variants_vesm_paper` counts their paper rows.

### Headline set and the paper's exact set

Decided 2026-11-16: the **headline metric uses every archive and every row** (109 pairs, Genebass-sourced). The
paper's Fig. 6b set is recovered from the data alone, without the paper's files: archives with `in_vesm_paper`
in the index, and within them the rows with `in_vesm_paper`. `check_genebass.py --paper-set` does exactly that and
must print `REPRODUCED`; the evaluation stage should offer the same switch.

### Caveats

- **pLDDT:** proteingym-base sets structure B-factors to 0 inside archives, so pLDDT is read from the original AFDB
  files.
- **Score layer sources:** AlphaMissense and popEVE are read from their own complete releases, not dbNSFP. The
  dbNSFP per-transcript fields left 17% of rows without AlphaMissense (17 proteins with none) and 54% without popEVE;
  see [scripts/clinical/layers/README.md](../../scripts/clinical/layers/README.md) and the coverage numbers in
  [layers_sources_progress.md](layers_sources_progress.md).
- **Licences:** popEVE is CC BY 4.0. AlphaMissense is CC BY-NC-SA 4.0. CADD and dbNSFP are for non-commercial use.
  PrimateAI-3D cannot be redistributed. The score columns therefore need a licence decision before any public
  release.

## Files and locations

| What | Where |
|---|---|
| Genebass tables + comparison | `$WORK/genebass_api/{genebass_api_hits.csv, genebass_api_betas.csv.gz, compare_vs_paper.md}` |
| Archives, no MSA (Genebass source) | `$WORK/genebass_api/build/splits/`, with layers in `splits_layers/` |
| **Final archives** (Genebass source, MSAs, layers, `in_vesm_paper` flags) | `$WORK/genebass/splits_layers/` (compressed, 6.6 GB), from rebuild job 31919891 |
| Same archives without layers, uncompressed | `$WORK/genebass/splits/` (63 GB; delete once the paper-set check passes) |
| MSAs (66 proteins, 40 GB) | `$WORK/genebass_msas/` |
| Earlier paper-sourced archives (superseded) | `datasets/clinical/genebass/splits/` (repo, git-ignored), `$WORK/genebass_nomsa_layers/` |

`WORK=/no_backup/diasfrazer/egalpern_pgb_clinical`.

Code, all in `scripts/clinical/`:

- `genebass_api.py`: Genebass tables and comparison.
- `build_genebass.py`: archives.
- `check_genebass.py`: Fig. 6b reproduction.
- `msa/`: MSA chain.
- `layers/`: score and mask layers.

The `.sbatch` wrappers run each step on genoa64.

## Open

- **PrimateAI-3D:** not for now (decided 2026-11-16); would need the licensed file from Illumina.
- **Hosting and licences** of the archives and score layers (see Caveats).

Done since: MSAs for the 5 new proteins (CST3, GPLD1, GPT, SHBG, SLC22A12; job 31810198); headline set and
`in_vesm_paper` flags; structures kept as AFDB v6 with fragments.
