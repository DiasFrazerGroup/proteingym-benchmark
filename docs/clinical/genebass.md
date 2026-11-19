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
| `score_alphamissense`, `score_alphamissense_source` | score layer | AlphaMissense (Zenodo 10.5281/zenodo.8208688) by UniProt + the row's SNV on the canonical transcript; for LPA (outdated AlphaMissense canonical sequence) the SNV on a non-canonical transcript with the same protein change; else by protein change |
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
- **MSA headers (decided 2026-11-17):** the archives carry the original ColabFold a3m, headers included. Each
  header holds the MMseqs2 hit annotation: name, bit score, identity, e-value, query start/end/length, target
  start/end/length. The alignment models need the bit scores to follow the EVE protocol (bit-score threshold
  selection, model branch `docs/clinical/models_long_proteins.md` §3.5-3.6). proteingym-base `MSA.dump` rewrites
  headers as `seq_<i>`, so `scripts/clinical/msa/restore_msa_headers.py` puts the original file back after the layer
  step, checking that every sequence is identical and in the same order and that headers are unique. Loaded with
  proteingym-base, `MSA.value` has no headers; `MSA.file` holds the raw a3m with them. Upstream issue to report:
  `MSA.dump` drops headers (like the pLDDT B-factors).

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
- **Score layer sources:** AlphaMissense and popEVE are read from their own complete releases, not dbNSFP, whose
  per-transcript fields left 17% of rows without AlphaMissense (17 proteins with none) and 54% without popEVE.
  Coverage now: AlphaMissense 100% of rows (66/66 proteins), popEVE 73.6% (61/66), CADD 99.9%. popEVE has no
  protein over 2,302 aa, so APOB, LRP2, PKD1, SPTA1 and STAB1 have none. Details and validation:
  [scripts/clinical/layers/README.md](../../scripts/clinical/layers/README.md),
  [layers_sources_progress.md](layers_sources_progress.md).
- **AlphaMissense vs the VESM paper:** identical on every row except the UGT1A family (O60656, P19224, P22309,
  P22310, Q9HAW7, Q9HAW9: in the shared exons the paper has UGT1A10's value for the SNV; ours is the archive's own
  protein) and LPA (the paper's values match neither AlphaMissense table; ours are AlphaMissense on the LPA transcript
  matching the current UniProt sequence, as dbNSFP had).
- **Licences:** popEVE is CC BY 4.0. AlphaMissense is CC BY-NC-SA 4.0. CADD and dbNSFP are for non-commercial use.
  PrimateAI-3D cannot be redistributed. The score columns therefore need a licence decision before any public
  release.

## Final dataset (<= 2,700 aa, filtered alignments)

Decided 2026-11-19 (user): the benchmark is capped at **2,700 aa**, and the archives carry **only the alignment rows
the models use** instead of the full ColabFold search. This is the final file set:

| | |
|---|---|
| Index | [`genebass_le2700.csv`](../../datasets/clinical/genebass/genebass_le2700.csv): 100 archives (pairs), 63 proteins, 6,066 rows; paper set (`in_vesm_paper`) 94 pairs |
| Per-protein alignment filter | [`genebass_le2700_msa_filter.tsv`](../../datasets/clinical/genebass/genebass_le2700_msa_filter.tsv) |
| Archives | `$WORK/genebass/final_le2700/` (100 archives, 868 MB; the 109 full archives in `splits_layers/` are 7.1 GB) |
| Groups | Lipid Metabolism 34, Liver Function 38, Renal Function 11, Other 17 |
| Dropped (> 2,700 aa) | APOB (4,563 aa, 5 pairs), LRP2 (4,655, 3), PKD1 (4,303, 1): 9 pairs, 1,725 rows |

- **Why 2,700 aa:** every protein up to 2,700 aa has exactly one full-length AlphaFold DB structure (checked on all
  100 archives), so the structure models need no fragments, and every alignment model run so far fits. The model
  side (branch `clinical/model-api`, `docs/clinical/models_long_proteins.md`) documents what each model would need
  beyond the cap.
- **Alignment = EVE-protocol filter of the ColabFold a3m**, the same code and parameters the alignment models use
  (`msa_filter.py`, vendored byte for byte from the model branch, git blob `d9bd2d4b…`; `filter_msas.py` refuses to
  run if the copy differs). Rows: fragments (> 50% gaps) and rows with non-standard residues removed, then the lowest
  bit-score threshold (MMseqs2 bits / L) giving L_cov >= 0.8 L and 10 L <= N <= 100k, relaxed to L_cov >= 0.7 L and
  N <= 200k, then the deepest alignment with L_cov >= 0.7 L. Each kept row is stored as its exact original a3m record
  (header with bit score, insertions), in the original order, query first. All columns are kept.
- **Result:** N from 123 to 99,164 rows (median 23,545; median N/L 37); a3m text 27.8 GB -> 1.9 GB. Rules: strict 44
  proteins, relaxed 4, low_depth 15. Index columns: `msa_n_sequences` (filtered N), `msa_n_sequences_colabfold`,
  `msa_bits_threshold`, `msa_N_over_L`, `msa_Lcov_frac`, `msa_filter_rule`, `msa_quality`.
- **Weak alignments (`msa_quality = low_depth`, N < 10 L; scores kept, flag for interpretation):** LPA (N 123,
  N/L 0.06: kringle repeats), PCSK9 2.0, GPLD1 2.2, CGNL1 2.8, PAPPA2 3.3, STAB1 4.0, B4GALNT3 4.3, MROH2A 5.2,
  MAP3K15 6.4, ANGPTL3 6.6, ALB 6.7, NPC1L1 8.2, RBM47 8.4, IGF1R 9.0, ABCA1 9.7. (EXOC3L4, flagged in an earlier
  no-bit-score run, keeps 13,224 rows with bit scores: ok.)
- **Score layers unchanged.** AlphaMissense and CADD cover all 63 proteins. **popEVE has no scores for SPTA1 and
  STAB1** (its release stops at 2,302 aa): those 3 pairs (322 rows) have NaN popEVE, so popEVE covers 97 of 100 pairs.
- **The model-side filter on these archives is a no-op** (idempotent: it keeps every row), so `*_filtered` models
  give the same predictions as on the full archives (RSALOR on APOA1, PCSK9, LDLR: max difference 4e-16), at a
  fraction of the memory (0.3–0.5 GB instead of ~13 GB, most of it the MSA load).

**Rebuild** (part of `scripts/clinical/msa/rebuild.sbatch`, after the header restore; alone:
`sbatch scripts/clinical/msa/filter_msas.sbatch`, 1.5 min, 71 GB peak on 8 workers):

```shell
$PY scripts/clinical/msa/filter_msas.py --index datasets/clinical/genebass/genebass_fig6b.csv \
  --archives $WORK/genebass/splits_layers --msa-dir $WORK/genebass_msas --out $WORK/genebass/final_le2700 \
  --out-index datasets/clinical/genebass/genebass_le2700.csv \
  --out-filter datasets/clinical/genebass/genebass_le2700_msa_filter.tsv --max-length 2700
$PY scripts/clinical/msa/check_final.py --index datasets/clinical/genebass/genebass_le2700.csv \
  --final $WORK/genebass/final_le2700 --source $WORK/genebass/splits_layers --msa-dir $WORK/genebass_msas \
  --filter-table datasets/clinical/genebass/genebass_le2700_msa_filter.tsv     # must print CHECK FINAL OK
```

`check_final.py` (job 32028275, `CHECK FINAL OK`): every non-MSA entry (manifest, assays, sequences, structure,
layers, slices) byte-identical to `splits_layers/`; one structure per archive; each MSA an ordered subset of the
original a3m (exact records); filter idempotent; filter table rows and sha256 match.

## Files and locations

| What | Where |
|---|---|
| **Final benchmark set** (<= 2,700 aa, filtered alignments) | `$WORK/genebass/final_le2700/` (868 MB) + `datasets/clinical/genebass/genebass_le2700.csv` |
| Genebass tables + comparison | `$WORK/genebass_api/{genebass_api_hits.csv, genebass_api_betas.csv.gz, compare_vs_paper.md}` |
| Archives, no MSA (Genebass source) | `$WORK/genebass_api/build/splits/`, with layers in `splits_layers/` |
| Full archives (109 pairs, full ColabFold MSAs with headers, layers, `in_vesm_paper` flags); source of the final set | `$WORK/genebass/splits_layers/` (7.1 GB); index `genebass_fig6b.csv` (its `archive` column now points to `splits_layers/`) |
| Same archives without layers, uncompressed | `$WORK/genebass/splits/` (63 GB; delete once the paper-set check passes) |
| MSAs (66 proteins, 40 GB) | `$WORK/genebass_msas/` |
| Earlier paper-sourced archives (superseded) | `datasets/clinical/genebass/splits/` (repo, git-ignored), `$WORK/genebass_nomsa_layers/` |

`WORK=/no_backup/diasfrazer/egalpern_pgb_clinical`.

Code, all in `scripts/clinical/`:

- `genebass_api.py`: Genebass tables and comparison.
- `build_genebass.py`: archives.
- `check_genebass.py`: Fig. 6b reproduction.
- `msa/`: MSA chain; `restore_msa_headers.py` (a3m headers back into the archives); `filter_msas.py`,
  `check_final.py` and the vendored `msa_filter.py` (final <= 2,700 aa set).
- `layers/`: score and mask layers.

The `.sbatch` wrappers run each step on genoa64.

## Open

- **PrimateAI-3D:** not for now (decided 2026-11-16); would need the licensed file from Illumina.
- **Hosting and licences** of the archives and score layers (see Caveats).

Done since: MSAs for the 5 new proteins (CST3, GPLD1, GPT, SHBG, SLC22A12; job 31810198); headline set and
`in_vesm_paper` flags; structures kept as AFDB v6 with fragments.
