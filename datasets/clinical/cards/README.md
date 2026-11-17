# Clinical benchmark cards as PG2 archives

One set of `.splits.pgdata` archives per clinical benchmark card, one archive per gene, so the unchanged zoo can
score them. The cards (variant sets, filters, metrics) are defined in the clinical-curation repo
(`scripts/16_benchmark_cards.py`, `benchmark_cards/`, `DECISIONS.md`); this folder holds the per-card archive
indexes. The archives themselves are generated (not committed) under
`WORK=/no_backup/diasfrazer/cpugh/pgb_clinical_cards/cards/<card>/splits/`.

## Cards

| Card | Variants | Target | Evaluation (as the card) |
|---|---|---|---|
| `clinvar_v0` | ClinVar germline missense, ≥2 stars, P/LP vs B/LB, genes with ≥10 of each | `label` (1 = P/LP) | pooled AUROC; per-gene mean AUROC over genes with both classes |
| `dnm_aptest_ge1` | popEVE ST5 `ap_test_set` de novo missense, genes with ≥1 case and ≥1 control | `dnm_case` (1 = SDD case) | as ClinVar |
| `dnm_aptest_ge2` | same, ≥2 and ≥2 (= DNM v1) | `dnm_case` | as ClinVar |
| `dnm_aptest_random300` | same, 300 random `ap_test_set` genes (seed 0; `dnm_aptest_random300_genes.tsv`). New: no card in clinical-curation yet | `dnm_case` | as ClinVar (most genes have one class: pooled AUROC is the usable metric) |
| `gnomad_af_clinvar_dnm_ge1` | gnomAD v4.1 exome PASS MANE-Select missense SNVs, non-UKB AF < 0.1%, AN ≥ 50% of max; genes of `clinvar_v0` ∪ `dnm_aptest_ge1` | `log10_af` (non-UKB) | Spearman ρ(score, log10 AF), pooled and per-gene mean over genes with ≥20 variants (lower = better) |
| `gnomad_af_clinvar_dnm_ge2` | same, genes of `clinvar_v0` ∪ `dnm_aptest_ge2` | `log10_af` | as above |
| `gnomad_af_random300` | same, the card's 300 random genes | `log10_af` | as above |

The reference cards `dnm_aptest_full` and `gnomad_af_full` (proteome scale) are not archived.

Every archive keeps all of the card's rows for that gene (cards select genes on all rows).

## Archives

`<card>/splits/<ENTRY>_<SUFFIX>.splits.pgdata`, suffix `CLINVAR_V0`, `DNM_APTEST` or `GNOMAD_V4_1_RARE`
(e.g. `SCN2A_HUMAN_DNM_APTEST`). Cards are separate sets: a gene in two cards has an archive in each.

- **Wild type:** Swiss-Prot 2026_03 canonical sequence.
- **Structures:** AlphaFold DB v6; proteins over 2,700 aa as AFDB fragments F1..Fn, each with its residue span
  (as the Genebass archives).
- **MSA:** ColabFold on UniRef100 2026_01, no post-alignment filtering (`scripts/clinical/msa`); Genebass MSAs are
  reused when their first record is the same wild type.
- **Split:** trivial 2-fold (zero-shot models ignore it).
- A protein without a usable AFDB model or MSA is archived without it; the index column `inputs_note` says why.

### Assay columns

| Column | Kind | Content |
|---|---|---|
| `sequence` | input | full-length mutated sequence |
| `label` / `dnm_case` / `log10_af` | **target** | see the card table |
| `mutant` | id | UniProt numbering, e.g. `R44W` |
| `mutant_source` | id | the change in the source's own protein numbering (RefSeq / MANE) |
| `variant_id` | id | GRCh38 `chrom-pos-ref-alt` (ClinVar: all SNVs causing the change, `;`-joined; DNMs: recovered from the GRCh37 source SNV, empty if not unique) |
| `mapping` | provenance | `am_snv_index` (SNV → UniProt protein change, AlphaMissense table) or `source_offset±k` |
| ClinVar: `hgvsp`, `source_transcript`, `clnsig`, `max_stars`, `n_records` | source | |
| DNM: `variant_id_grch37`, `source_protein`, `n_sdd`, `n_control` | source | popEVE ST5 coordinates are GRCh37 |
| gnomAD: `af`, `AC_non_ukb`, `AN_non_ukb`, `AC`, `AN`, `enst`, `mane_select` | source | |

**No model scores** are stored in the archives: they are benchmark data only (protein inputs, variants, targets).
The cards' reference scores (AlphaMissense, popEVE, CADD, ESM1b, VESM-3B, GPN-MSA) stay in clinical-curation and
in the working tables `$WORK/cards/<card>/variants.tsv.gz`, where `check_cards.py` uses them to reproduce the
card numbers. Cards measure performance on the variants scored by all six models; that subset is a property of the
reference predictors, not of the benchmark, so it is not flagged in the archives.

### Index files

`<card>.csv`: one row per archive (gene, UniProt, entry, length, rows, classes or singletons, complete rows,
mapping routes, structures, MSA depth, `inputs_note`, archive path).

## Rebuild

Scripts in `scripts/clinical/cards/` (settings in `env.sh`), all as sbatch jobs on genoa64, from the repo root:

```shell
bash scripts/clinical/cards/submit.sh prepare     # card rows -> UniProt-mapped rows, proteins, MSA queries
bash scripts/clinical/cards/submit.sh msa 150     # MSA search in chunks of 150 proteins, then build (afterok)
bash scripts/clinical/cards/submit.sh build       # build + check now (uses whatever MSAs exist)
```

`check_cards.py` compares every card with clinical-curation's `benchmark_cards/summary.tsv`: selection counts,
unmapped rows (≤3%), the six models' pooled and per-gene metrics recomputed from the archives, and the sequences.

## Licences

ClinVar is public domain and gnomAD summary data is released under CC0; the popEVE ST5 DNM table follows its
paper's terms (check before a public release). AlphaFold DB models are CC BY 4.0.
