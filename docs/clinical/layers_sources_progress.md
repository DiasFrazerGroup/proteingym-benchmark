# Layers from their own sources: status (branch clinical/layers-sources)

AlphaMissense and popEVE layers now come from their own complete releases instead of dbNSFP; CADD stays on dbNSFP.
Design and usage: [scripts/clinical/layers/README.md](../../scripts/clinical/layers/README.md).

## Code
- `sources.toml`: `[alphamissense]` (hg38 `table`, `protein_table` from aa_substitutions, `isoform_table` from
  isoforms_hg38), `[popeve]` (`zip`, `sequences`); `[snv_index]` removed (the AlphaMissense hg38 parquet is the SNV
  index); `[dbnsfp]` only for CADD.
- `build_snv_index.py` -> `$WORK/layers/alphamissense_{hg38,aa_substitutions,isoforms_hg38}.parquet`
  (jobs 31920801, 31921227, 31921954). `build_popeve_index.py` -> `$WORK/layers/popeve_sequences.parquet`
  (job 31921332; 18,332 proteins, max 2,302 aa). The old `am_hg38_snv_index.parquet` is untouched
  (`build_genebass.py` still defaults to it).
- `add_layers.py`: classes `AlphaMissense`, `PopEVE`, `Dbnsfp`; new column `score_alphamissense_source`
  (`snv` | `isoform_snv` | `protein` | `none`); RefSeq / sequence notes in column descriptions and
  `layers_summary.tsv`; `layers_am_snv_mismatch.tsv`.
- `compare_layers.py OLD NEW [--am-ref DIR]`: coverage / agreement report.

## Validation (109 no-MSA Genebass archives, 7,791 rows, 66 proteins; job 31922061)
Old = `$WORK/genebass_api/build/splits_layers` (dbNSFP), new = `$WORK/layers_v2_test/genebass_nomsa`; full report
in `$WORK/layers_v2_test/compare_nomsa.txt`.

| Layer | Rows old -> new | Archives with any | Proteins with any |
|---|---|---|---|
| AlphaMissense | 6,472 (83.1%) -> 7,791 (100%) | 83 -> 109 | 49 -> 66 |
| popEVE | 3,552 (45.6%) -> 5,734 (73.6%) | 61 -> 97 | 37 -> 61 |
| CADD | 7,784 -> 7,787 (99.9%) | 109 -> 109 | 66 -> 66 |

- AlphaMissense source per row: snv 7,651, isoform_snv 131 (LPA, transcript ENST00000316300.9), protein 9.
  0 rows whose SNV gives AlphaMissense another protein change than the row's.
- AlphaMissense new vs old (6,472 rows both): 5,714 identical, 665 differ by > 1e-3 (11 proteins). On every one of
  those proteins the new value equals the VESM paper's value; dbNSFP had taken another transcript's value at the
  same position number.
- AlphaMissense new vs VESM paper (pred = -pathogenicity, confirmed): 7,468 rows both, 6,977 identical. All 470
  differences are the UGT1A family (O60656, P19224, P22309, P22310, Q9HAW7, Q9HAW9; 411 rows: in the shared exons
  the paper used UGT1A10's (Q9HAW8) value for the SNV) and LPA (59 rows: the paper's values match neither the
  canonical (old 4,548-aa sequence) nor the isoform table; old dbNSFP = new isoform_snv values).
- Benchmark (mean -log10 Pearson P per archive): on the 80 archives scored by paper, old and new: paper 2.901,
  old 2.871, new 2.899. On the 105 archives with paper scores: paper 2.648, new 2.644. New on all 109: 2.695.
- popEVE new vs old (3,552 rows both): 3,445 identical; 107 differ, all P55017 (dbNSFP used NP_000330.3, a 1,030-aa
  isoform, at the same position number; new uses NP_001119580.2, identical to UniProt) and Q9NRM0 (dbNSFP values
  match no protein in popEVE v1.1; likely an older release). popEVE mean -log10 P on the 61 archives both score:
  old 2.356, new 2.474.
- popEVE still missing: P04114 APOB, P98164 LRP2, P98161 PKD1, P02549 SPTA1, Q9NY15 STAB1 (no popEVE protein,
  all > 2,302 aa); A0AV96 (6 rows), P07306 (3), Q8WTV0 (1): positions outside the aligned popEVE isoform.
  A6NES4 and Q8NB49 found by sequence (NP_001354436.1, NP_001340741.2; Swiss-Prot's NPs are not in popEVE).
- CADD still missing 4 rows (P04180 P430S, P55899 R234Q x2, P98161 E1614K): no CADD in dbNSFP for those SNVs.

## Final layered archives
- Job 31922391 (`--dependency=afterok:31919891`): `layers.sbatch $WORK/genebass/splits $WORK/genebass/splits_layers
  --force`. Rebuild 31919891 writes `splits_layers` with the old code first; this job overwrites it.
- Check after it ends: log `$WORK/logs/layers_final_31922391.log` ends with `LAYERS OK`; `layers_summary.tsv` there.

## Not done / notes
- `add_layers.zip_bytes` writes members with `ct or ZIP_DEFLATED`, so stored (ct = 0) members are deflated: layered
  archives are ~4x smaller than their inputs (6.6 vs 24 GB). Unchanged here (pre-existing behaviour).
- Nothing pushed or merged.
