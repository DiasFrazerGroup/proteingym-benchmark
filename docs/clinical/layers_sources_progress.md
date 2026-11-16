# Layers from their own sources: progress (branch clinical/layers-sources)

Goal: AlphaMissense and popEVE layers from their own complete releases instead of dbNSFP (which left 17% / 54% of
Genebass rows empty); CADD stays on dbNSFP.

## Done
- `sources.toml`: `[alphamissense]` (hg38 `table` + `protein_table` from aa_substitutions), `[popeve]` (zip +
  `sequences` index); `[snv_index]` removed (the AlphaMissense hg38 parquet is the SNV index); `[dbnsfp]` only CADD.
- `build_snv_index.py` builds both AlphaMissense parquets in `$WORK/layers/` (jobs 31920801, 31921227: OK).
  The old `am_hg38_snv_index.parquet` is untouched (build_genebass.py still defaults to it).
- `build_popeve_index.py` builds `$WORK/layers/popeve_sequences.parquet` (job 31921332: OK, 18,332 proteins, max
  length 2,302 aa).
- `add_layers.py`: `AlphaMissense` (by UniProt + SNV, else by protein change from aa_substitutions), `PopEVE`
  (Swiss-Prot DR RefSeq -> zip; else by identical / near-identical sequence), `Dbnsfp` (CADD only). New column
  `score_alphamissense_source`. RefSeq/sequence notes in descriptions and `layers_summary.tsv`.
- `compare_layers.py`: old vs new coverage/agreement + AlphaMissense vs reference predictions.
- Docs: layers README, genebass.md (columns, caveats), status.md.

## Tested
- First 109-archive run (job 31920906, before the aa_substitutions / sequence fallbacks): AM 83.1% -> 98.3%,
  popEVE 45.6% -> 70.6%; where values differ from dbNSFP, the new ones agree with the VESM paper (old dbNSFP picked
  other transcripts). Remaining paper differences: UGT1A family shared exons (paper used UGT1A10's value).
- Smoke test of fallbacks (LPA, MROH2A, ATP11C, ASGR1, APOB): OK.

## Jobs
- 31921444 layers_v2_test rerun -> $WORK/layers_v2_test/genebass_nomsa.

## Next steps
1. `compare_layers.py $WORK/genebass_api/build/splits_layers $WORK/layers_v2_test/genebass_nomsa --am-ref
   .../proteingym_clinical/genebass/predictions/paper__AlphaMissense`; record numbers here.
2. After rebuild 31919891 (REBUILD OK): `sbatch ... layers.sbatch $WORK/genebass/splits $WORK/genebass/splits_layers
   --force` (or `--dependency=afterok:31919891`).
