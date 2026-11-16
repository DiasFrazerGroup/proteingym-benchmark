# Layers from their own sources: progress (branch clinical/layers-sources)

Goal: AlphaMissense and popEVE layers from their own complete releases instead of dbNSFP (which left 17% / 54% of
Genebass rows empty); CADD stays on dbNSFP. Task brief: see the commits on this branch.

## Done
- `sources.toml`: new `[alphamissense]` (parquet `table` + `source`) and `[popeve]` (zip) sections; `[snv_index]`
  removed (the AlphaMissense parquet is the SNV index); `[dbnsfp]` now only for CADD.
- `build_snv_index.py`: writes `$WORK/layers/alphamissense_hg38.parquet` (adds `am_pathogenicity`, sorted by
  UniProt). Built by job 31920801 (OK: 71,697,556 SNVs, 19,117 proteins). The old `am_hg38_snv_index.parquet` is
  untouched (build_genebass.py and the running job 31919891 still use it).
- `add_layers.py`: classes `AlphaMissense` (by UniProt + SNV, else protein change; sequence rebuilt and aligned if it
  differs), `PopEVE` (zip read in place, UniProt -> RefSeq via Swiss-Prot DR lines, sequence check + alignment),
  `Dbnsfp` (CADD only, no transcript gating). New column `score_alphamissense_source`; RefSeq/sequence notes in column
  descriptions and `layers_summary.tsv`; `layers_am_snv_mismatch.tsv` written.
- `layers.sbatch` reads `[alphamissense] table`.

## Tested
- Smoke test locally on 4 archives (AATC, ALBU x3): 100% coverage of all layers ($WORK/layers_v2_test/smoke).

## Jobs
- 31920801 am_index: done.
- 31920906 layers_v2_test: 109 no-MSA Genebass archives -> $WORK/layers_v2_test/genebass_nomsa (submitted).

## Next steps
1. When 31920906 is done: `python <scratchpad>/validate_layers.py $WORK/genebass_api/build/splits_layers
   $WORK/layers_v2_test/genebass_nomsa` (old vs new coverage, AM new == old, AM == -paper pred, popEVE vs old).
   (Script copy: scripts/clinical/layers/validate_layers.py once committed.)
2. Update scripts/clinical/layers/README.md and docs/clinical/genebass.md (layer section, licences).
3. After job 31919891 ends with REBUILD OK: sbatch layers.sbatch $WORK/genebass/splits $WORK/genebass/splits_layers
   --force (or with --dependency=afterok:31919891).
