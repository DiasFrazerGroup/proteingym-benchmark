# rsalor_filtered, structuredca_filtered (alignment + one structure)

Status 2026-11-20: final, 63/63 proteins of the capped set. Predictions: `predictions/v2__rsalor_filtered/`,
`v2__structuredca_filtered/`.

## Changes

- Shared MSA filter ([msa_filter_shared.md](msa_filter_shared.md)).
- RSALOR: NaN for positions without RSA (the zoo scores them as the wild type).
- Zoo `model.py` copied verbatim; cards unchanged: RSALOR seqid weights 0.80, min_seqid 0.35; StructureDCA 0.9.

## Validation

| | T1 (APOA1, filter off vs zoo) | Filtered vs zoo, Spearman |
|---|---|---|
| RSALOR | max diff 0 | 0.92–0.999 (APOA1, PCSK9, LDLR, ABCA1) |
| StructureDCA | 1.2e-4 (solver numerics) | 0.79–0.99 |

Where the filter is strong (PCSK9: 100k → 1.4k rows), Spearman against the unfiltered zoo is lower, as expected.

The filter also makes StructureDCA feasible: the zoo timed out after more than 2 h on ABCA1; with the filter it takes
13 min.

## Genebass runs (CPU)

| | Total | Slowest protein | Unscored |
|---|---|---|---|
| RSALOR | 39 min | 2 min | 0 |
| StructureDCA | 303 min | 22 min | 0 |

Peak memory is at most 9 GB, mostly the MSA load. Capped set (own rows): StructureDCA 2.471, RSALOR 2.151.

The earlier `v2__rsalor` (T5 region run, W = 1,000) is superseded by `v2__rsalor_filtered`.

## Above 2,700 aa

- Both need **one full-length single-chain structure**. AFDB gives only fragments above 2,700 aa, so APOB, LRP2 and
  PKD1 fail (recorded as model limitations in `status.tsv`).
- Options, not done:
  - most-central AFDB fragment per position, which would be our addition, against the rule;
  - NaN (current).
- Compute is not the issue: T5 region runs took minutes at under 1 GB.
