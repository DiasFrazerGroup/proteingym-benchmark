# Clinical benchmark cards as PG2 archives

Status: 2026-11-17. Branch `clinical/cards-data` (from `clinical/genebass-data`). Dataset README:
[datasets/clinical/cards/README.md](../../datasets/clinical/cards/README.md).

## What and why

The provisional clinical benchmarks are defined in clinical-curation as **benchmark cards**
(`benchmark_cards/`, `scripts/16_benchmark_cards.py`). This turns each card into its own set of PG2 archives, one
`.splits.pgdata` per gene, so the zoo can score the missense variants on the protein. The archives hold
benchmark data only (wild type, structures, MSA, mutated sequences, target, provenance), **no model scores**,
following the Genebass archives.

Decisions (2026-11-17):
- **Cards:** every provisional card, plus a new random DNM set; the reference cards `dnm_aptest_full` and
  `gnomad_af_full` are not archived.
- **One set per card:** genes shared by several cards are archived in each.
- **New `dnm_aptest_random300`:** 300 genes drawn with seed 0 from all `ap_test_set` genes, mirroring
  `gnomad_af_random300`.
- **MSAs** for all card proteins with the Genebass recipe (ColabFold on UniRef100 2026_01, no filtering).
- **Random sets topped up (2026-11-20):** both random sets have exactly 300 genes, all proteins ≤ 2,700 aa.
  - Kept: the original draw's usable genes.
  - Replacements: drawn from the same pool with seed 1.
  - `dnm_aptest_random300`: 9 replaced (8 unmappable, TACC2 2,948 aa).
  - `gnomad_af_random300`: 2 replaced (MEF2A, RERGL).
  - Lists with origin: `datasets/clinical/cards/<card>_genes.tsv`.

## Results (final build with MSAs, 2026-11-19; random sets 2026-11-20)

| Card | Archives (genes) | Rows | Card rows unmapped | Target |
|---|---|---|---|---|
| `clinvar_v0` | 158 of 159 | 14,267 | 101 (0.70%) | `label` |
| `dnm_aptest_ge1` | 411 of 414 | 2,486 | 28 (1.11%) | `dnm_case` |
| `dnm_aptest_ge2` | 75 of 76 | 613 | 9 (1.45%) | `dnm_case` |
| `dnm_aptest_random300` | 300 (9 replaced) | 722 | 6 (0.82%) | `dnm_case` |
| `gnomad_af_clinvar_dnm_ge1` | 536 of 538 | 489,821 | 6,778 (1.37%) | `log10_af` |
| `gnomad_af_clinvar_dnm_ge2` | 221 of 222 | 261,647 | 2,330 (0.88%) | `log10_af` |
| `gnomad_af_random300` | 300 (2 replaced) | 127,261 | 1,866 (1.45%) | `log10_af` |

**Structures:** 1,092 proteins in all, every one with its AFDB v6 structure. RERGL, the one protein whose AFDB model
did not match Swiss-Prot 2026_03, left with the gnomAD random-set replacement.

**Every archive has its MSA:**
- 7 reused from Genebass;
- 1,077 searched by the chain in 8 jobs of 2–6 h each;
- 10 searched for the random-set replacements. Depths: median 199k sequences, range 222 to 2.97M. Archives (compressed): 104 GB over the seven
sets, from 3.6 GB (`dnm_aptest_ge2`) to 30 GB (`gnomad_af_clinvar_dnm_ge1`).

**Check (`check_cards.py`): CHECK OK.**
- **Selection:** the counts before mapping are identical to `summary.tsv` on all six existing cards (genes, rows,
  classes, complete rows, singletons).
- **Rows:** the archive rows are identical to the mapped rows.
- **Sequences:** every row's sequence is the wild type with exactly its `mutant`.
- **Loading:** sampled archives load with `Subsets.from_path`.
- **Metrics:** the six card models' metrics recomputed on the archived rows reproduce `clinvar_v0` exactly (max
  |diff| 4.8e-5). On the other cards they differ by ≤ 0.008 (DNM ≥2/≥2), ≤ 0.003 (DNM ≥1/≥1) and ≤ 0.001
  (gnomAD), from the unmapped rows.

## Mapping onto UniProt

- **Protein:** one per gene, the Swiss-Prot canonical accession most of the gene's SNVs map to in AlphaMissense,
  else Swiss-Prot by gene symbol.
- **By SNV:** rows are placed by their SNV through the AlphaMissense tables (99.7% of rows).
  - ClinVar and gnomAD use GRCh38.
  - **The popEVE DNMs are GRCh37**, so they go through `AlphaMissense_hg19`. Their GRCh38 `variant_id` is recovered
    from the hg38 index, and the source coordinate is kept as `variant_id_grch37`.
  - A first pass that treated the DNMs as GRCh38 lost 18–62% of their rows.
- **By source numbering:** otherwise, the source protein numbering is placed with a constant offset fitted per
  source isoform (≥90% of ≥10 rows matching).
  - This recovers genes missing from the AlphaMissense canonical set, such as COG4 and TTC8.

**Still unmapped** (listed in `$WORK/cards/<card>/unmapped_rows.tsv`): variants whose source isoform (MANE or
RefSeq) differs in exon structure from the UniProt canonical isoform, so that no single offset fits.
- Whole genes lost: CACNA1C (ClinVar, gnomAD), SHANK2, SLC4A11 and TTC8 (DNM ≥1/≥1), and MEF2A (gnomAD random300).
- Partly unmapped: CIC, RFX7 and others.
- Placing these needs the source protein sequences (MANE/RefSeq FASTA) and an alignment. That's a possible
  follow-up, or such genes could use the MANE protein as their wild type instead of UniProt canonical.

## Runs and files

- **MSA search:** the card proteins in chunks of 150 (`submit.sh msa`), each 64 CPUs / 400 GB on genoa64.
  - Jobs 31960150–31960164, 2026-11-18/19.
  - The per-user memory limit runs the chunks one after another.
- **Rebuild and check:** the per-card build jobs (31960166–31960175) ran after the search, followed by the check
  (31960176: CHECK OK).

Work dir: `/no_backup/diasfrazer/cpugh/pgb_clinical_cards`:
- logs and `logs/jobs.txt`;
- **final archives** in `cards/<card>/splits/`;
- MSAs in `msas/` (502 GB);
- staged per-protein inputs in `cards/build/proteins/` (copies of the MSAs; can be deleted, a rebuild recreates
  them).
