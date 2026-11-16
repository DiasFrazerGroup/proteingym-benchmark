# Annotation layers for dataset archives

[add_layers.py](add_layers.py) adds two layers to every assay row of any ProteinGym archive (`.splits.pgdata` or
`.pgarchive`; clinical or DMS). Both are extra `non_targets` columns. They are not model inputs and not targets: the
first is for comparison, the second for masking an analysis.

## Columns

| Column | Layer | Meaning | Source |
|---|---|---|---|
| `score_alphamissense` | scores | AlphaMissense, 0–1, higher = more pathogenic | AlphaMissense hg38 table (Zenodo 10.5281/zenodo.8208688) |
| `score_alphamissense_source` | scores | `snv`, `protein` or `none` (see below) | |
| `score_popeve` | scores | popEVE, more negative = more deleterious | popEVE v1.1 zip (UKBB 2025-03-12), per RefSeq protein |
| `score_cadd_phred` | scores | CADD PHRED, higher = more deleterious | dbNSFP 5.3.1a |
| `score_primateai3d` | scores | PrimateAI-3D, higher = more pathogenic | Illumina file, **only if `sources.toml` points to it** |
| `score_snv_source` | scores | SNVs used for CADD / PrimateAI-3D: `genomic`, `protein` or `none` (see below) | |
| `mask_plddt` | masks | pLDDT of the mutated residue (most central structure covering it; min over sites) | archive PDBs if they carry pLDDT, else AlphaFold DB v6 |
| `mask_plddt_source` | masks | `archive`, `afdb_v6` or `none` | |
| `mask_protein_length` | masks | length of the full canonical UniProt protein, not of the archive sequence | Swiss-Prot human 2026_03 |
| `uniprot_pos` | both | mutated position(s) on the canonical UniProt sequence | |

Each column's description in `manifest.lock` records its source and version.

Masking an analysis then reads, for example:

```python
df[(df.mask_plddt.astype(float) > 70) & (df.mask_protein_length.astype(int) >= 100)]
```

## How each row is resolved

- **Variant:** taken from the `mutant` column (`A12V`, or `A12V:K40E` for multi-mutants), numbered on the archive's
  wild type. If there is no `mutant` column, it is derived from the mutated sequence.
- **UniProt numbering:** the wild type is placed on the canonical UniProt sequence by identity, substring or local
  alignment (≥90% identity). DMS archives on a domain therefore still get UniProt positions and the full length.
- **pLDDT source:** proteingym-base rewrites structures when it builds an archive and sets every B-factor to
  0.00, so archive structures normally carry no pLDDT. PG2 archives are affected too. The pLDDT then comes from
  the AlphaFold DB v6 models of the UniProt entry (`[plddt] afdb_dir`). For PG2 archives, whose structures are
  their own AF2 predictions, this is a close proxy, not the same model.
- **pLDDT:** each structure is matched to the wild type, using the residue span in the metadata if present, otherwise
  by alignment. For each site, the value comes from the most central structure that covers it, which handles AFDB
  fragments. For multi-mutants it is the minimum over the mutated sites.
- **External scores:** each from its own complete release (`sources.toml`). Only single substitutions are
  scored; multi-mutants get empty values.
  - **AlphaMissense** (`[alphamissense]`): the official hg38 table (every missense SNV of the canonical transcript
    of every UniProt canonical isoform), turned once into a parquet sorted by UniProt (`build_snv_index.py`). Looked
    up by the archive's UniProt accession plus the row's GRCh38 SNV (`score_alphamissense_source = snv`), else by
    UniProt accession plus protein change (`protein`; mean over its SNVs, which agree because AlphaMissense is
    protein-level). AlphaMissense used UniProt 2021_02: its sequence is rebuilt from the wild-type letters of
    `protein_variant` and, if it differs from the current UniProt sequence, positions are mapped by alignment; the
    column description says which. `layers_am_snv_mismatch.tsv` lists rows whose SNV gives AlphaMissense a different
    protein change than the row's (UniProt numbering); the SNV value is kept.
  - **popEVE** (`[popeve]`): the v1.1 release zip, one CSV per RefSeq protein, read in place. UniProt → RefSeq comes
    from the `DR RefSeq` lines of the local Swiss-Prot file (no id-mapping service). Among the listed NP accessions
    found in the zip (by accession, any version), the one whose sequence (rebuilt from the wild-type letters of
    `mutant`) equals the UniProt canonical sequence is used, then the displayed isoform, then the one covering most
    positions; if it differs, positions are mapped by alignment. The NP used and how it matched are in the column
    description and in `layers_summary.tsv` (`popeve_refseq`, `popeve_sequence`).
  - **CADD** (dbNSFP, SNV-level): `score_snv_source` says which SNVs were used. `genomic`: the row's own GRCh38 SNV
    (`--snv-col`, default `variant_id`, e.g. `2-21041028-G-A`). `protein`: all SNVs causing the protein change, from
    the AlphaMissense table (canonical transcripts), averaged. `none`: no SNV found.
- **Why not dbNSFP for AlphaMissense and popEVE:** their dbNSFP fields are per transcript and had to be resolved by
  UniProt accession and position; on the 109 Genebass archives that left 17% of rows without AlphaMissense (17
  proteins with none) and 54% without popEVE. The values that were found agree with the sources above.

## Guarantees

- Rows and their order are unchanged, so `slices.json` stays valid.
- Each output archive is reloaded with `proteingym.base.Subsets` and its row counts checked.
- Inputs are never overwritten: `--out` must be a different folder.
- Existing layer columns are only replaced with `--force`.
- `layers_summary.tsv` in `--out` gives per-assay coverage of every column and the AlphaMissense / popEVE sequence
  matching.
- `compare_layers.py OLD NEW [--am-ref DIR]` compares two layered runs (coverage, agreement, AlphaMissense against
  reference predictions).

## Run (cluster)

```shell
sbatch -p genoa64 --qos=short -t 04:00:00 -c 4 --mem=48G --chdir=<repo> \
  scripts/clinical/layers/layers.sbatch IN_DIR OUT_DIR [--layers scores,masks] [--snv-col variant_id]
```

The first run builds the AlphaMissense parquet once (`build_snv_index.py`, ~2 min, path in `sources.toml`).

## Licences

AlphaMissense: CC BY-NC-SA 4.0 (non-commercial). popEVE: CC BY 4.0. CADD and dbNSFP: non-commercial use.
PrimateAI-3D: Illumina licence, not redistributable. Score columns therefore need a licence decision before any public
release of layered archives.

## Not yet available

**PrimateAI-3D** is distributed by Illumina after registration. dbNSFP's `PrimateAI_score` is the 2018 PrimateAI,
so it is not used under this name. Once the file is downloaded, set `[primateai3d] path` (and its column names) in
`sources.toml` and rerun with `--force`.
