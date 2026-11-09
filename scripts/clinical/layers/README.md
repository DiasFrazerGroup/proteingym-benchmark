# Annotation layers for dataset archives

[add_layers.py](add_layers.py) adds two layers to every assay row of any ProteinGym archive (`.splits.pgdata` or
`.pgarchive`; clinical or DMS). Both are extra `non_targets` columns. They are not model inputs and not targets: the
first is for comparison, the second for masking an analysis.

## Columns

| Column | Layer | Meaning | Source |
|---|---|---|---|
| `score_alphamissense` | scores | AlphaMissense, 0–1, higher = more pathogenic | dbNSFP 5.3.1a |
| `score_popeve` | scores | popEVE, more negative = more deleterious | dbNSFP 5.3.1a |
| `score_cadd_phred` | scores | CADD PHRED, higher = more deleterious | dbNSFP 5.3.1a |
| `score_primateai3d` | scores | PrimateAI-3D, higher = more pathogenic | Illumina file, **only if `sources.toml` points to it** |
| `score_snv_source` | scores | `genomic`, `protein` or `none` (see below) | |
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
- **External scores:** scored per SNV, then mapped back to the protein change. Only single substitutions are scored;
  multi-mutants get empty values.
  - `genomic`: the row's own GRCh38 SNV (`--snv-col`, default `variant_id`, e.g. `2-21041028-G-A`).
  - `protein`: all SNVs that cause the protein change, from the AlphaMissense hg38 table (canonical transcripts).
    Values are averaged over them; AlphaMissense and popEVE are protein-level, so only CADD varies.
  - `none`: no SNV found.
- **dbNSFP transcript choice:** dbNSFP fields are per transcript. The transcript is the one whose UniProt accession
  and protein position match, falling back to accession plus amino acids (for numbering offsets such as UGT1A5).

## Guarantees

- Rows and their order are unchanged, so `slices.json` stays valid.
- Each output archive is reloaded with `proteingym.base.Subsets` and its row counts checked.
- Inputs are never overwritten: `--out` must be a different folder.
- Existing layer columns are only replaced with `--force`.
- `layers_summary.tsv` in `--out` gives per-assay coverage of every column.

## Run (cluster)

```shell
sbatch -p genoa64 --qos=short -t 04:00:00 -c 4 --mem=48G --chdir=<repo> \
  scripts/clinical/layers/layers.sbatch IN_DIR OUT_DIR [--layers scores,masks] [--snv-col variant_id]
```

The first run builds the SNV index once (`build_snv_index.py`, parquet, path in `sources.toml`).

## Not yet available

**PrimateAI-3D** is distributed by Illumina after registration. dbNSFP's `PrimateAI_score` is the 2018 PrimateAI,
so it is not used under this name. Once the file is downloaded, set `[primateai3d] path` (and its column names) in
`sources.toml` and rerun with `--force`.
