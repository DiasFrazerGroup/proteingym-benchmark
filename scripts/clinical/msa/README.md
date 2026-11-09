# MSAs for the clinical datasets (PG2 protocol)

PG2 MSA protocol: **ColabFold MSA built on UniRef100 (January 2026 release) with no post-alignment filtering.**
(The PG2 `complete` archives also carry a jackhmmer MSA. We leave it out.)

## Recipe

[colabfold_search_uniref100.py](colabfold_search_uniref100.py) runs the monomer search of `colabfold_search`
(colabfold 1.5.5, `mmseqs_search_monomer`) with `filter=0`, against MMseqs2 DB of UniRef100 2026_01:

1. `search`: 3 profile iterations, `-e 0.1`, `-s 8`, `-a`, `--max-seqs 1000000`.
2. `align` with the final query profile: `-e 10`, `--max-accept 1000000`, `--alt-ali 10`.
3. `filterresult`: `--qid 0 --qsc -20 --diff 0 --max-seq-id 1.0` (the `filter=0` values).
4. `result2msa`: `--filter-msa 0`, a3m (`--msa-format-mode 6`).

The query comes first in each a3m, with its header renamed to the UniProt ID.

Departures from `colabfold_search`, all caused by searching UniRef100 directly:

- **No `expandaln`.** ColabFold expands UniRef30 cluster hits to their members. UniRef100 has no clusters to
  expand, so every member is searched directly.
- **No environmental database.** The protocol names UniRef100 only.
- **`--max-seqs 1000000` instead of ColabFold's 10,000.** ColabFold gets its depth from the cluster expansion,
  and PG2 ColabFold MSAs reach 1,000,000 sequences (median about 168k).

Tools:

- **MMseqs2 18:** the CRG EasyBuild build `18-git20261027-gompi-2025b-int64`, copied with its libraries to
  `$WORK/mmseqs18`. `/software` differs between genoa64 nodes (genoa64-06 has only v17), so the module path is
  not portable.
- **UniRef100:** `/db/uniprot/2026_01/uniref/uniref100/blast/db/uniref100.fasta`.

## Run (SLURM, independent of any interactive session)

```shell
bash scripts/clinical/msa/pipeline.sh      # from the repo root; job ids appended to $WORK/logs/jobs.txt
```

| Job | What | Resources (genoa64) |
|---|---|---|
| `msa_smoke` | Whole recipe on 2 proteins vs Swiss-Prot human; checks the a3m (query = wild type) | 8 CPU, 32 GB, ≤45 min |
| `u100_createdb` | UniRef100 FASTA (248 GB) → `$DB/uniref100` (no index) | 4 CPU, 32 GB, ≤10 h |
| `msa_search` | 61 Genebass wild types → `$WORK/genebass_msas/<UNIPROT>.a3m` (runs if smoke and createdb succeed) | 64 CPU, 400 GB, ≤2 days (predicted max cost €50) |
| `gb_rebuild` | `build_genebass.py --msa-dir` → `$WORK/genebass` (archives with MSAs), Fig. 6b check, index CSV copied into the repo | 8 CPU, 64 GB, ≤6 h |

Paths:

- `WORK=/no_backup/diasfrazer/egalpern_pgb_clinical`: logs in `logs/`, MSAs in `genebass_msas/`, archives in
  `genebass/splits/`.
- `DB=/no_backup/diasfrazer/egalpern_uniref100_2026_01_mmseqs`.

Each log ends with `SMOKE OK`, `CREATEDB OK`, `SEARCH OK` or `REBUILD OK` when that step succeeds.

Check the status with:

```shell
squeue -u $USER; cat /no_backup/diasfrazer/egalpern_pgb_clinical/logs/jobs.txt
tail -n 5 /no_backup/diasfrazer/egalpern_pgb_clinical/logs/*_<jobid>.log
```
