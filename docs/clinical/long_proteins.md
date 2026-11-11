# Zoo models on long human proteins: limitations, tests and lessons

Status: 2026-11-11. This complements the long-protein section of [clinical_benchmark.md](../clinical_benchmark.md).

The model API work (windowing, fragments, sparse scoring) is being done by the model-API side. This note collects
what we already know, what has been tested, and what other projects in the group taught us, so that work does not
start from zero.

## 1. The problem in numbers

Clinical genes are long. Two Genebass reference sets:

- **Fig. 6b set (103 gene–phenotype pairs, 61 proteins):**
  - 18 proteins are longer than 1,022 aa (the ESM-2 context).
  - 3 are longer than 2,700 aa, the length above which AlphaFold DB only provides fragments: APOB 4,563,
    LRP2 4,655 and PKD1 4,303.
- **All 332 pairs (185 proteins, 79,752 variants):**
  - The variant-weighted mean length is 1,811 aa.
  - The 52 proteins longer than 1,022 aa hold 57% of the variants.
  - 8 proteins are longer than 2,700 aa: NF1, PKHD1, PKD1, LRP1, APOB, LRP2, HERC1 and DST (7,570 aa). They hold
    20% of the variants.

Variants are sparse but rarely share a position: 79,752 variants sit at 62,453 distinct positions. Masked-marginal
cost scales with the number of positions, so sparsity saves less than one might hope.

## 2. Current zoo behaviour

Per model (from reading the code; see [clinical_benchmark.md](../clinical_benchmark.md#long-proteins) for line
references):

| Model | Long proteins today |
|---|---|
| `esm2` | **Fails** above 1,022 aa (hard cap in `can_model` and per instance). By default it scores the unmasked full sequence, one full forward pass per mutant. |
| `esmc` | Runs, because the cap is removed, but memory grows as L². `score_mutants` tokenises with `truncation=True`, so long inputs risk being silently truncated (inferred from the code, not observed). |
| `saprot` | No cap; trained on ≤1,024 tokens, so quality beyond that is untested. Needs structure-derived Foldseek 3Di tokens. |
| `esm-if1` | No cap. One full autoregressive pass per mutant, and the encoder is re-run each time. |
| `pottsmpnn` | One encoder pass per protein. Works to about 21k aa on a 24 GB GPU (measured, see §5). |
| `eve`, `evmutation(2)`, `gemme`, `rsalor`, `structuredca` | No cap, but they need an MSA matching the full wild type. Cost grows with L (EVmutation couplings ∝ L²·q²); see §4. |
| `blosum62` | Any length. |

Shared constraints, in evedesign and proteingym-base:

- **No sliding-window code exists** anywhere in evedesign, the zoo or proteingym-base.
- **Structures:** evedesign asserts a single chain and `first_index=1`, so residue numbers must equal wild-type
  positions.
- **proteingym-base zeroes the B-factor column when it builds an archive.** AlphaFold pLDDT is therefore lost inside
  every archive, PG2's own included. Any model or analysis that reads pLDDT from archive structures sees 0.00.
  - Found on 2026-11-09.
  - Our annotation layers take pLDDT from the original AFDB files instead (`scripts/clinical/layers`).

## 3. What has been implemented and tested (Genebass, 2026-11-06)

This is a prototype in the parallel folder `project/proteingym_clinical/`, kept aside for reference. Each windowed
model is a subclass of the zoo class, mounted into the unchanged zoo container (`models/clinical_window.py`,
`models/clinical_struct.py`). Scoring is masked marginal throughout, with a 1,022-residue context for every model.

| Model | Long-protein adaptation | Validation |
|---|---|---|
| `esm2_window` | One 1,022 window per mutated position, centred and shifted at the ends; `<cls>`/`<eos>` only at the true termini. The logic is ported from dark_energy (§5). | Short protein: equal to the zoo's `ESM2.score_mutants` (max diff 2.4e-7). Per-variant r against the dark-energy ESM2 masked scores: 0.9993 (L ≤ 1,022) and 0.99997 (L > 1,022). Identical benchmark numbers. |
| `esmc_window` | Same windows. | Short protein: equal to the zoo (4e-7). **APOB, 1,022 window vs full length: r = 0.921 (48 variants)**, so the window changes ESMC scores noticeably. |
| `saprot_clinical` | Residue and 3Di tokens cropped to the same window. | Short protein: exactly equal to the zoo. |
| `esm_if1_clinical` | Tiles of 1,022 residues every 511; the tile where the site is most central; score = ΣLL(mutant) − ΣLL(wild type), batched. | Equal to L·(zoo − zoo_wt), r = 1.000000. |
| `pottsmpnn` | Unchanged. | — |

Fragments: for proteins longer than 2,700 aa, each AFDB fragment was aligned to the benchmark sequence. Residues
differing between versions were threaded, and each position used the fragment where it is most central. 79,723 of
79,752 variants had a structure.

Throughput:

- Running one container start per model, instead of one per protein, gave identical outputs and cut the overhead
  from 37 s to 1.5 s per protein (`batch_train.py`).
- On one full H100, the esm2 and esmc windows (build, validation and scoring) took about 2.5 h. SaProt with PottsMPNN
  took 35 min; ESM-IF1 took 1 h 18 min.

Fig. 6b results, Pearson, mean −log10 P over 101 pairs (the paper's 103 minus the 2 UGT1A5 pairs not mapped at
the time). Recorded in the prototype folder; I didn't rerun them for this note.

| Model | Mean −log10 P |
|---|---|
| esm2_window | 2.15 |
| esmc_window | 1.99 |
| SaProt | 1.51 |
| PottsMPNN | 1.28 |
| ESM-IF1 | 1.13 |
| BLOSUM62 | 0.78 |

The windowed ESM2 reproduces the paper's ESM2-650M result (2.10) at the same level.

**Not tested yet:**

- alignment-based models on long proteins;
- ESMC at its native context;
- SaProt beyond 1,024 tokens without cropping;
- AF2-predicted structures (PG2 protocol) for long proteins.

## 4. MSAs for long proteins (PG2 protocol, built 2026-11-09)

The protocol is a ColabFold MSA on UniRef100 2026_01 with no post-alignment filtering. For the 61 Fig. 6b proteins:

- **Depth:** up to 1.16M sequences.
- **File sizes:** PKD1 3.9 GB, LRP2 3.6 GB, ABCG5 2.6 GB; 38 GB in total.
- **Implications for the alignment models:**
  - EVE: sequence weighting is O(depth²). The EVE card already subsamples to 200k (`max_msa_depth`).
  - EVmutation: memory grows as L²·q², which makes proteins over 4,000 aa impractical without restricting to a
    region.
  - proteingym-base loads each MSA into memory when it builds an archive, and the builder copies it into every
    archive of that protein.
- **Region policy (not decided):** whether variants outside the well-covered region of the MSA are dropped or
  reported as missing.

## 5. Lessons from other projects in the group

- **dark_energy, ESM windowing** (`dark_energy/scripts/esm1v/esm1v_core.py:129-181`):
  - masked mode uses one 1,022 window per position;
  - unmasked mode follows the ntranoslab (VESM) intervals: 1,022 residues overlapping by 511, keeping the window
    where the site is most central.

  This is the basis of `esm2_window`.
- **dark_energy, ESM-IF context-chain bug** (`dark_energy/scripts/esmif/esmif_core.py:214-244`, verified, not fixed):
  context chains are cut from the wrong offsets, so for proteins over 500 aa they can overlap the scored chunk and
  skip residues. Existing dark-energy ESM-IF scores for long proteins may be affected. Fix it before reusing that
  code.
- **blend_hamiltonians, PottsMPNN and AWSEM on long chains**
  (`blend_hamiltonians/prototype/docs/long_proteins.md`):
  - AFDB fragments are 1,400 aa long and start every 200 aa, so each interior residue is covered by up to 7.
  - dark_energy merges fragment scores with an unweighted mean over all fragments covering a position. Across
    fragments the scores disagree by up to 16.6 for PottsMPNN (mean 1.56).
  - A stitched full-length model is realistic locally but only approximate globally: consecutive fragments disagree
    on the arrangement of distant domains (RMSD 6–27 Å).
  - PottsMPNN scales to about 21k aa on 24 GB.
- **VESM masked vs unmasked** (`genebass/vesm`):
  - On Genebass, masked inference hurts VESM 650M (Fig. 6b 2.10 vs 2.41 unmasked) but helps ESM1v and ESM2.
  - On PG2 DMS, masked ≥ unmasked for all three, so the VESM deficit is specific to Genebass.
  - The scoring mode is therefore a model-card choice that matters, and both could be offered.
- **PoET with PG2 MSAs** (`rescaling_evolution/evo_benchmarks`): PG2 MSAs can contain the query twice. With a
  maximum similarity of 1.0, PoET can copy the wild type into its prompt, which uniformly penalises mutations. The
  fix is a sampling ensemble at 0.99/0.95.
- **AlphaFold limits:**
  - AFDB itself only provides single models up to 2,700 aa and fragments longer proteins.
  - Running AF2 with default settings (the PG2 structure protocol) on 4,300–4,700 aa chains is beyond typical GPU
    memory. I haven't tested this here.
  - So the PG2 protocol needs an explicit fragment rule for the long clinical proteins.

## 6. Recommendations for the model API

1. **Masked marginal with per-position windows** for every PLM, one window per unique position. For parity, use
   1,022 everywhere; offer the native context of ESMC/SaProt as a separate card.
2. **Fragments for structures:** one structure per fragment with its residue span in the manifest (our archives
   already store `span_start`/`span_end`). Use the most central fragment per position, not the dark-energy
   unweighted mean.
3. **pLDDT:** do not rely on archive B-factors; fix the proteingym-base structure writer, or pass pLDDT separately.
4. **MSAs:** decide the policy for positions outside the covered region and for depth caps for long proteins
   before running EVE/EVmutation on PKD1, LRP2 and APOB.
5. **Validate any window or fragment change** like §3: equal to the zoo on short proteins, then window vs full
   length where the model allows it.
