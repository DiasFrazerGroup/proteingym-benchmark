# evmutation2_filtered (zoo `evmutation2`: `evedesign.models.evmutation2:EVmutation2`, MSA-only small checkpoint)

Status 2026-11-20: final, 63/63 proteins of the capped set, 4,653 of 4,653 Genebass mutants scored.
Predictions: `predictions/v2__evmutation2_filtered/`.

## Changes

- Shared MSA filter ([msa_filter_shared.md](msa_filter_shared.md)). EVmutation2 scores every position, so no
  variant is NaN (0 of the Genebass rows).
- **Exact chunked encoder** (`chunked_attention.py`, card option `attn_chunk_bytes: 4000000000`; `null` = the
  unchanged code):
  - **Why.** evmutation2 0.0.2 (its vendored `alphafold3_pytorch`) has no chunk option, and several encoder steps
    allocate one tensor that grows as L³ or L² × channels. The unchanged model runs out of memory at **1,367 aa
    (IGF1R) on an 80 GB GPU**: TriangleAttention repeats the (heads, L, L) bias for each of the L rows and asks for
    one 38 GiB block.
  - **What.** Five modules get a patched `forward` that runs the same operations on blocks of rows written into one
    preallocated output:
    - TriangleAttention: per row, the shared bias, similarity and softmax;
    - OuterProductMean: (L, L, 32, 32);
    - Transition: pair SwiGLU 128 → 1,024 and the MSA transition;
    - TriangleMultiplication: the triangle einsum as the equivalent matmul, per output row;
    - MSAPairWeightedAveraging: per MSA sequence.
  - **Exact by construction.** Every step is independent per row, so the arithmetic per element is unchanged. A
    block is used only when one tensor would exceed the budget, so proteins below about 1,000 aa run the original
    code paths.
  - Inference only. Patches are installed at model construction; `install(None)` restores the original methods.
- The rest of the card is the zoo's: 16 encoder samples, 4 recycling steps, `encoder_max_num_msa` 2,048,
  16 decoding orders.
- **Runner:**
  - `validate_msa_filter.py` gets `--seed`, which seeds python/numpy/torch before build and score. EVmutation2 is
    stochastic: it samples MSA subsets and decoding orders.
  - It also records `gpu_peak_gb`.
  - `run_genebass.sbatch` passes `HP` overrides as `--hp`.

## Validation

- **T1** (APOA1, `msa_filter: false` vs the zoo card in the zoo image, seed 0, 500 singles + Genebass variants):
  max abs diff **0** on 535 of 535 variants, no NaN (jobs 32051467/8, RTX 2080 Ti). The wrapper, with the chunking
  installed (below budget at 267 aa), is the zoo model.
- **Chunked vs unchanged code** (PCSK9, 692 aa, seed 0, 400 singles + 64 Genebass variants, budget forced down to
  0.5 GB so every patched step runs in blocks):

  | Code version | Max abs diff vs unchanged | Score time | GPU peak |
  |---|---|---|---|
  | unchanged (`attn_chunk_bytes: null`) | — | 366 s | 13.6 GB |
  | chunked v1–v4 (each revision) | **0** (464 of 464 scores) | 366 s | 5.2 GB (v4) |

  Some Genebass proteins ran with earlier revisions, because a job reads the package at each protein's start. All
  revisions gave identical scores.
- **No released reference.** EVmutation2 is not in ProteinGym v1.3 or the popEVE release, so there is nothing to
  compare against beyond the zoo itself (T1).

## Genebass runs (GPU)

- ≤ 1,800 aa on a 3g.40gb slice; 2,040–2,570 aa on the 7g.80gb GPU (they run out of memory on 40 GB).
- Time is build (16 encoder samples × 4 recycles) plus scoring. It grows about as L².

  | Protein | L | Time (build + score) | GPU peak |
  |---|---|---|---|
  | APOA1 | 267 | 0.5 min | 2.1 GB |
  | IGF1R | 1,367 | 17 + 6 min | 19.3 GB |
  | MROH2A | 1,674 | 30 + 14 min | 28.9 GB |
  | PAPPA2 | 1,791 | 36 + 8 min | 32.9 GB |
  | LPA | 2,040 | 25 + 10 min | 42.7 GB (7g.80gb) |
  | ABCA1 | 2,261 | 37 + 26 min | 52.3 GB |
  | SPTA1 | 2,419 | 45 + 12 min | 59.8 GB |
  | STAB1 | 2,570 | 53 + 24 min | 67.4 GB |

- Total: 628 GPU-min (393 on 40 GB slices for the 59 proteins up to 1,791 aa, 235 on the 80 GB GPU for the 4
  largest). Peak host RAM is 64 GB (the 809k-row ColabFold MSA of STAB1 is loaded before the filter).
- Jobs 32028635, 32028636 (3g.40gb) and 32044587 (7g.80gb). The OOM records from 40 GB are kept in
  `$WORK/models/msa_filter/genebass/evm2_oom_40gb/`.
- Single run per protein, unseeded, as the zoo runs it.
- Capped set: 1.937 on the 5,750 variants scored by every model and 2.092 on all 6,066 variants (100 pairs).
  - That is below the other alignment models (EVmutation 2.208, GEMME 2.35, StructureDCA 2.47).
  - It is above CADD and the inverse-folding models.

## Above 2,700 aa

- **Memory is the limit.** The GPU peak grows about as L² (42.7 GB at 2,040 aa, 67.4 GB at 2,570 aa). Extrapolated
  as L² [estimate]:
  - about 85 GB at 2,900 aa, so the 80 GB GPU ends near 2,800 aa;
  - 190–220 GB for PKD1 (4,303 aa), APOB (4,563 aa) and LRP2 (4,655 aa).
- **Why chunking stops helping.** With the 4 GB budget no single block dominates any more, so what remains is
  probably the pair representation (L² × 128 channels) and the copies each encoder step keeps [not profiled].
- **Below the cap, chunking is enough:** nothing above 80 GB was needed.
- **What it would take:** a smaller chunk budget plus chunking of the pair-representation updates in place (more
  code in the vendored AF3 encoder), or a GPU with more memory. Time would be about 3–4 h per protein [estimate].
- **Alternative:** MSA regions (T5: W = 500 reproduces full length at r 0.93, 22–25 min per region). The user ruled
  this out (no MSA segments).
- **Current:** NaN for the three proteins (recorded as `missing` in `status.tsv`), like every model above the cap.
