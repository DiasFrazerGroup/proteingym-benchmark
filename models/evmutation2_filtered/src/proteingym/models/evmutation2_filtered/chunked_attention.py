"""Exact row chunking for EVmutation2's L^3 encoder steps (clinical: proteins up to 2,700 aa).

evmutation2 0.0.2 (vendored alphafold3_pytorch) materialises encoder tensors that grow as L^3, L^2 x 1,024 or
MSA depth x L x 512 and has no chunk option:
- TriangleAttention: the (h, L, L) bias is repeated for each of the L rows, then attention runs over (L*h, L, L)
  (one tensor ~41 GB at L = 1,367; the out-of-memory point of job 32024606);
- OuterProductMean: (L, L, 32, 32) outer products (~27 GB at L = 2,570);
- the pairwise Transition (SwiGLU, 128 -> 1,024 channels per pair, ~27 GB at L = 2,570) and TriangleMultiplication
  (128 -> 512 channel projections, ~13 GB); the MSA transition (64 -> 512 channels per MSA cell, ~11 GB);
- MSAPairWeightedAveraging: values, gates and output of size MSA depth x L x 256 (~5 GB each at L = 2,570).
All are independent per row i of the pair representation (Transition and the projections per element; the triangle
einsum per output row; the MSA steps per MSA sequence). The patched forwards below run the same operations on
blocks of rows written into one preallocated output (TriangleMultiplication's einsum as the equivalent matmul on
projections stored once in matmul layout): the same arithmetic per element, so scores agree with the unpatched code
up to float32 kernel noise (checked: docs/clinical/models_long_proteins.md). Inference only.

install(max_bytes): max_bytes = budget for one block's largest tensor; None restores the original methods (zoo
behaviour). Only proteins whose tensors exceed the budget are chunked.
"""
import os

import torch
from einops import einsum, rearrange, repeat
from evmutation2.alphafold3_pytorch import alphafold3 as af3

_ORIG = {"tri": af3.TriangleAttention.forward, "opm": af3.OuterProductMean.forward,
         "ff": af3.Transition.forward, "mul": af3.TriangleMultiplication.forward,
         "pwa": af3.MSAPairWeightedAveraging.forward}
_BUDGET = {"bytes": None}


def _rows(total_rows: int, bytes_per_row: int) -> int:
    return max(1, min(total_rows, int(_BUDGET["bytes"] // max(1, bytes_per_row))))


def _blocks(n: int, step: int, fn, dim: int = 1):
    """fn(s, e) for row blocks [s, e) along `dim`, written into one preallocated tensor (no list + cat copy)."""
    out = None
    for s in range(0, n, step):
        o = fn(s, min(n, s + step))
        if out is None:
            shape = list(o.shape)
            shape[dim] = n
            out = o.new_empty(shape)
        out.narrow(dim, s, o.shape[dim]).copy_(o)
        del o
    return out


def _tri_attn_forward(self, pairwise_repr, mask=None, **kwargs):
    b, n = pairwise_repr.shape[0], pairwise_repr.shape[1]
    heads = self.to_attn_bias[0].weight.shape[0]
    # per row: repeated bias + sim + softmax (float32), each heads x n x n
    per_row = 3 * heads * n * n * 4 * b
    if _BUDGET["bytes"] is None or per_row * n <= _BUDGET["bytes"]:
        return _ORIG["tri"](self, pairwise_repr, mask=mask, **kwargs)
    if self.need_transpose:
        pairwise_repr = rearrange(pairwise_repr, "b i j d -> b j i d")
    attn_bias = self.to_attn_bias(pairwise_repr)                     # b h n n, shared by every row
    def block(s, e):
        c = e - s
        ab = repeat(attn_bias, "b ... -> (b repeat) ...", repeat=c)   # same '(b repeat)' order as the original
        m = repeat(mask, "b ... -> (b repeat) ...", repeat=c) if mask is not None else None
        packed, unpack_one = af3.pack_one(pairwise_repr[:, s:e], "* n d")
        return unpack_one(self.attn(packed, mask=m, attn_bias=ab, **kwargs))

    out = _blocks(n, _rows(n, per_row), block)
    if self.need_transpose:
        out = rearrange(out, "b j i d -> b i j d")
    return self.dropout(out)


def _opm_forward(self, msa, *, mask=None, msa_mask=None):
    b, n = msa.shape[0], msa.shape[2]
    dh = self.to_hidden.weight.shape[0] // 2
    per_row = 3 * n * dh * dh * 4 * b                                 # outer product, mean, masked copy
    if _BUDGET["bytes"] is None or per_row * n <= _BUDGET["bytes"]:
        return _ORIG["opm"](self, msa, mask=mask, msa_mask=msa_mask)
    dtype = msa.dtype
    msa = self.norm(msa)
    a, bb = self.to_hidden(msa).chunk(2, dim=-1)
    if msa_mask is not None:
        mm = msa_mask.type(dtype)
        a = a * mm[..., None, None]
        bb = bb * mm[..., None, None]
        num_msa = msa_mask.type(dtype).sum(dim=-1).clamp(min=self.eps)   # (b,)
    else:
        num_msa = msa.shape[1]
    pair_mask = af3.to_pairwise_mask(mask) if mask is not None else None
    def block(s, e):
        op = einsum(a[:, :, s:e], bb, "b s i d, b s j e -> b i j d e")
        if msa_mask is not None:
            op = op / num_msa[:, None, None, None, None]
        else:
            op = op / num_msa
        op = rearrange(op, "... d e -> ... (d e)")
        if pair_mask is not None:
            op = op * pair_mask[:, s:e].type(dtype)[..., None]
        return self.to_pairwise_repr(op)

    return _blocks(n, _rows(n, per_row), block)


def _transition_forward(self, x):
    # per element (... d); chunk the leading rows of a pair representation (b n n d) only
    width = self.ff[0].weight.shape[0]                                # 2 x inner
    if _BUDGET["bytes"] is None or x.dim() != 4:
        return _ORIG["ff"](self, x)
    b, n = x.shape[0], x.shape[1]
    per_row = 2 * b * x.shape[2] * width * 4                          # projection + SwiGLU output
    if per_row * n <= _BUDGET["bytes"]:
        return _ORIG["ff"](self, x)
    return _blocks(n, _rows(n, per_row), lambda s, e: self.ff(x[:, s:e]))


def _tri_mul_forward(self, x, mask=None):
    b, n = x.shape[0], x.shape[1]
    width = self.left_right_proj[0].weight.shape[0]                   # 4 x hidden
    per_row = 2 * b * n * width * 4
    if _BUDGET["bytes"] is None or per_row * n <= _BUDGET["bytes"]:
        return _ORIG["mul"](self, x, mask=mask)
    step = _rows(n, per_row)
    hid = width // 4
    pm = rearrange(af3.to_pairwise_mask(mask), "... -> ... 1") if mask is not None else None
    outgoing = self.mix_einsum_eq.startswith("... i k d")
    # projections (per element), stored once in matmul layout (b, d, ., .):
    #   outgoing out[i,j] = sum_k left[i,k] right[j,k]:  A = left as (d,i,k), B = right as (d,j,k)  -> A @ B^T
    #   incoming out[i,j] = sum_k left[k,j] right[k,i]:  A = right as (d,i,k), B = left as (d,k,j) -> A @ B
    A = x.new_empty((b, hid, n, n))
    B = x.new_empty((b, hid, n, n))
    for s in range(0, n, step):
        e = min(n, s + step)
        l_, r_ = self.left_right_proj(x[:, s:e]).chunk(2, dim=-1)     # rows s:e of the (n, n) grid
        if pm is not None:
            l_, r_ = l_ * pm[:, s:e], r_ * pm[:, s:e]
        if outgoing:
            A[:, :, s:e] = l_.permute(0, 3, 1, 2)
            B[:, :, s:e] = r_.permute(0, 3, 1, 2)
        else:
            A[:, :, :, s:e] = r_.permute(0, 3, 2, 1)
            B[:, :, s:e] = l_.permute(0, 3, 1, 2)
        del l_, r_
    Bm = B.transpose(-1, -2) if outgoing else B

    def block(s, e):                                                  # output rows i in s:e
        o = torch.matmul(A[:, :, s:e], Bm).permute(0, 2, 3, 1)       # b i j d
        return self.to_out(self.to_out_norm(o)) * self.out_gate(x[:, s:e]).sigmoid()

    return _blocks(n, step, block)


def _pwa_forward(self, *, msa, pairwise_repr, mask=None):
    b, s_total, n = msa.shape[0], msa.shape[1], msa.shape[2]
    width = self.msa_to_values_and_gates[1].weight.shape[0]           # 2 x heads x dim_head
    per_row = 3 * b * n * width * 4                                   # values+gates, output, gated output
    if _BUDGET["bytes"] is None or per_row * s_total <= _BUDGET["bytes"]:
        return _ORIG["pwa"](self, msa=msa, pairwise_repr=pairwise_repr, mask=mask)
    bias = self.pairwise_repr_to_attn(pairwise_repr)                  # b h i j, shared by every MSA sequence
    if mask is not None:
        bias = bias.masked_fill(~rearrange(mask, "b j -> b 1 1 j"), af3.max_neg_value(bias))
    weights = bias.softmax(dim=-1)
    def block(s, e):                                                  # MSA sequences s:e
        values, gates = self.msa_to_values_and_gates(msa[:, s:e])
        return self.to_out(einsum(weights, values, "b h i j, b h s j d -> b h s i d") * gates.sigmoid())

    return _blocks(s_total, _rows(s_total, per_row), block)


def install(max_bytes: int | None):
    """Chunk rows when a tensor would exceed max_bytes; None = original evmutation2 code."""
    _BUDGET["bytes"] = None if max_bytes is None else int(max_bytes)
    if max_bytes is not None:   # allocator setting only (no numerical effect); read at the first CUDA allocation
        os.environ.setdefault("PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True")
    if max_bytes is None:
        af3.TriangleAttention.forward = _ORIG["tri"]
        af3.OuterProductMean.forward = _ORIG["opm"]
        af3.Transition.forward = _ORIG["ff"]
        af3.TriangleMultiplication.forward = _ORIG["mul"]
        af3.MSAPairWeightedAveraging.forward = _ORIG["pwa"]
    else:
        af3.TriangleAttention.forward = _tri_attn_forward
        af3.OuterProductMean.forward = _opm_forward
        af3.Transition.forward = _transition_forward
        af3.TriangleMultiplication.forward = _tri_mul_forward
        af3.MSAPairWeightedAveraging.forward = _pwa_forward
