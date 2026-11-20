"""EVmutation2 (evedesign EVmutation2) with the EVE-protocol MSA filter (msa_filter.py; docs/clinical/models_long_proteins.md §3.6).

The zoo class is unchanged; build() replaces the MSA by the filtered alignment, score() returns NaN for variants at
positions the model does not score. Filter parameters are popped from the card's hyper_parameters, the rest go to
the zoo class. `attn_chunk_bytes` (clinical, long proteins): run the encoder steps that grow as L^3 or L^2 x
channels on blocks of rows when a tensor would exceed this many bytes (chunked_attention.py; exact, same
arithmetic); null = the unchanged evmutation2 code.
"""
from evedesign.models.evmutation2 import EVmutation2

from . import chunked_attention
from .msa_filter import MSAFilterMixin


class EVmutation2Filtered(MSAFilterMixin, EVmutation2):
    name: str = "evmutation2_filtered"

    def __init__(self, attn_chunk_bytes: int | None = None, **kwargs):
        chunked_attention.install(attn_chunk_bytes)
        super().__init__(**self._pop_filter_params(kwargs))
