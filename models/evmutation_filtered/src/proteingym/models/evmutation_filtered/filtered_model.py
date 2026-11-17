"""EVmutation (evedesign EVcouplingsPLM, plmc) with the EVE-protocol MSA filter (msa_filter.py; docs/clinical/models_long_proteins.md §3.6).

The zoo class is unchanged; build() replaces the MSA by the filtered alignment, score() returns NaN for variants at
positions the model does not score. Filter parameters are popped from the card's hyper_parameters, the rest go to
the zoo class.
"""
from evedesign.models.evcouplings import EVcouplingsPLM

from .msa_filter import MSAFilterMixin


class EVmutationFiltered(MSAFilterMixin, EVcouplingsPLM):
    name: str = "evmutation_filtered"

    def __init__(self, **kwargs):
        super().__init__(**self._pop_filter_params(kwargs))
