"""ESM-2 with windowed masked-marginal scoring, for proteins of any length (see window.py)."""
from evedesign.models.esm2 import ESM2

from .window import WindowedMaskedMarginal


class ESM2Window(WindowedMaskedMarginal, ESM2):
    name: str = "ESM2Window"

    def __init__(self, window: int = 1022, special_tokens: str = "termini", **kwargs):
        super().__init__(**kwargs)
        self._init_window(window, special_tokens)
