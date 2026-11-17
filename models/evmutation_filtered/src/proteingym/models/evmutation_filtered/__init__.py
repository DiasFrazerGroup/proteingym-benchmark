# Keep the a3m bit scores when the shared runner converts the dataset (msa_filter.install_header_hook): the package
# is imported before __main__ runs, i.e. before evedesign.proteingym.dataset_to_evedesign is called.
from .msa_filter import install_header_hook

install_header_hook()
