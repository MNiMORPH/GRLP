# module will import as grlp (folder name)
# This sets which objects, functions, etc. are exposed
# e.g., 
# from .grlp import LongProfile
# I will just expose all of them until I start to have too many
# or want to restrict the UI
from .grlp import *
from .build_synthetic_network import *

# In addition, import and make visible the version number from _version.py
from ._version import __version__

# Public API for `from grlp import *` (import grlp still exposes everything as
# grlp.X). Curated so the star-import gives the classes and network-generation
# functions, not the imported third-party modules.
__all__ = [
    "LongProfile", "Network",
    "generate_random_network", "Shreve_Random_Network",
    "generate_x_domain", "generate_discharges", "generate_ssds",
    "generate_variable_widths", "generate_zs",
    "__version__",
]

