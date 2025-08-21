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

__all__ = ["__version__"]

