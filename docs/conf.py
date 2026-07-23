"""
Sphinx configuration for the GRLP documentation.

Docs are written in Markdown (MyST). API pages use autodoc against the installed
`grlp` package, so the package (and its numpy/scipy/matplotlib/networkx deps)
must be importable at build time — Read the Docs installs it via `pip install .`
(see ../.readthedocs.yaml); for a local build, install grlp into the same
environment first.
"""

import os
import sys

# Make the package importable for autodoc when GRLP is not pip-installed
# (e.g. a bare local `make html` from a source checkout).
sys.path.insert(0, os.path.abspath(".."))

# Render matplotlib figures without a display, in case any imported example or
# docstring touches pyplot during the build.
import matplotlib

matplotlib.use("Agg")

# -- Project information ------------------------------------------------------

project = "GRLP"
copyright = "2018–2026, Andrew D. Wickert and GRLP contributors"
author = "Andrew D. Wickert, Fergus McNab, Eric A. Barefoot"

try:
    from grlp import __version__ as release
except Exception:
    # Fall back to reading the single-source version file without importing the
    # full package (autodoc pages will be empty, but the build still succeeds).
    release = "unknown"
    _vfile = os.path.join(os.path.dirname(__file__), "..", "grlp", "_version.py")
    if os.path.exists(_vfile):
        with open(_vfile) as _f:
            for _line in _f:
                if _line.startswith("__version__"):
                    release = _line.split("=")[1].strip().strip('"').strip("'")
version = release

# -- General configuration ----------------------------------------------------

extensions = [
    "myst_parser",             # Markdown source
    "sphinx.ext.autodoc",      # pull docstrings from grlp
    "sphinx.ext.autosummary",  # summary tables
    "sphinx.ext.napoleon",     # NumPy / Google style docstrings
    "sphinx.ext.mathjax",      # render LaTeX math
    "sphinx.ext.viewcode",     # [source] links
    "sphinx.ext.intersphinx",  # cross-link to numpy/scipy docs
]

autosummary_generate = True
autodoc_member_order = "bysource"
autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "show-inheritance": True,
}

# Enable MyST features we use (dollar/`$$` math, aligned equations, etc.).
myst_enable_extensions = [
    "dollarmath",
    "amsmath",
    "colon_fence",
    "deflist",
]
myst_heading_anchors = 3

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
}

templates_path = ["_templates"]

# Files/dirs Sphinx should NOT try to build as pages:
#   - literature/     : working notes that link to copyrighted (untracked) PDFs
#   - network_solver.md: design note whose method names predate the solver
#                        extraction (assemble_by_walking -> solver.assemble);
#                        integrate after it is refreshed.
exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    "literature/*",
    "network_solver.md",
    "requirements.txt",
]

# -- HTML output --------------------------------------------------------------

html_theme = "furo"
html_title = f"GRLP {version}"
html_static_path = ["_static"]
