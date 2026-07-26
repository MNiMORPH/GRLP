"""
Sphinx configuration for the GRLP documentation.

Docs are written in Markdown (MyST). API pages use autodoc against the installed
`grlp` package, so the package (and its numpy/scipy/matplotlib/networkx deps)
must be importable at build time — Read the Docs installs it via `pip install .`
(see ../.readthedocs.yaml); for a local build, install grlp into the same
environment first.
"""

import glob
import os
import subprocess
import sys

# Make the package importable for autodoc when GRLP is not pip-installed
# (e.g. a bare local `make html` from a source checkout).
sys.path.insert(0, os.path.abspath(".."))

# Render matplotlib figures without a display, in case any imported example or
# docstring touches pyplot during the build.
import matplotlib

matplotlib.use("Agg")

# Copy the tutorial notebooks from examples/ into the docs source tree so myst-nb
# can render and execute them, keeping a single canonical copy in examples/. The
# copies (docs/tutorials/) are git-ignored.
import shutil

_here = os.path.abspath(os.path.dirname(__file__))
_tutorial_dst = os.path.join(_here, "tutorials")
os.makedirs(_tutorial_dst, exist_ok=True)
for _nb in ("example_1d.ipynb", "example_network.ipynb",
            "example_random_network.ipynb"):
    shutil.copyfile(os.path.join(_here, "..", "examples", _nb),
                    os.path.join(_tutorial_dst, _nb))

# Build the interactive browser demo (docs/interactive.md). `panel convert`
# compiles interactive_demo/grlp_panel.py into a standalone WebAssembly app that
# runs GRLP entirely in the reader's browser via Pyodide. The app and the wheels
# it loads (grlp, built fresh here; panel + bokeh, self-hosted from PyPI) live
# together in _static/interactive/ and are referenced by *relative* URLs, so the
# demo works wherever the docs are served with no external CDN. interactive.md
# embeds _static/interactive/grlp_panel.html in an <iframe>.
import re

import bokeh
import panel

_demo_src = os.path.join(_here, "..", "interactive_demo", "grlp_panel.py")
_demo_out = os.path.join(_here, "_static", "interactive")
os.makedirs(_demo_out, exist_ok=True)
for _old in glob.glob(os.path.join(_demo_out, "*.whl")):
    os.remove(_old)
# Fresh grlp wheel (matches the current source) + self-hosted panel/bokeh wheels.
subprocess.run([sys.executable, "-m", "pip", "wheel", os.path.join(_here, ".."),
                "--no-deps", "-w", _demo_out], check=True)
subprocess.run([sys.executable, "-m", "pip", "download",
                f"panel=={panel.__version__}", f"bokeh=={bokeh.__version__}",
                "--no-deps", "-d", _demo_out], check=True)
_grlp_whl = os.path.basename(glob.glob(os.path.join(_demo_out, "grlp-*.whl"))[0])
subprocess.run([sys.executable, "-m", "panel", "convert", _demo_src,
                "--to", "pyodide-worker", "--out", _demo_out,
                "--requirements", _grlp_whl, "numpy", "scipy", "networkx"],
               check=True, cwd=_demo_out)
# Point the panel + bokeh wheel URLs at the co-located copies (the holoviz CDN's
# bokeh wheel 403s; self-hosting also removes the run-time CDN dependency).
for _f in ("grlp_panel.js", "grlp_panel.html"):
    _fp = os.path.join(_demo_out, _f)
    if not os.path.exists(_fp):
        continue
    with open(_fp) as _fh:
        _txt = _fh.read()
    _txt = re.sub(r"https://cdn\.holoviz\.org/\S*?/(bokeh-[\d.]+-py3-none-any\.whl)",
                  r"\1", _txt)
    _txt = re.sub(r"https://cdn\.holoviz\.org/\S*?/(panel-[\d.]+-py3-none-any\.whl)",
                  r"\1", _txt)
    with open(_fp, "w") as _fh:
        _fh.write(_txt)

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
    "myst_nb",                 # Markdown source + executable notebooks
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

# Execute the tutorial notebooks at build time (fresh v3 outputs/plots), cached
# so unchanged notebooks are not re-run on every build.
nb_execution_mode = "cache"
nb_execution_timeout = 300

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
}

templates_path = ["_templates"]

# Files/dirs Sphinx should NOT try to build as pages:
#   - literature/ : working notes that link to copyrighted (untracked) PDFs
exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    "literature/*",
    "requirements.txt",
    # myst-nb writes executed copies of the tutorial notebooks here; they are a
    # build artifact, not source pages.
    "jupyter_execute",
    ".jupyter_cache",
]

# -- HTML output --------------------------------------------------------------

html_theme = "furo"
html_title = f"GRLP {version}"
html_static_path = ["_static"]
