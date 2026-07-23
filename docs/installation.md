# Installation

GRLP requires Python ≥ 3.9 and depends on `numpy`, `scipy`, `matplotlib`, and
`networkx`.

## From PyPI (release)

```sh
pip install grlp
```

If your system Python shields the core install from external packages, either
add `--break-system-packages` (usually fine, but packages could clash) or, more
robustly, install into a dedicated virtual environment or conda environment.

## From source (development install)

To run the latest unreleased code and have local edits take effect immediately:

```sh
# Clone the repository
gh repo clone MNiMORPH/GRLP
cd GRLP

# Editable install
pip install -e .
```

As above, you may need `--break-system-packages` or a dedicated environment.

You can also simply download the `grlp/` source and import it as a local module,
though managing the path of `grlp.py` relative to your driver script is less
convenient than an editable install.

## Verifying the install

```python
import grlp
print(grlp.__version__)
```

## Building these docs locally

The documentation is built with Sphinx (MyST Markdown). From a checkout with
`grlp` installed in the same environment:

```sh
pip install -r docs/requirements.txt
sphinx-build -b html docs docs/_build/html
```

Then open `docs/_build/html/index.html`.
