[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6968526.svg)](https://doi.org/10.5281/zenodo.6968526)

# gravel-river-long-profile
Long-profile evolution of gravel-bed rivers

`grlp.py` is module that contains the equations and solvers. `run_grlp.py` interfaces with these to, well, run a model based off of `grlp.py`. Comments inside `run_grlp.py` should help you learn how to write your own models. Please contact Andy Wickert for questions or assistance.


When you use **grlp**, please cite:

**Wickert, A. D. and T. F. Schildgen (2019), [Long-Profile Evolution of Transport-Limited Gravel-Bed Rivers](https://www.earth-surf-dynam.net/7/17/2019/esurf-7-17-2019.html), *Earth Surf. Dynam.*, *7*, 17–43, doi:10.5194/esurf-7-17-2019.**

Furthermore, if you use any of the features associated with the linearization or the plotting scripts noted in the folder [McNab_et_al_submitted_GRL](https://github.com/MNiMORPH/GRLP/tree/master/examples/McNab_et_al_submitted_GRL), please cite:

**M<sup>c</sup>Nab, F., T. F. Schildgen, J. M. Turowski, and A. D. Wickert (2022, submitted), Diverse responses of alluvial rivers to periodic environmental change.**


## Installation

### Via pip and PyPI

Releases will be sent to [PyPI](https://pypi.org/project/GRLP/).

To download and install the release version within your python system, use:

```
# Python 3
pip install grlp

# If you computer is shielding your core Python install from external
# packages, you have two options:
# First, you may simply ignore these blocks
# (Fine in my experience but could cause packages to clash)
pip install grlp --break-system-packages
# Second, you can build a separate environment for GRLP
```

### Locally with pip and incorporating ongoing code modifications

To install the unreleased code from this repository and/or to make changes to it locally and have this reflected immediately in how GRLP runs:

```
# Download the repository
gh repo clone awickert/GRLP

# Install it
# First, navigate to the root grlp directory. Then use:
pip install -e .

# As noted above, you may need to  --break-system-packages or create an
# environment for this to work.
```

You may always just download the `grlp` source from here and run it as a local (rather than system-wide installed) module.
But this can be inconvenient when needing to manage the directory of `grlp.py` relative to that of the driver `*.py` file that you are building to create your model run.

## Learning how to use GRLP

### 1-D river long profile

For a tutorial run the [Jupyter notebook contained within this package](https://github.com/awickert/GRLP/blob/master/example_1d.ipynb).
After installing Jupyter on your local machine, navigate to this directory in the terminal and type:
```
jupyter notebook
```
to launch it. Alternatively, a number of cloud-based services can help to host Jupyter notebooks.

### Network of 1-D river long profiles

* A [complete example](examples/run_network_5_segments_2020.py) for a river network, with plotting and descriptions
* An [example to set up synthetic river networks](https://github.com/awickert/GRLP/blob/master/examples/network_examples/network_examples.py
