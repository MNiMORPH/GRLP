[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3740658.svg)](https://doi.org/10.5281/zenodo.3740658)

# gravel-river-long-profile
Long-profile evolution of gravel-bed rivers

`grlp.py` is module that contains the equations and solvers. `run_grlp.py` interfaces with these to, well, run a model based off of `grlp.py`. Comments inside `run_grlp.py` should help you learn how to write your own models. Please contact Andy Wickert for questions or assistance.


## Sources to cite

### Base / always

When you use any version of **grlp**, please cite the following two sources:

**Wickert, A. D. and T. F. Schildgen (2019), [Long-Profile Evolution of Transport-Limited Gravel-Bed Rivers](https://www.earth-surf-dynam.net/7/17/2019/esurf-7-17-2019.html), *Earth Surf. Dynam.*, *7*, 17–43, doi:10.5194/esurf-7-17-2019.**

The version of the code that you used, which is found in **CITATION.cff**, automatically formatted by GitHub, and exposed in the web interface.

### Linearization and spectral response

Fergus M<sup>c</sup>Nab and Jens Turowski developed the mathematics and associated code to linearize the GRLP equation in order to take advantage of near-analytical solutions and therefore perform rapid tests with it. If you use any of the features associated with the linearization or the plotting scripts noted in the folder [M<sup>c</sup>Nab_et_al_GRL](https://github.com/MNiMORPH/GRLP/tree/master/examples/McNab_et_al_GRL), please cite:

**M<sup>c</sup>Nab, F., T. F. Schildgen, J. M. Turowski, and A. D. Wickert (2023), [Diverse responses of alluvial rivers to periodic environmental change](https://doi.org/10.1029/2023GL103075), *Geophys. Res. Lett.*, *50*, e2023GL103075, doi:10.1029/2023GL103075.**

### Network

Wickert designed and wrote the code for the networked version of GRLP, with assistance from M<sup>c</sup>Nab, including error finding and bug fixing. M<sup>c</sup>Nab in turn added functions to expanded network functionality and wrote and executed scientific code to understand river-network dynamics. Based on this, M<sup>c</sup>Nab led the study on river-network dynamics and wrote the associated paper. Therefore, if you use the river-network version of GRLP, please cite both of the following:

**Wickert, A. D., M<sup>c</sup>Nab, F., and Barefoot, E. (2025, August). [GRLP](https://doi.org/10.5281/zenodo.3740658) (Version 2.0.0). doi:10.5281/zenodo.3740658.**

**M<sup>c</sup>Nab, F., T. F. Schildgen, J. M. Turowski, and A. D. Wickert (2025), [Influence of network geometry on long-term morphodynamics of alluvial rivers](https://doi.org/10.5194/egusphere-2025-2468), *EGUsphere preprint*, doi:10.5194/egusphere-2025-2468.**


## Installation

### Via pip and PyPI

Releases will be sent to [PyPI](https://pypi.org/project/GRLP/).

To download and install the release version within your python system, use:

```sh
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

```sh
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

For a tutorial run the [Jupyter notebook contained within this package](https://github.com/awickert/GRLP/blob/master/examples/example_1d.ipynb).
After installing Jupyter on your local machine, navigate to the "examples" directory in the terminal and type:
```sh
jupyter notebook
```
to launch it. Alternatively, a number of cloud-based services can help to host Jupyter notebooks.

Alongside the Jupyter notebook is a file,
[run_grlp.py](https://github.com/MNiMORPH/GRLP/blob/master/examples/run_grlp.py), which replicates one of the figures from the Wickert and Schildgen (2019) article. It includes comments to describe how to set up a GRLP run, though the information is less extensive than that available in the Jupyter notebook.

Beyond these two, a set of examples are located in [the "one-dimensional" sub-directory of the "examples" folder](https://github.com/MNiMORPH/GRLP/tree/master/examples/one_dimensional).


### Network of 1-D river long profiles

A set of functional examples for river networks is avaialble within [the "network" sub-directory of the "examples" folder](https://github.com/MNiMORPH/GRLP/tree/master/examples/network).
