[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3740658.svg)](https://doi.org/10.5281/zenodo.3740658)
[![tests](https://github.com/MNiMORPH/GRLP/actions/workflows/tests.yml/badge.svg)](https://github.com/MNiMORPH/GRLP/actions/workflows/tests.yml)
[![Documentation Status](https://readthedocs.org/projects/grlp/badge/?version=latest)](https://grlp.readthedocs.io/en/latest/?badge=latest)

# gravel-river-long-profile

This model simulates the long-profile evolution of alluvial river valleys formed and evolved by gravel-bedded rivers.

Key assumptions are that:
* The river's bed and banks are formed of noncohesive gravel, and the rate and form of the river's morphological evolution is limited by its ability to move this gravel.
* The river always maintains a self-formed equilibrium width with mobile gravel banks (following [Parker, 1978](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/selfformed-straight-rivers-with-equilibrium-banks-and-mobile-bed-part-2-the-gravel-river/3AD9322C1939528ED73D409654E35E22)). This linearlizes the sediment-transport response to an increase in discharge. _Engineered or bedrock-walled rivers will require a fixed width, or one set by a different set of dynamics._
* Simulations may be based upon a single channel-forming discharge. (See [Blom et al., 2017](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2017JF004213).)
*  We are concerned with aggradation and incision over time scales affecting the full valley bottom (active fluvially worked surface), such that this is the area across which material must be added or removed for vertical change to take place.

The code-base structure, in short, is as follows:
* The **grlp** folder holds the core program
  * `grlp.py` contains the equations and the `LongProfile` and `Network` classes.
  * `solver.py` assembles and solves the network sparse-matrix system.
  * `build_synthetic_network.py` generates networks to run and test GRLP.
* The **examples** folder contains general examples (in the subfolders) as well as tutorial code for a one-dimensional model.
  * `run_grlp.py` contains comments intended to help you learn how to write your own models.
  * `example_1d.ipynb` is a Jupyter notebook containing more extensive tutorial information.
* The **preprocessing** folder holds an optional input-cleaning step that turns network geometry (e.g. from a DEM) into GRLP inputs; run it before GRLP, not as part of it.

## Documentation

Full documentation is hosted on Read the Docs: **<https://grlp.readthedocs.io>**. It covers installation, a quickstart, the theory and governing equations, the example scripts, and a complete API reference.

## Python API

GRLP's main entry points:
* **`grlp.LongProfile`** — a single gravel-bed river long profile. Set its geometry and forcing (`set_x`, `set_z`, `set_Q`, `set_B`; boundaries `set_Qs_input_upstream`/`set_S0` upstream and `set_z_bl`/`set_x_bl` downstream), then evolve with `evolve_threshold_width_river`.
* **`grlp.Network`** — a drainage network of long-profile segments joined at confluences, evolved together with `evolve_threshold_width_river_network`. Every solution is a network solution; a single profile is the trivial one-edge case.
* **`grlp.generate_random_network`** — build a random Shreve network to run without a DEM.

A minimal single-profile run:

```python
import grlp

lp = grlp.LongProfile()
lp.basic_constants()
lp.bedload_lumped_constants()
lp.set_hydrologic_constants()
lp.set_x(dx=1000, nx=251, x0=0)
lp.set_z(S0=7.2e-4, z1=0)
lp.set_Q(700.)
lp.set_B(10000.)
lp.set_niter(3)
lp.set_z_bl(0.)
lp.set_uplift_rate(0.)
lp.set_Qs_input_upstream(lp.k_Qs * lp.Q[0] * 7.2e-4**(7/6.))
lp.evolve_threshold_width_river(nt=10, dt=1e14)   # evolve toward steady state
```

Full class and function documentation is in the [API reference](https://grlp.readthedocs.io/en/latest/api.html).

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

**Wickert, A. D., M<sup>c</sup>Nab, F., and Barefoot, E. (2026). [GRLP](https://doi.org/10.5281/zenodo.3740658) (Version 2.1.0). doi:10.5281/zenodo.3740658.**

**M<sup>c</sup>Nab, F., T. F. Schildgen, J. M. Turowski, and A. D. Wickert (2025), [Influence of network geometry on long-term morphodynamics of alluvial rivers](https://doi.org/10.5194/esurf-13-1059-2025), *Earth Surf. Dynam.*, *13*, 1059–1092, doi:10.5194/esurf-13-1059-2025.**


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
gh repo clone MNiMORPH/GRLP

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

For more extensive network examples, including random network generation, see [the repository](https://doi.org/10.5281/zenodo.15524964) accompanying the paper [M<sup>c</sup>Nab et al. (2025, EGUsphere)](https://doi.org/10.5194/egusphere-2025-2468).
