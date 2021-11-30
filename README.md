[![DOI](https://zenodo.org/badge/159710833.svg)](https://zenodo.org/badge/latestdoi/159710833)

# gravel-river-long-profile
Long-profile evolution of gravel-bed rivers

`grlp.py` is module that contains the equations and solvers. `run_grlp.py` interfaces with these to, well, run a model based off of `grlp.py`. Comments inside `run_grlp.py` should help you learn how to write your own models. Please contact Andy Wickert for questions or assistance.


When you use **grlp**, please cite:

**Wickert, A. D. and T. F. Schildgen (2019), [Long-Profile Evolution of Transport-Limited Gravel-Bed Rivers](https://www.earth-surf-dynam.net/7/17/2019/esurf-7-17-2019.html), *Earth Surf. Dynam.*, *7*, 17–43, doi:10.5194/esurf-7-17-2019.**

## Learning how to use GRLP

### 1-D river long profile

For a tutorial run the [Jupyter notebook contained within this package](https://github.com/awickert/GRLP/blob/master/example_1d.ipynb). After installing Jupyter on your local machine, navigate to this directory in the terminal and type:
```
jupyter notebook
```
to launch it. Alternatively, a number of cloud-based services can help to host Jupyter notebooks.

### Network of 1-D river long profiles

* A [complete example](examples/run_network_5_segments_2020.py) for a river network, with plotting and descriptions
* An [example to set up synthetic river networks](https://github.com/awickert/GRLP/blob/master/examples/network_examples/network_examples.py
