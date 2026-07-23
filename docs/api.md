# API reference

Generated from the docstrings in the `grlp` package.

## `grlp.grlp`

The core module: the `LongProfile` and `Network` classes, the governing
equations, the semi-implicit solve, and the analysis methods (diffusivity,
gain/lag, Hack parameters, Strahler orders, Horton ratios).

```{eval-rst}
.. autoclass:: grlp.grlp.LongProfile
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: grlp.grlp.Network
   :members:
   :undoc-members:
   :show-inheritance:
```

## `grlp.solver`

The network solver: assembles and evolves one global sparse linear system over
all nodes of a network by walking the channel topology.

```{eval-rst}
.. automodule:: grlp.solver
   :members:
   :undoc-members:
```

## `grlp.build_synthetic_network`

Random-network generation (Shreve random binary trees) and the helpers that
populate a network with discharges, widths, and profiles.

```{eval-rst}
.. automodule:: grlp.build_synthetic_network
   :members:
   :undoc-members:
```
