# API reference

Generated from the docstrings in the `grlp` package.

## Which object do I use?

- **A single river long profile (1-D)** → **`grlp.LongProfile`**. The convenience
  wrapper: configure it and call `evolve_threshold_width_river`.
- **A river network** → **`grlp.Segment`** objects joined in a **`grlp.Network`**.
  Build the `Network` (directly or with `initialize` / `generate_random_network`)
  and evolve *it*.

A `Segment` is a pure network member — it holds the data and configuration but
**does not solve itself**. A `LongProfile` *composes* one `Segment` and the
one-edge `Network` that solves it, so the 1-D case reuses the exact same engine
as a full network. (Every solution is a network solution; 1-D is the trivial
one-edge case.)

## `grlp.grlp`

The core module: the `Segment` and `LongProfile` classes, the `Network` class,
the governing equations, the semi-implicit solve, and the analysis methods
(diffusivity, gain/lag, Hack parameters, Strahler orders, Horton ratios).

```{eval-rst}
.. autoclass:: grlp.grlp.Segment
   :members:
   :undoc-members:
   :show-inheritance:

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
