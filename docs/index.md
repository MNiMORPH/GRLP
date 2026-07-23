# GRLP: gravel-river long-profile evolution

**GRLP** simulates the long-profile evolution of alluvial river valleys formed
and reworked by gravel-bedded rivers, from single channels to full drainage
networks. It solves a physically based, transport-limited sediment-conservation
equation (Wickert & Schildgen, 2019) with a semi-implicit numerical scheme.

Key assumptions:

- The river's bed and banks are noncohesive gravel, and the rate and form of the
  river's evolution is limited by its ability to move that gravel
  (transport-limited).
- The channel maintains a self-formed, near-threshold **equilibrium width** with
  mobile gravel banks (Parker, 1978), which linearizes the transport response to
  changing discharge. Engineered or bedrock-walled rivers need a fixed width set
  by other dynamics.
- Evolution is driven by a single channel-forming (bank-full) discharge acting
  for an intermittency fraction of the time (Blom et al., 2017).
- Aggradation and incision act across the whole valley bottom — the area over
  which material must be added or removed for the bed to change elevation.

```{toctree}
:maxdepth: 2
:caption: Contents

installation
quickstart
theory
network_solver
examples
api
citing
references
```

```{toctree}
:maxdepth: 1
:caption: Tutorials

tutorials/example_1d
tutorials/example_network
tutorials/example_random_network
```

## Where things live

- **`grlp/`** — the package: `grlp.py` (the `LongProfile` and `Network` classes,
  equations, and analysis), `solver.py` (the network sparse-matrix solver), and
  `build_synthetic_network.py` (random-network generation).
- **`examples/`** — runnable, curated example scripts (see {doc}`examples`).
- **`preprocessing/`** — an optional input-cleaning step that turns network
  geometry (e.g. extracted from a DEM) into GRLP inputs; run before GRLP, not
  part of it.

## Indices

- {ref}`genindex`
- {ref}`modindex`
