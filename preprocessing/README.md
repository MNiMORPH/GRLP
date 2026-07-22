# `preprocessing/` — build GRLP network inputs from a DEM

This directory is a **step run before GRLP**, not part of the `grlp` package.
It turns a river network extracted from a DEM into a clean network that GRLP can
read. The cleaned network is meant to be **written to disk and committed**, so
the exact input to the model is inspectable and reproducible.

It is *not* installed by `pip install grlp` (only the `grlp/` package is). Work
from a repository clone.

## The three-stage pipeline

1. **DEM → network graph** — done in GRASS GIS (external), producing a node-link
   JSON (`grass-*.json`) whose edges carry `x`, `y`, `z`, `A` (drainage area) and
   `s` (upstream distance) sampled along each channel, and whose nodes carry a
   single-element `z`. *(The raw GRASS exports are not committed here.)*
2. **clean the graph** — *this directory*. Despike / smooth / enforce a monotone
   elevation envelope / coarsen the along-edge profiles.
3. **graph → GRLP `Network`** — GRLP reads the cleaned JSON (see the run scripts
   in `templates/network/`, which load the committed cleaned inputs).

## Files

- **`network_preprocessing.py`** — the cleaning library (pure NumPy / SciPy /
  NetworkX). Building blocks plus a `preprocess_network()` convenience:
  - despike + diffusive smooth: `clean_edge_values`
  - monotone envelope (pure NumPy): `monotone_envelope_all_edges`
  - smoothing: `symmetric_running_mean_all_edges`
  - coarsening: `downsample_edge_along_s`
  - node-link JSON I/O: `load_network_json`, `save_network_json`
- **`clean_network.py`** — a thin CLI driver that runs the default pipeline
  (monotone envelope → symmetric running mean → coarsen) end to end.

## Usage

```bash
# Reproduce the Whitewater input shipped in templates/network/
# (supply the raw GRASS export yourself; it is not committed).
python3 clean_network.py \
    grass-whitewater-3DEP-sand.json \
    whitewater-3DEP-sand-smoothed-coarsened.json \
    --outlet-node 0 --outlet-source-edge 20 0 --outlet-slope 4e-4 \
    --smooth-window 101 --keep-fraction 0.1
```

Or call the library directly:

```python
import network_preprocessing as npp
G = npp.load_network_json("grass-network.json")
H = npp.preprocess_network(G, node_attr="z", edge_attr="z",
                           smooth_window=101, keep_fraction=0.1)
npp.save_network_json(H, "network-cleaned.json")
```

Smoke tests live in `tests/test_network_preprocessing.py`.
