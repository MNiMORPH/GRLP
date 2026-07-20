# Step 1 design: Network-primary, single-thread as a special case

Design proposal responding to `unify-to-network.md` Step 1 (and folding in
Step 2, NetworkX). **Draft for review — no production code changed yet.**

## Goal

Make the network the primary abstraction and a single river a 1-segment
(1-edge) network, so there is one solver and one code path instead of the
current two. Keep a thin, backward-compatible `LongProfile(...)` façade for
1-D use.

## Current architecture (as-is)

Two parallel solve paths, with `LongProfile` doing double duty:

- **`LongProfile`** (~1300 lines) is *both* the single-thread model *and* a
  network segment. It holds the reach arrays (`x, z, Q, B`) and the padded
  `_ext` ghost arrays, has a standalone solver `evolve_threshold_width_river()`
  (via `build_matrices()`), *and* carries network-segment state
  (`upstream_segment_IDs` / `downstream_segment_IDs`, `network__*` methods).
- **`Network`** (~2000 lines) holds `list_of_LongProfile_objects`, stores
  topology as up/down-ID lists on each segment, assembles a block-diagonal
  matrix from each segment's block, links segments through the `_ext` ghosts
  (`update_z_ext_internal`, etc.), and solves the whole system in
  `evolve_threshold_width_river_network()`.

Consequences: duplicated solve logic (standalone vs `network__*`), `_ext` stored
in two shapes (single: ndarray; network: list-of-arrays), and topology encoded
as hand-managed ID lists.

## Enabling finding (measured, not assumed)

A **1-segment `Network` reproduces the standalone `LongProfile` to machine
precision** (max |Δz| = 1.4e-13 m; rel 5e-16 on a 20-node test). So
"single = 1-edge network" is already numerically exact — the façade can route
1-D runs through the network solver without changing results, and the existing
golden-master + analytical tests will catch any drift.

Also established in Step −1: the padded network solver is **correct** at
asymmetric / unequal-`dx` / multi-level junctions (conservation ~1e-14, ghosts
agree exactly). So the class/topology refactor here does **not** require
touching the solver internals; de-padding stays a separate Step 3.

## Target architecture

- **`Network` is primary** and owns: the topology (a NetworkX `DiGraph`, edges =
  river segments, nodes = junctions/endpoints — matching `v.stream.network` and
  what `perrot_net.py` already builds externally), plus the coupled
  block-diagonal solve. A single river is a 2-node, 1-edge graph.
- **`LongProfile` becomes the segment** — reach arrays + segment-local physics
  (coefficients, transport law, and the inherently single-segment analytical /
  spectral methods). Segments live as edge data on the graph.
- **Thin `LongProfile(...)` façade** preserves the current 1-D API
  (`set_x/z/Q/B`, `evolve_threshold_width_river`, the analytical/spectral
  methods). Under the hood a façade run builds/uses a 1-edge network and
  delegates the time stepping. The tutorial/examples keep working unchanged.
- **Nodes carry junction/boundary state** (outlet base level; confluence
  coupling); edges carry the reach arrays.

## Sequencing: fold NetworkX in now (recommended)

Do the NetworkX topology backend *with* this restructure rather than after, so
we don't rewrite the topology layer twice. This makes **`networkx` a hard
dependency** (add to `pyproject.toml`), which you've already approved.
`perrot_net.py` is the reference for the edge=segment / node=junction mapping.

## Backward compatibility

- `LongProfile(...)` + `evolve_threshold_width_river()` keep working (façade).
- `Network.initialize(...)` list-based API is retained (thin adapter that builds
  the `DiGraph`), so existing network examples keep working during the
  transition; the `DiGraph` becomes the internal source of truth.
- Deprecate, don't delete: where an internal method is superseded, warn and
  forward for one release.

## Phased plan (each phase gated by the test suite staying green)

- **B — NetworkX topology layer.** Introduce the `DiGraph` as the internal
  representation behind the current API; `initialize` builds it; topology
  queries (`find_up/downstream_IDs`, heads/mouths, Strahler, etc.) read the
  graph. No numerical change → golden master must stay bit-for-bit.
- **C — Network primary + `LongProfile` façade.** Route the 1-D path through the
  1-edge network; make `LongProfile` the segment class. Golden master + the
  1-seg≡standalone check guard it.
- **D — de-pad the solver (this is Step 3).** Separate, later; replace `_ext`
  with direct neighbor-reaching. Highest risk; not part of Step 1.

## Test strategy

- Add a **`test_single_equals_network.py`**: parametrized single-segment configs
  asserting `LongProfile.evolve_threshold_width_river` ≡ 1-edge `Network` to
  ~1e-12. This becomes the contract that lets the façade delegate.
- The existing **golden master** (`characterization_reference.npz`) is the
  invariant across B and C: setup code may change, recorded numbers may not.
- Run `pytest -m "not slow"` after every step; regenerate the golden master only
  on a *knowing* change (there should be none in B/C).

## Resolved decisions

1. **`networkx` is a hard dependency**, folded in with this step (add to
   `pyproject.toml`).
2. **The single-thread path delegates to the unified network solver**: the
   standalone `LongProfile` is the front end and, in the 1-segment special case,
   runs through the same block-diagonal solver. Backed by the measured
   machine-precision equivalence.
3. **Keep the names** `Network` / `LongProfile` for now; revisit a base-class
   split at Step 4 (FluvTree), when multiple segment *types* force it.
4. **Work directly on `master`**, committing granularly with the test suite
   green at each step. (No feature branch.)
