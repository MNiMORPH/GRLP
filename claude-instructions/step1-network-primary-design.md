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

## Phased plan (revised — see "Phase C findings" below)

- **B — NetworkX topology layer. DONE.** `Network.build_graph()` builds the
  `DiGraph` in `initialize`; `find_up/downstream_IDs` and the channel
  head/mouth lists (and, transitively, Strahler/Horton) now read the graph.
  Purely additive; golden master bit-for-bit. Guarded by `test_network_graph.py`
  and a per-head-`S0` ordering test.
- **C/D — de-pad the solver, which unifies by construction.** *Reordered:*
  replacing the padded `_ext` arrays with a single neighbor-walking matrix
  assembly is now the **vehicle** for unification, not a later step. One assembly
  ⇒ one BC discretization ⇒ single-segment (1-edge walk) and network agree by
  construction, and the `LongProfile` façade falls out as a 1-edge call. This
  dissolves the Phase C boundary-condition discrepancy rather than reconciling
  two soon-to-be-deleted implementations.

## Phase C findings (why the reorder)

Measured: a 1-segment `Network` equals the standalone `LongProfile` **only at
equilibrium**, not in the transient (for `niter > 1`). Two causes:

1. **Inert Picard loop** in `evolve_threshold_width_river_network`
   (`grlp.py` ~2754): it never wrote the iterate back into `z_ext[1:-1]`, so the
   semi-implicit coefficient never refined and `niter` was a no-op (the standalone
   writes `z_ext[1:-1] = spsolve(...)`). The author's own commented-out line at
   the segment-update block is the intended fix. **Lesson for the unified
   assembly: write the iterate into the interior before recomputing the
   coefficient (proper Picard), or `niter` does nothing.**
2. **Different upstream-BC discretization** (the real blocker). The standalone
   applies the upstream sediment-flux Neumann BC by *modifying the matrix*
   (`set_bcl_Neumann_LHS/RHS`, `grlp.py` ~489-571, called from `build_matrices`
   ~619). The network applies it via the *ghost node* `z_ext[0]` and does
   `pass` for head segments in
   `add_block_diagonal_matrix_upstream_boundary_conditions` (~1393). Both are
   valid; they give the same equilibrium but different transient fixed points
   (~3e-3). Reconciling ≈ deleting one implementation, which is what de-padding
   does anyway.

A partial Picard fix (cause 1) was written, verified equilibrium-preserving, then
**reverted** — it is throwaway against the de-pad rewrite; only its lesson is kept.

### Golden-safety of unifying the head Neumann BC (verified post-compact)

Two measured facts make the unification safe for the golden master:

1. **A de-padded single-segment assembly reproduces the standalone
   `build_matrices()` LHS+RHS to `0.0` absolute difference** (prototype
   `proto_depad_single.py`): computing the boundary gradients directly from the
   Neumann `S0` and the outlet `z_bl` — no stored `_ext` pad — and applying the
   standalone matrix-modification Neumann BC gives the identical tridiagonal
   system. So the single-segment interior + Neumann + Dirichlet stencil de-pads
   exactly.
2. **The network golden-master stores only steady-state values.**
   `build_network`/`run_topology_arrays` evolve `nt=1000, dt=3e10` to machine
   precision; every recorded array (`z_all`, `Qs_seg*`) is an equilibrium
   quantity. Because the two Neumann discretizations *agree at equilibrium*
   (the 1.4e-13 finding), applying the standalone matrix-mod Neumann BC at
   **all** channel heads — single-segment and network alike — cannot move any
   network golden value. Only the single-segment *transient* snapshots
   (`z_t0` at `nt=2`, `z_t1` at `nt=6`) move, and those are the ones slated to
   regenerate once confirmed more-correct.

Consequence: the junction stencil is re-expressed *faithfully* (same value
formulas, graph-walk indexing instead of block bookkeeping + `_ext` pads); the
head BC is *unified* to the matrix-mod form; nothing else changes.

## Test strategy

- **Hard invariants** the de-padded solver must satisfy (these validate a *new*
  discretization far better than bit-matching the old one — the point of the
  Step −1 suite): the analytical solutions (steady-state power law, uplift
  quadrature, linearized gain/lag) and sediment conservation, single-segment and
  network. These must hold to tolerance throughout.
- **Single-segment golden-master *transient* snapshots** (`z_t0`, `z_t1` in
  `characterization_reference.npz`) *will* change under the new assembly and be
  regenerated — but only after confirming the new transient is *more* correct
  (converges to the analytical / high-`niter` solution as the standalone does
  today). Equilibrium/`Q_s`/network golden values stay invariant.
- Run `pytest -m "not slow"` after each step.

## Resolved decisions

1. **`networkx` is a hard dependency**, folded in with this step (add to
   `pyproject.toml`).
2. **The single-thread path delegates to the unified network solver**: the
   standalone `LongProfile` is the front end and, in the 1-segment special case,
   runs through the same assembly. The equivalence is exact at equilibrium today;
   full transient equivalence comes from the single de-padded assembly (see
   "Phase C findings"), not from reconciling the two current solvers.
3. **Keep the names** `Network` / `LongProfile` for now; revisit a base-class
   split at Step 4 (FluvTree), when multiple segment *types* force it.
4. **Work directly on `master`**, committing granularly with the test suite
   green at each step. (No feature branch.)
