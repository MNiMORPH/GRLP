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
"single = 1-edge network" is numerically exact — the façade can route
1-D runs through the network solver without changing results.

> **Correction (see "Phase C findings, corrected").** That 1.4e-13 was measured
> with *uniform* discharge. With downstream-increasing `Q` the 1-segment network
> was off by ~8.9 m, due to two boundary ghost-discharge bugs — **now fixed**
> (commits `1deeea9`, `b5e0594`); post-fix the equivalence is genuinely ~1e-13
> for varying `Q` too. The lesson: "the existing golden-master + analytical
> tests will catch any drift" was *false* here — every fixture used uniform `Q`,
> so the bug was invisible. New varying-`Q` tests (`test_network_varying_Q.py`)
> close that gap.

Also established in Step −1: the padded network solver **conserves sediment** at
asymmetric / unequal-`dx` / multi-level junctions (~1e-14, ghosts agree exactly).
Conservation is exact, but note this does *not* imply the junction is
high-order: the confluence coupling is first-order (~0.8 m/junction, O(dx)) — a
conservative-but-first-order scheme. That is the de-pad's real target.

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
- **C/D — de-pad the solver.** *Rescoped (see "Phase C findings, corrected"
  below).* The single-segment vs. network discrepancy that originally motivated
  de-padding turned out to be **two boundary ghost-discharge bugs, now fixed**
  (commits `1deeea9`, `b5e0594`) — not a fundamental two-implementation
  conflict. With those fixed, a 1-segment network reproduces the standalone to
  ~1e-13, so single-segment *unification* no longer needs de-padding. What
  remains for the de-pad is the **confluence coupling**: an internal junction is
  still first-order (~0.8 m/junction, O(dx)) because the discharge is
  discontinuous there and the current stencil handles it with one-sided
  differences. Replacing the padded `_ext`/block machinery with a single
  neighbor-walking, **finite-volume flux-divergence** assembly is the vehicle to
  make the confluence second-order and conservative by construction (sum of
  tributary face-fluxes = outflow face-flux; discharge never differenced across
  the jump). That — not the boundary — is the real justification now.

## Phase C findings (corrected)

> The original Phase C write-up (two "causes" — an inert Picard loop and an
> irreconcilable boundary discretization — motivating de-padding as the
> unification vehicle) was **wrong on both counts**. It was measured on a set-up
> that conflated a real bug with intermittency and initial-condition
> mismatches. The corrected account, verified in detail, is below. The prior
> "Golden-safety" subsection (which asserted the two discretizations "agree at
> equilibrium") is superseded and removed.

**The 1-segment network did *not* equal the standalone even at equilibrium — it
was off by ~8.9 m at `dx = 1000` on a reach with downstream-increasing `Q`.**
Chased to ground (matched intermittency, matched IC, matched `S0`, grid
refinement, term-by-term matrix diff, reproduced identically on `v2.0.0` so not
a Phase B regression):

1. **The Picard loop is *active*, not inert.** The iterate converges across
   `niter` (0.12 → 3.5e-4 → 1e-6 per iteration). The earlier "inert Picard"
   claim was false; no Picard fix was needed.
2. **The discrepancy was two boundary ghost-discharge bugs** — the LHS stencil
   and Neumann matrix-modification are in fact *bit-identical* between the two
   solvers; the whole difference lived in one RHS term, the boundary `dQ/dx`.
   The network held the **ghost discharge constant** (`Q_ext[0] = Q[0]` at the
   head, `Q_ext[-1] = Q[-1]` at the mouth), collapsing the two-cell centered
   `dQ/dx` to a one-sided, first-order estimate; the standalone linearly
   extrapolates (`set_Q`: `2*Q[0]-Q[1]`, `2*Q[-1]-Q[-2]`), which is second-order.
   Neither boundary has a discharge discontinuity, so the linear ghost is simply
   correct. **Fixed** in `update_Q_ext_external_upstream` / `_downstream`
   (commits `1deeea9`, `b5e0594`). The error halved (first order) under grid
   refinement before the fix and quarters (second order) after; a single-segment
   network now reproduces the standalone to **~1e-13** given identical array
   inputs (`tests/test_network_varying_Q.py`).

Why it stayed hidden: every pre-existing network test used **uniform `Q` per
segment**, where `2*Q[0]-Q[1] == Q[0]` and `dQ/dx == 0`, so the bug is invisible
and no golden value moves. The earlier "1.4e-13 agreement" and the committed
"golden-safe" finding (`2c2bb60`) were measured in exactly that uniform-`Q`
regime — the *conclusion* (golden-safe) held, but only because the bug cannot
manifest there, not because the discretizations are equivalent in general.

Consequence for de-padding: the **boundary** is no longer a reason to de-pad —
it is fixed and the two solvers agree. The remaining first-order error is at the
**confluence**, where `Q` genuinely *is* discontinuous (each tributary carries
its own discharge; the reach below carries the sum), so a centered `dQ/dx`
straddles the jump and the current one-sided handling is first-order
(~0.8 m/junction, accumulating). That is the real target for the de-pad's
finite-volume flux-divergence form (see the C/D bullet above).

## Confluence handling: finite volume vs. finite difference (decision)

### The implied face flux of the current FD interior

The interior transport at node `i` (uniform grid, `dQ` term dropped for clarity)

```
T_i = (7/6) · C1_i / dx^2 · (2 z_i - z_{i-1} - z_{i+1})
    = [ F_{i+1/2} - F_{i-1/2} ] / dx ,     C1_i = C0 · |S_i|^{1/6} · Q_i / B_i
F_{i+1/2} = (7/6) · C1_i · (z_i - z_{i+1}) / dx
```

Two facts drive everything:

1. The `7/6` is the **tangent (Newton) linearization** of `Q_s = K Q S^{7/6}`
   (`d/dS S^{7/6} = 7/6 S^{1/6}`). `land_area_around_confluence` was built with
   the **secant** flux (frozen `|S|^{1/6}`, no `7/6`); the missing factor is the
   measured ~6/7 (0.85x) weakness at the confluence node.
2. `C1_i` is evaluated at the **node** (2-cell slope), so `F_{i+1/2}` is
   **double-valued** — nodes `i` and `i+1` compute the shared face's flux with
   different `C1`. The interior is therefore a **node-based, second-order but not
   strictly conservative** scheme. It does not drop into clean finite-volume
   form.

### Measured: naive finite volume is first-order here

A from-scratch single-valued-face-flux FV (secant flux, face-averaged `Q`,
node-centered cells; prototype `fv_prototype.py`) converges at **first order**
(order ~1.0, ~1.35 m at `dx = 1000` vs. the GRLP standalone, which is
second-order). So converting to FV is **not** a free accuracy win — a careless
FV downgrades the interior from second to first order. A second-order FV is
achievable but needs careful face reconstruction; it is real work, not a swap.

### Decision (two regimes)

- **Discharge-continuous junctions (single upstream segment)** use the standard
  second-order FD interior stencil with **neighbor-reaching** across the segment
  boundary — the confluence node is just a regular interior node. Exact
  (reduces to the single segment to ~1e-13); no special cell.
- **Genuine multi-tributary confluences (discharge discontinuity)** use a
  **flux-balance junction cell** — the FV idea (sum of tributary face-fluxes =
  outflow face-flux, cell area = `land_area`, no `dQ/dx` across the jump) — but
  with each face flux carrying the interior's **tangent (`7/6`) linearization**
  and constructed to reduce to the interior stencil in the single-tributary
  limit, so it stays second-order.
- A naive secant / face-averaged FV is **explicitly rejected** (measured
  first-order), as is a full model-wide FV conversion (it would change every
  interior number and move all golden values for no accuracy gain).

## Test strategy

- **Hard invariants** the de-padded solver must satisfy (these validate a *new*
  discretization far better than bit-matching the old one — the point of the
  Step −1 suite): the analytical solutions (steady-state power law, uplift
  quadrature, linearized gain/lag) and sediment conservation, single-segment and
  network. These must hold to tolerance throughout.
- **Single-segment golden values stay invariant.** The boundary fix already
  brought the single-segment network into ~1e-13 parity with the standalone, and
  a de-pad that preserves single-segment behavior leaves `z_t0`/`z_t1` and the
  single-segment characterization arrays untouched. (The boundary-ghost fixes
  themselves moved *no* golden value, since all fixtures use uniform `Q`.)
- **Network/confluence golden values *will* move** under a finite-volume
  confluence rewrite — even with uniform `Q` per segment, the discharge jumps at
  each junction, so the junction discretization changes the steady state by
  O(dx). Regenerate them only after confirming the new confluence is *more*
  correct against the hard invariants (analytical amplitude + exact conservation)
  and second-order under refinement.
- Run `pytest -m "not slow"` after each step.

## Resolved decisions

1. **`networkx` is a hard dependency**, folded in with this step (add to
   `pyproject.toml`).
2. **The single-thread path delegates to the unified network solver**: the
   standalone `LongProfile` is the front end and, in the 1-segment special case,
   runs through the same assembly. Equivalence is now **exact** (~1e-13, transient
   and steady) once the boundary ghost-discharge bugs are fixed — so delegation
   is safe today without any de-pad. The de-pad is about the *confluence*, not
   this single-segment equivalence.
3. **Keep the names** `Network` / `LongProfile` for now; revisit a base-class
   split at Step 4 (FluvTree), when multiple segment *types* force it.
4. **Work directly on `master`**, committing granularly with the test suite
   green at each step. (No feature branch.)
