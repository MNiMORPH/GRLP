# Solver extraction — design (for review before any code moves)

Goal: pull the **numerical engine** out of the `Network` class so the class is a
data structure + descriptors, and the solver is a separate component. Resolves
"solver core and plotting flourish in the same class" and makes "1D = trivial
N=1 network" literal. Behavior-preserving refactor; the 321-test suite guards
every step. Payoff is conceptual clarity, **not** correctness — so it is done
carefully, green at each step, not in a scramble.

## Guiding model (Andy's, adopted)

**Every solution is a network solution; a `LongProfile` is always solved as a
one-edge network.** There is exactly one solver interface — it takes a network —
and a single segment is the degenerate case, not a special path.

## Target shape

Three roles, but only **two** modules move apart (the data classes stay together
because they are mutually coupled):

- **`LongProfile`** (segment) — unchanged home. Holds the segment's *data*
  (`x, z, Q, B, A`, boundary scalars) **and** its *physics* (`basic_constants`,
  `build_LHS_coeff_C0`, the stencil coefficients `C`), analytical solutions, and
  the convenience solve-wrappers (`evolve_threshold_width_river`, `compute_Q_s`,
  `slope_area`) that wrap `self` in `Network([self])` — i.e. treat the segment as
  the one-edge network it is.
- **`Network`** (aggregate + descriptors) — unchanged home. A list of segments +
  topology (`build_graph`, `find_*_IDs`, channel-head/mouth lists), boundary
  setters, geometry setup (`get_z_lengths`, `compute_land_areas_around_confluences`),
  the diagnostic walk `compute_Q_s`, the analysis suite (`compute_*`,
  `find_hack_parameters`), and `plot()`. **These cohere** — they are all "things
  about a network you already have." Once the engine leaves, the class no longer
  mixes core with flourish.
- **`Solver`** (engine) — **NEW module `grlp/solver.py`**, module-level functions
  operating on a **duck-typed** network (they read `network.list_of_LongProfile_objects`,
  `network.graph`, `network.IDs`, `lp.C`, `lp.build_LHS_coeff_C0`, … but do **not**
  `import` the `Network`/`LongProfile` classes). This is what keeps the
  dependency one-way.

## Why `LongProfile` and `Network` stay in one module

`Network.initialize` **creates** `LongProfile` objects (network -> longprofile),
and `LongProfile.evolve` wraps itself in a `Network` (longprofile -> network).
That mutual reference is intrinsic to the model, not an accident. Splitting them
into separate modules forces a circular import (worked around only with
lazy/in-function imports — a smell for no gain). So they co-locate. The solver,
by contrast, needs *neither class by name* — only their duck-typed data — so it
extracts cleanly with a one-way dependency.

## Import graph (acyclic)

```
grlp.py (LongProfile, Network)  ->  solver.py  ->  numpy / scipy
```

`grlp.py` imports `solver` (Network's evolve calls it). `solver.py` imports
nothing from `grlp.py`; it operates on whatever network-shaped object it is
handed. No cycle.

## Move list

**-> `grlp/solver.py`** (become functions taking a network):
- `Network.assemble_by_walking(self, dt)` -> `solver.assemble(net, dt)`
- `Network._evolve_by_walking(self, nt, dt)` -> `solver.evolve(net, nt, dt)`
  (carries the Picard/`niter` loop)
- `Network.update_gravel_loss(self)` -> `solver.update_gravel_loss(net)`
  (called inside the evolve loop; a per-iteration source recompute)

**Stays on `Network`, becomes a thin caller:**
- `Network.evolve_threshold_width_river_network(self, nt, dt)` -> body becomes
  `solver.evolve(self, nt, dt)` (+ the existing `Qs_internal` bookkeeping).

**Unchanged** (Network): `initialize`, `build_graph`, `find_*_IDs`,
`create_list_of_channel_*`, `get_z_lengths`, `compute_land_areas_around_confluences`,
`set_niter`/`set_z_bl`/`set_x_bl`/`set_Qs_input_upstream`, `compute_Q_s`
(diagnostic), the analysis suite, `find_hack_parameters`, `plot`.

**Unchanged** (LongProfile): everything, including the `Network([self])` wrappers.

## Sequencing (green at each step; granular commits)

1. **In-place extraction.** Add `grlp/solver.py`. Move the three method bodies
   into `assemble(net, dt)` / `evolve(net, nt, dt)` / `update_gravel_loss(net)`.
   Point `Network.evolve_threshold_width_river_network` at `solver.evolve(self, …)`.
   Caller check (done): `_evolve_by_walking` and `update_gravel_loss` are
   internal-only -> delete cleanly. `assemble_by_walking` has **one external
   caller**, `tests/test_depad_walk.py` -> update that test to call
   `grlp.solver.assemble(net, DT)` (preferred: it tests the assembler, which now
   lives in the solver), *or* keep a one-line `Network.assemble_by_walking`
   delegator. Run full suite -> must be **321**.
2. **Confirm no import cycle**, `import grlp` clean, examples run. Commit.
3. **(Optional, later)** split the analysis suite into `grlp/analysis.py` the same
   way (functions on a duck-typed network, or leave as methods — lower priority;
   analysis-as-methods is fine). Not required to meet the goal.

## What this deliberately does NOT do

- **No `LongProfile` rename** — it is the published identity (Wickert & Schildgen
  2019); still resisted.
- **Does not split `LongProfile` from `Network`** — intrinsic coupling; see above.
- **Does not move analysis/plot off `Network`** — with the engine gone, they
  cohere as network descriptors/representations; moving them is optional polish,
  not part of resolving the core-vs-flourish concern.
- **No numerical change** — the goldens and 321 tests must be untouched; this is
  pure relocation.

## Risks / guards

- Behavior-preserving relocation; the 321-test suite (incl. characterization
  goldens, confluence conservation, single-segment parity) is the safety net at
  every commit.
- Main hazard is an overlooked external caller of a moved method or a subtle
  import ordering issue — caught by (a) grepping call sites before deletion and
  (b) `import grlp` + example runs after each step.
