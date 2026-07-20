# Resume: Step 1 solver de-pad (port the neighbor-walking assembler)

Pick-up point for the **de-pad port** into `grlp.py` (task: unify the solver by
replacing padded `_ext` arrays with graph-walking neighbor lookup). Read this
plus `step1-network-primary-design.md` (the design + the confluence decision) to
resume cold.

## Where we are (all committed)

Done and green (`pytest -m "not slow"` = 290 passing):

- **Boundary ghost-discharge bugs fixed** (channel head + river mouth): the
  network held `Q_ext` constant at the domain ends; now linearly extrapolated
  (`2*Q[0]-Q[1]`, `2*Q[-1]-Q[-2]`), matching the standalone. First-order →
  second-order. Guarded by `tests/test_network_varying_Q.py`.
- **Intermittency bugs fixed**: `Network.set_intermittency` overwrote the
  segment method with a float; the `Q_s_0 -> S0` conversion had a spurious `I`
  factor and wrong sign. Guarded by `tests/test_intermittency.py`.
- **Design doc corrected** (`step1-network-primary-design.md`): the old "inert
  Picard / irreconcilable scheme" Phase C story was wrong; replaced with the
  ghost-bug account. **Confluence decision recorded** there (two regimes; naive
  FV is first-order and rejected).
- **Safety-net golden**: `confluence_varying_Q` in the characterization suite
  pins current multi-tributary behavior (all other network fixtures are
  uniform-Q). It must **not** move during the easy-half de-pad.

Validated but **not yet in `grlp.py`**: the neighbor-walking assembler
reproduces `build_matrices` **bit-for-bit** (0.0) — see
`prototypes/walk_assemble.py`.

## The architecture (confirmed with Andy)

**Reuse GRLP's exact stencil formulas; walk the graph only for neighbor lookup.**
Do not re-derive the stencil or BCs — a from-scratch reimplementation was 0.125 m
off (`prototypes/neighbor_walk_v2.py`). The numbers going into
`C1`/`left`/`center`/`right`/`set_bcl_Neumann`/`set_bcr` are unchanged; only the
mechanism that supplies a node's upstream/downstream neighbor `z,Q,dx` changes
(from padded `z_ext`/`Q_ext` indexing to walking to the real neighbor node,
which at a single-upstream junction is the last node of the upstream segment).

## Next steps (the port)

- **2a — DONE** (commit `a52eeaf`). `Network.assemble_by_walking(dt)` builds the
  global LHS+RHS by walking; reproduces `build_matrices` bit-for-bit for a single
  segment (`tests/test_depad_walk.py`). Single-upstream junctions already walk
  across the boundary and get the interior stencil; multi-tributary confluences
  raise `NotImplementedError`. Additive; **not yet wired into the evolve loop**.
- **2b — 1-into-1 junction validated** (commit `5542635`). A 1-into-1 chain
  assembled by `assemble_by_walking` equals the single segment bit-for-bit
  (`tests/test_depad_walk.py::test_walk_1into1_chain_equals_single_segment`).
  The fix is proven at the assembly level; wiring into the evolve loop remains.
- **WIRING — DONE for single-segment + chains** (commit `e987159`).
  `evolve_threshold_width_river_network` dispatches to `_evolve_by_walking` when
  the network has no multi-tributary confluence (`any(len(upstream)>1)` is
  False). Proper Picard (RHS from `zold`, coefficient from the iterate).
  Verified: single-segment networks still match the standalone; a varying-Q
  1-into-1 chain now matches a single varying-Q segment to 0.0 (was ~0.8 m).
  `chain_uniform_B` did **not** move — with uniform Q the `dQ/dx` terms vanish,
  so the junction handling is irrelevant there (earlier "it will move" was
  wrong). Full suite 296 green.
- **2c — multi-tributary confluence (NEXT, the hard half).** The evolve branch
  still sends any network with a ≥2-tributary confluence down the old padded
  block-matrix path; the walker raises `NotImplementedError` on such a node.
  Two sub-options: (a) *delegate* — reproduce the current `land_area` block
  entries in the walker so `confluence_*`/`topology_*`/`confluence_varying_Q`
  goldens stay bit-identical, then flip the branch (banks the padless
  refactor but keeps the first-order junction); (b) *the real fix* — a
  second-order flux-balance junction cell (design doc "Confluence handling").
  **Findings (measured; do not re-walk):**
  1. *The current confluence is EXACT for uniform Q per segment* — the symmetric
     junction elevation is `600.0` at every `dx` (constant-slope profile, no
     discretization error). It is first-order **only** for within-segment-varying
     Q (the `confluence_varying_Q` case). So the standard confluences are fine;
     the "hard half" is a narrow edge case, and 2c-delegate keeps all standard
     confluence goldens exact.
  2. *Two re-derivations both fail conservation:* averaging interior stencils
     (~53% junction imbalance) and a one-sided flux-balance cell (~40%, junction
     z 388 vs 600). **Root cause:** each junction-adjacent face is shared between
     a face-based cell and a node-based interior neighbor (the tributary's last
     node, the confluence's downstream node) → the face flux is double-valued →
     mass leaks. A conserving confluence needs the *whole junction neighborhood*
     (confluence node + tributary last nodes + downstream node) to use
     **consistent shared face fluxes** on the junction faces.
  3. The current `land_area` code already arranges this two-sided handshake
     (`add_block_diagonal_matrix_upstream/downstream_boundary_conditions`). So
     **2c-delegate = reproduce those entries via the walk**, not re-derive them;
     re-derivation keeps hitting the two-sided-consistency wall. A clean
     second-order confluence (2c-b) is a separate, genuine numerical-design task
     (naive single-valued FV in `prototypes/fv_prototype.py` is first-order) —
     fresh session.
- **remove** the padded `_ext` machinery once 2c lands.
- **(later, hard half)** swap the delegated multi-tributary path for the
  second-order flux-balance junction cell (7/6 tangent flux, reduces to interior
  in the single-tributary limit). See the design doc's confluence section.

## Gotchas (each of these bit me)

1. **Validate against the array-Q convention, not `k_xQ`.** `make_long_profile`
   uses `set_Q(k_xQ, P_xQ)`, so its standalone evaluates the *exact power law* at
   the ghost. The network (and walker) use *linear* ghost extrapolation. They
   differ at O(dx^2). To check bit-for-bit, build the reference standalone with
   `set_Q(Q=array)` (see `walk_assemble.py`), or compare to a single-segment
   *network*.
2. **`set_z(z=...)` does NOT sync `z_ext[1:-1]`.** The evolve loop syncs it; a
   manual assembly at a chosen `z` must do `lp.z_ext[1:-1] = z` first, or
   `build_matrices` computes `C1` from stale interior z.
3. **`grlp` is a package**; the solver lives in the `grlp.grlp` submodule.
   Monkeypatching `spsolve` for diagnostics must target `grlp.grlp.spsolve`.

## Validated prototypes (in `prototypes/`)

- `walk_assemble.py` — neighbor-walk == `build_matrices`, 0.0. **Port from this.**
- `fv_prototype.py` — from-scratch FV; measured first-order (why naive FV is
  rejected).
- `neighbor_walk_v2.py` — junction no-op demo; also shows why not to re-roll BCs.
