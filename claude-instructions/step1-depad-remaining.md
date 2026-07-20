# Step 1 de-pad — remaining work (post-compaction handoff)

Authoritative status + instructions for finishing the de-pad. Supersedes the
stale body of `step1-depad-resume.md` (whose top banner is still correct).

## Status: de-pad functionally COMPLETE

`evolve_threshold_width_river_network` dispatches **all** networks to
`_evolve_by_walking`, which assembles the global matrix by walking the topology
(`assemble_by_walking`) — no `z_ext` used in the solve. Single-thread and
network are one solver. Done and committed (HEAD `7356725`):

- Boundary ghost-discharge bugs fixed (head + mouth); intermittency bugs fixed.
- 1-into-1 junction fixed (interior stencil across the boundary).
- **Multi-tributary confluence: conservative three-node junction cell** — each
  junction face carries one shared conductance `D` used by both adjacent nodes,
  so flux `D*(z_up-z_dn)` is single-valued → conserves to machine precision.
  Exact for uniform Q per segment (bit-matches old to 2.8e-13 → uniform-Q
  goldens preserved); first-order convergent. Fixes a real conservation bug
  (old `land_area` leaked ~0.85% for within-segment-varying Q).
- Dead padded evolve body deleted.
- Tests: `test_depad_walk.py`, `test_confluence_conservation.py` (single +
  nested varying-Q conservation + slow convergence-order), `test_network_varying_Q.py`,
  `test_intermittency.py`. Full suite **298 fast / 302 with slow**, green.
- Only `confluence_varying_Q` golden moved (regenerated to the conserving value).

## Remaining work, in order

### 1. Migrate the last `z_ext` consumers to `lp.z`  (prerequisite for deletion)
The walker does not maintain `z_ext`; a **compatibility sync** at the end of
`_evolve_by_walking` (in `grlp/grlp.py`, look for "Sync the padded z_ext") keeps
it fresh only for these consumers:
- `LongProfile.compute_Q_s` — computes slope from `z_ext[2:]-z_ext[:-2]`
  (2-cell). Rewrite to use `self.z` plus neighbor ghosts. **For a network node
  the upstream/downstream neighbor is across a segment boundary** — reuse the
  same neighbor-walk logic as `assemble_by_walking` (head: `2*z[0]-z[1]`-style
  ghost via S0; outlet: `z_bl`; junction: the real neighbor node). This is the
  one non-trivial rewrite; give it care and a test comparing to the current
  (sync-fed) `compute_Q_s` output before/after.
- `examples/network/*` scripts that read `z_ext` (NewNetwork_1segment.py,
  TestQ.py, whitewater_2m_net_*.py, netBLfall.py, ...): smoke-run them; switch
  any `z_ext` reads to `z`.

### 2. Drop the compatibility sync
Once #1 lands, delete the "Sync the padded z_ext" block at the end of
`_evolve_by_walking` and its `update_z_ext_internal` / `update_z_ext_external_upstream`
calls.

### 3. Delete the `_ext` machinery
Now safe to remove (verify no remaining refs first — some are referenced by
tests/examples, so update those too): `create_x_ext_lists`, `update_x_ext_*`,
`update_z_ext_*`, `update_Q_ext_*`, `update_dQ_ext_2cell`,
`network__compute_coefficient_time_varying`, `network__build_matrix_inner`,
`map_block_diagonal_matrix_blocks`,
`create_block_diagonal_matrix_with_internal_tridiagonals`,
`add_block_diagonal_matrix_upstream/downstream_boundary_conditions`,
`stack_RHS_vector`, the `z_ext`/`Q_ext` setup in `initialize`, and the padded
(`type(...) is list`) branches in `set_bcr_Dirichlet` / `set_bcl_Neumann_*`.
Run full suite (incl. slow) after.

### 4. Optional: second-order junction
The junction cell is first-order (correct + conserving). To raise it, the face
`D` uses a 1-cell face slope; a 2nd-order face reconstruction would lift it.
Not required.

### 5. Release-gated (each needs Andy's explicit go)
Lint; review the whole session diff; push. Then, only if asked: version bump
(`_version.py`, `CITATION.cff`, CHANGELOG roll), tag, release.

## Key facts / gotchas
- Validate walker vs the **array-Q** standalone (linear ghosts), not `k_xQ`
  (analytic ghost) — they differ at O(dx^2).
- Junction cell needs confluence segments >= 3 nodes, tributaries >= 2 (guarded
  with a clear `ValueError` in `assemble_by_walking`).
- Conservation is guaranteed by the **shared** `D` per junction face; that was
  the key. Stencil-averaging and one-sided cells leak.
- Everything is committed and push-ready; the solver change moved only the one
  intended golden.
