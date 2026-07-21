# Changelog

All notable changes to GRLP are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Entries for 2.0.0 and earlier summarize the existing
[GitHub Releases](https://github.com/MNiMORPH/GRLP/releases); follow each linked
version heading for the full notes.

## [Unreleased]

### Changed
- The networked solver is **de-padded**: it assembles its matrix by walking the
  topology to each node's real neighbour, instead of maintaining padded
  `z_ext`/`Q_ext` ghost arrays. Single-thread and network now share one solver
  (a single segment is a one-edge network). Results are unchanged to machine
  precision for single segments and for uniform-discharge networks; see the fix
  below for the one intended numerical change.
- `Network.compute_Q_s` (new): slope and sediment discharge at every network
  node, computed by walking the topology to each node's real neighbour (no
  `z_ext`), reusing the single-segment `LongProfile.compute_Q_s` relationship.
  Use it instead of calling `compute_Q_s` per segment on a network.
- Network boundary conditions no longer go through padded ghost arrays; they
  reuse the single-segment setters:
  - `Network.set_z_bl(z0)` is now the real base-level setter (the mouth-finding
    logic is folded in; it was previously a broken alias). It sets `z_bl` on the
    network's single outlet segment.
  - `Network.set_Qs_input_upstream(S0=None, Q_s_0=None)` replaces
    `update_z_ext_external_upstream`: it loops the channel heads and delegates to
    `LongProfile.set_Qs_input_upstream` for the `Q_s_0` case, or sets `S0`
    directly. Accepts a scalar or one value per head.
  - `LongProfile.set_Qs_input_upstream` now holds the upstream ghost-node
    elevation explicitly as `z_ghost_upstream`, writing the padded `z_ext[0]`
    only for a standalone segment (a network head has `z_ext = None`).
- A `Network` no longer maintains padded `x_ext` / `z_ext` / `Q_ext` arrays at
  all; the solver and its diagnostics walk the topology.

### Removed
- The network boundary methods `update_z_ext_external_upstream`,
  `update_z_ext_external_downstream`, and `update_x_ext_external_downstream`
  (**breaking**): use `set_Qs_input_upstream` and `set_z_bl`. Base level is the
  only downstream quantity the solver needs (its outlet ghost is
  `2*x[-1] - x[-2]`).
- The internal network padding machinery, now unused: `create_x_ext_lists`,
  `update_x_ext_*`, `update_dx_ext`, `update_dx_2cell`, `update_dx_ext_2cell`,
  `create_z_ext_lists`, `update_z_ext_internal`, `create_Q_ext_lists`,
  `update_Q_ext_*`, `update_dQ_ext_2cell`, and the dead padded block-matrix solve
  path.

### Fixed
- Networked solver, sediment conservation at confluences: the previous
  `land_area` junction discretization did **not** conserve sediment when
  discharge varied within a confluence's segments (~0.85% imbalance at steady
  state, and the junction elevation failed to converge under grid refinement).
  The de-padded solver uses a conservative three-node junction cell (a single
  shared sediment-flux conductance per junction face) that conserves to machine
  precision, is exact for uniform discharge per segment, and converges under
  refinement. Uniform-discharge networks are unchanged; only
  varying-discharge confluences move, to the conserving value. Guarded by
  `tests/test_confluence_conservation.py`.
- `Network.set_intermittency` now sets each segment's intermittency. It had
  assigned the value to the attribute holding the segment's bound
  `set_intermittency` *method* (`lp.set_intermittency = intermittency`) rather
  than calling it, so segment intermittency never changed and the method was
  overwritten with a float. Handles both the scalar and per-segment-list forms.
- Networked solver, `Q_s_0 → S0` boundary conversion: removed a spurious
  intermittency factor and corrected the sign. The network now performs this
  inversion through the single-segment `set_Qs_input_upstream`, so the two agree
  by construction. Intermittency scales the evolution rate (via `C0`), not the
  equilibrium slope, so it must not enter this geometric inversion; and the ghost
  placement `z_ghost_upstream = z[0] + S0·dx` requires a positive `S0` for a
  descending river. With
  `I ≠ 1` the old form gave the wrong slope by `I^(-6/7)` and the wrong sign, so a
  network driven by sediment supply diverged from the equivalent standalone by
  hundreds of metres. Latent: all tests and examples drive networks by `S0`.
- Networked solver, channel-head boundary: the upstream ghost discharge is now
  linearly extrapolated (`Q_ext[0] = 2*Q[0] - Q[1]`), matching the standalone
  single-segment solver, instead of held constant (`Q_ext[0] = Q[0]`). A channel
  head has no tributary junction and hence no discharge discontinuity, so the
  two-cell centered `dQ/dx` in the boundary sediment flux is well defined; the
  former constant value collapsed it to a one-sided, first-order estimate. On a
  reach with downstream-increasing discharge this biased the injected sediment
  flux and shifted the whole equilibrium profile by O(dx) — ~8.9 m at
  `dx = 1000` on the test domain — and the error only halved (first order) under
  grid refinement. It now converges second-order and tracks the analytical
  solution, matching the standalone. Networks with uniform discharge per segment
  (where `2*Q[0] - Q[1] == Q[0]`) are unaffected. Guarded by
  `tests/test_network_varying_Q.py`.
- Networked solver, river-mouth (outlet) boundary: the downstream ghost
  discharge is likewise linearly extrapolated (`Q_ext[-1] = 2*Q[-1] - Q[-2]`)
  rather than held constant, the sibling of the channel-head fix above. The
  single outlet has no tributary junction, so the same reasoning applies. With
  both boundary fixes a single-segment network now reproduces the standalone
  single-segment solver to machine precision given identical (array) inputs;
  `tests/test_network_varying_Q.py` asserts this parity and requires both fixes.

## [2.1.0] - 2026-07-20

### Added
- `LongProfile.analytical_threshold_width_uplift`: the correct semi-analytical
  steady-state long profile with uplift (or distributed base-level fall),
  replacing the never-working perturbation attempt.
- A comprehensive `pytest` test suite in `tests/` validating the model against
  analytical solutions (steady state, uplift, linearized gain/lag), checking
  sediment mass conservation, exercising the networked solver on non-trivial
  topologies, and pinning behavior with golden-master characterization tests.
  See `tests/README.md`.
- Continuous integration: a GitHub Actions workflow running the test suite
  across Python 3.9–3.13.

### Fixed
- `set_Sternberg_gravel_loss`: restored the per-kilometer → per-meter `/1000`
  unit conversion. Without it a genuine per-km abrasion coefficient overstated
  the sink by ~1000× and drove the bed elevation to `nan`; a per-km coefficient
  now yields Sternberg decay `Q_s = Q_s,0 * exp(-k x)`.
- `compute_e_folding_time`: corrected a method-name typo (`self.wavenumber` →
  `self.compute_wavenumber`) that made every call raise `AttributeError`.
- `set_Q`: now stores `P_xQ`, so the analytical solutions (which read
  `self.P_xQ`) stay consistent with the discharge field when `P_xQ` is
  overridden.

### Deprecated
- `LongProfile.analytical_threshold_width_perturbation`: incorrect (returns
  unphysical elevations) and retained only for the historical record. Use
  `analytical_threshold_width_uplift` instead.

## [2.0.0] - 2025-09-10
1D long-profiles and updated examples.

## [2.0.0-beta] - 2024-11-28
Version 2.0.0 beta.

## [2.0.0-alpha] - 2023-10-01
Version 2.0.0 alpha.

## [1.8.0] - 2023-07-04
Network tools.

## [1.7.0] - 2022-08-05
McNab 2022.

## [1.6.0] - 2021-12-11

## [1.5.0] - 2021-11-30
Networks and Jupyter.

## [1.4.1] - 2020-10-10
Check published version in Firefox.

## [1.4.0] - 2020-10-10
PyPI tools + v1.3.1 updates.

## [1.3.1] - 2020-10-08
Channel width and depth, and S0 sign.

## [1.3.0] - 2020-07-26
Valley width less important.

## [1.2.3] - 2020-04-10
PyPI Integration.

## [1.2.0] - 2020-04-09
Fully linked network.

## [1.1.0] - 2020-04-05
Ghost-node network.

## [1.0.0] - 2018-11-29
ESurf publication.

## [1.0.0-alpha] - 2018-05-22
ESurf-D Submission.

## [0.0] - 2018-05-06
ESurf-D initial release.

[Unreleased]: https://github.com/MNiMORPH/GRLP/compare/v2.1.0...HEAD
[2.1.0]: https://github.com/MNiMORPH/GRLP/compare/v2.0.0...v2.1.0
[2.0.0]: https://github.com/MNiMORPH/GRLP/releases/tag/v2.0.0
[2.0.0-beta]: https://github.com/MNiMORPH/GRLP/releases/tag/v2.0.0-beta
[2.0.0-alpha]: https://github.com/MNiMORPH/GRLP/releases/tag/v2.0.0-alpha
[1.8.0]: https://github.com/MNiMORPH/GRLP/releases/tag/v1.8.0
[1.7.0]: https://github.com/MNiMORPH/GRLP/releases/tag/v1.7.0
[1.6.0]: https://github.com/MNiMORPH/GRLP/releases/tag/v1.6.0
[1.5.0]: https://github.com/MNiMORPH/GRLP/releases/tag/v1.5.0
[1.4.1]: https://github.com/MNiMORPH/GRLP/releases/tag/v1.4.1
[1.4.0]: https://github.com/MNiMORPH/GRLP/releases/tag/v1.4.0
[1.3.1]: https://github.com/MNiMORPH/GRLP/releases/tag/v1.3.1
[1.3.0]: https://github.com/MNiMORPH/GRLP/releases/tag/v1.3.0
[1.2.3]: https://github.com/MNiMORPH/GRLP/releases/tag/v1.2.3
[1.2.0]: https://github.com/MNiMORPH/GRLP/releases/tag/v1.2.0
[1.1.0]: https://github.com/MNiMORPH/GRLP/releases/tag/v1.1.0
[1.0.0]: https://github.com/MNiMORPH/GRLP/releases/tag/v1.0.0
[1.0.0-alpha]: https://github.com/MNiMORPH/GRLP/releases/tag/v1.0.0-alpha
[0.0]: https://github.com/MNiMORPH/GRLP/releases/tag/v0.0
