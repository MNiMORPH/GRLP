# Changelog

All notable changes to GRLP are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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

[Unreleased]: https://github.com/MNiMORPH/GRLP/compare/v2.1.0...HEAD
[2.1.0]: https://github.com/MNiMORPH/GRLP/compare/v2.0.0...v2.1.0
