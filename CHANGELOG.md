# Changelog

All notable changes to GRLP are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Entries for 2.0.0 and earlier summarize the existing
[GitHub Releases](https://github.com/MNiMORPH/GRLP/releases); follow each linked
version heading for the full notes.

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
