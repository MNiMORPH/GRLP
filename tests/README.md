# GRLP test suite

Tests for the gravel-bed river long-profile model. They check both that the
model is internally self-consistent and that its numerical results reproduce
known analytical solutions of the transport-limited threshold-width equations.

## Running

```sh
pip install -e ".[test]"      # installs pytest
pytest                        # full suite (includes the slow benchmark)
pytest -m "not slow"          # fast suite only (default development loop)
pytest -m slow                # only the numerical spectral benchmark
```

## Layout

| File | What it checks | Reference |
|------|----------------|-----------|
| `conftest.py` | Shared `make_long_profile` factory and fixtures | — |
| `test_setup.py` | Setter internal consistency (grid, ghost nodes, power-law A/Q/B, boundary slope) before any time integration | — |
| `test_steady_state.py` | Numerical equilibrium vs. the analytical power-law bed `z = k_a x^β + c_a`; slope–area concavity | `analytical_threshold_width`; Wickert & Schildgen (2019) |
| `test_conservation.py` | Sediment discharge uniform at steady state (converging under grid refinement); reconciliation with intermittency | mass conservation |
| `test_uplift.py` | Numerical equilibrium **with uplift** converges (first order) to the semi-analytical solution | `analytical_threshold_width_uplift` |
| `test_network.py` | Confluence slopes, network-wide sediment conservation, uniform-discharge chain reduces to the exact linear profile | networked mass conservation |
| `test_spectral.py` | Analytical limits of the linearized gain/lag (quasi-static → 6/7 & 0; fast → 0; monotonicity, bounds) | `compute_z_gain` / `compute_z_lag`; McNab et al. (2023) |
| `test_spectral_benchmark.py` | **(slow)** Full numerical periodic forcing vs. the linearized gain/lag | as above; cf. `examples/McNab_et_al_GRL/Figure_S1_Benchmark.py` |
| `test_network_correctness.py` | Known-true physics on junction-stressing topologies (asymmetric `Q`, unequal `dx`, unequal lengths, multi-level) | network mass conservation; shared `network_helpers.NETWORK_TOPOLOGIES` |
| `test_single_segment_physics.py` | Irregular grid, distributed source/sink ≡ uplift, base-level fall, Sternberg loss | analytical; `set_source_sink_distributed`, `set_Sternberg_gravel_loss` |
| `test_api_and_metrics.py` | Response timescales, Strahler/Horton on a balanced tree, `Shreve_Random_Network` invariants | known metric values |
| `test_characterization.py` | Golden-master: current outputs of many prescribed-`B` configs and network topologies match a recorded reference to ~1e-9 | `characterization_reference.npz` |

## Known issues these tests surfaced (tracked as strict `xfail`)

Building this suite uncovered two live bugs, documented as `xfail` so they stay
visible and will flip to failing (`xpass`) the moment they are fixed:

- **`set_Sternberg_gravel_loss` units** — the finite-difference sink dropped the
  per-km → per-m `/1000` conversion that the original (commented) implementation
  had, so a genuine per-km loss overstates the sink ~1000× and drives `z` to
  `nan`. See `test_single_segment_physics.py::test_sternberg_gravel_loss_exponential_decay`.
- **`compute_e_folding_time` typo** — calls `self.wavenumber(n)`, but the method
  is `compute_wavenumber` → `AttributeError`. See
  `test_api_and_metrics.py::test_e_folding_time_formula`.

A diagnostic of the padded-array network solver at asymmetric / unequal-`dx` /
multi-level junctions found it **numerically correct** (conservation to ~1e-14,
two-sided junction ghosts agree exactly), so the planned de-padding refactor is
architectural cleanup, not a bug fix.

## Characterization (golden-master) tests

`test_characterization.py`, `characterization_configs.py`, and
`characterization_reference.npz` pin the *current* numerical behavior of the
model across many prescribed valley-width (`B`) configurations — uniform,
power-law, arbitrary `B(x)` array, with/without uplift and intermittency,
grain-size channel width/depth, and networks with per-segment widths — capturing
transient snapshots and equilibrium plus `B`-loaded quantities (diffusivity,
`Q_s`, slope). They exist to guard work that makes `B` a dynamic function of the
other variables: when that code is set to reproduce a prescribed `B`, these
confirm the results are unchanged. Deliberately `B`-loaded, because the
no-uplift equilibrium profile is itself `B`-independent.

To adopt new values after a *knowing* change, regenerate and review the diff:

```sh
python tests/generate_characterization_reference.py
git diff --stat        # an unexpected change here is the regression to catch
```

## Notes on the analytical solutions

* **No-uplift steady state** and the **linearized gain/lag** are closed form and
  are matched tightly (rel. err ~1e-4 and asymptotic limits).
* **Uplift steady state** has no elementary closed form for general exponents;
  `analytical_threshold_width_uplift` evaluates it by quadrature on a
  grid independent of the model, so the numerical model converges to it at first
  order rather than matching at a single resolution.
* `analytical_threshold_width_perturbation` is **deprecated and incorrect** (an
  early perturbation attempt that returns unphysical elevations); it is retained
  only for the historical record. Use `analytical_threshold_width_uplift`.
