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

## A known API sharp edge exercised here

`set_Q(P_xQ=...)` sets the discharge array but does **not** update `self.P_xQ`,
which the analytical solutions read. The test factory keeps them consistent; if
you change `P_xQ` yourself, set `lp.P_xQ` to match, or the analytical exponent
will disagree with the discharge field.
