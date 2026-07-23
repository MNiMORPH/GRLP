# Examples

Runnable example scripts live in the repository's `examples/` directory. They are
curated into groups; each script runs against the current (v3) API.

## Getting started

- **`examples/run_grlp.py`** — a commented driver that reproduces a figure from
  Wickert & Schildgen (2019); the gentlest introduction to setting up a run.
- **`examples/example_1d.ipynb`** — a Jupyter notebook with a more extensive
  1-D tutorial.

## One-dimensional (`examples/one_dimensional/`)

Single-profile examples that teach one feature each:

- **`Basic1D.py`**, **`Basic1D-var-dx.py`** — a minimal profile to steady state
  and a transient; the second uses a non-uniform grid.
- **`analytical_numerical.py`**, **`analytical_numerical_irregular_grid.py`**,
  **`analytical_numerical_uplift.py`**, **`analytical_concavity.py`** — compare
  the numerical solution against the analytical steady-state solution (§ Theory).
- **`base_level_fall_transience.py`**, **`sudden_bl_fall_U0.py`** — response to
  base-level fall.
- **`RioSantaCruz_set_S0_set_x_bl.py`** — force the upstream boundary by slope
  (`set_S0`) and move base level along a continental-shelf gradient
  (`set_x_bl` + `set_z_bl`), after Ruby et al. (2026).
- **`oscillatingQ_2panel.py`** — periodic discharge forcing (gain/lag behavior).
- **`test_uplift_equals_qs.py`** — a worked balance between uplift and supply.

## Networks (`examples/network/`)

- **`shreve_random_network.py`** — generate a random Shreve network, evolve it,
  and plot with `Network.plot()`.
- **`NewNetwork_1segment.py` … `NewNetwork_5segments.py`** — small,
  explicitly-constructed networks; good for understanding confluence handling.
- **`landslide5.py`**, **`landslideContinuous5.py`** — point-source sediment
  input (a single pulse, and continuous supply) representing landslides.
- **`netBLfall.py`**, **`netBLrise_animation.py`** — network response to
  base-level change; the second renders an animation.

## Spectral response (`examples/McNab_et_al_GRL/`)

Scripts reproducing figures from McNab et al. (2023) — gain/lag grids, profiles,
time series, equilibration times, and a benchmark against the analytical
solution. If you use these features, cite McNab et al. (2023) (see {doc}`citing`).

## Deprecated (`examples/deprecated/`)

Frozen scripts from earlier versions, kept for provenance. They are **not**
maintained against the v3 API and may not run; use the curated examples above.
