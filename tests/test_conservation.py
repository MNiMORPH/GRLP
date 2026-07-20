"""
Sediment mass conservation at steady state.

With no uplift and no sediment loss, the governing equation

    dz/dt = -1 / ((1 - lambda_p) * B) * dQ_s/dx

forces the bedload sediment discharge Q_s to be spatially uniform at steady
state (dz/dt = 0  =>  dQ_s/dx = 0), and equal to the supplied input Q_s_0.

Two subtleties the tests account for:

* ``compute_Q_s`` reports the long-term-averaged discharge, i.e. it includes
  the intermittency factor I.  The supplied ``Q_s_0`` is the during-flood
  value, so the two are reconciled by dividing the computed Q_s by I.
* The residual spatial variation in Q_s is a boundary-region discretization
  error that converges away as the grid is refined (first order in dx), which
  is itself tested here.
"""

import numpy as np

from conftest import make_long_profile, STEADY_NT, STEADY_DT


def _steady(**kwargs):
    lp = make_long_profile(**kwargs)
    lp.evolve_threshold_width_river(nt=STEADY_NT, dt=STEADY_DT)
    lp.compute_Q_s()
    return lp


def test_sediment_flux_uniform_at_steady_state():
    lp = _steady()
    spread = (lp.Q_s.max() - lp.Q_s.min()) / lp.Q_s.mean()
    # Uniform to a few percent on the canonical grid; see the convergence test.
    assert spread < 0.05


def test_sediment_flux_converges_with_resolution():
    # Halving dx should roughly halve the residual non-uniformity, confirming
    # the variation is discretization error rather than a conservation flaw.
    spreads = []
    for nx, dx in [(90, 1000.0), (180, 500.0), (360, 250.0)]:
        lp = _steady(nx=nx, dx=dx)
        spreads.append((lp.Q_s.max() - lp.Q_s.min()) / lp.Q_s.mean())
    # Monotonic decrease, and each refinement cuts the error by at least ~40%.
    assert spreads[1] < 0.65 * spreads[0]
    assert spreads[2] < 0.65 * spreads[1]


def test_sediment_flux_matches_supply():
    # During-flood flux (compute_Q_s divided by intermittency) matches the
    # imposed supply Q_s_0.
    lp = _steady()
    flood_flux = lp.Q_s / lp.intermittency
    np.testing.assert_allclose(flood_flux.mean(), lp.Q_s_0, rtol=0.03)


def test_sediment_flux_matches_supply_at_outlet_fine_grid():
    # At the downstream outlet on a fine grid the discretization error is
    # smallest, so the during-flood flux should recover Q_s_0 tightly.
    lp = _steady(nx=360, dx=250.0)
    flood_outlet = lp.Q_s[-1] / lp.intermittency
    assert not np.isnan(flood_outlet)
    assert abs(flood_outlet - lp.Q_s_0) / lp.Q_s_0 < 0.005


def test_intermittency_scales_reported_flux():
    # The reported (long-term) flux scales linearly with intermittency, while
    # the underlying during-flood transport is unchanged.
    lp_full = _steady(intermittency=1.0)
    lp_half = _steady(intermittency=0.5)
    ratio = lp_half.Q_s.mean() / lp_full.Q_s.mean()
    np.testing.assert_allclose(ratio, 0.5, rtol=1e-6)
