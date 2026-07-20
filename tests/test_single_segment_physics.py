"""
Single-segment physics beyond the steady-state and uplift cases: irregular
grids, distributed source/sink, base-level fall, and Sternberg gravel loss.
"""

import numpy as np
import pytest

import grlp


def _base(intermittency=1.0, x_ext=None, dx=1000.0, nx=90, x0=10000.0, S0=0.015):
    lp = grlp.LongProfile()
    lp.set_intermittency(intermittency)
    lp.basic_constants()
    lp.bedload_lumped_constants()
    lp.set_hydrologic_constants()
    if x_ext is not None:
        lp.set_x(x_ext=x_ext)
    else:
        lp.set_x(dx=dx, nx=nx, x0=x0)
    lp.set_z(S0=-S0, z1=0.0)
    lp.set_A(k_xA=1.0)
    lp.set_Q(k_xQ=1.43e-5, P_xQ=7 / 4.0 * 0.7)
    lp.P_xQ = 7 / 4.0 * 0.7
    lp.set_B(k_xB=25.0, P_xB=0.2)
    lp.set_niter(3)
    lp.set_z_bl(0.0)
    Qs0 = lp.k_Qs * lp.Q[0] * S0 ** (7 / 6.0)
    lp.set_Qs_input_upstream(Qs0)
    lp.set_uplift_rate(0.0)
    return lp, Qs0


# --------------------------------------------------------------------------- #
# Irregular (non-uniform dx) grid
# --------------------------------------------------------------------------- #

def test_irregular_grid_matches_analytical():
    # A deterministic, strongly non-uniform grid must still reach the same
    # analytical power-law equilibrium (the solution is a function of x, not of
    # the spacing).
    rng = np.random.default_rng(0)
    spacings = 500.0 + 1000.0 * rng.random(91)
    x_ext = np.concatenate(([8000.0], 8000.0 + np.cumsum(spacings)))
    lp, _ = _base(x_ext=x_ext)
    lp.evolve_threshold_width_river(nt=100, dt=1e13)
    z_analytical = lp.analytical_threshold_width()
    relief = lp.z.max() - lp.z.min()
    assert np.abs(lp.z - z_analytical).max() / relief < 1e-3


def test_irregular_grid_is_actually_irregular():
    # Guard the test above: confirm the spacing really is non-uniform.
    rng = np.random.default_rng(0)
    spacings = 500.0 + 1000.0 * rng.random(91)
    assert spacings.max() / spacings.min() > 2.0


# --------------------------------------------------------------------------- #
# Distributed source/sink equals uplift
# --------------------------------------------------------------------------- #

def test_distributed_source_equals_uplift():
    # Both uplift U and a distributed source ssd enter as dz/dt += (rate), so a
    # uniform source of magnitude R must produce exactly the profile that uplift
    # at rate R does.
    R = 5e-11
    lp_U, _ = _base()
    lp_U.set_uplift_rate(R)
    lp_U.evolve_threshold_width_river(nt=200, dt=1e12)

    lp_S, _ = _base()
    lp_S.set_uplift_rate(0.0)
    lp_S.set_source_sink_distributed(R)
    lp_S.evolve_threshold_width_river(nt=200, dt=1e12)

    np.testing.assert_allclose(lp_S.z, lp_U.z, atol=1e-9)


def test_distributed_sink_lowers_profile():
    # A negative source (sink) should lower the bed relative to no forcing.
    lp0, _ = _base()
    lp0.evolve_threshold_width_river(nt=200, dt=1e12)
    lpS, _ = _base()
    lpS.set_source_sink_distributed(-5e-11)
    lpS.evolve_threshold_width_river(nt=200, dt=1e12)
    assert lpS.z.mean() < lp0.z.mean()


# --------------------------------------------------------------------------- #
# Base-level fall
# --------------------------------------------------------------------------- #

def test_base_level_fall_reequilibrates():
    lp, _ = _base()
    lp.evolve_threshold_width_river(nt=100, dt=1e13)
    z_before = lp.z.copy()

    lp.set_z_bl(-50.0)
    lp.evolve_threshold_width_river(nt=200, dt=1e13)

    # New steady state matches the analytical power law again ...
    z_analytical = lp.analytical_threshold_width()
    relief = lp.z.max() - lp.z.min()
    assert np.abs(lp.z - z_analytical).max() / relief < 1e-3
    assert np.abs(lp.dz_dt).max() < 1e-15
    # ... and the whole profile has dropped by ~the base-level change.
    assert (z_before.mean() - lp.z.mean()) == pytest.approx(50.0, abs=1.0)


# --------------------------------------------------------------------------- #
# Sternberg gravel loss
# --------------------------------------------------------------------------- #

def test_sternberg_gravel_loss_exponential_decay():
    # A per-km abrasion coefficient makes bedload discharge decay downstream
    # following Sternberg's law, Qs = Qs0 * exp(-k * x).
    k_per_km = 0.02
    lp, Qs0 = _base(intermittency=1.0)
    lp.set_Sternberg_gravel_loss(gravel_fractional_loss_per_km=k_per_km)
    lp.evolve_threshold_width_river(nt=300, dt=1e13)
    lp.compute_Q_s()
    assert not np.isnan(lp.z).any()
    x_km = (lp.x - lp.x_ext[0]) / 1000.0
    predicted = Qs0 * np.exp(-k_per_km * x_km)
    assert np.abs(lp.Q_s - predicted).max() / Qs0 < 0.05
