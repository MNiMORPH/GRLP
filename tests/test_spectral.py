"""
Limiting behavior of the linearized (spectral) transient solution.

McNab et al. (2023) linearize the threshold-width equation to obtain
near-analytical gain and lag between a periodic forcing (in sediment or water
supply) and the valley response.  These are implemented in
``compute_z_gain`` / ``compute_z_lag`` (elevation) and ``compute_Qs_gain`` /
``compute_Qs_lag`` (sediment discharge).

Rather than run full numerical time series (that is the slow benchmark, kept
separate), these tests check the analytical series solution against its known
asymptotic limits:

* Quasi-static (period >> equilibration time): elevation gain -> 6/7
  everywhere and lag -> 0; sediment-discharge gain -> 1.
* Fast forcing (period << equilibration time): elevation gain -> 0.

They also check structural properties that must hold at every period: gain
increases monotonically with period, is bounded in [0, 6/7], and the elevation
lag is a positive fraction (< 1/4) of the forcing period across the physically
useful range.

The model uses uniform discharge and width, for which the initial constant
slope is already the equilibrium; no time integration is needed to evaluate the
linearized response.
"""

import numpy as np
import pytest

import grlp


SIX_SEVENTHS = 6.0 / 7.0


def make_uniform_profile():
    """Uniform-discharge, uniform-width profile at its (linear) equilibrium."""
    lp = grlp.LongProfile()
    lp.basic_constants()
    lp.bedload_lumped_constants()
    lp.set_hydrologic_constants()
    L, dx, Qw, Qs0 = 100e3, 1e3, 10.0, 0.001
    lp.set_x(x_ext=np.arange(0.0, L + dx, dx))
    lp.set_Q(Q=Qw)
    lp.set_B(B=98.12)
    lp.set_z(S0=(Qs0 / (lp.k_Qs * Qw)) ** (6.0 / 7.0))
    lp.set_Qs_input_upstream(Qs0)
    lp.set_z_bl(0.0)
    lp.set_uplift_rate(0.0)
    lp.set_niter(3)
    lp.compute_equilibration_time()
    return lp


@pytest.fixture(scope="module")
def uniform_profile():
    return make_uniform_profile()


# --------------------------------------------------------------------------- #
# Quasi-static (long-period) limits
# --------------------------------------------------------------------------- #

def test_z_gain_long_period_limit(uniform_profile):
    lp = uniform_profile
    gain = lp.compute_z_gain(1e4 * lp.equilibration_time)
    # Uniform 6/7 everywhere in the quasi-static limit.
    np.testing.assert_allclose(gain, SIX_SEVENTHS, atol=1e-4)


def test_z_lag_long_period_limit(uniform_profile):
    lp = uniform_profile
    P = 1e4 * lp.equilibration_time
    lag_frac = lp.compute_z_lag(P) / P
    # Response is in phase with the forcing.
    assert np.abs(lag_frac).max() < 1e-3


def test_Qs_gain_long_period_limit(uniform_profile):
    lp = uniform_profile
    gain = lp.compute_Qs_gain(1e4 * lp.equilibration_time, A_Qs=0.2)
    # Sediment discharge tracks its supply amplitude when forced slowly.
    np.testing.assert_allclose(gain, 1.0, atol=1e-4)


# --------------------------------------------------------------------------- #
# Fast-forcing (short-period) limit
# --------------------------------------------------------------------------- #

def test_z_gain_short_period_limit(uniform_profile):
    lp = uniform_profile
    gain = lp.compute_z_gain(1e-4 * lp.equilibration_time)
    # The valley cannot respond to very fast forcing.
    assert gain.max() < 1e-2


# --------------------------------------------------------------------------- #
# Structural properties at all periods
# --------------------------------------------------------------------------- #

def test_z_gain_monotonic_in_period(uniform_profile):
    lp = uniform_profile
    periods = np.logspace(-2, 2, 9) * lp.equilibration_time
    # Track a fixed interior point.
    idx = len(lp.x) // 2
    gains = np.array([lp.compute_z_gain(P)[idx] for P in periods])
    assert np.all(np.diff(gains) > 0)


def test_z_gain_bounded(uniform_profile):
    lp = uniform_profile
    for P in np.logspace(-2, 3, 6) * lp.equilibration_time:
        gain = lp.compute_z_gain(P)
        assert gain.min() >= -1e-9
        assert gain.max() <= SIX_SEVENTHS + 1e-6


def test_z_lag_positive_and_subquarter(uniform_profile):
    lp = uniform_profile
    # Across the physically useful range the elevation response lags the
    # forcing by a positive fraction of the period, below the diffusive
    # quarter-period bound.
    for P in np.logspace(-2, 2, 9) * lp.equilibration_time:
        lag_frac = lp.compute_z_lag(P)[0] / P
        assert 0.0 < lag_frac < 0.25
