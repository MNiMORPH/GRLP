"""
Slow benchmark: numerical periodic forcing vs. the linearized gain and lag.

This is the direct, expensive validation behind the fast limit checks in
``test_spectral.py``.  It reproduces the essence of
``examples/McNab_et_al_GRL/Figure_S1_Benchmark.py`` at reduced resolution:
force the upstream sediment supply sinusoidally, integrate the full nonlinear
model over several periods until it settles into a limit cycle, and measure the
amplitude ratio (gain) and phase lag of the elevation response.  These are
compared with the near-analytical ``compute_z_gain`` and ``compute_z_lag``.

Marked ``slow`` (a few seconds per forcing period); run with ``-m slow`` or
deselect with ``-m "not slow"``.
"""

from copy import deepcopy

import numpy as np
import pytest

import grlp


# Forcing amplitude (fractional) applied to the sediment supply.
FORCING_AMPLITUDE = 0.2
# Periods of settling + measurement, and time steps per period.
N_PERIODS = 8
N_MEASURE_PERIODS = 4
N_STEPS_PER_PERIOD = 200
# Interior slice, avoiding both boundaries (the outlet, where L - x -> 0, makes
# the gain normalization singular).
INTERIOR = slice(5, 90)


def _base_profile():
    """Uniform profile evolved to its steady (limit-cycle mean) state."""
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
    lp.evolve_threshold_width_river(nt=2000, dt=3.15e9)
    lp.compute_equilibration_time()
    lp._Qs0 = Qs0
    lp._Qw = Qw
    lp._L = L
    return lp


def _force_and_record(base, period):
    """
    Sinusoidally force the sediment supply and record the elevation time series
    over the measurement window (the last N_MEASURE_PERIODS periods).
    """
    lp = deepcopy(base)
    dt = period / N_STEPS_PER_PERIOD
    nt = N_PERIODS * N_STEPS_PER_PERIOD
    t = np.arange(nt) * dt
    scale = 1.0 + FORCING_AMPLITUDE * np.sin(2.0 * np.pi * t / period)
    z = np.zeros((nt, len(lp.x)))
    for j, s in enumerate(scale):
        lp.set_Qs_input_upstream(base._Qs0 * s)
        lp.evolve_threshold_width_river(nt=1, dt=dt)
        z[j, :] = lp.z
    i0 = (N_PERIODS - N_MEASURE_PERIODS) * N_STEPS_PER_PERIOD
    return t[i0:], z[i0:, :]


def _numerical_gain(base, period, t_win, z_win):
    S0 = (base._Qs0 / (base.k_Qs * base._Qw)) ** (6.0 / 7.0)
    amplitude = (z_win.max(axis=0) - z_win.min(axis=0)) / 2.0
    return amplitude / FORCING_AMPLITUDE / S0 / (base._L - base.x)


def _numerical_lag(period, t_win, z_win):
    # Least-squares projection onto the forcing frequency: z ~ a*sin - b*cos,
    # i.e. z ~ A_z * sin(w t - phi) with phi the phase lag.
    w = 2.0 * np.pi / period
    zc = z_win - z_win.mean(axis=0)
    a = 2.0 * np.mean(zc * np.sin(w * t_win)[:, None], axis=0)
    b = 2.0 * np.mean(zc * np.cos(w * t_win)[:, None], axis=0)
    phi = np.mod(np.arctan2(-b, a), 2.0 * np.pi)
    return phi / (2.0 * np.pi) * period


@pytest.fixture(scope="module")
def base_profile():
    return _base_profile()


@pytest.mark.slow
@pytest.mark.parametrize("period_over_Teq", [1.0, 3.0])
def test_numerical_gain_matches_linearized(base_profile, period_over_Teq):
    period = period_over_Teq * base_profile.equilibration_time
    t_win, z_win = _force_and_record(base_profile, period)
    num_gain = _numerical_gain(base_profile, period, t_win, z_win)
    ana_gain = base_profile.compute_z_gain(period)
    rel = np.abs(num_gain[INTERIOR] - ana_gain[INTERIOR]) / ana_gain[INTERIOR]
    assert rel.max() < 0.05


@pytest.mark.slow
@pytest.mark.parametrize("period_over_Teq", [1.0, 3.0])
def test_numerical_lag_matches_linearized(base_profile, period_over_Teq):
    period = period_over_Teq * base_profile.equilibration_time
    t_win, z_win = _force_and_record(base_profile, period)
    num_lag = _numerical_lag(period, t_win, z_win)
    ana_lag = base_profile.compute_z_lag(period)
    # Agreement to a small fraction of the forcing period.
    max_diff_frac = np.abs(num_lag[INTERIOR] - ana_lag[INTERIOR]).max() / period
    assert max_diff_frac < 0.03
