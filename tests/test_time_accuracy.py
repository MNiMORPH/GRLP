"""
Time-stepping accuracy of the threshold-width solver.

GRLP integrates the sediment-conservation equation with a semi-implicit scheme.
Its **time** accuracy is selectable (`set_time_integration`): backward Euler is
first-order (the default), BDF2 is second-order. This test verifies both orders
by self-convergence: it evolves the *same* transient to a fixed end time first
with a very small "reference" time step, then with a sequence of progressively
larger steps, and checks that the error against the reference

  (i) decreases monotonically as Δt shrinks, and
  (ii) does so at the expected order -- backward-Euler error halves each time Δt
       is halved (∝ Δt); BDF2 error quarters (∝ Δt²).

The transient is a step increase in uplift applied to a profile that starts at
steady state; the end time is one equilibration time (Paola et al., 1992), so the
valley is mid-adjustment and the time step genuinely matters (at steady state the
result is independent of Δt and the test would be vacuous).

See `docs/accuracy.md` and `claude-instructions/second-order-time-bdf2-design.md`.
The second-order case is the acceptance test for MNiMORPH/GRLP#16.
"""

import numpy as np

from conftest import make_long_profile, STEADY_NT, STEADY_DT


_YEAR = 3.15e7
_UPLIFT = 1.0e-4 / _YEAR            # 0.1 mm/yr step, in m/s -- drives the transient
_N_REF = 2048                       # reference (near-exact) number of steps
_N_STEPS = (8, 16, 32, 64, 128)     # test resolutions: Δt halves down the list


def _transient_final_z(T, n_steps, order=1, niter=3):
    """Bed elevation after evolving a steady-state profile through a step in
    uplift for total time ``T``, in ``n_steps`` equal steps, integrated with the
    given time-integration ``order`` (1 = backward Euler, 2 = BDF2)."""
    lp = make_long_profile(U=0.0, niter=niter)
    lp.evolve_threshold_width_river(nt=STEADY_NT, dt=STEADY_DT)    # to steady state
    lp.set_uplift_rate(_UPLIFT)
    lp.set_time_integration(order)
    lp.evolve_threshold_width_river(nt=n_steps, dt=T / n_steps)    # the transient
    return lp.z.copy()


def _equilibration_time():
    lp = make_long_profile()
    lp.evolve_threshold_width_river(nt=STEADY_NT, dt=STEADY_DT)
    lp.compute_equilibration_time()
    return lp.equilibration_time


def _convergence(order):
    """Self-convergence sweep at a given order; returns (dts, errs, table)."""
    T = _equilibration_time()
    z_ref = _transient_final_z(T, _N_REF, order=order)
    dts = np.array([T / n for n in _N_STEPS])
    errs = np.array([np.max(np.abs(_transient_final_z(T, n, order=order) - z_ref))
                     for n in _N_STEPS])
    table = "\n".join(
        "  n=%4d  dt=%9.1f yr  max_err=%.6f m" % (n, dt / _YEAR, e)
        for n, dt, e in zip(_N_STEPS, dts, errs))
    return dts, errs, table


def _rate(dts, errs):
    return np.polyfit(np.log(dts), np.log(errs), 1)[0]


def test_time_stepping_is_first_order():
    """Backward Euler (default): error falls at first order in Δt."""
    dts, errs, table = _convergence(order=1)
    assert np.all(np.diff(errs) < 0.0), "error not monotonic in Δt:\n" + table
    rate = _rate(dts, errs)
    assert 0.8 < rate < 1.3, (
        "expected first-order in time (rate ~1), measured %.3f:\n%s"
        % (rate, table))


def test_bdf2_is_second_order():
    """BDF2 (``set_time_integration(2)``): error falls at second order in Δt --
    the acceptance test for the second-order-in-time work (GRLP#16)."""
    dts, errs, table = _convergence(order=2)
    assert np.all(np.diff(errs) < 0.0), "error not monotonic in Δt:\n" + table
    rate = _rate(dts, errs)
    assert 1.7 < rate < 2.4, (
        "expected second-order in time (rate ~2), measured %.3f:\n%s"
        % (rate, table))
