"""
Time-stepping accuracy of the threshold-width solver.

GRLP integrates the sediment-conservation equation with a semi-implicit scheme
that is **first-order in time** (and second-order in space), so the temporal
truncation error scales like O(Δt).  This test verifies that scaling by
self-convergence: it evolves the *same* transient to a fixed end time first with
a very small "reference" time step, then with a sequence of progressively larger
steps, and checks that the error against the reference

  (i) decreases monotonically as Δt shrinks, and
  (ii) does so at first order -- the error roughly halves each time Δt is halved.

The transient is a step increase in uplift applied to a profile that starts at
steady state; the end time is one equilibration time (Paola et al., 1992), so the
valley is mid-adjustment and the time step genuinely matters (at steady state the
result is independent of Δt and the test would be vacuous).

Besides guarding the current behaviour, this is the acceptance-test scaffold for
the solver-accuracy enhancement (MNiMORPH/GRLP#16): a second-order-in-time scheme
would raise the measured convergence rate from ~1 toward ~2.
"""

import numpy as np

from conftest import make_long_profile, STEADY_NT, STEADY_DT


_YEAR = 3.15e7
_UPLIFT = 1.0e-4 / _YEAR            # 0.1 mm/yr step, in m/s -- drives the transient
_N_REF = 2048                       # reference (near-exact) number of steps
_N_STEPS = (8, 16, 32, 64, 128)     # test resolutions: Δt halves down the list


def _transient_final_z(T, n_steps, niter=3):
    """Bed elevation after evolving a steady-state profile through a step in
    uplift for total time ``T``, taken in ``n_steps`` equal time steps."""
    lp = make_long_profile(U=0.0, niter=niter)
    lp.evolve_threshold_width_river(nt=STEADY_NT, dt=STEADY_DT)    # to steady state
    lp.set_uplift_rate(_UPLIFT)
    lp.evolve_threshold_width_river(nt=n_steps, dt=T / n_steps)    # the transient
    return lp.z.copy()


def _equilibration_time():
    lp = make_long_profile()
    lp.evolve_threshold_width_river(nt=STEADY_NT, dt=STEADY_DT)
    lp.compute_equilibration_time()
    return lp.equilibration_time


def test_time_stepping_is_first_order():
    """Error against a fine-Δt reference falls at first order in Δt."""
    T = _equilibration_time()
    z_ref = _transient_final_z(T, _N_REF)

    dts = np.array([T / n for n in _N_STEPS])
    errs = np.array([np.max(np.abs(_transient_final_z(T, n) - z_ref))
                     for n in _N_STEPS])

    # Readable diagnostic if an assertion below fails.
    table = "\n".join(
        "  n=%4d  dt=%9.1f yr  max_err=%.6f m" % (n, dt / _YEAR, e)
        for n, dt, e in zip(_N_STEPS, dts, errs))

    # (i) Refining the step must reduce the error.
    assert np.all(np.diff(errs) < 0.0), "error not monotonic in Δt:\n" + table

    # (ii) First order: the log-log slope of error vs Δt is ~1 -- distinctly not
    # the ~2 that a second-order-in-time scheme would give.
    rate = np.polyfit(np.log(dts), np.log(errs), 1)[0]
    assert 0.8 < rate < 1.3, (
        "expected first-order in time (rate ~1), measured %.3f:\n%s"
        % (rate, table))
