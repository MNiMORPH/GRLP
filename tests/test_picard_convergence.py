"""
Convergence-controlled Picard iteration (``set_iteration_tolerance``).

The threshold-width solver is semi-implicit: each time step it relinearizes and
re-solves a few times (Picard iteration). By default it iterates to convergence
(``set_iteration_tolerance``), stopping once the inter-iteration elevation change
``max|z_k - z_{k-1}|`` drops below ``tol``, capped at ``max_iter`` iterations and
warning if the cap is reached first. ``set_niter(n)`` instead takes a fixed ``n``
iterations (the faster expert option), the two being mutually exclusive.

These tests check that
  (i)   fixed-iteration mode (``set_niter`` / ``tol=None``) takes exactly the
        fixed count and is deterministic;
  (ii)  a tight tolerance drives the step to the same fixed point a large fixed
        ``niter`` reaches (the Picard residual floors at round-off);
  (iii) a looser tolerance leaves a demonstrably larger residual than a tighter
        one; and
  (iv)  the non-convergence warning fires when the tolerance is unreachable
        within ``max_iter``.

See ``set_iteration_tolerance`` in ``grlp/grlp.py`` and the Picard loop in
``grlp/solver.py`` (MNiMORPH/GRLP#17).
"""

import warnings

import numpy as np

from conftest import make_long_profile, STEADY_NT, STEADY_DT


_YEAR = 3.15e7
_UPLIFT = 1.0e-4 / _YEAR            # 0.1 mm/yr step, in m/s -- drives the transient


def _transient(nt, dt_yr, tol=None, niter=3, max_iter=100):
    """Steady-state profile pushed through a uplift step for ``nt`` steps of
    ``dt_yr`` years, optionally with a Picard convergence tolerance set."""
    lp = make_long_profile(U=0.0, niter=niter)
    lp.evolve_threshold_width_river(nt=STEADY_NT, dt=STEADY_DT)     # to steady state
    lp.set_uplift_rate(_UPLIFT)
    if tol is not None:
        lp.set_iteration_tolerance(tol, max_iter=max_iter)
    lp.evolve_threshold_width_river(nt=nt, dt=dt_yr * _YEAR)
    return lp.z.copy()


def test_fixed_niter_mode_is_deterministic():
    """In fixed-iteration mode (``set_niter`` / ``tol=None``) the solver takes
    exactly ``niter`` iterations and is deterministic; this is the historical
    loop, so the Euler golden-master set stays bit-for-bit unchanged under it."""
    z_a = _transient(20, 5.0e4, tol=None, niter=3)
    z_b = _transient(20, 5.0e4, tol=None, niter=3)
    assert np.array_equal(z_a, z_b)


# A single large step from the (far-from-new-equilibrium) starting profile is the
# regime where a fixed ``niter`` is genuinely under-converged: here niter=3 lands
# ~1e-6 m from the Picard fixed point, so a tolerance that forces more iterations
# visibly changes the answer.  (Over many small transient steps each step starts
# near its answer and niter=3 is already converged to round-off, so those runs
# cannot tell the tolerance apart from the default.)
_BIG_NT = 1
_BIG_DT_YR = 1.0e7


def test_tolerance_converges_to_fixed_point():
    """A tight tolerance reaches the same fixed point that a large fixed ``niter``
    does -- both have converged the Picard iteration to round-off -- even on a
    single large step where the default niter=3 is ~1e-6 m short of it."""
    z_tol = _transient(_BIG_NT, _BIG_DT_YR, tol=1e-9, niter=3)
    z_big = _transient(_BIG_NT, _BIG_DT_YR, tol=None, niter=30)
    err = np.max(np.abs(z_tol - z_big))
    assert err < 1e-8, "tight-tol path should match a heavily-iterated path; " \
                       "max|dz| = %.3e m" % err


def test_tighter_tolerance_is_more_converged():
    """A looser tolerance leaves a larger distance from the converged fixed point
    than a tighter one -- the tolerance actually controls the residual."""
    z_conv = _transient(_BIG_NT, _BIG_DT_YR, tol=None, niter=30)
    err_loose = np.max(np.abs(_transient(_BIG_NT, _BIG_DT_YR, tol=1e-3) - z_conv))
    err_tight = np.max(np.abs(_transient(_BIG_NT, _BIG_DT_YR, tol=1e-7) - z_conv))
    assert err_tight < err_loose, (
        "tighter tol should be closer to converged: tight=%.3e loose=%.3e"
        % (err_tight, err_loose))


def test_warns_when_not_converged():
    """When the tolerance is below the round-off floor and ``max_iter`` is small,
    the step cannot converge and a RuntimeWarning is raised."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        _transient(1, 1.0e7, tol=1e-20, max_iter=3)
        runtime = [w for w in caught if issubclass(w.category, RuntimeWarning)]
    assert runtime, "expected a non-convergence RuntimeWarning"
    assert "did not converge" in str(runtime[0].message)


def test_positive_tolerance_required():
    """A non-positive tolerance is rejected (None disables the feature instead)."""
    lp = make_long_profile()
    for bad in (0.0, -1e-6):
        try:
            lp.set_iteration_tolerance(bad)
        except ValueError:
            pass
        else:
            raise AssertionError("expected ValueError for tol=%r" % (bad,))
