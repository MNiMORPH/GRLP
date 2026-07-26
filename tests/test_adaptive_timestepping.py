"""
Adaptive time stepping (MNiMORPH/GRLP#16, Phase 2).

``solver.evolve_adaptive`` advances the network by a total time ``T`` with
variable-step BDF2, sizing each step from a step-doubling error estimate so the
per-step error stays at or below the tolerance set by ``set_adaptive_timestep``.
It carries the BDF2 two-level history across steps (unlike ``solver.evolve``,
which re-bootstraps it each call) and advances with the finer two-half-step
solution (local extrapolation).

These tests check the two pieces:
 * the **variable-step BDF2 weights** stay second order on a non-uniform step
   sequence (the foundation the controller stands on); and
 * the **controller** actually controls -- achieved error falls as the tolerance
   tightens, is independent of the starting-step guess, and lands exactly on the
   requested end time.

See ``claude-instructions/adaptive-timestepping-design.md``.
"""

import warnings

import numpy as np

from conftest import make_long_profile, STEADY_NT, STEADY_DT
from grlp import solver


_YEAR = 3.15e7
_UPLIFT = 1.0e-4 / _YEAR


def _steady_bdf2(picard_tol=1e-10):
    """A steady-state profile with an uplift step applied, BDF2 selected, and the
    Picard solve converged (so the error estimate is not polluted by an
    under-converged nonlinear step)."""
    lp = make_long_profile(U=0.0)
    lp.evolve_threshold_width_river(nt=STEADY_NT, dt=STEADY_DT)
    lp.set_uplift_rate(_UPLIFT)
    lp.set_time_integration(2)
    lp.set_picard_tolerance(picard_tol, max_iter=50)
    return lp


def _equilibration_time():
    lp = make_long_profile()
    lp.evolve_threshold_width_river(nt=STEADY_NT, dt=STEADY_DT)
    lp.compute_equilibration_time()
    return lp.equilibration_time


def _fine_reference(T, n=4096):
    """The 'truth': the same transient integrated with many small fixed BDF2
    steps."""
    ref = _steady_bdf2()
    ref.evolve_threshold_width_river(nt=n, dt=T / n)
    return ref.z.copy()


# ----------------------------------------------------------------------------
# Variable-step BDF2 weights (the controller's foundation)
# ----------------------------------------------------------------------------
def _evolve_sequence(net, dts):
    """Step a network through a *prescribed* list of step sizes, carrying the
    BDF2 history and setting the variable-step weight ``omega = dt_n / dt_{n-1}``
    each step. Drives the real solver primitives (``_picard_step`` reads the
    variable-step weights assembled in ``solver.assemble``); the fixed-list
    stand-in for the controller, used to test the weights on their own."""
    segs = net.segments
    lengths = list(net.list_of_segment_lengths)
    starts = np.cumsum([0] + lengths)[:-1]
    picard_tol, max_iter, gravel_loss_active = solver._picard_config(net)
    bdf2 = getattr(net, "time_order", 1) == 2
    dt_prev = None
    for ti, dt in enumerate(dts):
        for seg in segs:
            seg.zold2 = None if ti == 0 else seg.zold
            seg.zold = seg.z.copy()
        net._bdf2_step = bdf2 and ti > 0
        net._bdf2_omega = 1.0 if dt_prev is None else dt / dt_prev
        solver._picard_step(net, dt, segs, starts, lengths,
                            gravel_loss_active, picard_tol, max_iter)
        net.t += dt
        dt_prev = dt


def _nonuniform_partition(T, N, seed=0):
    """``N`` positive step sizes summing to ``T`` and varying by up to ~2x."""
    rng = np.random.default_rng(seed)
    w = 1.0 + rng.uniform(0.0, 1.0, N)     # ratios in [1, 2)
    return w * (T / w.sum())


def test_variable_step_bdf2_is_second_order():
    """BDF2 stays second order on a genuinely *non-uniform* step sequence, in
    max(Δt) -- the variable-step weights (omega-dependent 3-level combination)
    are what preserve the order; with the uniform (3/2, 2, 1/2) weights the rate
    collapses toward one.  Prerequisite for step-doubling error control and for
    the step controller, which change Δt every step."""
    T = _equilibration_time()

    def final_z(N):
        dts = _nonuniform_partition(T, N)
        prof = _steady_bdf2(picard_tol=None)
        _evolve_sequence(prof._net(), dts)
        return prof.z.copy(), dts.max()

    z_ref, _ = final_z(4096)
    Ns = (16, 32, 64, 128, 256)
    maxdts, errs = [], []
    for N in Ns:
        z, mx = final_z(N)
        maxdts.append(mx)
        errs.append(np.max(np.abs(z - z_ref)))
    maxdts = np.array(maxdts)
    errs = np.array(errs)
    table = "\n".join("  N=%4d  max_dt=%9.1f yr  err=%.3e m" % (n, m / _YEAR, e)
                      for n, m, e in zip(Ns, maxdts, errs))
    assert np.all(np.diff(errs) < 0.0), "error not monotonic in max Δt:\n" + table
    rate = np.polyfit(np.log(maxdts), np.log(errs), 1)[0]
    assert 1.7 < rate < 2.4, (
        "variable-step BDF2 should stay ~2nd order, measured %.3f:\n%s"
        % (rate, table))


# ----------------------------------------------------------------------------
# The step controller
# ----------------------------------------------------------------------------
def _adaptive_final_z(tol, T, dt_init=None):
    lp = _steady_bdf2()
    lp.set_adaptive_timestep(tol)
    lp.evolve_threshold_width_river_adaptive(T, dt_init=dt_init)
    return lp


def test_adaptive_error_tracks_tolerance():
    """Tightening the tolerance monotonically reduces the achieved error against
    a fine fixed-step reference, and each run lands within a modest band of its
    tolerance (tol is a per-step local tolerance, so the path error is of
    comparable, not identical, magnitude)."""
    T = _equilibration_time()
    z_ref = _fine_reference(T)
    tols = [1e-2, 1e-3, 1e-4]
    errs = [np.max(np.abs(_adaptive_final_z(tol, T).z - z_ref)) for tol in tols]
    table = "\n".join("  tol=%.0e  err=%.3e m" % (t, e)
                      for t, e in zip(tols, errs))
    assert np.all(np.diff(errs) < 0.0), \
        "achieved error not monotonic in tol:\n" + table
    for tol, err in zip(tols, errs):
        assert 0.05 * tol < err < 20.0 * tol, \
            "achieved error out of band for tol=%.0e:\n%s" % (tol, table)


def test_adaptive_is_starting_step_independent():
    """The controller forgets its initial step guess: the same tolerance gives
    the same profile whether the first step is tiny or large."""
    T = _equilibration_time()
    z_small = _adaptive_final_z(1e-3, T, dt_init=T * 1e-4).z
    z_large = _adaptive_final_z(1e-3, T, dt_init=T * 0.1).z
    spread = np.max(np.abs(z_small - z_large))
    assert spread < 1e-4, \
        "adaptive result depends on the starting step (spread %.3e m)" % spread


def test_adaptive_reaches_target_time_exactly():
    """Adaptive stepping advances by exactly the requested total time (the last
    step is trimmed to land on it)."""
    lp = _steady_bdf2()
    t0 = lp.t
    T = _equilibration_time()
    lp.set_adaptive_timestep(1e-3)
    lp.evolve_threshold_width_river_adaptive(T)
    assert np.isclose(lp.t, t0 + T, rtol=1e-9, atol=0.0), \
        "final time %g != t0 + T = %g" % (lp.t, t0 + T)


def test_adaptive_step_grows_as_transient_settles():
    """Steps should be small early (fast adjustment right after the uplift step)
    and larger late (the profile is settling), so the mean late step exceeds the
    mean early step. Measured by counting error-estimate trials in each half."""
    T = _equilibration_time()
    calls = {"t": []}
    orig = solver._trial_step

    def counting(net, *a, **k):
        calls["t"].append(net.t)
        return orig(net, *a, **k)

    solver._trial_step = counting
    try:
        lp = _steady_bdf2()
        t0 = lp.t
        lp.set_adaptive_timestep(1e-3)
        lp.evolve_threshold_width_river_adaptive(T)
    finally:
        solver._trial_step = orig
    ts = np.array(calls["t"]) - t0
    early = np.sum(ts < T / 2)
    late = np.sum(ts >= T / 2)
    # More trials are spent in the fast early half than the settling late half.
    assert early > late, \
        "expected denser stepping early (fast transient): early=%d late=%d" \
        % (early, late)


def test_adaptive_requires_a_tolerance():
    """``evolve_adaptive`` without a configured tolerance is an error, not a
    silent no-op."""
    lp = _steady_bdf2()
    try:
        lp.evolve_threshold_width_river_adaptive(_equilibration_time())
    except ValueError:
        pass
    else:
        raise AssertionError("expected ValueError without set_adaptive_timestep")
