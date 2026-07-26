"""
Adaptive time stepping (MNiMORPH/GRLP#16, Phase 2), built up in increments.

The adaptive controller needs a stepping loop that carries the BDF2 two-level
history ``(z^n, z^{n-1})`` across steps -- unlike ``solver.evolve``, which
re-bootstraps that history on the first step of every call. ``solver.evolve_adaptive``
is that dedicated loop.

**Increment 1 (fixed step, no control yet).** These tests pin the scaffold: at a
fixed step with BDF2 requested, ``evolve_adaptive`` must reproduce ``evolve``
bit-for-bit and inherit its second-order-in-time convergence. Later increments
(embedded error estimate, step controller, variable-step weights) add tolerance
control and will bring their own acceptance tests; this equivalence guards the
history bookkeeping underneath them.

See ``claude-instructions/adaptive-timestepping-design.md``.
"""

import numpy as np

from conftest import make_long_profile, STEADY_NT, STEADY_DT
from grlp import solver


_YEAR = 3.15e7
_UPLIFT = 1.0e-4 / _YEAR
_N_REF = 2048
_N_STEPS = (8, 16, 32, 64, 128)


def _steady_bdf2():
    """A steady-state profile with a uplift step applied and BDF2 selected."""
    lp = make_long_profile(U=0.0)
    lp.evolve_threshold_width_river(nt=STEADY_NT, dt=STEADY_DT)
    lp.set_uplift_rate(_UPLIFT)
    lp.set_time_integration(2)
    return lp


def _final_z(loop, nt, dt):
    lp = _steady_bdf2()
    loop(lp._net(), nt, dt)
    return lp.z.copy()


def test_evolve_adaptive_fixed_step_matches_evolve_bitwise():
    """At a fixed step, the dedicated adaptive loop reproduces ``evolve``'s BDF2
    path exactly -- same bootstrap, same history, same solves."""
    for nt, yrs in [(1, 1.0e6), (5, 2.0e5), (20, 5.0e4)]:
        dt = yrs * _YEAR
        z_evolve = _final_z(solver.evolve, nt, dt)
        z_adaptive = _final_z(solver.evolve_adaptive, nt, dt)
        assert np.array_equal(z_evolve, z_adaptive), (
            "evolve_adaptive diverged from evolve at nt=%d, dt=%g yr "
            "(max|dz|=%.3e m)" % (nt, yrs, np.max(np.abs(z_evolve - z_adaptive))))


def test_evolve_adaptive_fixed_step_is_second_order():
    """The scaffold inherits BDF2's second-order-in-time convergence at fixed
    step (self-convergence rate ~2)."""
    lp = make_long_profile()
    lp.evolve_threshold_width_river(nt=STEADY_NT, dt=STEADY_DT)
    lp.compute_equilibration_time()
    T = lp.equilibration_time
    z_ref = _final_z(solver.evolve_adaptive, _N_REF, T / _N_REF)
    dts = np.array([T / n for n in _N_STEPS])
    errs = np.array([np.max(np.abs(_final_z(solver.evolve_adaptive, n, T / n) - z_ref))
                     for n in _N_STEPS])
    table = "\n".join("  n=%4d  err=%.3e m" % (n, e)
                      for n, e in zip(_N_STEPS, errs))
    assert np.all(np.diff(errs) < 0.0), "error not monotonic in Δt:\n" + table
    rate = np.polyfit(np.log(dts), np.log(errs), 1)[0]
    assert 1.7 < rate < 2.4, (
        "expected second-order (rate ~2), measured %.3f:\n%s" % (rate, table))


def _evolve_sequence(net, dts):
    """Step a network through a *prescribed* list of step sizes, carrying the
    BDF2 history and setting the variable-step weight ``omega = dt_n / dt_{n-1}``
    each step. Drives the real solver primitives (``_picard_step`` reads the
    variable-step weights assembled in ``solver.assemble``); it is the fixed-list
    stand-in for the step controller, used to test the weights on their own."""
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
    lp = make_long_profile()
    lp.evolve_threshold_width_river(nt=STEADY_NT, dt=STEADY_DT)
    lp.compute_equilibration_time()
    T = lp.equilibration_time

    def final_z(N):
        dts = _nonuniform_partition(T, N)
        prof = _steady_bdf2()
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
