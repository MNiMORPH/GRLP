"""
Does adaptive time stepping run faster than fixed-step BDF2 for GRLP?

This benchmark answers that quantitatively and reproducibly. It integrates a
representative transient two ways -- fixed-step BDF2 across a range of step
counts, and adaptive (step-doubling) BDF2 across a range of tolerances -- and
measures, for each run, the achieved accuracy against a fine reference and the
work done. Work is counted as the number of sparse **linear solves**
(``spsolve`` calls): a hardware-independent metric that captures every
Picard iteration, every step-doubling sub-step (one Δt step plus two Δt/2
steps), and every rejected step. Wall-clock time is reported alongside for
realism. Both schemes use the same converged Picard tolerance, so the
comparison isolates the *stepping scheme*, not the nonlinear solve.

Result (see ``main`` output): **adaptive stepping does not speed GRLP up.** At
matched accuracy it costs several times more linear solves than a well-chosen
fixed BDF2 step, for both a smooth (uplift-step) and a sharper (base-level-drop)
transient. GRLP's transients are diffusive and therefore smooth in time, so
uniform steps are already near the optimal placement, and step doubling's ~3x
per-step overhead is never recovered. Adaptive stepping's value for GRLP is
convenience (choose a tolerance instead of guessing a step) and automatic
resolution of unfamiliar timescales -- not speed.

Reproducible: deterministic, seeded only through the fixed canonical setup in
``tests/conftest.make_long_profile``. Run from the repository root::

    python benchmarks/adaptive_timestep_benchmark.py

Guarded (a small version) by ``tests/test_adaptive_benchmark.py``.
"""
import os
import sys
import time

import numpy as np

_here = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(_here, "..", "tests"))
from conftest import make_long_profile, STEADY_NT, STEADY_DT   # noqa: E402
import grlp.solver as solver                                    # noqa: E402

_YEAR = 3.15e7
_UPLIFT = 1.0e-4 / _YEAR             # 0.1 mm/yr uplift step (smooth transient)
_BL_DROP = -40.0                     # base-level drop [m] (sharper transient)
_PICARD_TOL = 1.0e-10                # both schemes: same converged nonlinear solve


# --- linear-solve counter (wraps solver.spsolve for the duration of a run) ---
class _SolveCounter:
    def __enter__(self):
        self._orig = solver.spsolve
        self.n = 0

        def counting(*args, **kwargs):
            self.n += 1
            return self._orig(*args, **kwargs)

        solver.spsolve = counting
        return self

    def __exit__(self, *exc):
        solver.spsolve = self._orig


def _steady(scenario):
    """Canonical profile at steady state, perturbed and set up for BDF2."""
    lp = make_long_profile(U=0.0)
    lp.evolve_threshold_width_river(nt=STEADY_NT, dt=STEADY_DT)
    if scenario == "uplift":
        lp.set_uplift_rate(_UPLIFT)
    elif scenario == "baselevel":
        lp.set_z_bl(_BL_DROP)
    else:
        raise ValueError("scenario must be 'uplift' or 'baselevel'")
    lp.set_time_integration(2)
    lp.set_iteration_tolerance(_PICARD_TOL, max_iter=50)
    return lp


def _equilibration_time():
    lp = make_long_profile()
    lp.evolve_threshold_width_river(nt=STEADY_NT, dt=STEADY_DT)
    lp.compute_equilibration_time()
    return lp.equilibration_time


def _reference(scenario, T, n=8192):
    ref = _steady(scenario)
    ref.evolve_threshold_width_river(nt=n, dt=T / n)
    return ref.z.copy()


def measure(scenario, T, z_ref, fixed_nts, adaptive_tols):
    """Return dict with, for fixed and adaptive: achieved error, linear-solve
    count, and wall time for each configuration."""
    fixed = {"err": [], "solves": [], "time": []}
    for nt in fixed_nts:
        lp = _steady(scenario)
        with _SolveCounter() as c:
            t0 = time.time()
            lp.evolve_threshold_width_river(nt=nt, dt=T / nt)
            wall = time.time() - t0
        fixed["err"].append(float(np.max(np.abs(lp.z - z_ref))))
        fixed["solves"].append(c.n)
        fixed["time"].append(wall)

    adaptive = {"err": [], "solves": [], "time": []}
    for tol in adaptive_tols:
        lp = _steady(scenario)
        lp.set_adaptive_timestep(tol)
        with _SolveCounter() as c:
            t0 = time.time()
            lp.evolve_threshold_width_river_adaptive(T)
            wall = time.time() - t0
        adaptive["err"].append(float(np.max(np.abs(lp.z - z_ref))))
        adaptive["solves"].append(c.n)
        adaptive["time"].append(wall)

    for d in (fixed, adaptive):
        for k in d:
            d[k] = np.array(d[k], float)
    return {"fixed": fixed, "adaptive": adaptive}


def solves_at_accuracy(err, solves, target_err):
    """Linear solves needed to reach ``target_err``, by log-log interpolation of
    the (error, solves) curve (error decreasing as solves increase)."""
    le, ls = np.log(err), np.log(solves)
    order = np.argsort(le)
    return float(np.exp(np.interp(np.log(target_err), le[order], ls[order])))


def cost_ratio(result, target_err):
    """Adaptive / fixed linear-solve cost to reach ``target_err`` (>1 = adaptive
    is more expensive)."""
    a = solves_at_accuracy(result["adaptive"]["err"], result["adaptive"]["solves"],
                           target_err)
    f = solves_at_accuracy(result["fixed"]["err"], result["fixed"]["solves"],
                           target_err)
    return a / f, a, f


def _print_table(name, result):
    print("\n===== %s =====" % name)
    print(" fixed BDF2:     nt      err[m]     solves   time[s]")
    for i in range(len(result["fixed"]["err"])):
        print("             %6d   %.3e   %6d   %.3f"
              % (result["fixed"]["_nt"][i], result["fixed"]["err"][i],
                 result["fixed"]["solves"][i], result["fixed"]["time"][i]))
    print(" adaptive:      tol      err[m]     solves   time[s]")
    for i in range(len(result["adaptive"]["err"])):
        print("            %.0e   %.3e   %6d   %.3f"
              % (result["adaptive"]["_tol"][i], result["adaptive"]["err"][i],
                 result["adaptive"]["solves"][i], result["adaptive"]["time"][i]))


def main():
    T = _equilibration_time()
    fixed_nts = [16, 32, 64, 128, 256]
    adaptive_tols = [1e-1, 1e-2, 1e-3, 1e-4]
    for scenario, label in [("uplift", "uplift step (smooth transient)"),
                            ("baselevel", "base-level drop (sharper transient)")]:
        z_ref = _reference(scenario, T)
        result = measure(scenario, T, z_ref, fixed_nts, adaptive_tols)
        result["fixed"]["_nt"] = fixed_nts
        result["adaptive"]["_tol"] = adaptive_tols
        _print_table(label, result)
        # Compare at an accuracy both schemes bracket.
        target = 2.0e-3
        ratio, a, f = cost_ratio(result, target)
        print("  -> to reach %.0e m: adaptive %.0f vs fixed %.0f linear solves "
              "(adaptive is %.1fx the cost)" % (target, a, f, ratio))
    print("\nVerdict: adaptive stepping costs several times more linear solves "
          "than fixed BDF2\nat matched accuracy for GRLP's smooth, diffusive "
          "transients. Its value is\nconvenience (a tolerance instead of a step "
          "guess), not speed.")


if __name__ == "__main__":
    main()
