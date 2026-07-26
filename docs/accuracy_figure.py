"""
Regenerate the time-stepping accuracy figure used in ``docs/accuracy.md``.

Runs GRLP's threshold-width solver on a single-profile uplift transient and
measures, as the time step lengthens, two complementary quantities:

  * **per-step (local) error** -- how far a single step lands from the truth,
    measured in a smoothly evolving part of the run; and
  * **final-answer (global) error** -- how far the whole run lands from the truth
    at a fixed end time (one equilibration time).

The reference ("truth") in each case is the *same* transient integrated with a
much smaller time step. The figure shows the current first-order-in-time
behaviour: per-step error grows like step^2, while the final-answer error -- the
one a modeller cares about -- grows like step (the small per-step errors
accumulate over the proportionally larger number of steps taken at small dt).

Reuses the canonical single-profile setup from the test suite
(``tests/conftest.make_long_profile``) so this figure and
``tests/test_time_accuracy.py`` describe the same model.  Requires matplotlib
(and pytest, for the shared fixture module).  Run from the repository root::

    python docs/accuracy_figure.py
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_here = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(_here, "..", "tests"))
from conftest import make_long_profile, STEADY_NT, STEADY_DT   # noqa: E402

YEAR = 3.15e7
UPLIFT = 1.0e-4 / YEAR      # 0.1 mm/yr step change in uplift, in m/s


def _steady():
    lp = make_long_profile(U=0.0)
    lp.evolve_threshold_width_river(nt=STEADY_NT, dt=STEADY_DT)
    return lp


def _mid_transient(t_mid, n=400):
    """A smoothly-evolving state, ``t_mid`` into the uplift transient."""
    lp = _steady()
    lp.set_uplift_rate(UPLIFT)
    lp.evolve_threshold_width_river(nt=n, dt=t_mid / n)
    return lp


def per_step_error(step_years):
    """Error from a single step of the given length, in a smooth part of run."""
    tau = step_years * YEAR
    one = _mid_transient(100e3 * YEAR)
    one.evolve_threshold_width_river(nt=1, dt=tau)
    ref = _mid_transient(100e3 * YEAR)
    ref.evolve_threshold_width_river(nt=256, dt=tau / 256)
    return np.max(np.abs(one.z - ref.z))


def final_answer_error(n_steps, T, z_ref, order=1):
    lp = _steady()
    lp.set_uplift_rate(UPLIFT)
    lp.set_time_integration(order)   # 1 = backward Euler, 2 = BDF2
    lp.evolve_threshold_width_river(nt=n_steps, dt=T / n_steps)
    return np.max(np.abs(lp.z - z_ref))


def main():
    # (A) per-step (local) error vs step length, in a smooth part of the run.
    stepsA = np.array([250, 500, 1000, 2000, 4000, 8000, 16000], float)
    locA = np.array([per_step_error(s) for s in stepsA])

    # (B) final-answer (global) error at one equilibration time vs step length.
    lp = _steady()
    lp.compute_equilibration_time()
    T = lp.equilibration_time
    ref = _steady()
    ref.set_uplift_rate(UPLIFT)
    ref.set_time_integration(2)          # BDF2 fine reference = the "truth"
    ref.evolve_threshold_width_river(nt=4096, dt=T / 4096)
    z_ref = ref.z.copy()
    nsB = [4, 8, 16, 32, 64, 128, 256]
    stepsB = np.array([T / n / YEAR for n in nsB])
    glob_be = np.array([final_answer_error(n, T, z_ref, order=1) for n in nsB])
    glob_bdf2 = np.array([final_answer_error(n, T, z_ref, order=2) for n in nsB])

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
    ax[0].loglog(stepsA, locA, "o-", color="C0")
    ax[0].loglog(stepsA, locA[0] * (stepsA / stepsA[0]) ** 2, "--", color="0.6",
                 label=r"$\propto$ step$^2$")
    ax[0].set_title("Per-step error\n(one step, smooth part of run)")
    ax[0].set_xlabel("time-step length [yr]")
    ax[0].set_ylabel("elevation error [m]")
    ax[0].legend()
    ax[0].grid(True, which="both", alpha=0.3)

    ax[1].loglog(stepsB, glob_be, "s-", color="C1",
                 label="backward Euler (default)")
    ax[1].loglog(stepsB, glob_bdf2, "^-", color="C2",
                 label="BDF2  (set_time_integration(2))")
    ax[1].loglog(stepsB, glob_be[-1] * (stepsB / stepsB[-1]) ** 1, "--",
                 color="0.7", label=r"$\propto$ step")
    ax[1].loglog(stepsB, glob_bdf2[0] * (stepsB / stepsB[0]) ** 2, ":",
                 color="0.7", label=r"$\propto$ step$^2$")
    ax[1].set_title("Final-answer error vs. time step\n(backward Euler vs. BDF2)")
    ax[1].set_xlabel("time-step length [yr]")
    ax[1].set_ylabel("elevation error [m]")
    ax[1].legend()
    ax[1].grid(True, which="both", alpha=0.3)

    fig.suptitle("GRLP time-stepping accuracy: first-order (default) vs. "
                 "second-order (BDF2)", y=1.02)
    fig.tight_layout()

    out = os.path.join(_here, "_static", "timestep_accuracy.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=110, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
