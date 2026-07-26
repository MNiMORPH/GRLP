"""
Guard the adaptive-vs-fixed cost finding (MNiMORPH/GRLP#16).

A small, fast version of ``benchmarks/adaptive_timestep_benchmark.py`` run in the
suite so the headline result stays reproducible and cannot silently rot:
**adaptive (step-doubling) BDF2 costs more linear solves than fixed BDF2 at
matched accuracy** for GRLP's smooth diffusive transients. If a future change
made adaptive competitive, this test would flag it (and should then be updated
with the new finding). The full sweep and both transients live in the benchmark
script; here we check the qualitative ordering on the smooth transient with a
coarse reference, to keep it quick.
"""
import os
import sys

import numpy as np

_BENCH = os.path.join(os.path.dirname(__file__), "..", "benchmarks")
sys.path.insert(0, os.path.abspath(_BENCH))
from adaptive_timestep_benchmark import (   # noqa: E402
    measure, cost_ratio, _reference, _equilibration_time)


def test_adaptive_costs_more_than_fixed_at_matched_accuracy():
    """At matched accuracy, adaptive uses clearly more linear solves than fixed
    BDF2 -- the benchmark's headline. Generous margin so it guards the finding,
    not a fragile exact ratio."""
    T = _equilibration_time()
    z_ref = _reference("uplift", T, n=1024)
    result = measure("uplift", T, z_ref,
                     fixed_nts=[16, 32, 64, 128],
                     adaptive_tols=[1e-2, 1e-3, 1e-4])
    # Both schemes converge (error falls as work grows).
    assert np.all(np.diff(result["fixed"]["err"]) < 0.0)
    assert np.all(np.diff(result["adaptive"]["err"]) < 0.0)
    ratio, a, f = cost_ratio(result, target_err=2.0e-3)
    assert ratio > 1.5, (
        "expected adaptive to cost more than fixed at matched accuracy; "
        "adaptive=%.0f fixed=%.0f solves (ratio %.2f)" % (a, f, ratio))
