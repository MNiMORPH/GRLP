"""
Numerical vs. analytical steady-state long profiles (no uplift).

For a transport-limited threshold-width gravel river with no uplift and no
sediment loss, Wickert & Schildgen (2019) give the equilibrium bed as a power
law in downstream distance,

    z(x) = k_a * x**beta + c_a,   beta = 1 - 6*P_xQ/7,

implemented in ``LongProfile.analytical_threshold_width``.  These tests evolve
the numerical model to steady state and confirm it reproduces that closed form
across a range of parameters, and that the run is genuinely at steady state
(dz/dt -> 0).
"""

import numpy as np
import pytest

from conftest import make_long_profile, STEADY_NT, STEADY_DT


# Parameter sets spanning slope, width exponent, discharge exponent,
# intermittency, and grid resolution.  Each is a kwargs dict for
# make_long_profile.
PARAM_SETS = [
    pytest.param({}, id="canonical"),
    pytest.param({"S0": 0.01, "P_xB": 0.0}, id="gentle-uniform-width"),
    pytest.param({"S0": 0.02, "P_xB": 0.4}, id="steep-widening"),
    pytest.param({"P_xQ": 1.0, "k_xQ": 1e-4}, id="P_xQ=1.0"),
    pytest.param({"P_xQ": 1.4, "k_xQ": 1e-6}, id="P_xQ=1.4"),
    pytest.param({"intermittency": 1.0}, id="no-intermittency"),
    pytest.param({"nx": 180, "dx": 500.0}, id="fine-grid"),
]


@pytest.fixture(params=PARAM_SETS)
def evolved_case(request):
    """A profile built from a parameter set and evolved to steady state."""
    lp = make_long_profile(**request.param)
    lp.evolve_threshold_width_river(nt=STEADY_NT, dt=STEADY_DT)
    return lp


def test_matches_analytical_power_law(evolved_case):
    lp = evolved_case
    z_analytical = lp.analytical_threshold_width()
    relief = lp.z.max() - lp.z.min()
    max_abs_err = np.abs(lp.z - z_analytical).max()
    # Relative to the total relief of the profile; the semi-implicit solver
    # reproduces the analytical power law to ~1e-4 on these grids.
    assert max_abs_err / relief < 1e-3


def test_reaches_steady_state(evolved_case):
    lp = evolved_case
    # dz/dt should be negligible compared with any plausible geomorphic rate.
    # A generous ceiling: << 1 mm per thousand years.
    assert np.abs(lp.dz_dt).max() < 1e-15


def test_endpoints_are_pinned(evolved_case):
    # The analytical solution is constructed to pass through the numerical
    # endpoints, so this is a self-consistency check on the boundary handling:
    # the interior fit is what test_matches_analytical_power_law verifies.
    lp = evolved_case
    z_analytical = lp.analytical_threshold_width()
    assert z_analytical[0] == pytest.approx(lp.z[0], rel=1e-9)
    assert z_analytical[-1] == pytest.approx(lp.z[-1], rel=1e-9)


def test_concavity_matches_discharge_exponent():
    # The bed slope is S = |dz/dx| ~ x**(P_a - 1), with P_a = 1 - 6*P_xQ/7,
    # and A ~ x**P_xA. A log S -- log A regression therefore recovers
    # concavity theta = (1 - P_a) / P_xA = (6*P_xQ/7) / P_xA.
    lp = make_long_profile()
    lp.evolve_threshold_width_river(nt=STEADY_NT, dt=STEADY_DT)
    lp.slope_area()
    theta_expected = (6.0 * lp.P_xQ / 7.0) / lp.P_xA
    assert lp.theta == pytest.approx(theta_expected, rel=1e-3)
    assert lp.thetaR2 > 0.999
