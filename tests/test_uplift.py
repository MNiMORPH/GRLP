"""
Numerical vs. analytical steady-state long profiles *with uplift*.

With uplift (or distributed base-level fall), the steady-state bed is given by
``LongProfile.analytical_threshold_width_uplift``: the long-term bedload
discharge grows downstream as uplifted material is exported,

    Q_s(x) = I * Q_s_0 + (1 - lambda_p) * U * \\int_{x0}^{x} B dx',

and the bed follows from the transport law and integration up from base level.

That reference is evaluated by quadrature on a grid independent of the model, so
the numerical solution does not match it to round-off at any single resolution;
instead it *converges* to it at first order as the grid is refined.  These tests
assert that convergence, which validates the numerical model and the analytical
solution against each other, and check the physically required trends.
"""

import numpy as np
import pytest

import grlp


def make_uplift_profile(nx, dx, U_mm_yr=10.0, P_xB=0.0, intermittency=1.0):
    lp = grlp.LongProfile()
    lp.set_intermittency(intermittency)
    lp.basic_constants()
    lp.bedload_lumped_constants()
    lp.set_hydrologic_constants()
    S0 = 0.01
    lp.set_x(dx=dx, nx=nx, x0=5e3)
    lp.set_z(S0=-S0, z1=0.0)
    lp.set_A(k_xA=1.0)
    lp.set_Q(k_xQ=1.433776163432246e-05, P_xQ=7 / 4.0 * 0.7)
    lp.set_B(k_xB=10.0, P_xB=P_xB)
    lp.set_uplift_rate(U_mm_yr * 1e-3 / 3.15e7)
    lp.set_niter(3)
    lp.set_z_bl(0.0)
    Qs0 = lp.k_Qs * lp.Q[0] * S0 ** (7 / 6.0)
    lp.set_Qs_input_upstream(Qs0)
    lp.evolve_threshold_width_river(nt=400, dt=1e13)
    return lp


def _max_rel_error(lp):
    z_analytical = lp.analytical_threshold_width_uplift()
    relief = lp.z.max() - lp.z.min()
    return np.abs(lp.z - z_analytical).max() / relief


UPLIFT_CASES = [
    pytest.param({}, id="uniform-width"),
    pytest.param({"P_xB": 0.2}, id="widening"),
    pytest.param({"intermittency": 0.7}, id="intermittent"),
    pytest.param({"U_mm_yr": 1.0}, id="slow-uplift"),
]


@pytest.mark.parametrize("case", UPLIFT_CASES)
def test_converges_to_analytical_uplift(case):
    # Halving dx should roughly halve the misfit (first-order scheme), and the
    # finest grid should be within ~1% of the analytical solution.
    errors = []
    for nx, dx in [(70, 1000.0), (140, 500.0), (280, 250.0)]:
        lp = make_uplift_profile(nx=nx, dx=dx, **case)
        errors.append(_max_rel_error(lp))
    assert errors[1] < 0.6 * errors[0]
    assert errors[2] < 0.6 * errors[1]
    assert errors[-1] < 0.01


def test_reaches_steady_state():
    lp = make_uplift_profile(nx=140, dx=500.0)
    # Everywhere the bed should be keeping pace with uplift: dz/dt << U.
    assert np.abs(lp.dz_dt).max() / lp.U < 1e-6


def test_uplift_steepens_relative_to_no_uplift():
    # Adding uplift increases sediment export downstream, requiring a steeper
    # (higher-relief) profile than the no-uplift case with the same supply.
    lp0 = make_uplift_profile(nx=140, dx=500.0, U_mm_yr=0.0)
    lpU = make_uplift_profile(nx=140, dx=500.0, U_mm_yr=10.0)
    assert (lpU.z.max() - lpU.z.min()) > (lp0.z.max() - lp0.z.min())


def test_no_uplift_limit_matches_flat_transport():
    # With U = 0 the uplift solution must reduce to the standard no-uplift
    # power-law solution. Anchor the closed-form solution to the same endpoints
    # (the uplift solution is pinned to exact base level and supply, whereas the
    # default closed form is pinned to the numerical bed's endpoints), so this
    # is a pure analytical-vs-analytical check of the quadrature.
    lp = make_uplift_profile(nx=140, dx=500.0, U_mm_yr=0.0)
    z_uplift = lp.analytical_threshold_width_uplift()
    z_base = lp.analytical_threshold_width(z0=z_uplift[0], z1=z_uplift[-1])
    relief = z_uplift.max() - z_uplift.min()
    assert np.abs(z_uplift - z_base).max() / relief < 1e-3


def test_missing_parameters_raise():
    # The analytical solution needs the power-law parameters; a profile built
    # by passing B directly (no k_xB/P_xB) should raise a clear error.
    lp = grlp.LongProfile()
    lp.basic_constants()
    lp.bedload_lumped_constants()
    lp.set_hydrologic_constants()
    lp.set_x(dx=1000.0, nx=50, x0=5e3)
    lp.set_z(S0=-0.01, z1=0.0)
    lp.set_A(k_xA=1.0)
    lp.set_Q(k_xQ=1.433776163432246e-05, P_xQ=7 / 4.0 * 0.7)
    lp.set_B(B=100.0)  # scalar width: k_xB / P_xB never stored
    lp.set_uplift_rate(1e-10)
    with pytest.raises(ValueError):
        lp.analytical_threshold_width_uplift()


def test_deprecated_perturbation_warns():
    lp = make_uplift_profile(nx=50, dx=1000.0)
    lp.analytical_threshold_width()  # sets k_a, P_a used by the old method
    with pytest.warns(DeprecationWarning):
        lp.analytical_threshold_width_perturbation()


def test_uplift_defaults_to_zero():
    # A freshly constructed LongProfile has zero uplift, so U is a usable
    # number, not an unset attribute.
    assert grlp.LongProfile().U == 0.0


def test_evolve_without_set_uplift_rate():
    # Uplift is optional: a profile must evolve without ever calling
    # set_uplift_rate. The solver reads lp.U, so the constructor default
    # (U = 0) has to be in place -- otherwise this raises
    # AttributeError: 'LongProfile' object has no attribute 'U'.
    lp = grlp.LongProfile()
    lp.basic_constants()
    lp.bedload_lumped_constants()
    lp.set_hydrologic_constants()
    lp.set_x(dx=1000.0, nx=50, x0=5e3)
    lp.set_z(S0=-0.01, z1=0.0)
    lp.set_A(k_xA=1.0)
    lp.set_Q(k_xQ=1.433776163432246e-05, P_xQ=7 / 4.0 * 0.7)
    lp.set_B(k_xB=10.0, P_xB=0.0)
    lp.set_niter(3)
    lp.set_z_bl(0.0)
    Qs0 = lp.k_Qs * lp.Q[0] * 0.01 ** (7 / 6.0)
    lp.set_Qs_input_upstream(Qs0)
    # No set_uplift_rate call.
    lp.evolve_threshold_width_river(nt=10, dt=1e13)
    assert np.all(np.isfinite(lp.z))
