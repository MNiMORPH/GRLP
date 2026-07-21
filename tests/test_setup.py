"""
Unit tests for the ``LongProfile`` setter methods.

These check that the model is internally self-consistent *before* any time
integration: that the domain, elevations, drainage area, discharge, and width
arrays have the right shapes and satisfy the relationships (linear grid, linear
initial bed, power laws) that the setters promise.  The extended ("_ext")
arrays, which carry the ghost nodes used for boundary conditions, are checked
against the interior arrays.
"""

import numpy as np
import pytest

from conftest import make_long_profile
import grlp


# --------------------------------------------------------------------------- #
# Domain: set_x
# --------------------------------------------------------------------------- #

def test_set_x_regular_grid_shapes_and_spacing():
    lp = grlp.LongProfile()
    lp.set_x(dx=1000.0, nx=90, x0=10000.0)

    assert lp.nx == 90
    assert len(lp.x) == 90
    assert lp.x[0] == pytest.approx(10000.0)
    np.testing.assert_allclose(np.diff(lp.x), 1000.0)


def test_set_x_ext_brackets_interior_by_one_cell():
    lp = grlp.LongProfile()
    lp.set_x(dx=1000.0, nx=90, x0=10000.0)

    # x_ext adds exactly one ghost node on each side.
    assert len(lp.x_ext) == lp.nx + 2
    np.testing.assert_allclose(lp.x_ext[1:-1], lp.x)
    assert lp.x_ext[0] == pytest.approx(lp.x[0] - 1000.0)
    assert lp.x_ext[-1] == pytest.approx(lp.x[-1] + 1000.0)


def test_set_x_length_and_dx_ext():
    lp = grlp.LongProfile()
    lp.set_x(dx=1000.0, nx=90, x0=10000.0)

    np.testing.assert_allclose(lp.dx_ext, 1000.0)
    assert lp.L == pytest.approx(lp.x_ext[-1] - lp.x_ext[0])


def test_set_x_from_x_ext_recovers_interior():
    x_ext = np.arange(0.0, 101000.0, 1000.0)
    lp = grlp.LongProfile()
    lp.set_x(x_ext=x_ext)

    np.testing.assert_allclose(lp.x, x_ext[1:-1])
    assert lp.nx == len(x_ext) - 2


# --------------------------------------------------------------------------- #
# Initial bed: set_z
# --------------------------------------------------------------------------- #

def test_set_z_linear_slope_and_downstream_boundary():
    S0 = 1.5e-2
    lp = grlp.LongProfile()
    lp.set_x(dx=1000.0, nx=90, x0=10000.0)
    lp.set_z(S0=-S0, z1=0.0)

    # Bed descends downstream at the prescribed slope.
    np.testing.assert_allclose(np.diff(lp.z) / np.diff(lp.x), -S0)
    # z1 is the elevation at the downstream end.
    assert lp.z[-1] == pytest.approx(0.0)


def test_set_z_is_linear_from_S0():
    S0 = 1.5e-2
    lp = grlp.LongProfile()
    lp.set_x(dx=1000.0, nx=90, x0=10000.0)
    lp.set_z(S0=-S0, z1=0.0)

    # set_z(S0) lays down a straight line of slope -S0. z_ext is no longer
    # stored; the upstream ghost is set later, by set_Qs_input_upstream.
    np.testing.assert_allclose(np.diff(lp.z) / np.diff(lp.x), -S0)


# --------------------------------------------------------------------------- #
# Power-law transfer functions: set_A, set_Q, set_B
# --------------------------------------------------------------------------- #

def test_set_A_power_law():
    lp = grlp.LongProfile()
    lp.set_hydrologic_constants()
    lp.set_x(dx=1000.0, nx=90, x0=10000.0)
    lp.set_A(k_xA=1.0)

    np.testing.assert_allclose(lp.A, 1.0 * lp.x ** lp.P_xA)
    np.testing.assert_allclose(lp.A_ext[1:-1], lp.A)


def test_set_Q_power_law():
    lp = grlp.LongProfile()
    lp.set_hydrologic_constants()
    lp.set_x(dx=1000.0, nx=90, x0=10000.0)
    lp.set_A(k_xA=1.0)
    k_xQ, P_xQ = 1.43e-5, 7 / 4.0 * 0.7
    lp.set_Q(k_xQ=k_xQ, P_xQ=P_xQ)

    np.testing.assert_allclose(lp.Q, k_xQ * lp.x ** P_xQ)


def test_set_B_power_law():
    lp = grlp.LongProfile()
    lp.set_hydrologic_constants()
    lp.set_x(dx=1000.0, nx=90, x0=10000.0)
    lp.set_A(k_xA=1.0)
    lp.set_B(k_xB=25.0, P_xB=0.2)

    np.testing.assert_allclose(lp.B, 25.0 * lp.x ** 0.2)


def test_set_B_scalar_broadcasts():
    lp = grlp.LongProfile()
    lp.set_x(dx=1000.0, nx=90, x0=10000.0)
    lp.set_B(B=100.0)

    assert lp.B.shape == lp.x.shape
    np.testing.assert_allclose(lp.B, 100.0)


# --------------------------------------------------------------------------- #
# Hydrologic constants and sediment-supply boundary condition
# --------------------------------------------------------------------------- #

def test_hydrologic_constants_default_consistency():
    lp = grlp.LongProfile()
    lp.set_hydrologic_constants()
    # P_xQ is the product of the inverse-Hack and area--discharge exponents.
    assert lp.P_xQ == pytest.approx(lp.P_xA * lp.P_AQ)


def test_Qs_input_sets_consistent_boundary_slope():
    # The boundary slope S0 recovered from Q_s_0 must reproduce Q_s_0 through
    # the transport law k_Qs * Q * S0**(7/6).
    lp = make_long_profile()
    Qs0_recovered = lp.k_Qs * lp.Q[0] * lp.S0 ** (7 / 6.0)
    assert Qs0_recovered == pytest.approx(lp.Q_s_0, rel=1e-12)


def test_Qs_input_sets_upstream_ghost_node():
    # set_Qs_input_upstream lifts the upstream ghost node so its slope to the
    # first interior node equals the transport slope S0.
    lp = make_long_profile()
    ghost_slope = (lp.z_ghost_upstream - lp.z[0]) / lp.dx[0]
    assert ghost_slope == pytest.approx(lp.S0, rel=1e-12)


def test_intermittency_is_stored():
    lp = grlp.LongProfile()
    lp.set_intermittency(0.3)
    assert lp.intermittency == pytest.approx(0.3)
