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

from conftest import make_long_profile, STEADY_NT, STEADY_DT
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


def test_set_x_ghost_positions_bracket_interior():
    lp = grlp.LongProfile()
    lp.set_x(dx=1000.0, nx=90, x0=10000.0)

    # The boundary ghost-node positions sit exactly one cell outside the
    # interior grid (the padded x_ext array is no longer stored).
    assert lp.x_ghost_upstream == pytest.approx(lp.x[0] - 1000.0)
    assert lp.x_ghost_downstream == pytest.approx(lp.x[-1] + 1000.0)


def test_set_x_length_and_uniform_spacing():
    lp = grlp.LongProfile()
    lp.set_x(dx=1000.0, nx=90, x0=10000.0)

    np.testing.assert_allclose(lp.dx, 1000.0)
    # L spans the two boundary ghost nodes.
    assert lp.L == pytest.approx(lp.x_ghost_downstream - lp.x_ghost_upstream)


def test_set_x_from_x_ext_recovers_interior():
    x_ext = np.arange(0.0, 101000.0, 1000.0)
    lp = grlp.LongProfile()
    lp.set_x(x_ext=x_ext)

    np.testing.assert_allclose(lp.x, x_ext[1:-1])
    assert lp.nx == len(x_ext) - 2


def test_set_x_from_x_alone_builds_usable_grid():
    # Passing x alone is a documented option: it must build a usable grid with
    # linearly extrapolated ghost positions, not fall through to sys.exit.
    x = np.arange(10000.0, 10000.0 + 1000.0 * 50, 1000.0)
    lp = grlp.LongProfile()
    lp.set_x(x=x, verbose=False)

    assert lp.nx == 50
    np.testing.assert_allclose(lp.x, x)
    assert lp.x_ghost_upstream == pytest.approx(x[0] - 1000.0)
    assert lp.x_ghost_downstream == pytest.approx(x[-1] + 1000.0)
    assert lp.L == pytest.approx(lp.x_ghost_downstream - lp.x_ghost_upstream)


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
    # The boundary ghost drainage areas follow the same power law (the padded
    # A_ext array is no longer stored).
    assert lp.A_ghost_upstream == pytest.approx(1.0 * lp.x_ghost_upstream ** lp.P_xA)
    assert lp.A_ghost_downstream == pytest.approx(
        1.0 * lp.x_ghost_downstream ** lp.P_xA)


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


def test_set_S0_recovers_slope_and_matches_Qs_input():
    # Forcing by slope must reproduce forcing by the equivalent sediment supply:
    # set_S0 derives that supply (intermittency-free inverse; the helper runs
    # at I=0.8) and delegates to set_Qs_input_upstream.
    ref = make_long_profile(S0=0.015)
    alt = make_long_profile(S0=0.015)
    alt.set_S0(0.015)
    assert alt.S0 == pytest.approx(0.015, rel=1e-12)
    assert alt.Q_s_0 == pytest.approx(ref.Q_s_0, rel=1e-12)
    ref.evolve_threshold_width_river(nt=STEADY_NT, dt=STEADY_DT)
    alt.evolve_threshold_width_river(nt=STEADY_NT, dt=STEADY_DT)
    np.testing.assert_allclose(alt.z, ref.z, rtol=0, atol=1e-12)


def test_set_S0_boundary_slope_responds_to_discharge():
    # Under set_S0 the boundary is stored as a sediment supply, so raising the
    # water discharge lowers the boundary slope -- a physical forcing drives the
    # boundary, rather than S0 being pinned independent of Q. S0 ~ Q**(-6/7).
    lp = make_long_profile()
    lp.set_S0(0.015)
    S0_before = lp.S0
    lp.set_Q(Q=2.0 * lp.Q)
    assert lp.S0 < S0_before
    assert lp.S0 == pytest.approx(S0_before * 2.0 ** (-6 / 7.0), rel=1e-12)


def test_intermittency_is_stored():
    lp = grlp.LongProfile()
    lp.set_intermittency(0.3)
    assert lp.intermittency == pytest.approx(0.3)


# --------------------------------------------------------------------------- #
# Base level as a point: set_bl
# --------------------------------------------------------------------------- #

def test_set_bl_sets_both_coordinates():
    # set_bl places base level at the point (x, z) -- position and elevation
    # together -- equivalent to calling set_x_bl and set_z_bl.
    lp = make_long_profile()
    lp.set_bl(x=123456.0, z=7.0)
    assert lp.x_ghost_downstream == pytest.approx(123456.0)
    assert lp.x_bl == pytest.approx(123456.0)
    assert lp.z_bl == pytest.approx(7.0)


def test_set_bl_moves_one_coordinate_at_a_time():
    # Passing only one of x, z moves base level in that coordinate alone.
    lp = make_long_profile()
    lp.set_bl(x=5000.0, z=3.0)
    lp.set_bl(z=9.0)                     # elevation only
    assert lp.z_bl == pytest.approx(9.0)
    assert lp.x_ghost_downstream == pytest.approx(5000.0)
    lp.set_bl(x=8000.0)                  # position only
    assert lp.x_ghost_downstream == pytest.approx(8000.0)
    assert lp.z_bl == pytest.approx(9.0)


def test_set_bl_matches_individual_setters():
    a = make_long_profile()
    b = make_long_profile()
    a.set_bl(x=4321.0, z=3.0)
    b.set_x_bl(4321.0)
    b.set_z_bl(3.0)
    assert a.x_ghost_downstream == b.x_ghost_downstream
    assert a.z_bl == b.z_bl
