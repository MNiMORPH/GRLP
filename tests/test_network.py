"""
River-network steady-state behavior and analytical checks.

The networked solver couples single-segment long profiles across confluences.
At steady state with no uplift and no loss, the physics reduces to statements
that can be checked exactly:

* Every channel-head segment carries the slope S0 imposed as its boundary
  condition.
* Bedload sediment discharge is conserved through the whole network: a trunk
  transports the sum of the during-flood supplies of the tributaries above it
  (Q_s = k_Qs * Q * S**(7/6), summed).  In a symmetric confluence the summed
  supply grows in step with discharge, so the trunk slope also equals S0.
* A degenerate "chain" network with uniform discharge must reproduce the exact
  single-segment linear profile z = S0 * (x_bl - x).

These are all verified below to near machine precision, which exercises the
block-matrix assembly, the internal (confluence) boundary conditions, and the
external upstream/downstream boundary conditions together.
"""

import numpy as np
import pytest

import grlp


# Steady state is reached to machine precision on these small domains with far
# fewer, larger steps than a naive integration would use (see the tuning in the
# test-development notes): ~5 s wall time.
NET_NT = 2000
NET_DT = 1.0e9

S0 = 0.015
QH = 5.0
B = 100.0
DX = 2000.0
NSEG = 4


def _uniform_segment_arrays(n, value):
    return value * np.ones(n)


def make_confluence_network():
    """
    Symmetric Y network: two heads (IDs 0, 1), each discharge QH, joining a
    trunk (ID 2) of discharge 2*QH.  Evolved to steady state.
    """
    x = [
        DX * np.arange(1, NSEG + 1, dtype=float),          # head 0
        DX * np.arange(1, NSEG + 1, dtype=float),          # head 1
        DX * np.arange(NSEG + 1, 2 * NSEG + 1, dtype=float),  # trunk
    ]
    net = grlp.Network()
    net.initialize(
        x_bl=DX * (2 * NSEG + 1),
        z_bl=0.0,
        S0=[S0, S0],  # one per channel head; scalar branch is not wired up
        Q_s_0=None,
        upstream_segment_IDs=[[], [], [0, 1]],
        downstream_segment_IDs=[[2], [2], []],
        x=x,
        z=[np.zeros(NSEG) for _ in range(3)],
        Q=[
            _uniform_segment_arrays(NSEG, QH),
            _uniform_segment_arrays(NSEG, QH),
            _uniform_segment_arrays(NSEG, 2 * QH),
        ],
        B=[_uniform_segment_arrays(NSEG, B) for _ in range(3)],
    )
    net.set_niter(3)
    net.get_z_lengths()
    net.evolve_threshold_width_river_network(nt=NET_NT, dt=NET_DT)
    return net


def make_chain_network():
    """
    Two segments in series with uniform discharge QH throughout: head (ID 0) ->
    downstream (ID 1).  Should collapse to a single uniform-slope profile.
    """
    x = [
        DX * np.arange(1, NSEG + 1, dtype=float),
        DX * np.arange(NSEG + 1, 2 * NSEG + 1, dtype=float),
    ]
    net = grlp.Network()
    net.initialize(
        x_bl=DX * (2 * NSEG + 1),
        z_bl=0.0,
        S0=[S0],
        Q_s_0=None,
        upstream_segment_IDs=[[], [0]],
        downstream_segment_IDs=[[1], []],
        x=x,
        z=[np.zeros(NSEG), np.zeros(NSEG)],
        Q=[_uniform_segment_arrays(NSEG, QH), _uniform_segment_arrays(NSEG, QH)],
        B=[_uniform_segment_arrays(NSEG, B), _uniform_segment_arrays(NSEG, B)],
    )
    net.set_niter(3)
    net.get_z_lengths()
    net.evolve_threshold_width_river_network(nt=NET_NT, dt=NET_DT)
    return net


def _during_flood_Qs(lp):
    """During-flood bedload discharge on the interior faces of a segment."""
    S = np.abs(np.diff(lp.z) / np.diff(lp.x))
    Q_mid = (lp.Q[:-1] + lp.Q[1:]) / 2.0
    return lp.k_Qs * Q_mid * S ** (7 / 6.0)


@pytest.fixture(scope="module")
def confluence_net():
    return make_confluence_network()


@pytest.fixture(scope="module")
def chain_net():
    return make_chain_network()


def test_head_segments_have_imposed_slope(confluence_net):
    for ID in confluence_net.list_of_channel_head_segment_IDs:
        lp = confluence_net.list_of_LongProfile_objects[ID]
        slope = -np.diff(lp.z) / np.diff(lp.x)
        np.testing.assert_allclose(slope, S0, rtol=1e-8)


def test_symmetric_trunk_slope_equals_S0(confluence_net):
    # Summed supply grows in proportion to discharge, so the trunk slope is
    # unchanged from the heads in a symmetric confluence.
    trunk = confluence_net.list_of_LongProfile_objects[2]
    slope = -np.diff(trunk.z) / np.diff(trunk.x)
    np.testing.assert_allclose(slope, S0, rtol=1e-8)


def test_network_sediment_conservation(confluence_net):
    net = confluence_net
    k_Qs = net.list_of_LongProfile_objects[0].k_Qs
    Qs_supply = sum(
        k_Qs * net.list_of_LongProfile_objects[ID].Q[0] * S0 ** (7 / 6.0)
        for ID in net.list_of_channel_head_segment_IDs
    )
    # Each head carries its own share; the trunk carries the full sum.
    for ID in [0, 1]:
        lp = net.list_of_LongProfile_objects[ID]
        np.testing.assert_allclose(_during_flood_Qs(lp), Qs_supply / 2.0, rtol=1e-6)
    trunk = net.list_of_LongProfile_objects[2]
    np.testing.assert_allclose(_during_flood_Qs(trunk), Qs_supply, rtol=1e-6)


def test_uniform_chain_reduces_to_linear_profile(chain_net):
    net = chain_net
    x_all = np.hstack([lp.x for lp in net.list_of_LongProfile_objects])
    z_all = np.hstack([lp.z for lp in net.list_of_LongProfile_objects])
    x_bl = DX * (2 * NSEG + 1)
    # z_bl = 0 at x_bl, uniform slope S0 upstream.
    z_linear = S0 * (x_bl - x_all)
    assert np.abs(z_all - z_linear).max() < 1e-6


def test_uniform_chain_has_uniform_slope(chain_net):
    net = chain_net
    x_all = np.hstack([lp.x for lp in net.list_of_LongProfile_objects])
    z_all = np.hstack([lp.z for lp in net.list_of_LongProfile_objects])
    slope = -np.diff(z_all) / np.diff(x_all)
    np.testing.assert_allclose(slope, S0, rtol=1e-8)
