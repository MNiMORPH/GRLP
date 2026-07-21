"""
Correctness of the networked solver on non-trivial topologies.

Unlike ``test_network.py`` (small symmetric baseline cases), these assert
*known-true physics* on the network configurations that stress the junction
coupling the hardest: asymmetric discharge, unequal grid spacing across a
junction, unequal tributary lengths, and multi-level (Strahler-3) trees.

They also served as a diagnostic for the planned solver refactor (replacing the
padded ``_ext`` arrays with direct neighbor-reaching): they confirm the current
padded implementation is numerically correct at these junctions -- sediment is
conserved to machine precision and the two-sided junction ghost nodes agree
exactly -- so that work is an architectural cleanup, not a bug fix.

At steady state with no uplift/loss, bedload discharge is conserved through the
whole network: a segment carries the summed during-flood supply of every head
above it, and each channel head carries the imposed slope ``S0``.
"""

import numpy as np
import pytest

import grlp
from network_helpers import (
    NETWORK_TOPOLOGIES,
    build_network,
    during_flood_Qs,
    expected_during_flood_Qs,
    heads_above,
)

S0 = 0.015


def test_per_head_S0_maps_to_channel_head_order():
    # Distinct boundary slopes per head must map to the heads in
    # list_of_channel_head_segment_IDs order. All other network tests use a
    # uniform S0, so this is the guard on channel-head ordering (which the
    # topology-graph migration must preserve, since it sets sediment supply).
    d = 2000.0
    net = grlp.Network()
    net.initialize(
        x_bl=d * 9, z_bl=0.0, S0=[0.01, 0.02], Q_s_0=None,
        upstream_segment_IDs=[[], [], [0, 1]],
        downstream_segment_IDs=[[2], [2], []],
        x=[d * np.arange(1, 5.0), d * np.arange(1, 5.0), d * np.arange(5, 9.0)],
        z=[np.zeros(4)] * 3,
        Q=[5 * np.ones(4), 5 * np.ones(4), 10 * np.ones(4)],
        B=[100.0 * np.ones(4)] * 3,
    )
    net.set_niter(3)
    net.get_z_lengths()
    net.evolve_threshold_width_river_network(nt=1000, dt=3e10)
    assert net.list_of_channel_head_segment_IDs == [0, 1]
    slopes = [
        float(np.mean(-np.diff(net.list_of_LongProfile_objects[i].z)
                      / np.diff(net.list_of_LongProfile_objects[i].x)))
        for i in (0, 1)
    ]
    assert slopes[0] == pytest.approx(0.01, rel=1e-6)
    assert slopes[1] == pytest.approx(0.02, rel=1e-6)


@pytest.fixture(scope="module", params=list(NETWORK_TOPOLOGIES),
                ids=list(NETWORK_TOPOLOGIES))
def net_case(request):
    spec = NETWORK_TOPOLOGIES[request.param]
    net = build_network(spec["x"], spec["Q"], spec["up"], spec["down"],
                        spec["x_bl"], S0=S0)
    return net, spec


def test_reaches_steady_state(net_case):
    net, spec = net_case
    for lp in net.list_of_LongProfile_objects:
        assert np.abs(lp.dz_dt).max() < 1e-13


def test_sediment_conserved_through_network(net_case):
    net, spec = net_case
    up = spec["up"]
    for lp in net.list_of_LongProfile_objects:
        expected = expected_during_flood_Qs(net, lp.ID, up, S0)
        np.testing.assert_allclose(during_flood_Qs(lp), expected, rtol=1e-8)


def test_channel_heads_carry_imposed_slope(net_case):
    net, spec = net_case
    up = spec["up"]
    for lp in net.list_of_LongProfile_objects:
        if len(up[lp.ID]) == 0:
            slope = -np.diff(lp.z) / np.diff(lp.x)
            np.testing.assert_allclose(slope, S0, rtol=1e-8)


def test_sediment_conserved_across_each_junction(net_case):
    # The physical invariant (replacing the old padded-ghost-consistency check,
    # which tested a now-bypassed implementation detail): at every confluence the
    # summed during-flood sediment supply of the tributaries equals the outflow.
    # The de-padded walking solver's three-node junction cell conserves this to
    # machine precision.
    net, spec = net_case
    segs = net.list_of_LongProfile_objects
    for ds in segs:
        if len(ds.upstream_segment_IDs) > 1:
            influx = sum(during_flood_Qs(segs[t])[-1]
                         for t in ds.upstream_segment_IDs)
            outflux = during_flood_Qs(ds)[0]
            np.testing.assert_allclose(outflux, influx, rtol=1e-9)


def test_elevation_decreases_downstream_within_each_segment(net_case):
    net, spec = net_case
    for lp in net.list_of_LongProfile_objects:
        assert np.all(np.diff(lp.z) < 0)


def test_junction_elevation_drop_matches_slope(net_case):
    # Across each junction the elevation drop from an upstream segment's last
    # node to the downstream segment's first node is consistent with a positive
    # transport slope (no spurious step or reversal at the confluence).
    net, spec = net_case
    segs = net.list_of_LongProfile_objects
    for lp in segs:
        for ds_id in lp.downstream_segment_IDs:
            ds = segs[ds_id]
            gap = ds.x[0] - lp.x[-1]
            junction_slope = (lp.z[-1] - ds.z[0]) / gap
            assert junction_slope > 0


def test_network_base_level_x_position_honored():
    # Network.initialize now applies x_bl (it was previously accepted but
    # ignored). A base level moved seaward in x, at fixed z_bl, must shift the
    # mouth ghost and change the solution; the walker reads x_ghost_downstream
    # at the outlet. Networks that pass x_bl at the mirror (x[-1] + dx) are
    # unchanged -- this exercises a genuine move.
    x = [np.arange(0.0, 30000.0, 1000.0)]
    Q = [10.0 * np.ones(len(x[0]))]
    up, down = [[]], [[]]
    mouth_end = x[0][-1]

    mirror = build_network(x, Q, up, down, x_bl=mouth_end + 1000.0)
    moved = build_network(x, Q, up, down, x_bl=mouth_end + 6000.0)

    seg_mirror = mirror.list_of_LongProfile_objects[0]
    seg_moved = moved.list_of_LongProfile_objects[0]
    # x_bl reached the mouth's ghost position
    assert seg_moved.x_ghost_downstream == pytest.approx(mouth_end + 6000.0)
    # base level farther downstream at fixed z_bl -> gentler drop -> higher bed
    assert seg_moved.z[-1] > seg_mirror.z[-1] + 1.0


def test_network_gravel_loss_tracks_through_confluence():
    # Sternberg gravel loss is a purely local sink, but because compute_Q_s
    # walks the topology, the abrasion accumulates along the whole downstream
    # flow path -- through the confluence -- with no per-path bookkeeping. A
    # no-loss network conserves sediment (Q_s ~ constant down the trunk); with
    # loss, Q_s decays monotonically, and the trunk inherits an already-abraded
    # supply from its tributaries (loss upstream of the junction carried across).
    spec = NETWORK_TOPOLOGIES["symmetric_confluence"]
    k_per_km = 0.05

    def built(gravel):
        net = build_network(spec["x"], spec["Q"], spec["up"], spec["down"],
                            spec["x_bl"], evolve=False)
        if gravel:
            for lp in net.list_of_LongProfile_objects:
                lp.gravel_fractional_loss_per_km = k_per_km
        net.evolve_threshold_width_river_network(nt=1000, dt=3.0e10)
        net.compute_Q_s()
        return net

    net0 = built(gravel=False)
    netg = built(gravel=True)
    mouth = net0.list_of_channel_mouth_segment_IDs[0]
    trunk0 = net0.list_of_LongProfile_objects[mouth]
    trunkg = netg.list_of_LongProfile_objects[mouth]

    # No loss: sediment conserved along the trunk (roughly constant Q_s).
    assert trunk0.Q_s.std() / trunk0.Q_s.mean() < 0.05
    # With loss: Q_s decreases monotonically down the trunk (ongoing abrasion).
    assert np.all(np.diff(trunkg.Q_s) < 0)
    # Through the confluence: the trunk head already carries less than the no-loss
    # supply (tributaries abraded upstream), and the deficit grows toward the mouth.
    assert trunkg.Q_s[0] < trunk0.Q_s[0]
    assert (trunk0.Q_s[-1] - trunkg.Q_s[-1]) > (trunk0.Q_s[0] - trunkg.Q_s[0])
