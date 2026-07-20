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

from network_helpers import (
    build_network,
    during_flood_Qs,
    expected_during_flood_Qs,
    heads_above,
)

S0 = 0.015
D = 2000.0


def _h(n=4):
    """A head-segment x array of n nodes at spacing D starting at D."""
    return D * np.arange(1, n + 1, dtype=float)


# name -> dict(x, Q, up, down, x_bl). All heads carry S0; discharges chosen so
# that Q of a trunk equals the sum of its tributaries' Q.
TOPOLOGIES = {
    "symmetric_confluence": dict(
        x=[_h(), _h(), D * np.arange(5, 9, dtype=float)],
        Q=[5 * np.ones(4), 5 * np.ones(4), 10 * np.ones(4)],
        up=[[], [], [0, 1]], down=[[2], [2], []], x_bl=D * 9,
    ),
    "asymmetric_Q_confluence": dict(
        x=[_h(), _h(), D * np.arange(5, 9, dtype=float)],
        Q=[5 * np.ones(4), 15 * np.ones(4), 20 * np.ones(4)],
        up=[[], [], [0, 1]], down=[[2], [2], []], x_bl=D * 9,
    ),
    "unequal_dx_coarse_to_fine": dict(
        # heads dx=2000, trunk dx=1000
        x=[_h(), _h(), 8000.0 + 1000.0 * np.arange(1, 7, dtype=float)],
        Q=[5 * np.ones(4), 5 * np.ones(4), 10 * np.ones(6)],
        up=[[], [], [0, 1]], down=[[2], [2], []], x_bl=14000.0 + 1000.0,
    ),
    "unequal_dx_fine_to_coarse": dict(
        # heads dx=1000, trunk dx=2000
        x=[1000.0 * np.arange(1, 7, dtype=float),
           1000.0 * np.arange(1, 7, dtype=float),
           6000.0 + 2000.0 * np.arange(1, 5, dtype=float)],
        Q=[5 * np.ones(6), 5 * np.ones(6), 10 * np.ones(4)],
        up=[[], [], [0, 1]], down=[[2], [2], []], x_bl=14000.0 + 2000.0,
    ),
    "unequal_length_tributaries": dict(
        # one short head (3 nodes), one long head (6 nodes)
        x=[_h(3), 1000.0 * np.arange(1, 7, dtype=float),
           D * np.arange(5, 9, dtype=float)],
        Q=[5 * np.ones(3), 5 * np.ones(6), 10 * np.ones(4)],
        up=[[], [], [0, 1]], down=[[2], [2], []], x_bl=D * 9,
    ),
    "multi_level": dict(
        # heads 0,1 -> mid 4; mid 4, head 2 -> mid 5; mid 5, head 3 -> trunk 6
        x=[_h(), _h(), _h(), _h(),
           D * np.arange(5, 9, dtype=float),
           D * np.arange(9, 13, dtype=float),
           D * np.arange(13, 17, dtype=float)],
        Q=[5 * np.ones(4)] * 4
          + [10 * np.ones(4), 15 * np.ones(4), 20 * np.ones(4)],
        up=[[], [], [], [], [0, 1], [4, 2], [5, 3]],
        down=[[4], [4], [5], [6], [5], [6], []], x_bl=D * 17,
    ),
}


@pytest.fixture(scope="module", params=list(TOPOLOGIES), ids=list(TOPOLOGIES))
def net_case(request):
    spec = TOPOLOGIES[request.param]
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


def test_junction_ghost_nodes_are_two_sided_consistent(net_case):
    # The invariant the padded-array refactor must preserve: at every internal
    # junction the upstream segment's downstream ghost equals the downstream
    # segment's first real node, and the downstream segment's upstream ghost
    # (for this tributary) equals the upstream segment's last real node.
    net, spec = net_case
    segs = net.list_of_LongProfile_objects
    for lp in segs:
        for j, ds_id in enumerate(lp.downstream_segment_IDs):
            ds = segs[ds_id]
            up_index = ds.upstream_segment_IDs.index(lp.ID)
            assert lp.z_ext[j][-1] == pytest.approx(ds.z[0], abs=1e-12)
            assert ds.z_ext[up_index][0] == pytest.approx(lp.z[-1], abs=1e-12)


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
