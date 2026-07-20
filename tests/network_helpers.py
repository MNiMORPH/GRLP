"""
Helpers for building and checking arbitrary GRLP river networks in tests.

A network is specified by parallel lists (one entry per segment):
  x    : list of 1-D node-position arrays
  Q    : list of discharge arrays
  up   : list of upstream-segment-ID lists
  down : list of downstream-segment-ID lists
and a single downstream base-level position ``x_bl``.

The single sediment-supply slope ``S0`` is applied at every channel head, so the
during-flood bedload supply of a head is ``k_Qs * Q_head[0] * S0**(7/6)`` and,
with no loss, is conserved downstream; a segment carries the sum of the supplies
of all heads above it (see :func:`expected_during_flood_Qs`).
"""

import numpy as np

import grlp

# Large steps drive these small domains to steady state (machine precision) in
# few iterations; sized for the largest (multi-level) network used in the tests.
DEFAULT_NT = 1000
DEFAULT_DT = 3.0e10


def build_network(x, Q, up, down, x_bl, S0=0.015, B=100.0,
                  nt=DEFAULT_NT, dt=DEFAULT_DT, evolve=True):
    """Build (and by default evolve to steady state) a network."""
    n_heads = sum(1 for i in range(len(x)) if len(up[i]) == 0)
    net = grlp.Network()
    net.initialize(
        x_bl=x_bl,
        z_bl=0.0,
        S0=[S0] * n_heads,
        Q_s_0=None,
        upstream_segment_IDs=up,
        downstream_segment_IDs=down,
        x=x,
        z=[np.zeros(len(xi)) for xi in x],
        Q=Q,
        B=[B * np.ones(len(xi)) for xi in x],
    )
    net.set_niter(3)
    net.get_z_lengths()
    if evolve:
        net.evolve_threshold_width_river_network(nt=nt, dt=dt)
    return net


def heads_above(seg_id, up):
    """Return the channel-head segment IDs at or upstream of ``seg_id``."""
    if len(up[seg_id]) == 0:
        return [seg_id]
    result = []
    for u in up[seg_id]:
        result += heads_above(u, up)
    return result


def during_flood_Qs(lp):
    """During-flood bedload discharge on the interior faces of a segment."""
    S = np.abs(np.diff(lp.z) / np.diff(lp.x))
    Q_mid = (lp.Q[:-1] + lp.Q[1:]) / 2.0
    return lp.k_Qs * Q_mid * S ** (7 / 6.0)


def expected_during_flood_Qs(net, seg_id, up, S0):
    """Sum of the during-flood supplies of all heads above ``seg_id``."""
    segs = net.list_of_LongProfile_objects
    return sum(
        segs[h].k_Qs * segs[h].Q[0] * S0 ** (7 / 6.0)
        for h in heads_above(seg_id, up)
    )


# --------------------------------------------------------------------------- #
# Shared catalog of test network topologies.
#
# Used both by the correctness tests (assert physics) and by the
# characterization suite (capture golden-master outputs), so that every topology
# whose behavior is pinned is also checked for correctness.
# --------------------------------------------------------------------------- #

_D = 2000.0


def _head(n=4):
    """A head-segment x array of n nodes at spacing _D starting at _D."""
    return _D * np.arange(1, n + 1, dtype=float)


NETWORK_TOPOLOGIES = {
    "symmetric_confluence": dict(
        x=[_head(), _head(), _D * np.arange(5, 9, dtype=float)],
        Q=[5 * np.ones(4), 5 * np.ones(4), 10 * np.ones(4)],
        up=[[], [], [0, 1]], down=[[2], [2], []], x_bl=_D * 9,
    ),
    "asymmetric_Q_confluence": dict(
        x=[_head(), _head(), _D * np.arange(5, 9, dtype=float)],
        Q=[5 * np.ones(4), 15 * np.ones(4), 20 * np.ones(4)],
        up=[[], [], [0, 1]], down=[[2], [2], []], x_bl=_D * 9,
    ),
    "unequal_dx_coarse_to_fine": dict(
        x=[_head(), _head(), 8000.0 + 1000.0 * np.arange(1, 7, dtype=float)],
        Q=[5 * np.ones(4), 5 * np.ones(4), 10 * np.ones(6)],
        up=[[], [], [0, 1]], down=[[2], [2], []], x_bl=14000.0 + 1000.0,
    ),
    "unequal_dx_fine_to_coarse": dict(
        x=[1000.0 * np.arange(1, 7, dtype=float),
           1000.0 * np.arange(1, 7, dtype=float),
           6000.0 + 2000.0 * np.arange(1, 5, dtype=float)],
        Q=[5 * np.ones(6), 5 * np.ones(6), 10 * np.ones(4)],
        up=[[], [], [0, 1]], down=[[2], [2], []], x_bl=14000.0 + 2000.0,
    ),
    "unequal_length_tributaries": dict(
        x=[_head(3), 1000.0 * np.arange(1, 7, dtype=float),
           _D * np.arange(5, 9, dtype=float)],
        Q=[5 * np.ones(3), 5 * np.ones(6), 10 * np.ones(4)],
        up=[[], [], [0, 1]], down=[[2], [2], []], x_bl=_D * 9,
    ),
    "multi_level": dict(
        x=[_head(), _head(), _head(), _head(),
           _D * np.arange(5, 9, dtype=float),
           _D * np.arange(9, 13, dtype=float),
           _D * np.arange(13, 17, dtype=float)],
        Q=[5 * np.ones(4)] * 4
          + [10 * np.ones(4), 15 * np.ones(4), 20 * np.ones(4)],
        up=[[], [], [], [], [0, 1], [4, 2], [5, 3]],
        down=[[4], [4], [5], [6], [5], [6], []], x_bl=_D * 17,
    ),
    "balanced_tree": dict(
        x=[_head(), _head(), _head(), _head(),
           _D * np.arange(5, 9, dtype=float),
           _D * np.arange(5, 9, dtype=float),
           _D * np.arange(9, 13, dtype=float)],
        Q=[5 * np.ones(4)] * 4
          + [10 * np.ones(4), 10 * np.ones(4), 20 * np.ones(4)],
        up=[[], [], [], [], [0, 1], [2, 3], [4, 5]],
        down=[[4], [4], [5], [5], [6], [6], []], x_bl=_D * 13,
    ),
}


def run_topology_arrays(spec, S0=0.015):
    """
    Build+evolve a topology from NETWORK_TOPOLOGIES and return golden-master
    arrays: concatenated bed elevations and per-segment during-flood Q_s.
    """
    net = build_network(spec["x"], spec["Q"], spec["up"], spec["down"],
                        spec["x_bl"], S0=S0)
    out = {"z_all": np.hstack([lp.z for lp in net.list_of_LongProfile_objects])}
    for lp in net.list_of_LongProfile_objects:
        out["Qs_seg%d" % lp.ID] = during_flood_Qs(lp)
    return out
