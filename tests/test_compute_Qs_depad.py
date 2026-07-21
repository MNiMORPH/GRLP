"""
De-padded slope / sediment-discharge diagnostics (Step 1 de-pad).

``compute_Q_s`` no longer reads a maintained padded ``z_ext`` array. Two paths:

* ``LongProfile.compute_Q_s`` -- a lone segment, reconstructing its head (``S0``)
  and outlet (``z_bl``) ghosts from ``self.z``;
* ``Network.compute_Q_s`` -- walks the topology to each node's real neighbor and
  sets ``S`` / ``Q_s`` on every segment.

For a *single-segment* network the two must agree to machine precision (same
physics, same ghosts), given the same elevations. And calling the single-segment
method on a *networked* segment -- whose across-junction neighbor it cannot see
-- must fail loudly rather than silently return a wrong slope.

These guards do not depend on ``z_ext`` being synced, so they remain valid once
the compatibility sync and the ``_ext`` machinery are removed.
"""

import numpy as np
import pytest

import grlp


S0 = 1.5e-2
DT = 1e12
NT = 800


def _standalone(x, Q, B, dx, niter=5):
    lp = grlp.LongProfile()
    lp.set_intermittency(1.0)
    lp.basic_constants()
    lp.bedload_lumped_constants()
    lp.set_hydrologic_constants()
    lp.set_x(dx=dx, nx=len(x), x0=x[0])
    lp.set_z(S0=-S0)
    lp.set_A(k_xA=1.0)
    lp.set_Q(Q=Q.copy())
    lp.set_B(B=B.copy())
    lp.set_uplift_rate(0.0)
    lp.set_niter(niter)
    lp.set_Qs_input_upstream(lp.k_Qs * lp.Q[0] * S0 ** (7 / 6.0))
    lp.set_z_bl(0.0)
    return lp


def _single_segment_network(x, Q, B, dx, niter=5):
    net = grlp.Network()
    net.initialize(
        x_bl=x[-1] + dx, z_bl=0.0, S0=[S0], Q_s_0=None,
        upstream_segment_IDs=[[]], downstream_segment_IDs=[[]],
        x=[x.copy()], z=[np.zeros_like(x)], Q=[Q.copy()], B=[B.copy()],
    )
    net.list_of_LongProfile_objects[0].set_intermittency(1.0)
    net.set_niter(niter)
    net.get_z_lengths()
    return net


def _confluence_network():
    d = 1000.0
    xt = d * np.arange(1, 6.0)
    xtr = d * np.arange(6, 11.0)
    Qt = 2.0 + 0.1 * np.arange(5)
    Qtrunk = 4.0 + 0.2 * np.arange(5)
    B = 100.0 * np.ones(5)
    net = grlp.Network()
    net.initialize(
        x_bl=xtr[-1] + d, z_bl=0.0, S0=[S0, S0], Q_s_0=None,
        upstream_segment_IDs=[[], [], [0, 1]],
        downstream_segment_IDs=[[2], [2], []],
        x=[xt.copy(), xt.copy(), xtr.copy()],
        z=[np.zeros(5), np.zeros(5), np.zeros(5)],
        Q=[Qt.copy(), Qt.copy(), Qtrunk.copy()],
        B=[B.copy(), B.copy(), B.copy()],
    )
    net.set_niter(5)
    net.get_z_lengths()
    return net


def test_network_compute_Qs_matches_standalone_single_segment():
    """Walk-based Network.compute_Q_s reproduces the standalone diagnostic for a
    single-segment network. Elevations are copied across first, so this isolates
    the slope/flux computation (ghost reconstruction + neighbor walk) from any
    O(dx**2) drift between the two solvers."""
    dx = 1000.0
    x = dx * np.arange(10, 100.0)
    Q = 2.0 + 0.05 * np.arange(len(x))      # downstream-increasing (ghosts bite)
    B = 100.0 * np.ones(len(x))

    net = _single_segment_network(x, Q, B, dx)
    net.evolve_threshold_width_river_network(nt=NT, dt=DT)
    seg = net.list_of_LongProfile_objects[0]

    a = _standalone(x, Q, B, dx)
    a.set_z(z=seg.z.copy())                 # identical elevations
    a.compute_Q_s()
    net.compute_Q_s()

    np.testing.assert_allclose(seg.S, a.S, rtol=0, atol=1e-12)
    np.testing.assert_allclose(seg.Q_s, a.Q_s, rtol=0, atol=1e-12)


def test_longprofile_compute_Qs_raises_on_networked_segment():
    """A networked segment cannot see its across-junction neighbor; the
    single-segment method must refuse rather than return a wrong slope, and
    point the caller at Network.compute_Q_s."""
    net = _confluence_network()
    net.evolve_threshold_width_river_network(nt=200, dt=1e10)
    for seg in net.list_of_LongProfile_objects:
        with pytest.raises(ValueError, match="Network.compute_Q_s"):
            seg.compute_Q_s()


def test_network_compute_Qs_populates_every_segment():
    """Network.compute_Q_s sets S and Q_s on every segment, finite and with
    downstream sediment transport (Q_s > 0 for a descending network)."""
    net = _confluence_network()
    net.evolve_threshold_width_river_network(nt=1000, dt=1e10)
    net.compute_Q_s()
    for seg in net.list_of_LongProfile_objects:
        assert seg.S is not None and seg.Q_s is not None
        assert seg.S.shape == seg.z.shape
        assert np.all(np.isfinite(seg.S))
        assert np.all(np.isfinite(seg.Q_s))
        assert np.all(seg.Q_s > 0)
