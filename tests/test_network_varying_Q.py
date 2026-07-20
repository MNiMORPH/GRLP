"""
A single-segment network must match the standalone solver and the analytical
solution when discharge varies downstream.

This closes a test gap that hid a first-order accuracy bug at channel heads.
Every other network test uses *uniform* discharge per segment, where the
upstream ghost discharge ``2*Q[0] - Q[1]`` equals ``Q[0]`` and the boundary
``dQ/dx`` term vanishes -- so the way the ghost discharge is extrapolated makes
no difference and the bug is invisible.

With a downstream-increasing ``Q`` the choice matters.  The channel head has no
tributary junction and hence no discharge discontinuity, so the two-cell
centered ``dQ/dx`` in the boundary sediment flux is well defined and
second-order.  The former constant ghost (``Q_ext[0] = Q[0]``) collapsed it to a
one-sided, first-order estimate, biasing the injected flux and the whole
equilibrium profile by O(dx) -- ~8.9 m at ``dx = 1000`` on this domain, which is
what this test would catch.

A single-segment network (one channel head straight to base level, no
confluence) is the same physics as a standalone ``LongProfile``; the two must
agree, and both must reproduce ``analytical_threshold_width``.
"""

import numpy as np

import grlp


def _single_segment_network(x, Q, B, S0, dx, niter=4):
    net = grlp.Network()
    net.initialize(
        x_bl=x[-1] + dx, z_bl=0.0, S0=[S0], Q_s_0=None,
        upstream_segment_IDs=[[]], downstream_segment_IDs=[[]],
        x=[x.copy()], z=[np.zeros_like(x)], Q=[Q.copy()], B=[B.copy()],
    )
    # Match the standalone reference's intermittency (steady state is in fact
    # intermittency-independent, but keep it explicit and matched).
    net.list_of_LongProfile_objects[0].set_intermittency(1.0)
    net.set_niter(niter)
    net.get_z_lengths()
    return net


def test_single_segment_network_matches_standalone_and_analytical(
        long_profile_factory):
    S0 = 1.5e-2
    a = long_profile_factory(intermittency=1.0)  # varying Q (k_xQ, P_xQ)
    x = a.x.copy()
    dx = x[1] - x[0]
    assert a.Q[1] != a.Q[0], "test must use downstream-varying Q to bite"
    a.set_z(z=np.zeros_like(x))
    a.set_z_bl(0.0)

    net = _single_segment_network(x, a.Q, a.B, S0, dx)

    a.evolve_threshold_width_river(nt=800, dt=1e12)
    net.evolve_threshold_width_river_network(nt=800, dt=1e12)

    zN = net.list_of_LongProfile_objects[0].z
    za = a.analytical_threshold_width()

    # Second-order boundary flux tracks the analytical (and the standalone) to
    # well under a metre; the bug gave ~8.9 m.  The tolerance sits above the
    # remaining first-order *downstream* (outlet) ghost residual (~0.2 m at this
    # dx) and should tighten once that ghost is likewise made second-order.
    assert np.abs(zN - a.z).max() < 0.5
    assert np.abs(zN - za).max() < 0.5
