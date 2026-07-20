"""
Sediment conservation at a multi-tributary confluence with discharge varying
*within* each segment.

Every other network conservation test uses uniform Q per segment, where the
land_area confluence discretization is exact. With Q varying within the
segments the current confluence does NOT conserve sediment at steady state
(~0.85% imbalance, converged -- not a transient) and the junction elevation does
not converge under grid refinement. This was a real bug (the padded land_area confluence leaked ~0.85%);
the de-padded walking solver's conservative three-node junction cell fixes it.
"""
import numpy as np
import grlp


def test_varying_Q_confluence_conserves_sediment():
    nseg = 30
    dx = 1000.0
    xt = dx * np.arange(1, nseg + 1.0)
    xtr = dx * np.arange(nseg + 1, 2 * nseg + 1.0)
    Qt = 2.0 + 0.05 * np.arange(nseg)          # varying within each tributary
    Qtr = 4.0 + 0.10 * np.arange(nseg)
    net = grlp.Network()
    net.initialize(
        x_bl=dx * (2 * nseg + 1), z_bl=0.0, S0=[1.5e-2, 1.5e-2], Q_s_0=None,
        upstream_segment_IDs=[[], [], [0, 1]],
        downstream_segment_IDs=[[2], [2], []],
        x=[xt.copy(), xt.copy(), xtr.copy()], z=[np.zeros(nseg)] * 3,
        Q=[Qt.copy(), Qt.copy(), Qtr.copy()], B=[100.0 * np.ones(nseg)] * 3,
    )
    for lp in net.list_of_LongProfile_objects:
        lp.set_intermittency(1.0)
    net.set_niter(4)
    net.get_z_lengths()
    net.evolve_threshold_width_river_network(nt=3000, dt=1e12)

    seg = net.list_of_LongProfile_objects
    kQs = seg[0].k_Qs

    def qs_last(lp, first=False):
        S = -np.diff(lp.z) / np.diff(lp.x)
        Qf = (lp.Q[:-1] + lp.Q[1:]) / 2.0
        flux = kQs * Qf * S ** (7 / 6.0)
        return flux[0] if first else flux[-1]

    influx = qs_last(seg[0]) + qs_last(seg[1])
    outflux = qs_last(seg[2], first=True)
    assert abs(outflux - influx) / influx < 1e-6


def test_nested_confluences_conserve_with_varying_Q():
    """The three-node junction cell must compose: sediment conserved at every
    junction of a nested (Strahler-3) network with discharge varying within each
    segment. Guards that the de-padded confluence works beyond a single symmetric
    junction -- nested tributaries whose own inflow is itself a confluence."""
    ns = 20
    dx = 1000.0

    def X(k):
        return dx * np.arange(k * ns + 1, (k + 1) * ns + 1.0)

    def Qx(base, slope, k):
        return base + slope * (X(k) - X(k)[0]) / (ns * dx)

    # heads 0,1 -> mid 2 ; head 3 + mid 2 -> trunk 4 -> outlet
    x = [X(0), X(1), X(2), X(3), X(4)]
    Q = [Qx(2.0, 1.0, 0), Qx(2.0, 1.0, 1), Qx(4.0, 2.0, 2),
         Qx(4.0, 2.0, 3), Qx(8.0, 4.0, 4)]
    B = [100.0 * np.ones(ns)] * 5
    net = grlp.Network()
    net.initialize(
        x_bl=dx * (5 * ns + 1), z_bl=0.0, S0=[1.5e-2, 1.5e-2, 1.5e-2], Q_s_0=None,
        upstream_segment_IDs=[[], [], [0, 1], [], [2, 3]],
        downstream_segment_IDs=[[2], [2], [4], [4], []],
        x=[xi.copy() for xi in x], z=[np.zeros(ns)] * 5,
        Q=[qi.copy() for qi in Q], B=[bi.copy() for bi in B],
    )
    for lp in net.list_of_LongProfile_objects:
        lp.set_intermittency(1.0)
    net.set_niter(5)
    net.get_z_lengths()
    net.evolve_threshold_width_river_network(nt=4000, dt=1e12)

    segs = net.list_of_LongProfile_objects
    kQs = segs[0].k_Qs

    def qs(lp, first=False):
        S = -np.diff(lp.z) / np.diff(lp.x)
        Qf = (lp.Q[:-1] + lp.Q[1:]) / 2.0
        f = kQs * Qf * S ** (7 / 6.0)
        return f[0] if first else f[-1]

    for c in segs:                        # every confluence conserves
        if len(c.upstream_segment_IDs) > 1:
            influx = sum(qs(segs[t]) for t in c.upstream_segment_IDs)
            np.testing.assert_allclose(qs(c, first=True), influx, rtol=1e-9)
