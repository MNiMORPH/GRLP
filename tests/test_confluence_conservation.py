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
import pytest
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


@pytest.mark.slow
def test_confluence_convergence_order():
    """The confluence must be CONVERGENT under grid refinement (the old
    land_area scheme was not -- its junction elevation drifted). With a fixed
    Q(x) field, the junction elevation converges at first order (~1). This guards
    against a regression that breaks convergence; a fixed field (not
    resolution-dependent Q) is essential or the test refines a different problem.
    """
    def junction_z(nseg):
        dx = 40000.0 / nseg
        xt = dx * np.arange(1, nseg + 1.0)
        xtr = dx * np.arange(nseg + 1, 2 * nseg + 1.0)
        Qt = 2.0 + xt / 40000.0             # fixed Q(x), independent of nseg
        Qtr = 4.0 + 2 * xtr / 40000.0
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
        net.set_niter(5)
        net.get_z_lengths()
        net.evolve_threshold_width_river_network(nt=4000, dt=1e12)
        return net.list_of_LongProfile_objects[2].z[0]

    z = {n: junction_z(n) for n in (10, 20, 40, 80)}
    d1, d2, d3 = z[20] - z[10], z[40] - z[20], z[80] - z[40]
    orders = [np.log2(abs(d1 / d2)), np.log2(abs(d2 / d3))]
    # converges (differences shrink); at least first order.
    assert min(orders) > 0.8, "junction elevation not converging: orders=%r" % orders


def test_confluence_conserves_with_uplift():
    """With uplift, the confluence cell is not a pure conserver -- it *generates*
    sediment over its land area, exactly (1-lambda)*U*A. This is the quantitative
    guard on the junction source balance: it fails (short by 6/7) if the cell
    uses the linearized Jacobian conductance, and holds only with the
    semi-implicit sediment-flux coefficient (conductance*dz = Q_s). Complements
    the supply-only conservation tests above."""
    nseg = 30
    dx = 1000.0
    xt = dx * np.arange(1, nseg + 1.0)
    xtr = dx * np.arange(nseg + 1, 2 * nseg + 1.0)
    U = 5.0e-4 / 3.15e7                          # ~0.5 mm/yr, in m/s
    net = grlp.Network()
    net.initialize(
        x_bl=dx * (2 * nseg + 1), z_bl=0.0, S0=[1.5e-2, 1.5e-2], Q_s_0=None,
        upstream_segment_IDs=[[], [], [0, 1]],
        downstream_segment_IDs=[[2], [2], []],
        x=[xt.copy(), xt.copy(), xtr.copy()], z=[np.zeros(nseg)] * 3,
        Q=[10.0 * np.ones(nseg), 10.0 * np.ones(nseg), 20.0 * np.ones(nseg)],
        B=[100.0 * np.ones(nseg)] * 3,
    )
    for lp in net.list_of_LongProfile_objects:
        lp.set_intermittency(1.0)
        lp.set_uplift_rate(U)                    # initialize hardcodes U=0
    net.set_niter(4)
    net.get_z_lengths()
    net.evolve_threshold_width_river_network(nt=6000, dt=1e12)

    seg = net.list_of_LongProfile_objects
    kQs = seg[0].k_Qs

    conf = seg[2]
    # Balance over the confluence node's own control volume (land area A): the
    # flux out its downstream face, minus the flux entering across the junction-
    # crossing faces (each tributary's terminal node -> the confluence node),
    # equals the uplift sediment generated in that cell, (1-lambda)*U*A. The
    # tributaries' *interior* faces would span three source-bearing nodes, so the
    # crossing faces are the correct control surface.
    S_out = -(conf.z[1] - conf.z[0]) / (conf.x[1] - conf.x[0])
    outflux = kQs * (conf.Q[0] + conf.Q[1]) / 2.0 * S_out ** (7 / 6.0)
    influx = sum(
        kQs * seg[t].Q[-1]
        * ((seg[t].z[-1] - conf.z[0]) / (conf.x[0] - seg[t].x[-1])) ** (7 / 6.0)
        for t in (0, 1))
    generation = (1 - conf.lambda_p) * U * conf.land_area_around_confluence
    assert (outflux - influx) == pytest.approx(generation, rel=1e-6)


def test_three_way_confluence_conserves():
    """Three tributaries meeting at one node (an N-way junction, N > 2). The
    is_confluence assembly loops over all upstream tributaries, so it handles any
    number, but every other confluence test uses binary junctions. Sediment must
    conserve across the 3-way junction (supply-driven, no uplift)."""
    nseg = 25
    dx = 1000.0
    xt = dx * np.arange(1, nseg + 1.0)
    xtr = dx * np.arange(nseg + 1, 2 * nseg + 1.0)
    net = grlp.Network()
    net.initialize(
        x_bl=dx * (2 * nseg + 1), z_bl=0.0, S0=[1.5e-2, 1.5e-2, 1.5e-2],
        Q_s_0=None,
        upstream_segment_IDs=[[], [], [], [0, 1, 2]],
        downstream_segment_IDs=[[3], [3], [3], []],
        x=[xt.copy(), xt.copy(), xt.copy(), xtr.copy()],
        z=[np.zeros(nseg)] * 4,
        Q=[5.0 * np.ones(nseg)] * 3 + [15.0 * np.ones(nseg)],
        B=[100.0 * np.ones(nseg)] * 4,
    )
    for lp in net.list_of_LongProfile_objects:
        lp.set_intermittency(1.0)
    net.set_niter(4)
    net.get_z_lengths()
    net.evolve_threshold_width_river_network(nt=3000, dt=1e12)

    seg = net.list_of_LongProfile_objects
    kQs = seg[0].k_Qs

    def qs(lp, first=False):
        S = -np.diff(lp.z) / np.diff(lp.x)
        Qf = (lp.Q[:-1] + lp.Q[1:]) / 2.0
        f = kQs * Qf * S ** (7 / 6.0)
        return f[0] if first else f[-1]

    influx = sum(qs(seg[t]) for t in (0, 1, 2))
    outflux = qs(seg[3], first=True)
    assert outflux == pytest.approx(influx, rel=1e-6)
