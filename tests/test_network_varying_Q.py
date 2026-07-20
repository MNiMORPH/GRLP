"""
A single-segment network must match the standalone solver and the analytical
solution when discharge varies downstream.

This closes a test gap that hid first-order accuracy bugs at the domain
boundaries.  Every other network test uses *uniform* discharge per segment,
where the extrapolated ghost discharge ``2*Q[0] - Q[1]`` equals ``Q[0]`` (and
likewise at the outlet) and the boundary ``dQ/dx`` term vanishes -- so how the
ghost discharge is set makes no difference and the bug is invisible.

With a downstream-increasing ``Q`` it matters.  Neither a channel head nor the
single river mouth has a tributary junction, hence no discharge discontinuity,
so the two-cell centered ``dQ/dx`` in the boundary sediment flux is well defined
and second-order.  The former constant ghosts (``Q_ext[0] = Q[0]`` upstream and
``Q_ext[-1] = Q[-1]`` downstream) collapsed it to one-sided, first-order
estimates, biasing the injected flux and the equilibrium profile by O(dx) --
~8.9 m at ``dx = 1000`` on this domain for the upstream boundary alone.

A single-segment network (one channel head straight to base level, no
confluence) is the same physics as a standalone ``LongProfile``, so:

* it must reproduce the analytical steady state (to the model's discretization
  error), and
* given *identical* inputs -- an array ``Q`` (so the standalone also extrapolates
  its ghosts rather than evaluating an analytic ``k_xQ * x**P_xQ``) -- it must
  reproduce the standalone to machine precision.
"""

import numpy as np
import pytest

import grlp


S0 = 1.5e-2
DT = 1e12
NT = 1000


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


def _standalone_from_arrays(x, Q, B, dx, niter=5):
    """Standalone LongProfile fed *array* Q/B, so it extrapolates its ghost
    discharges (set_Q: 2*Q[0]-Q[1], 2*Q[-1]-Q[-2]) exactly as the network must
    -- the apples-to-apples reference for machine-precision parity."""
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


def test_single_segment_network_matches_analytical(long_profile_factory):
    """Physical correctness: the varying-Q single-segment network tracks the
    analytical steady state.  Bites hard without the head fix (~8.9 m)."""
    a = long_profile_factory(intermittency=1.0)  # varying Q (k_xQ, P_xQ)
    x = a.x.copy()
    dx = x[1] - x[0]
    assert a.Q[1] != a.Q[0], "test must use downstream-varying Q to bite"
    a.set_z(z=np.zeros_like(x))
    a.set_z_bl(0.0)

    net = _single_segment_network(x, a.Q, a.B, dx)
    a.evolve_threshold_width_river(nt=NT, dt=DT)
    net.evolve_threshold_width_river_network(nt=NT, dt=DT)
    zN = net.list_of_LongProfile_objects[0].z
    za = a.analytical_threshold_width()  # a is evolved -> correct amplitude

    # The network extrapolates its ghost discharges (it has only array Q, not
    # the analytic k_xQ*x**P_xQ form), so it tracks both the evolved standalone
    # and the analytical to O(dx**2) -- well under a metre here; the pre-fix bug
    # shifted the amplitude by ~8.9 m.
    assert np.abs(zN - a.z).max() < 0.5
    assert np.abs(zN - za).max() < 0.5


def test_single_segment_network_identical_to_standalone(long_profile_factory):
    """Exact equivalence: with matched array inputs the single-segment network
    reproduces the standalone to machine precision.  Requires *both* boundary
    ghost fixes (upstream and downstream); either one missing leaves an O(dx)
    residual and fails this assertion."""
    ref = long_profile_factory(intermittency=1.0)
    x = ref.x.copy()
    dx = x[1] - x[0]
    Q = ref.Q.copy()
    B = ref.B.copy()

    a = _standalone_from_arrays(x, Q, B, dx)
    a.set_z(z=np.zeros_like(x))
    a.set_z_bl(0.0)
    net = _single_segment_network(x, Q, B, dx)

    a.evolve_threshold_width_river(nt=NT, dt=DT)
    net.evolve_threshold_width_river_network(nt=NT, dt=DT)
    zN = net.list_of_LongProfile_objects[0].z

    np.testing.assert_allclose(zN, a.z, rtol=0, atol=1e-9)


def test_channel_head_ghost_discharge_is_linear_extrapolation(
        long_profile_factory):
    """Direct guard on the upstream fix: with varying Q the channel-head ghost
    discharge must be the linear extrapolation 2*Q[0]-Q[1], not the old constant
    Q[0]."""
    ref = long_profile_factory(intermittency=1.0)
    x = ref.x.copy()
    dx = x[1] - x[0]
    net = _single_segment_network(x, ref.Q, ref.B, dx)
    seg = net.list_of_LongProfile_objects[0]
    assert seg.Q_ext[0][0] == pytest.approx(2 * seg.Q[0] - seg.Q[1])
    # Distinct from the former constant value (only because Q varies here).
    assert seg.Q_ext[0][0] != pytest.approx(seg.Q[0])


def test_outlet_ghost_discharge_is_linear_extrapolation(long_profile_factory):
    """Direct guard on the downstream fix: the river-mouth ghost discharge must
    be 2*Q[-1]-Q[-2], not the old constant Q[-1]."""
    ref = long_profile_factory(intermittency=1.0)
    x = ref.x.copy()
    dx = x[1] - x[0]
    net = _single_segment_network(x, ref.Q, ref.B, dx)
    seg = net.list_of_LongProfile_objects[0]
    assert seg.Q_ext[0][-1] == pytest.approx(2 * seg.Q[-1] - seg.Q[-2])
    assert seg.Q_ext[0][-1] != pytest.approx(seg.Q[-1])


@pytest.mark.slow
def test_single_segment_network_second_order_convergence(long_profile_factory):
    """The boundary sediment flux must be second-order: the error against the
    analytical steady state falls ~4x per halving of dx.  The pre-fix constant
    ghost was first-order (error only halved), so this fails without the fix.
    """
    x0 = 10000.0
    L = 88000.0
    errs = []
    for nx in (45, 89, 177):
        dx = L / (nx - 1)
        a = long_profile_factory(intermittency=1.0, nx=nx, dx=dx, x0=x0, niter=5)
        x = a.x.copy()
        a.set_z(z=np.zeros_like(x))
        a.set_z_bl(0.0)
        net = _single_segment_network(x, a.Q, a.B, dx)
        a.evolve_threshold_width_river(nt=1200, dt=DT)
        net.evolve_threshold_width_river_network(nt=1200, dt=DT)
        zN = net.list_of_LongProfile_objects[0].z
        za = a.analytical_threshold_width()
        errs.append(np.abs(zN - za).max())
    orders = [np.log2(errs[i] / errs[i + 1]) for i in range(len(errs) - 1)]
    assert min(orders) > 1.8, "expected second-order; got orders=%r (errs=%r)" % (
        orders, errs)


def test_ghost_discharge_linear_on_all_heads_of_a_confluence():
    """The upstream fix must apply to every channel head, not just a lone
    segment: build a two-tributary confluence with downstream-increasing Q in
    each tributary and check both heads get the linear ghost discharge."""
    d = 1000.0
    xt = d * np.arange(1, 6.0)          # two 5-node tributaries
    xtr = d * np.arange(6, 11.0)        # 5-node trunk below the junction
    Qt = 2.0 + 0.1 * np.arange(5)       # increasing discharge along a tributary
    Qtrunk = 4.0 + 0.2 * np.arange(5)   # trunk carries ~the summed, still rising
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
    heads = net.list_of_channel_head_segment_IDs
    assert len(heads) == 2
    for head in heads:
        seg = net.list_of_LongProfile_objects[head]
        for Q_ext_array in seg.Q_ext:
            assert Q_ext_array[0] == pytest.approx(2 * seg.Q[0] - seg.Q[1])
            assert Q_ext_array[0] != pytest.approx(seg.Q[0])
