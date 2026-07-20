"""
Sediment conservation at a multi-tributary confluence with discharge varying
*within* each segment.

Every other network conservation test uses uniform Q per segment, where the
land_area confluence discretization is exact. With Q varying within the
segments the current confluence does NOT conserve sediment at steady state
(~0.85% imbalance, converged -- not a transient) and the junction elevation does
not converge under grid refinement. This is a real bug in the confluence
handling, surfaced during the Step 1 de-pad; the test is xfail until the
confluence is fixed to conserve (see claude-instructions/step1-depad-resume.md).
"""
import numpy as np
import pytest

import grlp


@pytest.mark.xfail(reason="multi-tributary confluence does not conserve sediment "
                          "for within-segment-varying Q (known bug; de-pad 2c)",
                   strict=True)
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
