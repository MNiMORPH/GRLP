"""
Parity harness for unifying the single-segment solver with the network walker.

Goal: a single segment run as a *one-edge network* must reproduce the standalone
``LongProfile`` solver bit-for-bit, across the range of single-segment
configurations. Once that holds, ``LongProfile.evolve_threshold_width_river`` can
become a thin wrapper around a one-edge ``Network`` and the padded single-segment
path can be retired.

Both solvers are fed the *same array* ``Q`` so both extrapolate their ghost
discharges linearly; the standalone ``set_Q(k_xQ, P_xQ)`` analytic-power-law path
would differ at O(dx**2) and is deliberately avoided here.

Known gap, deliberately not covered: Sternberg gravel loss, which the walker does
not yet apply in its Picard loop.
"""

import numpy as np
import pytest

import grlp


S0 = 1.5e-2
DX = 1000.0
NX = 40
X0 = 10000.0
NT = 600
DT = 1e12


def _x():
    return X0 + DX * np.arange(NX)


def _standalone(Q, B, U=0.0, I=1.0):
    lp = grlp.LongProfile()
    lp.set_intermittency(I)
    lp.basic_constants()
    lp.bedload_lumped_constants()
    lp.set_hydrologic_constants()
    lp.set_x(dx=DX, nx=NX, x0=X0)
    lp.set_z(S0=-S0)
    lp.set_A(k_xA=1.0)
    lp.set_Q(Q=Q.copy())
    lp.set_B(B=B.copy())
    lp.set_uplift_rate(U)
    lp.set_niter(3)
    lp.set_Qs_input_upstream(lp.k_Qs * lp.Q[0] * S0 ** (7 / 6.0))
    lp.set_z_bl(0.0)
    lp.set_z(z=np.zeros(NX))
    return lp


def _one_edge_network(Q, B, U=0.0, I=1.0):
    x = _x()
    net = grlp.Network()
    net.initialize(
        x_bl=x[-1] + DX, z_bl=0.0, S0=[S0], Q_s_0=None,
        upstream_segment_IDs=[[]], downstream_segment_IDs=[[]],
        x=[x.copy()], z=[np.zeros(NX)], Q=[Q.copy()], B=[B.copy()],
    )
    seg = net.list_of_LongProfile_objects[0]
    seg.set_intermittency(I)
    seg.set_uplift_rate(U)   # initialize hardcodes 0; set it here for parity
    net.set_niter(3)
    net.get_z_lengths()
    return net


def _uniform_B(x):
    return 100.0 * np.ones_like(x)


def _powerlaw_B(x):
    return 25.0 * x ** 0.2


def _arbitrary_B(x):
    return 40.0 + 30.0 * (x - x[0]) / (x[-1] - x[0])


CONFIGS = {
    "uniform_B":    dict(B=_uniform_B),
    "powerlaw_B":   dict(B=_powerlaw_B),
    "arbitrary_B":  dict(B=_arbitrary_B),
    "intermittency": dict(B=_uniform_B, I=0.6),
    "uplift":       dict(B=_uniform_B, U=10.0e-3 / 3.15e7),
}


@pytest.mark.parametrize("name", list(CONFIGS))
def test_one_edge_network_matches_standalone(name):
    cfg = CONFIGS[name]
    x = _x()
    Q = 2.0 + 0.05 * np.arange(NX)          # downstream-increasing: ghosts bite
    B = cfg["B"](x)
    U = cfg.get("U", 0.0)
    I = cfg.get("I", 1.0)

    a = _standalone(Q, B, U=U, I=I)
    net = _one_edge_network(Q, B, U=U, I=I)

    a.evolve_threshold_width_river(nt=NT, dt=DT)
    net.evolve_threshold_width_river_network(nt=NT, dt=DT)
    zN = net.list_of_LongProfile_objects[0].z

    np.testing.assert_allclose(zN, a.z, rtol=0, atol=1e-9)
