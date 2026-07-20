"""
Flow intermittency ``I`` in the networked solver.

Intermittency is the fraction of time the channel is geomorphically active.  It
enters the model only through the lumped coefficient ``C0`` (i.e. it multiplies
the time step), so it sets the *rate* of landscape evolution, not the
equilibrium form:

* the steady-state profile is independent of ``I``;
* halving ``I`` and doubling ``dt`` give identical transients (equal ``I*dt``);
* ``I`` must **not** appear in the ``Q_s_0 -> S0`` boundary conversion, which is
  a purely geometric inversion of the transport law.

These tests would have caught two bugs:

1. ``Network.set_intermittency`` assigned the value to the segment's *method*
   attribute instead of calling it, so segment intermittency never changed (and
   the method was clobbered).
2. ``update_z_ext_external_upstream`` divided by ``I`` and used the wrong sign
   when converting ``Q_s_0`` to ``S0``, so a network driven by sediment supply
   diverged from the equivalent standalone once ``I != 1``.
"""

import numpy as np
import pytest

import grlp
from network_helpers import NETWORK_TOPOLOGIES, build_network


S0 = 1.5e-2
DT = 1e12


def _confluence(evolve=False):
    spec = NETWORK_TOPOLOGIES["symmetric_confluence"]
    return build_network(spec["x"], spec["Q"], spec["up"], spec["down"],
                         spec["x_bl"], evolve=evolve)


def _single_segment_network(x, Q, B, dx, I, niter=5):
    net = grlp.Network()
    net.initialize(
        x_bl=x[-1] + dx, z_bl=0.0, S0=[S0], Q_s_0=None,
        upstream_segment_IDs=[[]], downstream_segment_IDs=[[]],
        x=[x.copy()], z=[np.zeros_like(x)], Q=[Q.copy()], B=[B.copy()],
    )
    net.set_intermittency(I)
    net.set_niter(niter)
    net.get_z_lengths()
    return net


# --- Network.set_intermittency actually sets the segments -------------------

def test_network_set_intermittency_scalar_takes_effect():
    net = _confluence()
    net.set_intermittency(0.6)
    for seg in net.list_of_LongProfile_objects:
        assert seg.intermittency == 0.6
        # The method must survive (the bug overwrote it with a float).
        assert callable(seg.set_intermittency)


def test_network_set_intermittency_per_segment_list():
    net = _confluence()
    net.set_intermittency([0.5, 0.7, 0.9])
    got = [seg.intermittency for seg in net.list_of_LongProfile_objects]
    assert got == [0.5, 0.7, 0.9]


# --- Q_s_0 -> S0 conversion: no intermittency, correct sign ----------------

@pytest.mark.parametrize("I", [1.0, 0.8, 0.5])
def test_Qs0_to_S0_omits_intermittency_and_is_positive(long_profile_factory, I):
    a = long_profile_factory(intermittency=I)
    x = a.x.copy()
    dx = x[1] - x[0]
    Qs0 = a.k_Qs * a.Q[0] * S0 ** (7 / 6.0)
    a.set_Qs_input_upstream(Qs0)  # standalone reference: no I, +sign

    net = grlp.Network()
    net.initialize(
        x_bl=x[-1] + dx, z_bl=0.0, S0=None, Q_s_0=float(Qs0),
        upstream_segment_IDs=[[]], downstream_segment_IDs=[[]],
        x=[x.copy()], z=[np.zeros_like(x)], Q=[a.Q.copy()], B=[a.B.copy()],
    )
    net.set_intermittency(I)
    net.update_z_ext_external_upstream(Q_s_0=float(Qs0))
    seg = net.list_of_LongProfile_objects[0]

    assert seg.S0 == pytest.approx(a.S0)   # equals standalone, independent of I
    assert seg.S0 > 0                       # correct sign for a descending river


def test_network_driven_by_Qs0_matches_standalone(long_profile_factory):
    """A network driven by sediment supply reproduces the standalone driven by
    the same supply (to the known ghost-extrapolation residual), for I != 1."""
    I = 0.5
    a = long_profile_factory(intermittency=I)
    x = a.x.copy()
    dx = x[1] - x[0]
    Qs0 = a.k_Qs * a.Q[0] * S0 ** (7 / 6.0)
    a.set_Qs_input_upstream(Qs0)
    a.set_z(z=np.zeros_like(x))
    a.set_z_bl(0.0)

    net = grlp.Network()
    net.initialize(
        x_bl=x[-1] + dx, z_bl=0.0, S0=None, Q_s_0=float(Qs0),
        upstream_segment_IDs=[[]], downstream_segment_IDs=[[]],
        x=[x.copy()], z=[np.zeros_like(x)], Q=[a.Q.copy()], B=[a.B.copy()],
    )
    net.set_intermittency(I)
    net.set_niter(5)
    net.get_z_lengths()

    a.evolve_threshold_width_river(nt=500, dt=DT)
    net.evolve_threshold_width_river_network(nt=500, dt=DT)
    zN = net.list_of_LongProfile_objects[0].z
    # The bug gave hundreds of metres (wrong sign/magnitude); the residual here
    # is only the array-vs-analytic ghost-extrapolation difference (~0.2 m).
    assert np.abs(zN - a.z).max() < 0.5


# --- physics: rate vs. form -------------------------------------------------

def test_steady_state_independent_of_intermittency(long_profile_factory):
    ref = long_profile_factory()
    x = ref.x.copy()
    dx = x[1] - x[0]
    z_hi = _single_segment_network(x, ref.Q, ref.B, dx, I=1.0)
    z_lo = _single_segment_network(x, ref.Q, ref.B, dx, I=0.4)
    z_hi.evolve_threshold_width_river_network(nt=500, dt=DT)
    z_lo.evolve_threshold_width_river_network(nt=500, dt=DT)
    a = z_hi.list_of_LongProfile_objects[0].z
    b = z_lo.list_of_LongProfile_objects[0].z
    np.testing.assert_allclose(a, b, rtol=0, atol=1e-9)


def test_evolution_rate_scales_with_intermittency(long_profile_factory):
    """I enters only through C0 (~ I*dt), so halving I and doubling dt yields an
    identical transient."""
    ref = long_profile_factory()
    x = ref.x.copy()
    dx = x[1] - x[0]
    fast = _single_segment_network(x, ref.Q, ref.B, dx, I=1.0)
    slow = _single_segment_network(x, ref.Q, ref.B, dx, I=0.5)
    fast.evolve_threshold_width_river_network(nt=3, dt=5e11)
    slow.evolve_threshold_width_river_network(nt=3, dt=1e12)  # 2x dt, 0.5x I
    a = fast.list_of_LongProfile_objects[0].z
    b = slow.list_of_LongProfile_objects[0].z
    np.testing.assert_allclose(a, b, rtol=0, atol=1e-9)
