"""
Shared fixtures and helpers for the GRLP test suite.

The central helper is :func:`make_long_profile`, which builds a single-segment
:class:`grlp.LongProfile` with the same canonical parameter set used in
``examples/run_grlp.py`` (a reproduction of Fig. 2 of Wickert & Schildgen,
2019).  Individual tests override whichever parameters they need to probe.

Parameter conventions follow the source code:

* ``S0``    -- upstream (boundary) transport slope [-], entered as a positive
              magnitude; the bed descends downstream.
* ``P_xQ``  -- discharge--distance exponent.  The default equals
              ``P_xA * P_AQ = 7/4 * 0.7 = 49/40``, so it is consistent with the
              hydrologic constants set by ``set_hydrologic_constants``, which is
              what ``analytical_threshold_width`` reads via ``self.P_xQ``.
"""

import numpy as np
import pytest

import grlp


# A time step large enough that a few dozen steps reach steady state for the
# canonical domain, established empirically (dz/dt falls to ~machine zero).
STEADY_DT = 1.0e13
STEADY_NT = 50


def make_long_profile(
    S0=1.5e-2,
    P_xB=0.2,
    nx=90,
    dx=1000.0,
    x0=10000.0,
    k_xA=1.0,
    k_xQ=1.43e-5,
    P_xQ=7 / 4.0 * 0.7,
    k_xB=25.0,
    intermittency=0.8,
    U=0.0,
    niter=3,
    z_bl=0.0,
):
    """
    Construct and fully initialize a single-segment ``LongProfile``.

    The returned object has x, z, A, Q, B, the sediment-supply boundary
    condition (from ``S0``), and base level all set, but has **not** yet been
    evolved through time.  Call ``lp.evolve_threshold_width_river(nt, dt)`` to
    advance it.
    """
    lp = grlp.LongProfile()
    lp.set_intermittency(intermittency)
    lp.basic_constants()
    lp.bedload_lumped_constants()
    lp.set_hydrologic_constants()
    lp.set_x(dx=dx, nx=nx, x0=x0)
    lp.set_z(S0=-S0)
    lp.set_A(k_xA=k_xA)
    lp.set_Q(k_xQ=k_xQ, P_xQ=P_xQ)
    # set_Q sets the Q array but does not update self.P_xQ, which
    # analytical_threshold_width reads.  Keep them consistent so the analytical
    # power-law exponent matches the discharge field actually used.
    lp.P_xQ = P_xQ
    lp.set_B(k_xB=k_xB, P_xB=P_xB)
    lp.set_uplift_rate(U)
    lp.set_niter(niter)
    # Sediment supply consistent with the boundary slope S0.
    Qs0 = lp.k_Qs * lp.Q[0] * S0 ** (7 / 6.0)
    lp.set_Qs_input_upstream(Qs0)
    lp.set_z_bl(z_bl)
    return lp


@pytest.fixture
def long_profile_factory():
    """Expose :func:`make_long_profile` as a fixture for parametrized use."""
    return make_long_profile


@pytest.fixture
def steady_profile():
    """
    A canonical single-segment profile evolved to steady state (no uplift).

    Convergence is verified by the caller where relevant; here we simply
    integrate long enough that dz/dt is negligible for the default domain.
    """
    lp = make_long_profile()
    lp.evolve_threshold_width_river(nt=STEADY_NT, dt=STEADY_DT)
    return lp
