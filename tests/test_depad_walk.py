"""
Step 1 de-pad: the neighbor-walking global assembler must reproduce the
single-segment solver bit-for-bit.

``Network.assemble_by_walking`` builds the global LHS/RHS by walking the topology
to each node's real neighbor instead of reading the padded ``z_ext``/``Q_ext``
ghosts. For a single-segment network it must be identical to
``LongProfile.build_matrices`` -- validated against an *array-Q* standalone,
which uses the same linear ghost-discharge convention (2*Q[0]-Q[1],
2*Q[-1]-Q[-2]) as the network (not the analytic k_xQ ghost of make_long_profile).
"""

import numpy as np

import grlp


DT = 1e12


def _array_Q_standalone_matrices(x, Q, B, S0, zc):
    """build_matrices output for a standalone fed array Q/B, assembled at zc."""
    lp = grlp.LongProfile()
    lp.set_intermittency(1.0)
    lp.basic_constants()
    lp.bedload_lumped_constants()
    lp.set_hydrologic_constants()
    lp.set_x(dx=x[1] - x[0], nx=len(x), x0=x[0])
    lp.set_z(S0=-S0)
    lp.set_A(k_xA=1.0)
    lp.set_Q(Q=Q.copy())
    lp.set_B(B=B.copy())
    lp.set_uplift_rate(0.0)
    lp.set_niter(5)
    lp.set_Qs_input_upstream(lp.k_Qs * lp.Q[0] * S0 ** (7 / 6.0))
    lp.set_z(z=zc.copy())
    lp.set_z_bl(0.0)
    lp.z_ext[1:-1] = zc            # evolve normally syncs this; do it by hand
    lp.build_LHS_coeff_C0(dt=DT)
    lp.zold = zc.copy()
    lp.build_matrices()
    return lp.LHSmatrix.toarray(), lp.RHS.copy()


def test_walk_reproduces_build_matrices_single_segment(long_profile_factory):
    ref = long_profile_factory(intermittency=1.0)     # varying Q (k_xQ, P_xQ)
    x = ref.x.copy(); Q = ref.Q.copy(); B = ref.B.copy()
    n = len(x); dx = x[1] - x[0]; S0 = 1.5e-2
    zc = np.linspace(30.0, 1.0, n)                     # arbitrary smooth profile

    Mref, Rref = _array_Q_standalone_matrices(x, Q, B, S0, zc)

    net = grlp.Network()
    net.initialize(
        x_bl=x[-1] + dx, z_bl=0.0, S0=[S0], Q_s_0=None,
        upstream_segment_IDs=[[]], downstream_segment_IDs=[[]],
        x=[x.copy()], z=[zc.copy()], Q=[Q.copy()], B=[B.copy()],
    )
    net.list_of_LongProfile_objects[0].set_intermittency(1.0)
    net.set_niter(5)
    net.get_z_lengths()

    Mw, Rw = net.assemble_by_walking(DT)

    np.testing.assert_allclose(Mw.toarray(), Mref, rtol=0, atol=1e-9)
    np.testing.assert_allclose(Rw, Rref, rtol=0, atol=1e-9)


def test_walk_1into1_chain_equals_single_segment(long_profile_factory):
    """A 1-into-1 chain assembled by the walker equals the single segment: the
    confluence node's neighbor walk crosses the segment boundary and gets the
    ordinary interior stencil, so splitting a reach is an exact no-op (vs. the
    ~0.8 m/junction error of the land_area code). This is the de-pad's
    discharge-continuous junction fix."""
    ref = long_profile_factory(intermittency=1.0)
    x = ref.x.copy(); Q = ref.Q.copy(); B = ref.B.copy()
    n = len(x); dx = x[1] - x[0]; S0 = 1.5e-2
    zc = np.linspace(30.0, 1.0, n)

    def _walk(net):
        for lp in net.list_of_LongProfile_objects:
            lp.set_intermittency(1.0)
        net.set_niter(5); net.get_z_lengths()
        return net.assemble_by_walking(DT)

    single = grlp.Network()
    single.initialize(
        x_bl=x[-1] + dx, z_bl=0.0, S0=[S0], Q_s_0=None,
        upstream_segment_IDs=[[]], downstream_segment_IDs=[[]],
        x=[x.copy()], z=[zc.copy()], Q=[Q.copy()], B=[B.copy()])
    Ms, Rs = _walk(single)

    k = n // 2
    chain = grlp.Network()
    chain.initialize(
        x_bl=x[-1] + dx, z_bl=0.0, S0=[S0], Q_s_0=None,
        upstream_segment_IDs=[[], [0]], downstream_segment_IDs=[[1], []],
        x=[x[:k].copy(), x[k:].copy()], z=[zc[:k].copy(), zc[k:].copy()],
        Q=[Q[:k].copy(), Q[k:].copy()], B=[B[:k].copy(), B[k:].copy()])
    Mc, Rc = _walk(chain)

    np.testing.assert_allclose(Mc.toarray(), Ms.toarray(), rtol=0, atol=1e-9)
    np.testing.assert_allclose(Rc, Rs, rtol=0, atol=1e-9)
