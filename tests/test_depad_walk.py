"""
The neighbor-walking global assembler: structural invariants.

``Network.assemble_by_walking`` builds the global LHS/RHS by walking the topology
to each node's real neighbor. Splitting a reach into a 1-into-1 chain must be an
exact no-op: the confluence node's neighbor walk crosses the segment boundary and
gets the ordinary interior stencil (vs. the ~0.8 m/junction error of the old
``land_area`` code). This is the de-pad's discharge-continuous junction fix.

(The walker's bit-for-bit agreement with the former padded single-segment solver
is now covered end-to-end by the characterization goldens and
``test_single_segment_unification``; that padded assembler -- ``build_matrices``
-- has since been removed, so there is no longer a separate reference to compare
the assembled matrix against.)
"""

import numpy as np

import grlp


DT = 1e12


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
