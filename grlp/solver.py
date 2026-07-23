"""
The GRLP network solver, extracted from the Network class.

These functions evolve a river-network long-profile model by walking the topology
and assembling a single global sparse system. They operate on a *duck-typed*
network -- any object exposing ``segments``,
``list_of_segment_lengths``, ``niter``, ``t``, and segments carrying the
per-segment physics -- so this module imports neither the ``Network`` nor the
``Segment`` class. That keeps the dependency one-way (grlp.py -> solver.py)
and lets a lone ``Segment`` be solved as the one-edge network it is: one
solver path for every case.
"""

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve


def assemble(net, dt):
    """
    De-padded global assembly (Step 1 de-pad).

    Build the global LHS matrix and RHS by *walking the topology* to each
    node's upstream/downstream neighbor, instead of reading the padded
    ``z_ext``/``Q_ext`` ghost arrays. The per-node stencil coefficients come
    from :meth:`Segment.build_LHS_coeff_C0` -- only neighbor lookup
    changes -- so for a single segment this reproduces the standalone solver
    bit-for-bit. Channel heads apply the sediment-flux Neumann boundary
    condition, the outlet the base-level Dirichlet condition; at a
    single-upstream junction the confluence node reaches across the segment
    boundary and so gets the ordinary interior stencil (fixing the
    first-order ``land_area`` junction handling in the 1-into-1 case).

    Multi-tributary confluences are not yet handled here (they will delegate
    to the existing junction code, or later a flux-balance cell); a node with
    more than one upstream segment raises ``NotImplementedError``.

    Returns ``(LHSmatrix, RHS)`` over the global node vector, ordered by
    segment as in ``list_of_segment_lengths``. Additive: not yet wired into
    the evolve loop.
    """
    segs = net.segments
    lengths = list(net.list_of_segment_lengths)
    starts = np.cumsum([0] + lengths)[:-1]
    n = int(np.sum(lengths))
    # The three-node junction cell reaches to the confluence's second
    # interior node and to each tributary's second-to-last node, so segments
    # adjacent to a multi-tributary confluence must be long enough. Fail
    # clearly rather than with an IndexError.
    for seg in segs:
        if len(seg.upstream_segment_IDs) > 1:
            if lengths[seg.ID] < 3:
                raise ValueError(
                    "Walking solver: confluence segment %d needs >= 3 nodes "
                    "(has %d)." % (seg.ID, lengths[seg.ID]))
            for ID in seg.upstream_segment_IDs:
                if lengths[ID] < 2:
                    raise ValueError(
                        "Walking solver: tributary segment %d into confluence "
                        "%d needs >= 2 nodes (has %d)."
                        % (ID, seg.ID, lengths[ID]))
    for seg in segs:
        seg.build_LHS_coeff_C0(dt=dt)
    rows = []
    cols = []
    vals = []
    RHS = np.zeros(n)
    for seg in segs:
        s = seg.ID
        offset = starts[s]
        L = lengths[s]
        # per-node source term (RHS additions: ssd + fining/subsidence + U)
        src = (np.asarray(seg.ssd)
               + np.asarray(seg.downstream_fining_subsidence_equivalent)
               + np.asarray(seg.U)) * dt
        src = np.broadcast_to(src, (L,))
        # RHS uses the start-of-step elevation (zold) during Picard
        # iteration; the coefficient (C1) uses the current iterate seg.z.
        # When zold is unset (static assembly, e.g. tests) it equals seg.z.
        z_rhs = getattr(seg, "zold", None)
        if z_rhs is None or np.size(z_rhs) != L:
            z_rhs = seg.z
        for i in range(L):
            # g: this node's index in the flattened global node vector
            # (segment offset + local index i)
            g = offset + i
            # ===== multi-tributary junction: shared-flux three-node cell ====
            # Conservation is by construction: each junction FACE carries one
            # shared conductance used identically by both adjacent nodes, so
            # the sediment flux conductance*(z_up - z_down) is single-valued.
            # conductance/land_area matches the interior coupling magnitude,
            # so a junction-adjacent node is consistent with the ordinary
            # interior stencil on its other face.
            def _face_conductance(z_up, z_down, Q, x_up, x_down, C0):
                # Semi-implicit sediment-flux coefficient: conductance * dz
                # = k_Qs * I * Q * (S/sinuosity)**(7/6) = Q_s at Picard
                # convergence (S**(1/6) from the current iterate times one
                # implicit power of S). This is the flux coefficient itself,
                # NOT the linearized (7/6) Jacobian dQ_s/dS: using the
                # Jacobian converges the junction cell to the *linearized*
                # balance, which under-applies distributed sources (uplift,
                # ssd, gravel loss) by exactly 6/7 at every confluence.
                L_face = x_down - x_up
                return C0 * Q \
                       * (np.abs(z_up - z_down) / L_face) ** (1 / 6.) \
                       / L_face
            is_confluence = (i == 0 and len(seg.upstream_segment_IDs) > 1)
            down_is_confluence = (
                i == L - 1 and seg.downstream_segment_IDs
                and len(segs[seg.downstream_segment_IDs[0]]
                        .upstream_segment_IDs) > 1)
            up_is_confluence = (i == 1 and len(seg.upstream_segment_IDs) > 1)
            if is_confluence:
                A_confluence = seg.land_area_around_confluence
                conductance_down = _face_conductance(
                    seg.z[0], seg.z[1], 0.5 * (seg.Q[0] + seg.Q[1]),
                    seg.x[0], seg.x[1], seg.C0)
                conductance_sum = conductance_down
                rows.append(g)
                cols.append(g + 1)
                vals.append(-conductance_down / A_confluence)
                for ID in seg.upstream_segment_IDs:
                    upseg = segs[ID]
                    upseg_g = starts[ID] + lengths[ID] - 1
                    conductance_upseg = _face_conductance(
                        upseg.z[-1], seg.z[0], upseg.Q[-1],
                        upseg.x[-1], seg.x[0], seg.C0)
                    conductance_sum += conductance_upseg
                    rows.append(g)
                    cols.append(upseg_g)
                    vals.append(-conductance_upseg / A_confluence)
                rows.append(g)
                cols.append(g)
                vals.append(1. + conductance_sum / A_confluence)
                RHS[g] = z_rhs[i] + src[i]
                continue
            if down_is_confluence:
                downseg = segs[seg.downstream_segment_IDs[0]]
                downseg_g = starts[downseg.ID]
                land_area = seg.B[-1] * 0.5 * ((seg.x[-1] - seg.x[-2])
                                              + (downseg.x[0] - seg.x[-1]))
                conductance_downseg = _face_conductance(
                    seg.z[-1], downseg.z[0], seg.Q[-1],
                    seg.x[-1], downseg.x[0], downseg.C0)  # shared with confluence
                conductance_up = _face_conductance(
                    seg.z[-2], seg.z[-1], 0.5 * (seg.Q[-2] + seg.Q[-1]),
                    seg.x[-2], seg.x[-1], seg.C0)
                rows.append(g)
                cols.append(g - 1)
                vals.append(-conductance_up / land_area)
                rows.append(g)
                cols.append(downseg_g)
                vals.append(-conductance_downseg / land_area)
                rows.append(g)
                cols.append(g)
                vals.append(
                    1. + (conductance_up + conductance_downseg) / land_area)
                RHS[g] = z_rhs[i] + src[i]
                continue
            if up_is_confluence:
                land_area = seg.B[1] * 0.5 * ((seg.x[1] - seg.x[0])
                                             + (seg.x[2] - seg.x[1]))
                conductance_up = _face_conductance(
                    seg.z[0], seg.z[1], 0.5 * (seg.Q[0] + seg.Q[1]),
                    seg.x[0], seg.x[1], seg.C0)  # shared with confluence
                conductance_down = _face_conductance(
                    seg.z[1], seg.z[2], 0.5 * (seg.Q[1] + seg.Q[2]),
                    seg.x[1], seg.x[2], seg.C0)
                rows.append(g)
                cols.append(g - 1)
                vals.append(-conductance_up / land_area)
                rows.append(g)
                cols.append(g + 1)
                vals.append(-conductance_down / land_area)
                rows.append(g)
                cols.append(g)
                vals.append(
                    1. + (conductance_up + conductance_down) / land_area)
                RHS[g] = z_rhs[i] + src[i]
                continue
            # --- upstream neighbor (or head ghost) ---
            if i > 0:
                up_g = g - 1
                z_up = seg.z[i - 1]
                x_up = seg.x[i - 1]
                Q_up = seg.Q[i - 1]
                is_head = False
            elif len(seg.upstream_segment_IDs) == 0:
                is_head = True
                up_g = None
                x_up = 2 * seg.x[0] - seg.x[1]
                z_up = seg.z[0] + seg.S0 * (seg.x[0] - x_up)
                if seg.Q_ghost_upstream is not None:
                    Q_up = seg.Q_ghost_upstream
                else:
                    Q_up = 2 * seg.Q[0] - seg.Q[1]
            elif len(seg.upstream_segment_IDs) == 1:
                is_head = False
                upseg = segs[seg.upstream_segment_IDs[0]]
                up_g = starts[upseg.ID] + lengths[upseg.ID] - 1
                z_up = upseg.z[-1]
                x_up = upseg.x[-1]
                Q_up = upseg.Q[-1]
            else:
                raise NotImplementedError(
                    "assemble_by_walking: multi-tributary confluence "
                    "(segment %d has %d upstream segments) not yet handled"
                    % (s, len(seg.upstream_segment_IDs)))
            # --- downstream neighbor (or outlet ghost) ---
            if i < L - 1:
                down_g = g + 1
                z_down = seg.z[i + 1]
                x_down = seg.x[i + 1]
                Q_down = seg.Q[i + 1]
                is_outlet = False
            elif len(seg.downstream_segment_IDs) == 0:
                is_outlet = True
                down_g = None
                # Base-level node: elevation z_bl at position
                # x_ghost_downstream (settable via set_x_bl to move the
                # mouth in x; defaults to one cell beyond the last node).
                if seg.x_ghost_downstream is not None:
                    x_down = seg.x_ghost_downstream
                else:
                    x_down = 2 * seg.x[-1] - seg.x[-2]
                z_down = seg.z_bl
                if seg.Q_ghost_downstream is not None:
                    Q_down = seg.Q_ghost_downstream
                else:
                    Q_down = 2 * seg.Q[-1] - seg.Q[-2]
            else:
                is_outlet = False
                downseg = segs[seg.downstream_segment_IDs[0]]
                down_g = starts[downseg.ID]
                z_down = downseg.z[0]
                x_down = downseg.x[0]
                Q_down = downseg.Q[0]
            # --- stencil (coefficients from build_LHS_coeff_C0) ---
            dx_up = seg.x[i] - x_up
            dx_down = x_down - seg.x[i]
            dx_2cell = x_down - x_up
            dQ_2cell = Q_down - Q_up
            S = np.abs(z_down - z_up) / dx_2cell
            C1 = seg.C0 * S ** (1 / 6.) * seg.Q[i] / seg.B[i]
            center = -C1 / dx_2cell * (7 / 3. * (-1 / dx_up - 1 / dx_down)) + 1.
            left = -C1 / dx_2cell * (7 / 3. / dx_up - dQ_2cell / seg.Q[i] / dx_2cell)
            right = -C1 / dx_2cell * (7 / 3. / dx_down + dQ_2cell / seg.Q[i] / dx_2cell)
            rhs_g = z_rhs[i] + src[i]
            if is_head:                       # set_bcl_Neumann
                right = -C1 / dx_2cell * 7 / 3. * (1 / dx_up + 1 / dx_down)
                rhs_g += seg.S0 * C1 * (7 / 3. / dx_up - dQ_2cell / seg.Q[i] / dx_2cell)
            if is_outlet:                     # set_bcr_Dirichlet
                rhs_g += z_down * C1 / dx_2cell * (
                    7 / 3. * (1 / (seg.x[-1] - seg.x[-2])
                              + 1 / (x_down - seg.x[-1])) / 2.
                    + dQ_2cell / seg.Q[i] / dx_2cell)
            rows.append(g)
            cols.append(g)
            vals.append(center)
            if up_g is not None:
                rows.append(g)
                cols.append(up_g)
                vals.append(left)
            if down_g is not None:
                rows.append(g)
                cols.append(down_g)
                vals.append(right)
            RHS[g] = rhs_g
    LHSmatrix = sparse.csr_matrix((vals, (rows, cols)), shape=(n, n))
    return LHSmatrix, RHS

def evolve(net, nt, dt):
    """
    Time-step the network through the de-padded walking assembler
    (:meth:`assemble_by_walking`). Used for networks with no multi-tributary
    confluence (single segments and 1-into-1 chains); the walker handles
    those exactly and fixes the first-order junction. Proper Picard: the RHS
    is frozen at the start-of-step elevation (``zold``) while the coefficient
    relinearizes on the current iterate each iteration.
    """
    net.dt = dt
    segs = net.segments
    lengths = list(net.list_of_segment_lengths)
    starts = np.cumsum([0] + lengths)[:-1]
    # Sternberg gravel loss enters as a distributed sink that depends on the
    # (evolving) sediment discharge, so it must be relinearized each Picard
    # iteration. Skip the recompute entirely when no segment sets it.
    gravel_loss_active = any(
        seg.gravel_fractional_loss_per_km is not None for seg in segs)
    for ti in range(int(nt)):
        for seg in segs:
            seg.zold = seg.z.copy()
        for _ in range(int(net.niter)):
            if gravel_loss_active:
                update_gravel_loss(net)
            LHS, RHS = assemble(net, dt)
            out = spsolve(sparse.csr_matrix(LHS), RHS)
            for seg in segs:
                s = seg.ID
                seg.z = out[starts[s]:starts[s] + lengths[s]]
        net.t += dt
        for seg in segs:
            seg.t = net.t
            seg.dz_dt = (seg.z - seg.zold) / dt

def update_gravel_loss(net):
    """
    Recompute the Sternberg gravel-abrasion sink on every segment that sets
    it, from the current profile, so it stays consistent within the
    semi-implicit Picard iteration (the sink depends on Q_s, which depends
    on the evolving z). compute_Q_s walks the topology, so the sediment
    discharge -- and hence the abrasion -- accumulates through the network:
    a grain abrades along its whole downstream path, across confluences.
    """
    net.compute_Q_s()
    for seg in net.segments:
        if seg.gravel_fractional_loss_per_km is not None:
            seg.downstream_fining_subsidence_equivalent = \
                - seg.gravel_fractional_loss_per_km / 1000. * seg.Q_s \
                / ( (1 - seg.lambda_p) * seg.B )

