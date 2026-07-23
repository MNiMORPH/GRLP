"""
The GRLP network solver, extracted from the Network class.

These functions evolve a river-network long-profile model by walking the topology
and assembling a single global sparse system. They operate on a *duck-typed*
network -- any object exposing ``list_of_LongProfile_objects``,
``list_of_segment_lengths``, ``niter``, ``t``, and segments carrying the
per-segment physics -- so this module imports neither the ``Network`` nor the
``LongProfile`` class. That keeps the dependency one-way (grlp.py -> solver.py)
and lets a lone ``LongProfile`` be solved as the one-edge network it is: one
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
    from :meth:`LongProfile.build_LHS_coeff_C0` -- only neighbor lookup
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
    segs = net.list_of_LongProfile_objects
    lengths = list(net.list_of_segment_lengths)
    starts = np.cumsum([0] + lengths)[:-1]
    n = int(np.sum(lengths))
    # The three-node junction cell reaches to the confluence's second
    # interior node and to each tributary's second-to-last node, so segments
    # adjacent to a multi-tributary confluence must be long enough. Fail
    # clearly rather than with an IndexError.
    for lp in segs:
        if len(lp.upstream_segment_IDs) > 1:
            if lengths[lp.ID] < 3:
                raise ValueError(
                    "Walking solver: confluence segment %d needs >= 3 nodes "
                    "(has %d)." % (lp.ID, lengths[lp.ID]))
            for ID in lp.upstream_segment_IDs:
                if lengths[ID] < 2:
                    raise ValueError(
                        "Walking solver: tributary segment %d into confluence "
                        "%d needs >= 2 nodes (has %d)."
                        % (ID, lp.ID, lengths[ID]))
    for lp in segs:
        lp.build_LHS_coeff_C0(dt=dt)
    rows = []
    cols = []
    vals = []
    RHS = np.zeros(n)
    for lp in segs:
        s = lp.ID
        offset = starts[s]
        L = lengths[s]
        # per-node source term (RHS additions: ssd + fining/subsidence + U)
        src = (np.asarray(lp.ssd)
               + np.asarray(lp.downstream_fining_subsidence_equivalent)
               + np.asarray(lp.U)) * dt
        src = np.broadcast_to(src, (L,))
        # RHS uses the start-of-step elevation (zold) during Picard
        # iteration; the coefficient (C1) uses the current iterate lp.z.
        # When zold is unset (static assembly, e.g. tests) it equals lp.z.
        z_rhs = getattr(lp, "zold", None)
        if z_rhs is None or np.size(z_rhs) != L:
            z_rhs = lp.z
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
            is_confluence = (i == 0 and len(lp.upstream_segment_IDs) > 1)
            down_is_confluence = (
                i == L - 1 and lp.downstream_segment_IDs
                and len(segs[lp.downstream_segment_IDs[0]]
                        .upstream_segment_IDs) > 1)
            up_is_confluence = (i == 1 and len(lp.upstream_segment_IDs) > 1)
            if is_confluence:
                A_confluence = lp.land_area_around_confluence
                conductance_down = _face_conductance(
                    lp.z[0], lp.z[1], 0.5 * (lp.Q[0] + lp.Q[1]),
                    lp.x[0], lp.x[1], lp.C0)
                conductance_sum = conductance_down
                rows.append(g)
                cols.append(g + 1)
                vals.append(-conductance_down / A_confluence)
                for ID in lp.upstream_segment_IDs:
                    upseg = segs[ID]
                    upseg_g = starts[ID] + lengths[ID] - 1
                    conductance_upseg = _face_conductance(
                        upseg.z[-1], lp.z[0], upseg.Q[-1],
                        upseg.x[-1], lp.x[0], lp.C0)
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
                downseg = segs[lp.downstream_segment_IDs[0]]
                downseg_g = starts[downseg.ID]
                land_area = lp.B[-1] * 0.5 * ((lp.x[-1] - lp.x[-2])
                                              + (downseg.x[0] - lp.x[-1]))
                conductance_downseg = _face_conductance(
                    lp.z[-1], downseg.z[0], lp.Q[-1],
                    lp.x[-1], downseg.x[0], downseg.C0)  # shared with confluence
                conductance_up = _face_conductance(
                    lp.z[-2], lp.z[-1], 0.5 * (lp.Q[-2] + lp.Q[-1]),
                    lp.x[-2], lp.x[-1], lp.C0)
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
                land_area = lp.B[1] * 0.5 * ((lp.x[1] - lp.x[0])
                                             + (lp.x[2] - lp.x[1]))
                conductance_up = _face_conductance(
                    lp.z[0], lp.z[1], 0.5 * (lp.Q[0] + lp.Q[1]),
                    lp.x[0], lp.x[1], lp.C0)  # shared with confluence
                conductance_down = _face_conductance(
                    lp.z[1], lp.z[2], 0.5 * (lp.Q[1] + lp.Q[2]),
                    lp.x[1], lp.x[2], lp.C0)
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
                z_up = lp.z[i - 1]
                x_up = lp.x[i - 1]
                Q_up = lp.Q[i - 1]
                is_head = False
            elif len(lp.upstream_segment_IDs) == 0:
                is_head = True
                up_g = None
                x_up = 2 * lp.x[0] - lp.x[1]
                z_up = lp.z[0] + lp.S0 * (lp.x[0] - x_up)
                if lp.Q_ghost_upstream is not None:
                    Q_up = lp.Q_ghost_upstream
                else:
                    Q_up = 2 * lp.Q[0] - lp.Q[1]
            elif len(lp.upstream_segment_IDs) == 1:
                is_head = False
                upseg = segs[lp.upstream_segment_IDs[0]]
                up_g = starts[upseg.ID] + lengths[upseg.ID] - 1
                z_up = upseg.z[-1]
                x_up = upseg.x[-1]
                Q_up = upseg.Q[-1]
            else:
                raise NotImplementedError(
                    "assemble_by_walking: multi-tributary confluence "
                    "(segment %d has %d upstream segments) not yet handled"
                    % (s, len(lp.upstream_segment_IDs)))
            # --- downstream neighbor (or outlet ghost) ---
            if i < L - 1:
                down_g = g + 1
                z_down = lp.z[i + 1]
                x_down = lp.x[i + 1]
                Q_down = lp.Q[i + 1]
                is_outlet = False
            elif len(lp.downstream_segment_IDs) == 0:
                is_outlet = True
                down_g = None
                # Base-level node: elevation z_bl at position
                # x_ghost_downstream (settable via set_x_bl to move the
                # mouth in x; defaults to one cell beyond the last node).
                if lp.x_ghost_downstream is not None:
                    x_down = lp.x_ghost_downstream
                else:
                    x_down = 2 * lp.x[-1] - lp.x[-2]
                z_down = lp.z_bl
                if lp.Q_ghost_downstream is not None:
                    Q_down = lp.Q_ghost_downstream
                else:
                    Q_down = 2 * lp.Q[-1] - lp.Q[-2]
            else:
                is_outlet = False
                downseg = segs[lp.downstream_segment_IDs[0]]
                down_g = starts[downseg.ID]
                z_down = downseg.z[0]
                x_down = downseg.x[0]
                Q_down = downseg.Q[0]
            # --- stencil (coefficients from build_LHS_coeff_C0) ---
            dx_up = lp.x[i] - x_up
            dx_down = x_down - lp.x[i]
            dx_2cell = x_down - x_up
            dQ_2cell = Q_down - Q_up
            S = np.abs(z_down - z_up) / dx_2cell
            C1 = lp.C0 * S ** (1 / 6.) * lp.Q[i] / lp.B[i]
            center = -C1 / dx_2cell * (7 / 3. * (-1 / dx_up - 1 / dx_down)) + 1.
            left = -C1 / dx_2cell * (7 / 3. / dx_up - dQ_2cell / lp.Q[i] / dx_2cell)
            right = -C1 / dx_2cell * (7 / 3. / dx_down + dQ_2cell / lp.Q[i] / dx_2cell)
            rhs_g = z_rhs[i] + src[i]
            if is_head:                       # set_bcl_Neumann
                right = -C1 / dx_2cell * 7 / 3. * (1 / dx_up + 1 / dx_down)
                rhs_g += lp.S0 * C1 * (7 / 3. / dx_up - dQ_2cell / lp.Q[i] / dx_2cell)
            if is_outlet:                     # set_bcr_Dirichlet
                rhs_g += z_down * C1 / dx_2cell * (
                    7 / 3. * (1 / (lp.x[-1] - lp.x[-2])
                              + 1 / (x_down - lp.x[-1])) / 2.
                    + dQ_2cell / lp.Q[i] / dx_2cell)
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
    segs = net.list_of_LongProfile_objects
    lengths = list(net.list_of_segment_lengths)
    starts = np.cumsum([0] + lengths)[:-1]
    # Sternberg gravel loss enters as a distributed sink that depends on the
    # (evolving) sediment discharge, so it must be relinearized each Picard
    # iteration. Skip the recompute entirely when no segment sets it.
    gravel_loss_active = any(
        lp.gravel_fractional_loss_per_km is not None for lp in segs)
    for ti in range(int(nt)):
        for lp in segs:
            lp.zold = lp.z.copy()
        for _ in range(int(net.niter)):
            if gravel_loss_active:
                update_gravel_loss(net)
            LHS, RHS = assemble(net, dt)
            out = spsolve(sparse.csr_matrix(LHS), RHS)
            for lp in segs:
                s = lp.ID
                lp.z = out[starts[s]:starts[s] + lengths[s]]
        net.t += dt
        for lp in segs:
            lp.t = net.t
            lp.dz_dt = (lp.z - lp.zold) / dt

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
    for lp in net.list_of_LongProfile_objects:
        if lp.gravel_fractional_loss_per_km is not None:
            lp.downstream_fining_subsidence_equivalent = \
                - lp.gravel_fractional_loss_per_km / 1000. * lp.Q_s \
                / ( (1 - lp.lambda_p) * lp.B )

