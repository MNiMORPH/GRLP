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

import warnings

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
    # Volume-first time integration: the solver conserves stored sediment
    # *volume* V, not bed elevation z (see claude-instructions/
    # second-order-time-bdf2-design.md). We assemble the elevation-based system
    # below, then row-scale it by the storage Jacobian J = dV/dz and swap the
    # RHS history into V-space. For the rectangular default (constant B) V is
    # linear in z, J is constant, and this reduces to an exact row-scaling that
    # leaves the solution unchanged to machine precision; a dynamic-width valley
    # (issue #19) makes it a genuine change that keeps BDF2 second-order.
    Jstore = np.ones(n)      # dV/dz at the current Picard iterate
    Vhist = np.zeros(n)      # time-level history in V-space (BE or BDF2)
    Vcorr = np.zeros(n)      # nonlinear-map linearization correction (0 if linear)
    z_hist = np.zeros(n)     # the elevation history added to RHS (to subtract off)
    for seg in segs:
        s = seg.ID
        offset = starts[s]
        L = lengths[s]
        # per-node source term (RHS additions: ssd + fining/subsidence + U)
        src = (np.asarray(seg.ssd)
               + np.asarray(seg.downstream_fining_subsidence_equivalent)
               + np.asarray(seg.U)) * dt
        src = np.broadcast_to(src, (L,))
        # Time-derivative discretization (see claude-instructions/
        # second-order-time-bdf2-design.md). Backward Euler (default) uses the
        # start-of-step elevation z^n (seg.zold) with a unit time-diagonal.
        # BDF2 (net.time_order == 2, once a two-level history exists) uses the
        # three-level derivative (3 z^{n+1} - 4 z^n + z^{n-1})/(2 dt): the
        # time-diagonal becomes 3/2 and the RHS history becomes (4 z^n - z^{n-1})/2.
        # The dt * operator (C1) and dt * source terms are identical either way.
        # RHS/coefficient split during Picard: the history uses the frozen
        # start-of-step levels; the coefficient (C1) uses the current iterate
        # seg.z. When zold is unset (static assembly, e.g. tests) it equals seg.z.
        zold = getattr(seg, "zold", None)
        if zold is None or np.size(zold) != L:
            zold = seg.z
        zold2 = getattr(seg, "zold2", None)
        use_bdf2 = (getattr(net, "_bdf2_step", False)
                    and zold2 is not None and np.size(zold2) == L)
        if use_bdf2:
            # Variable-step BDF2: the three-level derivative on a non-uniform
            # history. With omega = dt_n / dt_{n-1} the weights are
            # a = (1+2 omega)/(1+omega), b = 1+omega, c = omega^2/(1+omega), so the
            # time term is [a z^{n+1} - b z^n + c z^{n-1}]/dt_n; the coefficient a
            # is the time-diagonal and (b z^n - c z^{n-1}) is the RHS history. For
            # uniform steps (omega = 1) these are (3/2, 2, 1/2), reproducing the
            # fixed-step BDF2 bit-for-bit. Uniform stepping (evolve) leaves
            # net._bdf2_omega = 1; adaptive stepping varies it each step.
            w = getattr(net, "_bdf2_omega", 1.0)
            bdf2_b = 1.0 + w
            bdf2_c = w * w / (1.0 + w)
            z_rhs = bdf2_b * zold - bdf2_c * zold2
            time_diag = (1.0 + 2.0 * w) / (1.0 + w)
        else:
            z_rhs = zold
            time_diag = 1.0
        # Volume-first pieces for this segment (applied as a transform after the
        # matrix is assembled). J = dV/dz at the iterate; the history is carried
        # in V-space (V^n for BE, (4 V^n - V^{n-1})/2 for BDF2); the correction
        # is the nonlinear-map linearization term (identically zero when V is
        # linear in z, i.e. the rectangular constant-B default).
        sl = slice(offset, offset + L)
        z_hist[sl] = z_rhs
        Jstore[sl] = seg.storage_jacobian(seg.z)
        Vold = seg.storage_volume(zold)
        if use_bdf2:
            Vhist[sl] = bdf2_b * Vold - bdf2_c * seg.storage_volume(zold2)
        else:
            Vhist[sl] = Vold
        Vcorr[sl] = time_diag * (Jstore[sl] * seg.z - seg.storage_volume(seg.z))
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
                vals.append(time_diag + conductance_sum / A_confluence)
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
                    time_diag + (conductance_up + conductance_downseg) / land_area)
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
                    time_diag + (conductance_up + conductance_down) / land_area)
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
            center = -C1 / dx_2cell * (7 / 3. * (-1 / dx_up - 1 / dx_down)) + time_diag
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
    # Volume-first transform: row-scale the elevation system by J = dV/dz and
    # swap the RHS history into V-space. Each row i of the elevation system is
    # `(time_diag + flux_stencil) z = z_hist_i + (src+bcs)_i`; scaling the whole
    # row by J_i turns the time diagonal into J_i*time_diag (= d/dz of the V
    # storage term) and un-normalizes the flux, while the RHS becomes
    # `Vhist_i + Vcorr_i + J_i*(src+bcs)_i` with `(src+bcs)_i = RHS_i - z_hist_i`.
    # Rectangular default: J constant, Vhist = J*z_hist, Vcorr = 0 -> an exact
    # row-scaling (solution unchanged to machine precision).
    rows = np.asarray(rows)
    vals = np.asarray(vals) * Jstore[rows]
    LHSmatrix = sparse.csr_matrix((vals, (rows, cols)), shape=(n, n))
    RHS = Vhist + Vcorr + Jstore * (RHS - z_hist)
    return LHSmatrix, RHS

def _picard_config(net):
    """Read the semi-implicit iteration controls off the network. Returns
    ``(picard_tol, max_iter, gravel_loss_active)``. In the default fixed-niter
    mode ``picard_tol`` is None and ``max_iter`` is ``net.niter``; with a
    tolerance set (``set_picard_tolerance``) ``max_iter`` is the safety cap."""
    picard_tol = getattr(net, "picard_tol", None)
    if picard_tol is None:
        max_iter = int(net.niter)
    else:
        max_iter = int(getattr(net, "picard_max_iter", 100))
    # Sternberg gravel loss enters as a distributed sink that depends on the
    # (evolving) sediment discharge, so it must be relinearized each Picard
    # iteration. Skip the recompute entirely when no segment sets it.
    gravel_loss_active = any(
        seg.gravel_fractional_loss_per_km is not None for seg in net.segments)
    return picard_tol, max_iter, gravel_loss_active


def _picard_step(net, dt, segs, starts, lengths,
                 gravel_loss_active, picard_tol, max_iter):
    """Run the semi-implicit (Picard) iteration for one step at the network's
    current history/BDF2 state (``seg.zold``/``seg.zold2``/``net._bdf2_step``
    already set by the caller), updating ``seg.z`` in place. The RHS is frozen at
    the start-of-step history while the coefficient relinearizes on the current
    iterate each iteration. Iterates a fixed ``max_iter`` times in fixed-niter
    mode; in tolerance mode it stops once the inter-iteration elevation change
    ``max|z_k - z_{k-1}|`` falls below ``picard_tol`` and warns if the cap is
    hit first. Shared by ``evolve`` and ``evolve_adaptive``."""
    converged = picard_tol is None   # fixed-niter mode: nothing to check
    change = 0.0
    for k in range(max_iter):
        if gravel_loss_active:
            update_gravel_loss(net)
        LHS, RHS = assemble(net, dt)
        out = spsolve(sparse.csr_matrix(LHS), RHS)
        change = 0.0
        for seg in segs:
            s = seg.ID
            znew = out[starts[s]:starts[s] + lengths[s]]
            if picard_tol is not None:
                change = max(change, np.max(np.abs(znew - seg.z)))
            seg.z = znew
        if picard_tol is not None and change < picard_tol:
            converged = True
            break
    if not converged:
        warnings.warn(
            "Picard iteration did not converge to tol=%g in %d iterations "
            "(last change %g m) at t=%g; result may be under-converged."
            % (picard_tol, max_iter, change, net.t),
            RuntimeWarning)


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
    bdf2_requested = getattr(net, "time_order", 1) == 2
    picard_tol, max_iter, gravel_loss_active = _picard_config(net)
    segs = net.segments
    lengths = list(net.list_of_segment_lengths)
    starts = np.cumsum([0] + lengths)[:-1]
    for ti in range(int(nt)):
        for seg in segs:
            # Advance the time-level history for BDF2: z^{n-1} <- z^n, z^n <- z.
            # The first step of an evolve() call has no z^{n-1}, so BDF2
            # bootstraps with one backward-Euler step there.
            seg.zold2 = None if ti == 0 else seg.zold
            seg.zold = seg.z.copy()
        net._bdf2_step = bdf2_requested and ti > 0
        _picard_step(net, dt, segs, starts, lengths,
                     gravel_loss_active, picard_tol, max_iter)
        net.t += dt
        for seg in segs:
            seg.t = net.t
            seg.dz_dt = (seg.z - seg.zold) / dt


def evolve_adaptive(net, nt, dt):
    """
    Adaptive-time-step evolution loop (MNiMORPH/GRLP#16, Phase 2).

    Unlike :func:`evolve`, which re-bootstraps the BDF2 history on the first step
    of *every* call, this is a **dedicated loop** that carries the two-level
    history ``(z^n, z^{n-1})`` across all steps -- the prerequisite for both a
    per-step local-error estimate and step-size control (see
    ``claude-instructions/adaptive-timestepping-design.md``).

    **Increment 1 (this version): fixed step, no control yet.** It takes ``nt``
    steps of size ``dt``, bootstrapping the first step with backward Euler and
    using BDF2 thereafter, so with ``set_time_integration(2)`` it reproduces
    :func:`evolve` bit-for-bit at fixed ``dt``. The embedded error estimator, the
    step controller, and variable-step BDF2 weights land in later increments;
    the signature will grow a tolerance then.
    """
    net.dt = dt
    bdf2_requested = getattr(net, "time_order", 1) == 2
    picard_tol, max_iter, gravel_loss_active = _picard_config(net)
    segs = net.segments
    lengths = list(net.list_of_segment_lengths)
    starts = np.cumsum([0] + lengths)[:-1]
    dt_prev = dt
    for ti in range(int(nt)):
        for seg in segs:
            seg.zold2 = None if ti == 0 else seg.zold
            seg.zold = seg.z.copy()
        net._bdf2_step = bdf2_requested and ti > 0
        # Variable-step BDF2 needs omega = dt_n / dt_{n-1}; here dt is fixed so
        # omega = 1 every step (identical to evolve). The controller will vary it.
        net._bdf2_omega = dt / dt_prev
        _picard_step(net, dt, segs, starts, lengths,
                     gravel_loss_active, picard_tol, max_iter)
        net.t += dt
        dt_prev = dt
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

