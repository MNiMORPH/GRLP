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
        # second-order-time-bdf2-design.md). Backward Euler (time_order == 1) uses
        # the start-of-step elevation z^n (seg.zold) with a unit time-diagonal.
        # BDF2 (net.time_order == 2, the default, once a two-level history exists)
        # uses the
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
        # Channel-deposit fraction as a per-node array (1 by default).  The
        # gravel transport and storage both act over the effective width
        # f_ch * B (the space the channel gravel actually occupies), so the
        # flux conductance is divided by f_ch * B below.  This breaks the
        # f_ch cancellation the volume-first row-scaling would otherwise cause
        # (Jstore carries f_ch*B; dividing the flux by f_ch*B leaves the flux
        # f_ch-free while storage keeps it, so f_ch couples).
        f_ch_seg = np.broadcast_to(np.asarray(seg.f_ch, dtype=float), (L,))
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
                # Effective confluence-cell area with the channel-deposit
                # fraction: each contributing cell carries its own f_ch (the
                # gravel occupies the effective width f_ch*B).  Recomputed here
                # from the current f_ch (not the frozen land_area_around_
                # confluence, which is set once from B) -- it mirrors that sum
                # with f_ch on each cell and reduces to it exactly at f_ch = 1.
                A_confluence = (0.5 * (seg.x[1] - seg.x[0])
                                * f_ch_seg[0] * seg.B[0])
                for _ID in seg.upstream_segment_IDs:
                    _up = segs[_ID]
                    _up_fch = np.broadcast_to(
                        np.asarray(_up.f_ch, dtype=float), _up.B.shape)
                    A_confluence += (0.5 * (seg.x[0] - _up.x[-1])
                                     * _up_fch[-1] * _up.B[-1])
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
                land_area = f_ch_seg[-1] * seg.B[-1] * 0.5 * (
                    (seg.x[-1] - seg.x[-2]) + (downseg.x[0] - seg.x[-1]))
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
                land_area = f_ch_seg[1] * seg.B[1] * 0.5 * (
                    (seg.x[1] - seg.x[0]) + (seg.x[2] - seg.x[1]))
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
            C1 = seg.C0 * S ** (1 / 6.) * seg.Q[i] / (f_ch_seg[i] * seg.B[i])
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
    tolerance set (``set_iteration_tolerance``) ``max_iter`` is the safety cap."""
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
                 gravel_loss_active, picard_tol, max_iter, in_picard=False):
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
        if in_picard:
            # Tight coupling: recompute the valley geometry from this iterate so
            # the next assemble sees storage (f_ch*B) consistent with the current
            # z.  update_valley is idempotent (it resets to the frozen Bold/Hold),
            # so repeating it across iterations does not accumulate.
            for seg in segs:
                seg.dz_dt = (seg.z - seg.zold) / dt
                seg.update_valley(dt)
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
    in_picard = getattr(net, "_valley_in_picard", False)
    if in_picard and bdf2_requested:
        warnings.warn(
            "In-Picard valley coupling makes the storage term nonlinear in z, "
            "so BDF2 is no longer strictly second order; consider "
            "set_time_integration(1) (backward Euler).", RuntimeWarning)
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
            # Freeze the start-of-step valley geometry; the in-Picard coupling
            # resets to it on every iterate (update_valley is idempotent in it).
            seg.Bold = None if seg.B is None else seg.B.copy()
            seg.Hold = None if seg.H_valley is None else seg.H_valley.copy()
        net._bdf2_step = bdf2_requested and ti > 0
        _picard_step(net, dt, segs, starts, lengths,
                     gravel_loss_active, picard_tol, max_iter, in_picard)
        net.t += dt
        for seg in segs:
            seg.t = net.t
            seg.dz_dt = (seg.z - seg.zold) / dt
            # Between-step coupling (default): advance the valley geometry once,
            # after the solve (no-op unless valley dynamics are on).  The
            # in-Picard coupling instead advances it each iteration, inside
            # _picard_step, so skip the hand-off here.
            if not in_picard:
                seg.update_valley(dt)


def _trial_step(net, z_curr, z_prev, dt, dt_prev, use_bdf2,
                segs, starts, lengths, gravel_loss_active, picard_tol, max_iter):
    """One step-doubling trial from the two-level state ``(z_curr, z_prev)``.

    Takes one step of ``dt`` (``z_big``) and two steps of ``dt/2`` (``z_small``),
    both from the same start, and returns ``(z_small, est)`` with
    ``est = max|z_big - z_small|``. ``z_small`` is the finer (advanced) solution
    (local extrapolation); ``est`` estimates its error up to a constant factor.
    With ``use_bdf2=False`` all three sub-steps are backward Euler (the bootstrap,
    which has no ``z_prev``); otherwise they are variable-step BDF2, the first
    half-step reaching back to the ``dt_prev``-spaced ``z_prev`` (omega < 1)."""
    def _solve(dt_sub, omega):
        net._bdf2_step = use_bdf2
        net._bdf2_omega = omega
        _picard_step(net, dt_sub, segs, starts, lengths,
                     gravel_loss_active, picard_tol, max_iter)
        return [seg.z.copy() for seg in segs]

    def _set(zold, zold2):
        for i, seg in enumerate(segs):
            seg.zold2 = None if zold2 is None else zold2[i]
            seg.zold = zold[i]
            seg.z = zold[i].copy()

    hist2 = z_prev if use_bdf2 else None
    # z_big: one step of dt
    _set(z_curr, hist2)
    z_big = _solve(dt, (dt / dt_prev) if use_bdf2 else 1.0)
    # z_small: two steps of dt/2 (second half is uniform, omega = 1)
    _set(z_curr, hist2)
    z_half = _solve(dt / 2.0, (0.5 * dt / dt_prev) if use_bdf2 else 1.0)
    _set(z_half, z_curr if use_bdf2 else None)
    z_small = _solve(dt / 2.0, 1.0)
    est = max(float(np.max(np.abs(b - s))) for b, s in zip(z_big, z_small))
    return z_small, est


def evolve_adaptive(net, T, dt_init=None):
    """
    Adaptive-time-step integration (MNiMORPH/GRLP#16, Phase 2).

    Advance the network by a total time ``T`` [s], choosing each step size so the
    per-step step-doubling error estimate stays at or below the tolerance set by
    ``set_adaptive_timestep``. Unlike :func:`evolve`, this carries the BDF2
    two-level history across steps and resizes the step as the transient speeds up
    or settles: small steps where the profile changes fast, large where it is
    smooth.

    Method: variable-step BDF2 (second order), bootstrapped with one
    backward-Euler step; each step is estimated by step doubling (one step of Δt
    vs two of Δt/2 from the same state) and the network is advanced with the finer
    two-half-step solution (local extrapolation). Step size is chosen by an
    I-controller, ``Δt_new = Δt · safety · (tol/est)^{1/(p+1)}`` (p = 2, or 1 for
    the backward-Euler bootstrap), clamped by a max growth/shrink ratio and by
    ``dt_min``/``dt_max``; a step whose estimate exceeds ``tol`` is rejected and
    retried smaller. ``dt_init`` is only the starting guess (the bootstrap is
    itself error-controlled); it defaults to ``adaptive_dt_init`` or a small
    fraction of ``T``.

    Turn on ``set_iteration_tolerance`` for a converged per-step solve; the estimate
    and step both assume the nonlinear iteration has converged.
    """
    tol = getattr(net, "adaptive_tol", None)
    if tol is None:
        raise ValueError("adaptive stepping needs a tolerance; call "
                         "set_adaptive_timestep(tol) before evolve_adaptive")
    if getattr(net, "_valley_in_picard", False):
        warnings.warn(
            "In-Picard valley coupling is not supported by the adaptive solver; "
            "using the between-step coupling.", RuntimeWarning)
    dt_min = getattr(net, "adaptive_dt_min", 0.0)
    dt_max = getattr(net, "adaptive_dt_max", np.inf)
    safety = getattr(net, "adaptive_safety", 0.9)
    max_grow = getattr(net, "adaptive_max_grow", 5.0)
    max_shrink = getattr(net, "adaptive_max_shrink", 0.1)
    if dt_init is None:
        dt_init = getattr(net, "adaptive_dt_init", None)
    picard_tol, max_iter, gravel_loss_active = _picard_config(net)
    segs = net.segments
    lengths = list(net.list_of_segment_lengths)
    starts = np.cumsum([0] + lengths)[:-1]

    t_end = net.t + T
    eps = 1e-9 * max(abs(T), 1.0)

    def _propose(dt_used, est, p):
        # I-controller: grow/shrink toward holding est == tol; est == 0 => grow max.
        if est <= 0.0:
            factor = max_grow
        else:
            factor = safety * (tol / est) ** (1.0 / (p + 1.0))
        factor = min(max_grow, max(max_shrink, factor))
        return min(dt_max, max(dt_min, dt_used * factor))

    def _advance(z_curr, z_prev, dt_used):
        # commit the accepted step: publish z onto the segments, update t/dz_dt.
        net.t += dt_used
        net.dt = dt_used
        for i, seg in enumerate(segs):
            seg.zold = z_prev[i]
            seg.z = z_curr[i]
            seg.t = net.t
            seg.dz_dt = (z_curr[i] - z_prev[i]) / dt_used
            # Advance the valley geometry between steps (no-op unless valley
            # dynamics are switched on).  Freeze Bold/Hold first so the
            # idempotent reset in update_valley is a no-op here.
            seg.Bold = None if seg.B is None else seg.B.copy()
            seg.Hold = None if seg.H_valley is None else seg.H_valley.copy()
            seg.update_valley(dt_used)

    z_curr = [seg.z.copy() for seg in segs]
    # ---- bootstrap: one error-controlled backward-Euler step (no z_prev) ----
    dt = dt_init if dt_init else min(dt_max, max(dt_min, 1e-2 * T))
    dt = min(dt, t_end - net.t)
    z_prev = z_curr
    dt_prev = None
    first = True
    while net.t < t_end - eps:
        dt = min(max(dt, dt_min), t_end - net.t)
        use_bdf2 = not first
        p = 1 if first else 2
        rejects = 0
        while True:
            z_small, est = _trial_step(
                net, z_curr, z_prev, dt, dt_prev, use_bdf2,
                segs, starts, lengths, gravel_loss_active, picard_tol, max_iter)
            at_floor = dt <= dt_min * (1.0 + 1e-12)
            if est <= tol or at_floor or rejects >= 50:
                break
            # reject: shrink and retry the same step
            shrink = max(max_shrink, safety * (tol / est) ** (1.0 / (p + 1.0)))
            dt = max(dt_min, dt * min(shrink, 0.9))
            rejects += 1
        if est > tol and at_floor:
            warnings.warn(
                "Adaptive stepping hit dt_min=%g with est=%g > tol=%g at t=%g; "
                "step is under-resolved." % (dt_min, est, tol, net.t),
                RuntimeWarning)
        z_prev, z_curr = z_curr, z_small
        _advance(z_curr, z_prev, dt)
        dt_prev = dt
        dt = _propose(dt, est, p)
        first = False


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

