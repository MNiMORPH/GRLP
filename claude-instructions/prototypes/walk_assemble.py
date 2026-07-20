"""
Validate a neighbor-walking assembler against LongProfile.build_matrices,
bit-for-bit, for a single segment. Neighbor lookup is by global index (up[g],
dn[g]); ghosts are reconstructed from the BCs exactly as GRLP does. If this
matches build_matrices to 0.0, the same walk (with up[g] pointing across a
segment boundary at a single-upstream junction) fixes 1-into-1 by construction.
"""
import sys
sys.path.insert(0, "/home/awickert/models/GRLP/tests")
import numpy as np
from scipy.sparse import spdiags
from conftest import make_long_profile


def walk_assemble(z, x, Q, B, up, dn, C0, S0, z_bl, ssd, dfse, U, dt):
    n = len(z)
    left = np.zeros(n); center = np.zeros(n); right = np.zeros(n)
    bcl = 0.0; bcr = 0.0
    for g in range(n):
        # --- neighbors (ghost the boundary sides exactly as GRLP set_x/set_Q) ---
        if up[g] is None:                       # channel head
            dx_up = x[g + 0] - (2 * x[g] - x[g + 1]) if n > 1 else 0.0
            # GRLP: x_ext[0] = 2*x[0]-x[1]; ghost z via Neumann; ghost Q linear
            x_up = 2 * x[0] - x[1]; z_up = z[0] + S0 * (x[0] - x_up)
            Q_up = 2 * Q[0] - Q[1]
        else:
            x_up = x[up[g]]; z_up = z[up[g]]; Q_up = Q[up[g]]
        if dn[g] is None:                        # outlet
            x_dn = 2 * x[-1] - x[-2]; z_dn = z_bl; Q_dn = 2 * Q[-1] - Q[-2]
        else:
            x_dn = x[dn[g]]; z_dn = z[dn[g]]; Q_dn = Q[dn[g]]
        dxL = x[g] - x_up; dxR = x_dn - x[g]; dx2 = x_dn - x_up
        dQ2 = Q_dn - Q_up
        S = abs(z_dn - z_up) / dx2
        C1 = C0 * S ** (1 / 6.0) * Q[g] / B[g]
        left[g] = -C1 / dx2 * (7 / 3.0 / dxL - dQ2 / Q[g] / dx2)
        center[g] = -C1 / dx2 * (7 / 3.0 * (-1 / dxL - 1 / dxR)) + 1.0
        right[g] = -C1 / dx2 * (7 / 3.0 / dxR + dQ2 / Q[g] / dx2)
        if up[g] is None:                        # set_bcl_Neumann
            right[g] = -C1 / dx2 * 7 / 3.0 * (1 / dxL + 1 / dxR)
            bcl = S0 * C1 * (7 / 3.0 / dxL - dQ2 / Q[g] / dx2)
        if dn[g] is None:                        # set_bcr_Dirichlet (z_bl)
            bcr = z_bl * C1 / dx2 * (7 / 3.0 * (1 / (x[-1] - x[-2]) + 1 / (x_dn - x[-1])) / 2.0
                                     + dQ2 / Q[g] / dx2)
    left = np.roll(left, -1); right = np.roll(right, 1)
    LHS = spdiags(np.vstack((left, center, right)), [-1, 0, 1], n, n, format="csr")
    RHS = np.hstack((bcl + z[0], z[1:-1], bcr + z[-1])) + ssd * dt + dfse * dt + U * dt
    return LHS.toarray(), RHS


# reference: rebuild with ARRAY Q so the standalone uses linear ghosts
# (2*Q[0]-Q[1], 2*Q[-1]-Q[-2]) -- the same convention the network/walk uses.
import grlp
_ref = make_long_profile(nx=20, intermittency=1.0)
x = _ref.x.copy(); Q = _ref.Q.copy(); B = _ref.B.copy(); n = 20
lp = grlp.LongProfile(); lp.set_intermittency(1.0)
lp.basic_constants(); lp.bedload_lumped_constants(); lp.set_hydrologic_constants()
lp.set_x(dx=x[1] - x[0], nx=n, x0=x[0]); lp.set_z(S0=-1.5e-2)
lp.set_A(k_xA=1.0); lp.set_Q(Q=Q.copy()); lp.set_B(B=B.copy())
lp.set_uplift_rate(0.0); lp.set_niter(5)
lp.set_Qs_input_upstream(lp.k_Qs * lp.Q[0] * 1.5e-2 ** (7 / 6.0)); lp.set_z_bl(0.0)
zc = np.linspace(30.0, 1.0, n)
lp.set_z(z=zc.copy()); lp.set_z_bl(0.0)
lp.z_ext[1:-1] = zc            # sync interior ghosts (evolve normally does this)
lp.build_LHS_coeff_C0(dt=1e12); lp.zold = zc.copy(); lp.build_matrices()
Mref = lp.LHSmatrix.toarray(); Rref = lp.RHS.copy()

up = [None] + list(range(0, n - 1)); dn = list(range(1, n)) + [None]
Mw, Rw = walk_assemble(zc, x, Q, B, up, dn, lp.C0, lp.S0, 0.0,
                       lp.ssd, lp.downstream_fining_subsidence_equivalent, lp.U, 1e12)
print("max|M_walk - M_buildmatrices| =", np.abs(Mw - Mref).max())
print("max|R_walk - R_buildmatrices| =", np.abs(Rw - Rref).max())
