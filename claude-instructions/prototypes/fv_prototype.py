"""
Finite-volume prototype for the GRLP transport law, over a node/face graph so
that a single segment, a 1-into-1 chain, and a real confluence all run through
the identical code. Purpose: validate that

  (a) FV single-segment reproduces the analytical steady state (2nd order),
  (b) splitting a segment (1-into-1) is an exact no-op,
  (c) a 2-tributary confluence conserves sediment and gives each reach the
      correct threshold-width slope with the right amplitude.

Governing law (Exner + threshold-width transport):
  (1-lp) B dz/dt = -dQs/dx,   Qs = K Q S^{7/6},  S = -dz/dx,  K = k_Qs I / sin^{7/6}

FV form: integrate over a cell of plan-view area A = B*dx; convert dQs/dx to a
balance of sediment fluxes through the cell's faces:
  (1-lp) A_i dz_i/dt = sum(Qs through upstream faces) - Qs(downstream face)

Face flux (implicit, Picard-frozen |S|^{1/6}):
  Qs_face = K * Q_face * |S_face^old|^{1/6} * (z_up - z_dn)/dx_face
The SAME formula is used on every face, so influx == outflux by construction
(conservative), and a face's Q is single-valued -> no dQ/dx across a junction.
"""
import sys
sys.path.insert(0, "/home/awickert/models/GRLP/tests")
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
from conftest import make_long_profile

# --- physical constants pulled from GRLP so we match it ---
_lp = make_long_profile(intermittency=1.0)
k_Qs = _lp.k_Qs
lambda_p = _lp.lambda_p
sinuosity = _lp.sinuosity
I = 1.0
K = k_Qs * I / sinuosity ** (7 / 6.0)   # Qs = K Q S^{7/6}


class FVGraph:
    """Nodes carry (x, Q, B, z); faces connect an upstream node to a downstream
    node.  Channel-head nodes receive an imposed sediment supply Qs0; the outlet
    node is Dirichlet at z_bl."""
    def __init__(self):
        self.x = []; self.Q = []; self.B = []
        self.faces = []            # (up_node, dn_node, Qface, dxface)
        self.supply = {}           # head_node -> Qs0
        self.outlet = None; self.z_bl = 0.0

    def add_node(self, x, Q, B):
        self.x.append(x); self.Q.append(Q); self.B.append(B)
        return len(self.x) - 1

    def add_face(self, up, dn, Qface):
        dxf = abs(self.x[dn] - self.x[up])
        self.faces.append((up, dn, Qface, dxf))

    def cell_dx(self, i):
        # half-distance to each neighbor across the faces touching node i
        halves = []
        for (u, d, Qf, dxf) in self.faces:
            if u == i or d == i:
                halves.append(dxf / 2.0)
        return sum(halves) if halves else 1.0

    def solve(self, nt=4000, dt=1e12, niter=6):
        n = len(self.x)
        z = np.zeros(n)
        z[self.outlet] = self.z_bl
        A = np.array([self.B[i] * self.cell_dx(i) for i in range(n)])
        for step in range(nt):
            zold = z.copy()
            for _ in range(niter):
                rows = []; cols = []; vals = []
                rhs = (1 - lambda_p) * A / dt * zold
                diag = (1 - lambda_p) * A / dt
                for (u, d, Qf, dxf) in self.faces:
                    S = abs(z[u] - z[d]) / dxf
                    Dface = K * Qf * S ** (1 / 6.0) / dxf   # secant conductance
                    # flux u->d = Dface*(z_u - z_d): leaves u, enters d
                    diag_u = Dface; diag_d = Dface
                    rows += [u, u, d, d]
                    cols += [u, d, d, u]
                    vals += [Dface, -Dface, Dface, -Dface]
                for i in range(n):
                    rows.append(i); cols.append(i); vals.append(diag[i])
                M = sparse.csr_matrix((vals, (rows, cols)), shape=(n, n))
                # channel-head supply: constant influx on the RHS
                for h, Qs0 in self.supply.items():
                    rhs[h] += Qs0
                # outlet Dirichlet
                o = self.outlet
                M = M.tolil()
                M.rows[o] = [o]; M.data[o] = [1.0]; rhs[o] = self.z_bl
                z = spsolve(M.tocsr(), rhs)
        return z


def analytical(x, Q, S0, x_bl_slope_node=0):
    """Threshold-width steady state fit to boundary slope S0 at the head,
    integrated using the same nodal Q (so it is the discrete-consistent target)."""
    Qs0 = k_Qs * Q[0] * S0 ** (7 / 6.0)
    S = (Qs0 / (k_Qs * Q)) ** (6 / 7.0)          # local slope from Qs0
    # integrate slope from outlet (z=0) upstream
    z = np.zeros_like(x)
    for i in range(len(x) - 2, -1, -1):
        z[i] = z[i + 1] + S[i] * (x[i + 1] - x[i])
    return z


# ============================================================
# (a) single segment vs analytical, and convergence order
# ============================================================
print("=== (a) FV single-segment vs analytical ===")
S0 = 1.5e-2
prev = None
for nx in (45, 89, 177):
    lp = make_long_profile(nx=nx, dx=88000.0 / (nx - 1), x0=10000.0, intermittency=1.0)
    x = lp.x; Q = lp.Q; B = lp.B
    g = FVGraph()
    for i in range(nx):
        g.add_node(x[i], Q[i], B[i])
    for i in range(nx - 1):
        g.add_face(i, i + 1, 0.5 * (Q[i] + Q[i + 1]))
    g.supply[0] = k_Qs * Q[0] * S0 ** (7 / 6.0)
    g.outlet = nx - 1; g.z_bl = 0.0
    z = g.solve()
    za = analytical(x, Q, S0)
    err = np.abs(z - za).max()
    order = np.log2(prev / err) if prev else float("nan")
    print(f"  nx={nx:4d} dx={88000.0/(nx-1):7.1f}  |z_FV - analytical|max={err:.5f}  order={order:.3f}")
    prev = err
