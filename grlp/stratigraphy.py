"""
Valley-fill stratigraphic recorder for GRLP.

Records the channel-deposit fraction ``f_ch`` as a function of elevation and
downstream distance during a run -- a synthetic valley-fill section.  Each
downstream node keeps a ``(z, f_ch)`` polyline that grows on aggradation and is
truncated on incision (an unconformity).  The polyline is compressed with
vertical-metric Douglas-Peucker line simplification: a vertex is kept only if
dropping it and linearly interpolating across the gap would change ``f_ch`` by
more than a tolerance.  So sharp inflections (facies contacts) are held, gentle
ones are dropped, and the stored size is set by the composition tolerance, not
by the time step (it converges as ``dt -> 0``).

The recorder is a passive diagnostic: it reads ``z`` and ``f_ch`` each step and
never feeds back into the solve.

CAVEAT (see set_valley_coupling): with the default between-step valley coupling,
``f_ch`` lags one step, which -- amplified by the overbank positive feedback --
inflates transient features (a spurious sub-steady ``f_ch`` minimum just above a
surface) at coarse ``dt``.  For quantitative transient stratigraphy use the
in-Picard coupling or a finer ``dt``; the gross section is robust either way.
"""

import numpy as np


def douglas_peucker_vertical(z, f, tol):
    """
    Douglas-Peucker simplification of a ``(z, f)`` polyline using the *vertical*
    (``f``) deviation from the linear interpolant, so the tolerance ``tol`` is a
    bound on the ``f`` reconstruction error.  ``z`` must be monotonic.  Returns
    the sorted indices of the vertices to keep (always including the endpoints).
    """
    z = np.asarray(z, float)
    f = np.asarray(f, float)
    n = len(z)
    if n < 3:
        return list(range(n))
    keep = np.zeros(n, bool)
    keep[0] = keep[-1] = True
    stack = [(0, n - 1)]
    while stack:
        i0, i1 = stack.pop()
        if i1 <= i0 + 1:
            continue
        dz = z[i1] - z[i0]
        if dz == 0:
            f_lin = np.full(i1 - i0 + 1, f[i0])
        else:
            f_lin = f[i0] + (f[i1] - f[i0]) * (z[i0:i1 + 1] - z[i0]) / dz
        dev = np.abs(f[i0:i1 + 1] - f_lin)
        dev[0] = dev[-1] = 0.0
        k = int(np.argmax(dev))
        if dev[k] > tol:
            keep[i0 + k] = True
            stack.append((i0, i0 + k))
            stack.append((i0 + k, i1))
    return list(np.nonzero(keep)[0])


class StratRecorder:
    """
    Per-node ``(z, f_ch)`` stratigraphic columns with on-the-fly compression.

    Parameters
    ----------
    nx : int
        Number of downstream nodes.
    tol : float
        Maximum ``f_ch`` reconstruction error kept by the simplification.
    max_raw : int
        Compress a column (Douglas-Peucker) once its stored vertices exceed this,
        so peak memory is bounded (and independent of ``dt``).  ``None`` keeps
        every sample until queried (exact, but memory grows with the step count).
    record_thickness : bool
        Also accumulate the exact channel-deposit thickness ``sum(f_ch * dz)`` per
        node at full resolution (a cheap scalar per node), for a volume-closure
        diagnostic against the compressed section.
    """

    def __init__(self, nx, tol=0.02, max_raw=4096, record_thickness=False):
        self.nx = int(nx)
        self.tol = float(tol)
        self.max_raw = max_raw
        self.record_thickness = bool(record_thickness)
        self.z = [[] for _ in range(self.nx)]
        self.f = [[] for _ in range(self.nx)]
        # Exact channel-deposit thickness sum(f_ch*dz), accumulated on aggradation
        # and drawn down on scour (volume-closure diagnostic).
        self._exact_thickness = np.zeros(self.nx) if self.record_thickness else None

    def record(self, z, f_ch):
        """Record one step's surface: append where the bed rose, scour where it
        fell.  ``z`` is the bed elevation array; ``f_ch`` a scalar or array."""
        z = np.asarray(z, float)
        f_ch = np.broadcast_to(np.asarray(f_ch, float), z.shape)
        for i in range(self.nx):
            zi, fi = self.z[i], self.f[i]
            # Scour: remove anything at or above the new bed (unconformity).
            while zi and zi[-1] >= z[i]:
                if self._exact_thickness is not None and len(zi) >= 2:
                    self._exact_thickness[i] -= fi[-1] * (zi[-1] - zi[-2])
                zi.pop(); fi.pop()
            if self._exact_thickness is not None and zi:
                self._exact_thickness[i] += f_ch[i] * (z[i] - zi[-1])
            zi.append(float(z[i]))
            fi.append(float(f_ch[i]))
            if self.max_raw is not None and len(zi) > self.max_raw:
                idx = douglas_peucker_vertical(zi, fi, self.tol)
                self.z[i] = [zi[k] for k in idx]
                self.f[i] = [fi[k] for k in idx]

    def columns(self):
        """Return the compressed columns as a list of ``(z, f_ch)`` arrays."""
        out = []
        for i in range(self.nx):
            idx = douglas_peucker_vertical(self.z[i], self.f[i], self.tol)
            out.append((np.asarray(self.z[i])[idx], np.asarray(self.f[i])[idx]))
        return out

    def section(self, x, ns=200):
        """
        Reconstruct ``f_ch`` on a structured, boundary-conforming mesh for
        plotting.  Each column is resampled onto ``ns`` levels between its basal
        (onset) and final surfaces, so the mesh follows the wedge boundaries
        exactly (no triangulation).  Returns ``(Xmesh, Zmesh, Fmesh)``, each of
        shape ``(ns, nx)``; ``Fmesh`` is ``nan`` where a column has no deposit.
        """
        x = np.asarray(x, float)
        s = np.linspace(0.0, 1.0, ns)
        Xm = np.broadcast_to(x, (ns, self.nx)).copy()
        Zm = np.full((ns, self.nx), np.nan)
        Fm = np.full((ns, self.nx), np.nan)
        for i, (zc, fc) in enumerate(self.columns()):
            if len(zc) < 2 or zc[-1] <= zc[0]:
                continue
            zlev = zc[0] + s * (zc[-1] - zc[0])
            Zm[:, i] = zlev
            Fm[:, i] = np.interp(zlev, zc, fc)
        return Xm, Zm, Fm

    def plot_section(self, x, ax=None, cmap='viridis', ns=200):
        """Plot the reconstructed ``f_ch(x, z)`` section (Gouraud-shaded; vector
        in PDF/SVG).  Returns the ``(ax, QuadMesh)``."""
        from matplotlib import pyplot as plt
        if ax is None:
            _, ax = plt.subplots(figsize=(9, 5))
        Xm, Zm, Fm = self.section(x, ns=ns)
        pc = ax.pcolormesh(Xm, Zm, Fm, shading='gouraud', cmap=cmap,
                           vmin=0, vmax=1)
        ax.set_xlabel('Downstream distance')
        ax.set_ylabel('Elevation')
        return ax, pc

    def channel_deposit_thickness(self):
        """Channel-deposit-equivalent thickness ``integral(f_ch dz)`` per node,
        from the compressed section (trapezoidal)."""
        out = np.zeros(self.nx)
        for i, (zc, fc) in enumerate(self.columns()):
            if len(zc) >= 2:
                out[i] = np.trapezoid(fc, zc) if hasattr(np, 'trapezoid') \
                         else np.trapz(fc, zc)
        return out

    def volume_closure(self):
        """Volume-closure diagnostic: the per-node channel-deposit thickness from
        the compressed section vs the exact accumulated ``sum(f_ch dz)``.  Returns
        ``(compressed, exact, max_abs_difference)``; requires ``record_thickness``.
        The two should agree to within ~``tol`` times the column thickness."""
        if self._exact_thickness is None:
            raise ValueError("volume_closure requires record_thickness=True")
        compressed = self.channel_deposit_thickness()
        exact = self._exact_thickness.copy()
        return compressed, exact, float(np.max(np.abs(compressed - exact)))
