"""
Valley-fill stratigraphic recorder for GRLP.

Records one or more per-layer deposit *attributes* -- the channel-deposit
fraction ``f_ch``, the deposit ``age`` (the model time at which each layer was
laid down), and any further attribute added later -- as a function of elevation
and downstream distance during a run.  Each downstream node keeps a ``z``
(elevation) polyline that grows on aggradation and is truncated on incision (an
unconformity), with a parallel value array for every attribute.

Each attribute is compressed with the *same* machinery: vertical-metric
Douglas-Peucker line simplification, which keeps a vertex only if dropping it and
linearly interpolating across the gap would change that attribute by more than
its tolerance.  A layer boundary is kept if *any* attribute needs it (the union),
so all attributes share one set of layer boundaries.  Storage is set by the
attribute tolerances, not by the time step (it converges as ``dt -> 0``).

The recorder is a passive diagnostic: it reads the bed elevation and the
attribute values each step and never feeds back into the solve.

CAVEAT (see set_valley_coupling): with the between-step valley coupling ``f_ch``
lags one step, which -- amplified by the overbank feedback -- inflates transient
features at coarse ``dt``.  Use the in-Picard coupling (the default once valley
dynamics are active) or a finer ``dt`` for quantitative transient stratigraphy.
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
    Per-node stratigraphic columns with on-the-fly, multi-attribute compression.

    Parameters
    ----------
    nx : int
        Number of downstream nodes.
    tol : float
        Douglas-Peucker tolerance for ``f_ch`` (the channel-deposit fraction).
    max_raw : int
        Compress a column once its stored vertices exceed this, bounding peak
        memory (and making it ``dt``-independent).  ``None`` keeps every sample
        until queried (exact, but memory grows with the step count).
    record_thickness : bool
        Also accumulate the exact channel-deposit thickness ``sum(f_ch * dz)`` per
        node at full resolution, for a volume-closure diagnostic.
    attribute_tol : dict, optional
        Additional attributes to record beyond ``f_ch``, mapping name -> tolerance.
        A tolerance of ``None`` means "relative": 1% of that attribute's current
        span (convenient for ``age``, whose absolute scale is the run duration).
        Example: ``{'age': None}``.
    """

    def __init__(self, nx, tol=0.02, max_raw=4096, record_thickness=False,
                 attribute_tol=None):
        self.nx = int(nx)
        self.max_raw = max_raw
        self.record_thickness = bool(record_thickness)
        # Attribute tolerances; f_ch is always present and first.
        self.tol = {'f_ch': float(tol)}
        if attribute_tol:
            self.tol.update(attribute_tol)
        self.attributes = tuple(self.tol.keys())
        self.z = [[] for _ in range(self.nx)]
        self.data = {a: [[] for _ in range(self.nx)] for a in self.attributes}
        # Exact channel-deposit thickness sum(f_ch*dz), accumulated on aggradation
        # and drawn down on scour (volume-closure diagnostic).
        self._exact_thickness = np.zeros(self.nx) if self.record_thickness else None

    # -- backward-compatible view of the channel-deposit fraction ---------- #
    @property
    def f(self):
        """The raw ``f_ch`` value lists (one per node)."""
        return self.data['f_ch']

    def record(self, z, f_ch=None, **values):
        """Record one step's surface: append where the bed rose, scour where it
        fell.  ``z`` is the bed elevation array; ``f_ch`` (positional or keyword)
        and any further keyword are attributes (``age=...``, ``C_10Be=...``),
        each a scalar or array."""
        if f_ch is not None:
            values['f_ch'] = f_ch
        z = np.asarray(z, float)
        vb = {a: np.broadcast_to(np.asarray(values[a], float), z.shape)
              for a in self.attributes}
        for i in range(self.nx):
            zi = self.z[i]
            di = [self.data[a][i] for a in self.attributes]
            # Scour: remove anything at or above the new bed (unconformity).
            while zi and zi[-1] >= z[i]:
                if self._exact_thickness is not None and len(zi) >= 2:
                    self._exact_thickness[i] -= \
                        self.data['f_ch'][i][-1] * (zi[-1] - zi[-2])
                zi.pop()
                for d in di:
                    d.pop()
            if self._exact_thickness is not None and zi:
                self._exact_thickness[i] += vb['f_ch'][i] * (z[i] - zi[-1])
            zi.append(float(z[i]))
            for a, d in zip(self.attributes, di):
                d.append(float(vb[a][i]))
            if self.max_raw is not None and len(zi) > self.max_raw:
                self._compress(i)

    def _tolerance(self, attribute, f):
        """Absolute tolerance for ``attribute`` given its current values ``f``;
        a stored tolerance of ``None`` means 1% of the value span."""
        tol = self.tol[attribute]
        if tol is not None:
            return tol
        span = float(np.ptp(f)) if len(f) else 0.0
        return 0.01 * span if span > 0 else np.inf

    def _kept_indices(self, i):
        """Layer-boundary indices for node ``i``: the union of each attribute's
        Douglas-Peucker vertices, so every attribute is within its tolerance."""
        z = np.asarray(self.z[i], float)
        keep = set()
        for a in self.attributes:
            f = np.asarray(self.data[a][i], float)
            keep.update(douglas_peucker_vertical(z, f, self._tolerance(a, f)))
        return sorted(keep)

    def _compress(self, i):
        idx = self._kept_indices(i)
        self.z[i] = [self.z[i][k] for k in idx]
        for a in self.attributes:
            self.data[a][i] = [self.data[a][i][k] for k in idx]

    def columns(self, attribute='f_ch'):
        """Return the compressed columns as a list of ``(z, value)`` arrays for
        the named attribute, on the shared (union) layer boundaries."""
        out = []
        for i in range(self.nx):
            idx = self._kept_indices(i)
            out.append((np.asarray(self.z[i])[idx],
                        np.asarray(self.data[attribute][i])[idx]))
        return out

    def section(self, x, attribute='f_ch', ns=200):
        """
        Reconstruct an attribute on a structured, boundary-conforming mesh for
        plotting.  Each column is resampled onto ``ns`` levels between its basal
        (onset) and final surfaces, so the mesh follows the wedge boundaries
        exactly (no triangulation).  Returns ``(Xmesh, Zmesh, Vmesh)``, each of
        shape ``(ns, nx)``; ``Vmesh`` is ``nan`` where a column has no deposit.
        """
        x = np.asarray(x, float)
        s = np.linspace(0.0, 1.0, ns)
        Xm = np.broadcast_to(x, (ns, self.nx)).copy()
        Zm = np.full((ns, self.nx), np.nan)
        Vm = np.full((ns, self.nx), np.nan)
        for i, (zc, vc) in enumerate(zip(*self._columns_z_and(attribute))):
            if len(zc) < 2 or zc[-1] <= zc[0]:
                continue
            zlev = zc[0] + s * (zc[-1] - zc[0])
            Zm[:, i] = zlev
            Vm[:, i] = np.interp(zlev, zc, vc)
        return Xm, Zm, Vm

    def _columns_z_and(self, attribute):
        cols = self.columns(attribute)
        return [c[0] for c in cols], [c[1] for c in cols]

    def plot_section(self, x, attribute='f_ch', ax=None, cmap='viridis', ns=200,
                     vmin=None, vmax=None):
        """Plot the reconstructed ``attribute(x, z)`` section (Gouraud-shaded;
        vector in PDF/SVG).  ``f_ch`` defaults to a 0-1 colour range; other
        attributes autoscale unless ``vmin``/``vmax`` are given.  Returns the
        ``(ax, QuadMesh)``."""
        from matplotlib import pyplot as plt
        if ax is None:
            _, ax = plt.subplots(figsize=(9, 5))
        Xm, Zm, Vm = self.section(x, attribute=attribute, ns=ns)
        if attribute == 'f_ch' and vmin is None and vmax is None:
            vmin, vmax = 0.0, 1.0
        pc = ax.pcolormesh(Xm, Zm, Vm, shading='gouraud', cmap=cmap,
                           vmin=vmin, vmax=vmax)
        ax.set_xlabel('Downstream distance')
        ax.set_ylabel('Elevation')
        return ax, pc

    def channel_deposit_thickness(self):
        """Channel-deposit-equivalent thickness ``integral(f_ch dz)`` per node,
        from the compressed section (trapezoidal)."""
        out = np.zeros(self.nx)
        for i, (zc, fc) in enumerate(self.columns('f_ch')):
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
