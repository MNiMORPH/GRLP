"""
Clean an extracted river-network graph into a GRLP-ready input.

This is a **preprocessing step**, run before GRLP, not part of the model: it
takes a river network already extracted from a DEM (a NetworkX node-link graph,
e.g. a GRASS GIS export) and cleans the along-edge elevation profiles --
despiking, smoothing, enforcing a monotone envelope, and optionally coarsening
-- to produce the network that GRLP then reads. The cleaned graph is meant to be
written to disk and committed, so the exact model input is inspectable and
reproducible (see clean_network.py for a driver).

Pure NumPy / SciPy / NetworkX; no scikit-learn. Operates on a NetworkX graph in
which each node carries a single-element ``z`` list and each edge carries arrays
of ``x``, ``y``, ``z``, ``A`` (drainage area) and ``s`` (upstream distance) sampled
along its length.
"""

import numpy as np
from scipy.interpolate import interp1d
import json
import networkx as nx
import copy
from networkx.readwrite import json_graph


#############
# GRAPH I/O #
#############

def load_network_json(path):
    """Load a river-network graph from a node-link JSON file."""
    with open(path) as f:
        data = json.load(f)
    return json_graph.node_link_graph(data)


def _to_jsonable(obj):
    """
    Use this for export.
    Recursively convert objects into JSON-serializable forms:
    - numpy arrays -> lists
    - numpy scalars -> Python scalars
    - containers (dict/list/tuple) -> converted elementwise
    """

    # --- NumPy arrays ---
    if isinstance(obj, np.ndarray):
        # convert to list, then recurse (in case nested arrays, objects)
        return _to_jsonable(obj.tolist())

    # --- NumPy scalar types (float32, int64, etc.) ---
    if isinstance(obj, np.generic):
        return obj.item()  # returns a plain Python scalar

    # --- Containers: list/tuple ---
    if isinstance(obj, (list, tuple)):
        return type(obj)(_to_jsonable(x) for x in obj)

    # --- Dicts ---
    if isinstance(obj, dict):
        # Ensure values are converted; keys should already be JSON-OK
        return {k: _to_jsonable(v) for k, v in obj.items()}

    # Everything else is returned as-is; must already be JSON-serializable
    return obj


def make_json_safe_graph(G):
    """
    Use this for export
    Return a deep-copied version of G where:
    - all numpy arrays are converted to lists
    - all numpy scalar types are converted to Python scalars
    - node, edge, and graph attributes are processed
    """

    # deep copy so we don't mutate the original
    H = copy.deepcopy(G)

    # graph-level attributes
    for k in list(H.graph.keys()):
        H.graph[k] = _to_jsonable(H.graph[k])

    # node attributes
    for n, data in H.nodes(data=True):
        for k in list(data.keys()):
            data[k] = _to_jsonable(data[k])

    # edge attributes
    for u, v, data in H.edges(data=True):
        for k in list(data.keys()):
            data[k] = _to_jsonable(data[k])

    return H


def save_network_json(G, path):
    """Write a river-network graph to a node-link JSON file (numpy-safe)."""
    Hout = make_json_safe_graph(G)
    data = json_graph.node_link_data(Hout)
    with open(path, "w") as f:
        json.dump(data, f, indent=2)


#############
# UTILITIES #
#############

def bfs_upward(G, start):
    """
    Breadth-first search going upwards.
    Use this to update overall distance upstream from all outlets
    in a single sweep through the network, relying on those network components
    closer to the outlet to be updated first
    """
    R = G.reverse(copy=False)  # just a view, no data duplication
    for node in nx.bfs_tree(R, start):
        yield node


###############################
# 1D PROFILE OPERATIONS       #
###############################

def despike_1d(y, window=2, k=3.5):
    """
    Remove local outliers from a 1D array using a sliding median + MAD.

    y: 1D numpy array
    window: radius for local window (actual window size = 2*window + 1)
    k: threshold in units of local MAD
    """
    y = np.asarray(y, dtype=float)
    n = len(y)
    if n < 3:
        return y.copy()

    x = y.copy()
    is_outlier = np.zeros(n, dtype=bool)

    for i in range(1, n - 1):  # typically keep endpoints as-is here
        lo = max(0, i - window)
        hi = min(n, i + window + 1)
        neighborhood = y[lo:hi]

        median = np.median(neighborhood)
        mad = np.median(np.abs(neighborhood - median))
        if mad == 0:
            continue

        z = 0.6745 * (y[i] - median) / mad
        if np.abs(z) > k:
            is_outlier[i] = True

    # Replace outliers with local median of non-outlier neighbors
    for i in np.where(is_outlier)[0]:
        lo = max(0, i - window)
        hi = min(n, i + window + 1)
        neighborhood = x[lo:hi]
        good = ~is_outlier[lo:hi]
        if np.any(good):
            x[i] = np.median(neighborhood[good])
        else:
            x[i] = np.median(neighborhood)  # fallback

    return x


def smooth_1d_fixed_ends(y, left_val, right_val, alpha=0.35, n_iters=15):
    """
    Diffusive smoothing of a 1D signal with fixed endpoints (Dirichlet BC).

    y: 1D numpy array (original data along the edge)
    left_val, right_val: fixed boundary values from the nodes
    alpha: smoothing strength per iteration (0 < alpha < 1)
    n_iters: number of iterations
    """
    y = np.asarray(y, dtype=float)
    n = len(y)
    if n == 0:
        return y.copy()
    if n == 1:
        # Degenerate edge; just pick something between the two nodes
        return np.array([(left_val + right_val) / 2.0])

    x = y.copy()
    x[0] = left_val
    x[-1] = right_val

    for _ in range(n_iters):
        x_new = x.copy()
        # update only interior points
        for i in range(1, n - 1):
            neighbor_mean = 0.5 * (x[i - 1] + x[i + 1])
            x_new[i] = (1 - alpha) * x[i] + alpha * neighbor_mean
        # re-enforce fixed endpoints each iteration
        x_new[0] = left_val
        x_new[-1] = right_val
        x = x_new

    return x


def moving_average(y, k):
    if k is None or k <= 1:
        return np.asarray(y, dtype=float)
    k = int(k)
    kernel = np.ones(k) / k
    return np.convolve(y, kernel, mode="same")


# DOWNSTREAM CLAMPING -- MY IDEA WITH THE CHAT IMPLEMENTATION

def running_min_downstream(values):
    """
    Enforce a monotonically non-increasing profile by walking forward
    and clamping everything to the lowest value seen so far.

    values: 1D array-like
    returns: 1D numpy array (same length), monotone non-increasing
    """
    vals = np.asarray(values, dtype=float)
    if vals.size == 0:
        return vals.copy()

    out = np.empty_like(vals)
    current_min = vals[0]
    out[0] = current_min

    for i in range(1, len(vals)):
        if vals[i] < current_min:
            current_min = vals[i]
        out[i] = current_min

    return out


# Should not need
def running_max_upstream(values):
    """
    Enforce a monotonically non-decreasing profile by walking forward
    and clamping everything to the highest value seen so far.
    """
    vals = np.asarray(values, dtype=float)
    if vals.size == 0:
        return vals.copy()

    out = np.empty_like(vals)
    current_max = vals[0]
    out[0] = current_max

    for i in range(1, len(vals)):
        if vals[i] > current_max:
            current_max = vals[i]
        out[i] = current_max

    return out


def monotone_envelope_edge(values, parent_val, child_val, enforce_child=False):
    """
    Apply a running-min / running-max envelope along one edge.

    values      : 1D array of edge samples, along parent->child
    parent_val  : scalar elevation at parent node
    child_val   : scalar elevation at child node
    enforce_child : if True, do a final affine rescale so last sample == child_val

    Returns a 1D numpy array, monotone between parent_val and child_val.
    """
    vals = np.asarray(values, dtype=float)
    n = len(vals)
    if n == 0:
        return vals.copy()
    if n == 1:
        return np.array([(parent_val + child_val) / 2.0], dtype=float)

    # Initialize: enforce parent at first sample
    vals = vals.copy()
    vals[0] = parent_val

    downhill = parent_val > child_val  # elevation decreases along the edge?

    if downhill:
        # running min (non-increasing)
        out = running_min_downstream(vals)
    else:
        # running max (non-decreasing)
        out = running_max_upstream(vals)

    if enforce_child:
        # Affine rescale so out[-1] exactly equals child_val
        # This preserves monotonicity.
        eps = 1e-12
        a = (child_val - parent_val) / (out[-1] - out[0] + eps)
        b = parent_val - a * out[0]
        out = a * out + b

    return out


def symmetric_running_mean_edge(values, window=51, parent_val=None, child_val=None):
    """
    Symmetric (centered) running mean on a 1D array.

    - window: odd integer (we'll force it to be odd).
    - Near boundaries, the window is truncated to available indices.
    - If parent_val / child_val are given, we clamp the first and last
      samples back to those exact values after smoothing.
    """
    vals = np.asarray(values, dtype=float)
    n = len(vals)
    if n == 0 or window <= 1:
        # optionally clamp endpoints
        out = vals.copy()
        if parent_val is not None and n > 0:
            out[0] = parent_val
        if child_val is not None and n > 1:
            out[-1] = child_val
        return out

    # ensure window is odd
    if window % 2 == 0:
        window += 1
    half = window // 2

    out = np.empty_like(vals)

    for i in range(n):
        left  = max(0, i - half)
        right = min(n, i + half + 1)  # slice is [left:right)
        out[i] = vals[left:right].mean()

    # optionally re-attach exact node values at ends
    if parent_val is not None:
        out[0] = parent_val
    if child_val is not None and n > 1:
        out[-1] = child_val

    return out


##############################
# ALL-EDGES (GRAPH) WRAPPERS #
##############################

def clean_edge_values(
    G,
    node_attr="value",
    edge_attr="values",
    despike_window=2,
    despike_k=3.5,
    alpha=0.35,
    n_iters=15
):
    """
    For each edge in G:
      - take its values array
      - despike locally
      - smooth with fixed endpoints given by node_attr of endpoints
    """
    for u, v, data in G.edges(data=True):
        if edge_attr not in data:
            continue

        edge_vals = np.asarray(data[edge_attr], dtype=float)
        if edge_vals.size == 0:
            continue

        # node values (single-item lists)
        left_val = float(G.nodes[u][node_attr][0])
        right_val = float(G.nodes[v][node_attr][0])

        # 1) despike along the edge
        de_noised = despike_1d(edge_vals, window=despike_window, k=despike_k)

        # 2) smooth with fixed endpoint values from nodes
        smoothed = smooth_1d_fixed_ends(de_noised, left_val, right_val,
                                        alpha=alpha, n_iters=n_iters)

        # write back
        data[edge_attr] = smoothed.tolist()

    return G


def monotone_envelope_all_edges(
    G,
    node_attr="value",
    edge_attr="values",
    enforce_child=False,
):
    """
    For each directed edge (u->v), replace its 'values' with a
    strongly-smoothed monotone envelope between node elevations.
    """
    for u, v, data in G.edges(data=True):
        if edge_attr not in data:
            continue

        vals = data[edge_attr]
        if not vals:
            continue

        parent_val = float(G.nodes[u][node_attr][0])
        child_val  = float(G.nodes[v][node_attr][0])

        smoothed = monotone_envelope_edge(
            vals,
            parent_val,
            child_val,
            enforce_child=enforce_child,
        )

        data[edge_attr] = smoothed.tolist()

    return G


def symmetric_running_mean_all_edges(
    G,
    node_attr="value",
    edge_attr="values",
    window=51,
    clamp_endpoints=True,
):
    """
    Apply a symmetric running mean to every edge's values.

    Run this AFTER you've already enforced monotonicity
    (e.g. with monotone_envelope_all_edges).
    """
    for u, v, data in G.edges(data=True):
        if edge_attr not in data:
            continue

        vals = data[edge_attr]
        if not vals:
            continue

        parent_val = float(G.nodes[u][node_attr][0])
        child_val  = float(G.nodes[v][node_attr][0])

        smoothed = symmetric_running_mean_edge(
            vals,
            window=window,
            parent_val=parent_val if clamp_endpoints else None,
            child_val=child_val  if clamp_endpoints else None,
        )

        data[edge_attr] = smoothed.tolist()

    return G


def downsample_edge_along_s(
    G,
    keep_fraction=0.1,
    attrs=("x", "y", "z", "A", "s"),
):
    """
    Return a NEW graph where each edge has had (x,y,z,A,s) resampled
    along s so that:

      - s is evenly spaced between min(s) and max(s)
      - endpoints are excluded (inset by one spacing each side)
      - ~ keep_fraction of original points remain
      - x, y, z, A are interpolated on s_new
      - Original graph G is left untouched
    """

    # Make a deep copy so we don't modify original
    H = copy.deepcopy(G)

    x_attr, y_attr, z_attr, A_attr, s_attr = attrs

    for u, v, data in H.edges(data=True):

        # Retrieve original arrays
        s = np.asarray(data.get(s_attr, None))
        x = np.asarray(data.get(x_attr, None))
        y = np.asarray(data.get(y_attr, None))
        z = np.asarray(data.get(z_attr, None))
        A = np.asarray(data.get(A_attr, None))

        # Skip edges missing any attribute
        if any(arr is None for arr in (s, x, y, z, A)):
            continue

        # All arrays must have same length
        n = len(s)
        if not (len(x) == len(y) == len(z) == len(A) == n):
            continue

        # Need >=3 points to create interior
        if n < 3:
            continue

        # Sort by s
        idx = np.argsort(s)
        s = s[idx]
        x = x[idx]
        y = y[idx]
        z = z[idx]
        A = A[idx]

        # How many interior samples?
        n_new = int(round(keep_fraction * n))
        n_new = max(1, n_new)

        # Even spacing in s, excluding endpoints
        s_min, s_max = s[0], s[-1]
        ds = (s_max - s_min) / (n_new + 1)

        # Interior points only
        s_new = s_min + ds * np.arange(1, n_new + 1)
        s_new = s_new[::-1]  # Maintain original descending direction

        # Linear interpolators
        fx = interp1d(s, x, kind="linear")
        fy = interp1d(s, y, kind="linear")
        fz = interp1d(s, z, kind="linear")
        fA = interp1d(s, A, kind="linear")

        # Interpolate
        x_new = fx(s_new)
        y_new = fy(s_new)
        z_new = fz(s_new)
        A_new = fA(s_new)

        # Store back on the NEW graph
        data[s_attr] = s_new.tolist()
        data[x_attr] = x_new.tolist()
        data[y_attr] = y_new.tolist()
        data[z_attr] = z_new.tolist()
        data[A_attr] = A_new.tolist()

    return H


############
# PIPELINE #
############

def preprocess_network(
    G,
    node_attr="z",
    edge_attr="z",
    enforce_child=False,
    smooth_window=101,
    clamp_endpoints=True,
    coarsen=True,
    keep_fraction=0.1,
    coarsen_attrs=("x", "y", "z", "A", "s"),
):
    """
    Standard cleaning pipeline for a DEM-extracted river network: enforce a
    monotone elevation envelope along each edge, smooth it with a symmetric
    running mean, then (optionally) coarsen the along-edge sampling.

    Returns a NEW graph; G is not modified. For an alternative despike +
    diffusive-smooth pass, call clean_edge_values directly instead.
    """
    H = copy.deepcopy(G)
    monotone_envelope_all_edges(H, node_attr=node_attr, edge_attr=edge_attr,
                                enforce_child=enforce_child)
    symmetric_running_mean_all_edges(H, node_attr=node_attr, edge_attr=edge_attr,
                                     window=smooth_window,
                                     clamp_endpoints=clamp_endpoints)
    if coarsen:
        H = downsample_edge_along_s(H, keep_fraction=keep_fraction,
                                    attrs=coarsen_attrs)
    return H
