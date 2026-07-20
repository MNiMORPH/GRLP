import numpy as np
from scipy.interpolate import interp1d
import json
import networkx as nx
import copy
from networkx.readwrite import json_graph
from matplotlib import pyplot as plt

#########################
# LOAD JSON AS NETWORKX #
#########################

homedir = '/home/awickert/Dropbox/Papers/InProgress/WhitewaterHapp/GIS/'

with open(homedir+"grass-whitewater-3DEP-sand.json") as f:
    data = json.load(f)

# Convert JSON structure back to a graph
G = json_graph.node_link_graph(data)

# IMPORTANT
# Convert outlet to real value to not have nan propagate up network
# Use a reasonable sand-bed slope and 2 m cell size
G.nodes[0]['z'] = [ G.edges[20,0]['z'][-1] - 4E-4 ]

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

###############################################
# SMOOTH AND DESPIKE VALUES ON NETWORKX EDGES #
###############################################

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

import networkx as nx

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

H0 = monotone_envelope_all_edges(G, enforce_child=False, node_attr='z', edge_attr='z')



# NOW SMOOTH THESE

import numpy as np
from sklearn.isotonic import IsotonicRegression


def moving_average(y, k):
    if k is None or k <= 1:
        return np.asarray(y, dtype=float)
    k = int(k)
    kernel = np.ones(k) / k
    return np.convolve(y, kernel, mode="same")


def smooth_monotone_steps(
    y_edge,
    parent_val,
    child_val,
    window=51,           # big window => strong smoothing
    extra_bins=None      # optional further coarsening
):
    """
    Take a monotone edge profile between parent_val and child_val and
    make it smoother but still monotone.

    y_edge    : 1D array, already monotone along parent->child
    parent_val: scalar at start node
    child_val : scalar at end node
    window    : moving-average window (odd integer, e.g. 31, 51, 101)
    extra_bins: if not None, additionally bin before isotonic for even more smoothing
    """
    y_edge = np.asarray(y_edge, dtype=float)
    n = len(y_edge)
    if n == 0:
        return y_edge.copy()
    if n == 1:
        return np.array([(parent_val + child_val) / 2.0], dtype=float)

    # Decide direction from node values (elevation)
    downhill = parent_val > child_val   # elevation should decrease?

    # Flip downhill to uphill so we can always use increasing=True
    if downhill:
        y = -y_edge
        p = -parent_val
        c = -child_val
    else:
        y = y_edge.copy()
        p = parent_val
        c = child_val

    # 1) Strong moving average → kills big steps
    y_ma = moving_average(y, window)

    # 2) Optional coarse binning for *even stronger* smoothing
    t = np.linspace(0.0, 1.0, n)
    X_fit, y_fit = t, y_ma

    if extra_bins is not None and extra_bins < n:
        bins = np.linspace(0.0, 1.0, extra_bins + 1)
        idx = np.digitize(t, bins) - 1  # 0..extra_bins-1

        Xb, yb = [], []
        for b in range(extra_bins):
            mask = idx == b
            if not np.any(mask):
                continue
            Xb.append(t[mask].mean())
            yb.append(y_ma[mask].mean())
        X_fit = np.array(Xb)
        y_fit = np.array(yb)

    # Sort just in case
    order = np.argsort(X_fit)
    X_fit = X_fit[order]
    y_fit = y_fit[order]

    # 3) Isotonic regression to restore exact monotonicity
    ir = IsotonicRegression(increasing=True, out_of_bounds="clip")
    ir.fit(X_fit, y_fit)
    y_iso = ir.predict(t)

    # 4) Rescale to hit node values exactly
    eps = 1e-12
    if np.allclose(y_iso[-1], y_iso[0], atol=1e-8):
        # If isotonic is flat, just use a straight line between nodes
        y_rescaled = np.linspace(p, c, n)
    else:
        a = (c - p) / (y_iso[-1] - y_iso[0] + eps)
        b = p - a * y_iso[0]
        y_rescaled = a * y_iso + b

    # 5) Flip back if this was a downhill edge
    if downhill:
        y_final = -y_rescaled
    else:
        y_final = y_rescaled

    return y_final

def smooth_steps_all_edges(
    G,
    node_attr="value",
    edge_attr="values",
    window=51,
    bin_ratio=0.05,
    min_bins=20,
    max_bins=80,
):
    """
    Second-stage smoother: takes already-monotone edge values and
    rounds off large steps into smoother, still-monotone curves.
    """
    for u, v, data in G.edges(data=True):
        if edge_attr not in data:
            continue

        vals = np.asarray(data[edge_attr], dtype=float)
        n = vals.size
        if n == 0:
            continue

        parent_val = float(G.nodes[u][node_attr][0])
        child_val  = float(G.nodes[v][node_attr][0])

        # choose number of bins proportional to n
        extra_bins = int(bin_ratio * n)
        extra_bins = max(min_bins, extra_bins)
        extra_bins = min(max_bins, extra_bins, n - 1)

        smoothed = smooth_monotone_steps(
            vals,
            parent_val,
            child_val,
            window=window,
            extra_bins=extra_bins if extra_bins >= 2 else None,
        )

        data[edge_attr] = smoothed.tolist()

    return G

H1 = monotone_envelope_all_edges(H0, node_attr="z", edge_attr="z",
                                enforce_child=False)


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

H2 = symmetric_running_mean_all_edges(H1,
                                     node_attr="z",
                                     edge_attr="z",
                                     window=101,          # try 51–151
                                     clamp_endpoints=True)

# It works!


# FINALLY, coarsen.

def downsample_edge_along_s_new(
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
        s_new = s_new[::-1] # Maintain original descending direction
        #print("HI!")
        #print(s_new)

        # Linear interpolators
        fx = interp1d(s, x, kind="linear")
        fy = interp1d(s, y, kind="linear")
        fz = interp1d(s, z, kind="linear")
        fA = interp1d(s, A, kind="linear")

        # Interpolate
        print( len(s_new) )
        x_new = fx(s_new)
        print( len(x_new) )
        y_new = fy(s_new)
        print( len(y_new) )
        z_new = fz(s_new)
        print( len(z_new) )
        A_new = fA(s_new)
        print( len(A_new) )

        # Store back on the NEW graph
        data[s_attr] = s_new.tolist()
        data[x_attr] = x_new.tolist()
        data[y_attr] = y_new.tolist()
        data[z_attr] = z_new.tolist()
        data[A_attr] = A_new.tolist()

    return H

H = downsample_edge_along_s_new(H2)


#######################
# PLOTTING EVERYTHING #
#######################

net_toplot = H

for n in bfs_upward(net_toplot, 0):
    edges = net_toplot.in_edges(n)
    print(edges)


#plt.ion()
plt.figure()
for n in bfs_upward(net_toplot, 0):
    edges = net_toplot.in_edges(n)
    for parent, child in edges:
        plt.plot( net_toplot.edges[parent,child]['s'], net_toplot.edges[parent,child]['z'], 'k-', linewidth=3, alpha=1 )
        plt.plot( net_toplot.nodes[parent]['s'], net_toplot.nodes[parent]['z'], 'ko', alpha=1 )
plt.xlabel('Upstream distance')
plt.ylabel('Elevation') # hard-coded for now
plt.show()



##########
# EXPORT #
##########

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

# IMPORTANT -- maybe. Maybe it would be fine without.
# Maybe I should just give the outlet nodes real elevations.
# But is that possible?
# Return the outlet node to nan after finishing the interpolation
#H.nodes[0]['z'] = [ np.nan ]

# KEEP AS REAL VALUE FOR VALID RIVER MOUTH.
# FIX IN EXPORT CODE.

# Create JSON-safe graph for export
# This converts numpy arrays to lists
outjson = 'whitewater-3DEP-sand-smoothed-coarsened.json'
Hout = make_json_safe_graph(H)

# export H, the JSON-safe data file
if outjson is not None:
    print("Exporting JSON of NetworkX river-network object.")
    data = json_graph.node_link_data(Hout)
    with open(outjson, "w") as f:
        json.dump(data, f, indent=2)
    print("Export complete.")

