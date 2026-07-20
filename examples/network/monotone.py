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

import numpy as np
from sklearn.isotonic import IsotonicRegression

import numpy as np
from sklearn.isotonic import IsotonicRegression

def monotone_smooth_edge_sklearn(
    y_edge,
    left_val,
    right_val,
    node_weight=20.0,
    edge_weight=1.0,
    n_knots=None,
    increasing=True,
    enforce_endpoints=True,
):
    """
    Smooth a 1D noisy edge signal into a monotone curve using sklearn's
    IsotonicRegression.

    y_edge      : 1D array of samples along the edge
    left_val    : scalar value at start node
    right_val   : scalar value at end node
    node_weight : relative weight of node values vs edge samples
    edge_weight : weight of each edge sample
    n_knots     : if not None, subsample/average to this many points before fitting
    increasing  : True => non-decreasing, False => non-increasing
    enforce_endpoints : if True, rescale so first/last exactly match node values
    """
    y_edge = np.asarray(y_edge, dtype=float)
    n = len(y_edge)
    if n == 0:
        return y_edge.copy()
    if n == 1:
        return np.array([(left_val + right_val) / 2.0])

    # Parameter along edge (0..1); replace with cumulative distance if you have it
    t_edge = np.linspace(0.0, 1.0, n)

    # Build training set: node anchors + ALL edge samples
    X = np.concatenate(([0.0], t_edge, [1.0]))
    y = np.concatenate(([left_val], y_edge, [right_val]))

    w = np.full_like(y, edge_weight, dtype=float)
    w[0] = node_weight
    w[-1] = node_weight

    # Optional: reduce to n_knots by binning in t
    if n_knots is not None and n_knots < X.size:
        bins = np.linspace(0.0, 1.0, n_knots + 1)
        idx = np.digitize(X, bins) - 1  # 0..n_knots-1

        Xb, yb, wb = [], [], []
        for b in range(n_knots):
            mask = idx == b
            if not np.any(mask):
                continue
            w_sum = w[mask].sum()
            if w_sum == 0:
                continue
            Xb.append(X[mask].mean())
            yb.append((w[mask] * y[mask]).sum() / w_sum)
            wb.append(w_sum)

        X_fit = np.array(Xb)
        y_fit = np.array(yb)
        w_fit = np.array(wb)
    else:
        X_fit, y_fit, w_fit = X, y, w

    # Ensure X_fit is sorted (IsotonicRegression expects sorted X)
    order = np.argsort(X_fit)
    X_fit = X_fit[order]
    y_fit = y_fit[order]
    w_fit = w_fit[order]

    ir = IsotonicRegression(
        increasing=increasing,
        out_of_bounds="clip",
    )
    ir.fit(X_fit, y_fit, sample_weight=w_fit)

    y_smooth = ir.predict(t_edge)

    if enforce_endpoints:
        # Affine rescaling so that endpoints match exactly
        eps = 1e-12
        a = (right_val - left_val) / (y_smooth[-1] - y_smooth[0] + eps)
        b = left_val - a * y_smooth[0]
        y_smooth = a * y_smooth + b

    return y_smooth

def smooth_all_edges_monotone(
    G,
    node_attr="value",
    edge_attr="values",
    node_weight=20.0,
    edge_weight=1.0,
    n_knots=None,
    force_increasing=True,
):
    """
    Replace each edge's 'values' array with a smooth, monotone version
    using isotonic regression anchored to node values.
    """
    for u, v, data in G.edges(data=True):
        if edge_attr not in data:
            continue

        y_edge = np.asarray(data[edge_attr], dtype=float)
        if y_edge.size == 0:
            continue

        left_val = float(G.nodes[u][node_attr][0])
        right_val = float(G.nodes[v][node_attr][0])

        # Decide direction
        if force_increasing:
            inc = True
        else:
            inc = left_val <= right_val  # auto choose inc/dec based on endpoints

        y_smooth = monotone_smooth_edge_sklearn(
            y_edge,
            left_val,
            right_val,
            node_weight=node_weight,
            edge_weight=edge_weight,
            n_knots=n_knots,
            increasing=inc,
            enforce_endpoints=True,
        )

        data[edge_attr] = y_smooth.tolist()

    return G

"""
H = smooth_all_edges_monotone(
    G,
    node_attr="z",
    edge_attr="z",
    node_weight=100.0,   # not crazy huge
    edge_weight=1.0,
    n_knots=20,        # or None to keep full resolution
    force_increasing=False,
)
"""






import numpy as np
from sklearn.isotonic import IsotonicRegression

def monotone_smooth_edge_elevation(
    y_edge,
    parent_val,
    child_val,
    pre_smooth_window=None,
):
    """
    Smooth a 1D edge signal into a monotone curve between parent_val and child_val.

    - If child_val >= parent_val: enforce non-decreasing along the edge.
    - If child_val <  parent_val: enforce non-increasing along the edge.
    - Shape comes from isotonic regression on edge data (optionally pre-smoothed).
    - Then we rescale so endpoints match parent_val and child_val exactly.
    """
    y_edge = np.asarray(y_edge, dtype=float)
    n = len(y_edge)
    if n == 0:
        return y_edge.copy()
    if n == 1:
        # trivial: just interpolate between nodes
        return np.array([(parent_val + child_val) / 2.0], dtype=float)

    # Decide direction by node elevations
    increasing = child_val >= parent_val

    # Optional pre-smoothing (moving average) to kill high-frequency noise
    y_for_fit = y_edge.copy()
    if pre_smooth_window is not None and pre_smooth_window > 1:
        k = pre_smooth_window
        kernel = np.ones(k) / k
        y_for_fit = np.convolve(y_for_fit, kernel, mode="same")

    # Parameter along edge (0..1) – can be distance if you have it
    t = np.linspace(0.0, 1.0, n)

    # If edge should be decreasing, flip the data so isotonic can stay "increasing"
    flipped = False
    if not increasing:
        flipped = True
        y_for_fit = -y_for_fit
        parent_val, child_val = -parent_val, -child_val   # flip endpoints too

    # Isotonic regression on edge data ONLY (no endpoints here)
    ir = IsotonicRegression(increasing=True, out_of_bounds="clip")
    ir.fit(t, y_for_fit)
    y_iso = ir.predict(t)

    # If fit is (numerically) constant, fall back to straight line between endpoints
    if np.allclose(y_iso, y_iso[0], atol=1e-8):
        y_line = np.linspace(parent_val, child_val, n)
    else:
        # Affine rescale isotonic shape to hit parent/child values exactly
        eps = 1e-12
        a = (child_val - parent_val) / (y_iso[-1] - y_iso[0] + eps)
        b = parent_val - a * y_iso[0]
        y_line = a * y_iso + b

    # If we flipped for a decreasing edge, flip back
    if flipped:
        y_line = -y_line

    return y_line


import networkx as nx

def smooth_all_edges_monotone_elevation(
    G,
    node_attr="value",
    edge_attr="values",
    pre_smooth_window=None,
):
    """
    For each directed edge (parent -> child), replace its 'values' list with a
    smooth monotone profile between the parent and child node elevations.
    """
    for u, v, data in G.edges(data=True):
        if edge_attr not in data:
            continue

        y_edge = np.asarray(data[edge_attr], dtype=float)
        if y_edge.size == 0:
            continue

        parent_val = float(G.nodes[u][node_attr][0])
        child_val  = float(G.nodes[v][node_attr][0])

        y_smooth = monotone_smooth_edge_elevation(
            y_edge,
            parent_val,
            child_val,
            pre_smooth_window=pre_smooth_window,
        )

        data[edge_attr] = y_smooth.tolist()

    return G


#H = smooth_all_edges_monotone_elevation(G, pre_smooth_window=1E10)







import numpy as np
from sklearn.isotonic import IsotonicRegression


def moving_average(y, k):
    if k is None or k <= 1:
        return y
    k = int(k)
    kernel = np.ones(k) / k
    return np.convolve(y, kernel, mode="same")


def super_smooth_edge_monotone(
    y_edge,
    parent_val,
    child_val,
    bin_ratio=0.05,      # fraction of samples used as bins (5%)
    min_bins=20,
    max_bins=100,
    pre_ma_window=15,    # strong pre-smoothing
    post_ma_window=None  # optional extra post-smoothing
):
    """
    VERY strong monotone smoothing between parent_val and child_val.

    - Pre-smooth with a wide moving average.
    - Bin along the edge, average per bin.
    - Run isotonic regression on bin-averaged data.
    - Evaluate on full resolution.
    - Rescale so endpoints match node values exactly.
    - Handles both uphill and downhill edges.
    """
    y_edge = np.asarray(y_edge, dtype=float)
    n = len(y_edge)
    if n == 0:
        return y_edge.copy()
    if n == 1:
        return np.array([(parent_val + child_val) / 2.0], dtype=float)

    # Decide direction from node values (in elevation)
    increasing = child_val >= parent_val

    # Pre-smoothing: kill high-frequency spikes hard
    y_smooth = moving_average(y_edge, pre_ma_window)

    # Parameter along edge (0..1)
    t = np.linspace(0.0, 1.0, n)

    # If edge should be decreasing, flip signs so we can still do "increasing" fit
    flipped = False
    if not increasing:
        flipped = True
        y_smooth = -y_smooth
        parent_val, child_val = -parent_val, -child_val

    # ---- BINNING / COARSENING ----
    # Choose number of bins proportional to n, but clipped
    n_bins = int(bin_ratio * n)
    n_bins = max(min_bins, n_bins)
    n_bins = min(max_bins, n_bins)
    n_bins = min(n_bins, n)   # can't have more bins than samples

    # Bin edges in t
    bins = np.linspace(0.0, 1.0, n_bins + 1)
    idx = np.digitize(t, bins) - 1  # bin index 0..n_bins-1

    Xb, yb = [], []
    for b in range(n_bins):
        mask = idx == b
        if not np.any(mask):
            continue
        Xb.append(t[mask].mean())
        yb.append(y_smooth[mask].mean())

    Xb = np.array(Xb)
    yb = np.array(yb)

    # Guard: if for some reason binning collapses too much, fall back to straight line
    if len(Xb) < 2:
        line = np.linspace(parent_val, child_val, n)
        return -line if flipped else line

    # ---- ISOTONIC ON BINNED DATA ----
    order = np.argsort(Xb)
    Xb = Xb[order]
    yb = yb[order]

    ir = IsotonicRegression(increasing=True, out_of_bounds="clip")
    ir.fit(Xb, yb)

    # Predict on full grid
    y_iso = ir.predict(t)

    # Optional extra smoothing of isotonic output (still monotone because kernel >= 0)
    y_iso = moving_average(y_iso, post_ma_window or 1)

    # If isotonic is essentially flat, just use a straight line
    if np.allclose(y_iso, y_iso[0], atol=1e-8):
        y_line = np.linspace(parent_val, child_val, n)
    else:
        # Rescale isotonic "shape" to match node endpoints exactly
        eps = 1e-12
        a = (child_val - parent_val) / (y_iso[-1] - y_iso[0] + eps)
        b = parent_val - a * y_iso[0]
        y_line = a * y_iso + b

    # Flip back for downhill edges
    if flipped:
        y_line = -y_line

    return y_line


def super_smooth_all_edges(
    G,
    node_attr="value",
    edge_attr="values",
    bin_ratio=0.05,
    min_bins=20,
    max_bins=100,
    pre_ma_window=15,
    post_ma_window=None
):
    """
    Very strong monotone smoothing for every edge in G.
    """
    for u, v, data in G.edges(data=True):
        if edge_attr not in data:
            continue

        y_edge = np.asarray(data[edge_attr], dtype=float)
        if y_edge.size == 0:
            continue

        parent_val = float(G.nodes[u][node_attr][0])
        child_val  = float(G.nodes[v][node_attr][0])

        y_smooth = super_smooth_edge_monotone(
            y_edge,
            parent_val,
            child_val,
            bin_ratio=bin_ratio,
            min_bins=min_bins,
            max_bins=max_bins,
            pre_ma_window=pre_ma_window,
            post_ma_window=post_ma_window,
        )

        data[edge_attr] = y_smooth.tolist()

    return G


H = super_smooth_all_edges(
    G,
    node_attr="z",
    edge_attr="z",
    bin_ratio=0.3,   # 3% of samples per edge
    min_bins=20,
    max_bins=60,
    pre_ma_window=14, # heavy smoothing
    post_ma_window=7  # light final polish
)






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


