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


net_toplot = clean_edge_values(G, 'z', 'z', despike_window=200, despike_k=2., alpha=0.35, n_iters=200)

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

###################
# COARSEN NETWORK #
###################

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

H = downsample_edge_along_s_new(H0)

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

