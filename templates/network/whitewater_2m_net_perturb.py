#! /usr/bin/python3

# Basic libraries
import numpy as np
from matplotlib import pyplot as plt
import importlib

# Libraries for NetworkX
import json
import networkx as nx
from networkx.readwrite import json_graph
from collections import deque

# Libraries for GRLP
import grlp

#plt.ion()
plt.ioff()


#%load_ext autoreload
#%autoreload

importlib.reload(grlp)
#del net


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

def subgraph_with_parents(G, node):
    # 1. All recursive parents of `node`
    parents = nx.ancestors(G, node)   # set of nodes

    # 2. Induced subgraph on parents + node itself
    nodes_to_keep = parents | {node}
    H = G.subgraph(nodes_to_keep).copy()   # copy() to get an independent graph

    return H

def subgraph_from_edge_child_side(G, u, v):
    """
    Subgraph containing:
    - the edge (u, v)
    - node v
    - all ancestors (parents, grandparents, ...) of v
    """
    # all recursive parents of v
    parents = nx.ancestors(G, v)

    nodes_to_keep = parents | {u, v}
    H = G.subgraph(nodes_to_keep).copy()
    return H

def converging_subgraph_for_edge(G, u, v, mark_attr="is_final"):
    """
    Return a subgraph H consisting of:
      - nodes u and v
      - all ancestors of u
      - all ancestors of v
      - all edges among those nodes

    Additionally, the edge (u, v) in H is marked with edge attribute `mark_attr = True`.
    """

    # all ancestors of u and v
    anc_u = nx.ancestors(G, u)
    anc_v = nx.ancestors(G, v)

    nodes_to_keep = anc_u | anc_v | {u, v}

    # induced subgraph (copy so we can safely modify attributes)
    H = G.subgraph(nodes_to_keep).copy()

    # ensure the final edge is present and mark it
    if H.has_edge(u, v):
        H.edges[u, v][mark_attr] = True
    else:
        # In case something weird happens, you could also re-add it from G:
        # H.add_edge(u, v, **G.edges[u, v])
        raise ValueError(f"Final edge ({u}, {v}) not found in subgraph")

    return H

def walk_one_branch_up(G, start, choose_parent=None):
    """
    Walk upward from `start` following ONE parent at each step.
    
    Parameters
    ----------
    G : DiGraph
        Directed graph with parent → child edges.
    start : node
        Starting node.
    choose_parent : function or None
        Optional selector: choose_parent(list_of_parents) -> parent
        If None, the first parent is chosen.
    
    Yields
    ------
    nodes on the path upward, including the start.
    """
    node = start
    yield node

    while True:
        parents = list(G.predecessors(node))

        if not parents:
            break  # reached the root

        # choose which parent to follow
        if choose_parent is None:
            parent = parents[0]  # default: take first
        else:
            parent = choose_parent(parents)

        yield parent

        node = parent

def one_branch_subtree(G, start, choose_parent=None):
    # 1. Get the node sequence along the branch (child → ... → root)
    path_nodes = list(walk_one_branch_up(G, start, choose_parent))
    
    # 2. Build a fresh DiGraph
    H = nx.DiGraph()
    
    # 3. Add nodes with attributes
    for n in path_nodes:
        H.add_node(n, **G.nodes[n])
    
    # 4. Add edges along the path, preserving attributes
    #    Remember: in G, direction is parent -> child
    #    path_nodes is [child, parent, grandparent, ...]
    for child, parent in zip(path_nodes[:-1], path_nodes[1:]):
        # edge in G is (parent, child)
        if G.has_edge(parent, child):
            H.add_edge(parent, child, **G.edges[parent, child])
        else:
            # only if the original graph might not be a pure tree
            raise ValueError(f"No edge {parent} -> {child} in G")
    
    return H

def one_then_full_ancestor_subgraph(G, start, choose_parent=None):
    """
    From `start`:
      - choose ONE of its parents (first step up)
      - above that parent, take the FULL ancestor tree (both parents, etc.)

    Returns:
      H : DiGraph
          subgraph induced by {start} + chosen parent + all its ancestors
      chosen_parent : node
          the parent selected at the first step
    """

    parents = list(G.predecessors(start))
    if not parents:
        raise ValueError(f"{start!r} has no parents")

    if choose_parent is None:
        chosen_parent = parents[0]
    else:
        chosen_parent = choose_parent(parents)

    # all ancestors of the chosen parent (both parents at every level above)
    anc = nx.ancestors(G, chosen_parent)

    nodes_to_keep = anc | {chosen_parent, start}
    H = G.subgraph(nodes_to_keep).copy()

    return H, chosen_parent

######################
# STANDARD VARIABLES #
######################

# Change later?
dt = 3.15E7 # 1 year
_B = 200 # uniform


#########################################
# IMPORT NETWORK: LOAD JSON AS NETWORKX #
#########################################

with open("/home/awickert/Dropbox/Papers/InProgress/WhitewaterHapp/GIS/whitewater-3DEP-sand-smoothed-coarsened.json") as f:
    data = json.load(f)

"""
with open("/home/awickert/Dropbox/Papers/InProgress/WhitewaterHapp/GIS/whitewater-sand-smoothed-coarsened.json") as f:
    data = json.load(f)

# See if the smoothing propagated the nan upwards
with open("/home/awickert/Dropbox/Papers/InProgress/WhitewaterHapp/GIS/grass-whitewater-sand.json") as f:
    data2 = json.load(f)
"""

# Convert JSON structure back to a graph
G = json_graph.node_link_graph(data)
"""
G2 = json_graph.node_link_graph(data2)
"""

#####################
# USE WHOLE NETWORK #
#####################

outlet_node = 0

H = G

#H, bottom_node = one_then_full_ancestor_subgraph(G, outlet_node, None)

"""
# Plot to test that it worked -- yes.
plt.figure()
#for n in bfs_upward(H, 1258):
for n in H:
    edges = H.out_edges(n)
    for parent, child in edges:
        print (parent, child)
        plt.plot( H.edges[parent,child]['s'], H.edges[parent,child]['z'], 'k-', linewidth=3, alpha=1 )
        plt.plot( H.nodes[parent]['s'], H.nodes[parent]['z'], 'ko', alpha=1 )
    plt.plot( H.nodes[outlet_node]['s'], H.nodes[outlet_node]['z'], 'ko', alpha=1 )
plt.xlabel('Upstream distance')
plt.ylabel('Elevation') # hard-coded for now
plt.show()
"""

##################################################
# SET UP FOR GRLP: NUMBER AND LENGTH OF SEGMENTS #
##################################################

# Number of segments
nseg = len(H.edges)

# Number of elements
numel = []
for u, v, data in H.edges(data=True):
    arr = data.get("x", None)
    if arr is None:
        # no attribute "x" on this edge, skip or append 0
        continue
    # assume arr is a list (or numpy array), take its length
    # Add one because of the upstream cell that is on the node
    # But that will be taken up by the downstream edge (channel segment)
    numel.append(len(arr)+1)

# We can probably use the nodes above for the upstream segment IDs
# These will then be empty for the upstream-most segments

##########################################################
# SEGMENT IDs: USE THE NODE VALUES UPSTREAM TO SET THESE #
##########################################################

"""
#if node != bottom_node:
raw_upstream_segment_IDs = []
for node in H.nodes():
    parents = list(H.predecessors(node))
    raw_upstream_segment_IDs.append(parents)   # [] if no parents
upstream_segment_IDs = raw_upstream_segment_IDs

raw_downstream_segment_IDs = []
for node in H.nodes():
    child = list(H.successors(node))
    # Outlet has no downstream node in the subset
    # Taken care of when we use H
    #if child[0] == outlet_node:
    #    child = []
    raw_downstream_segment_IDs.append(child)
downstream_segment_IDs = raw_downstream_segment_IDs

# Actually, I think we want downstream IDs to be the node itself
# This means that our segment names are set by their upstream ID
# Probably not, actually. This should be the next downstream
# toseg, fromseg... therefore, seg is implicitly its own number
raw_downstream_segment_IDs = []
for node in H.nodes():
    #child = list(H.successors(node))
    # Outlet has no downstream node in the subset
    # Taken care of when we use H
    #if child[0] == outlet_node:
    #    child = []
    if node == outlet_node:
        raw_downstream_segment_IDs.append([])
    else:
        raw_downstream_segment_IDs.append(node)
downstream_segment_IDs = raw_downstream_segment_IDs


for node in H.nodes():
    print(node)

# This is probably the real way to do it
# MAYBE NOT: I JUST NEED TO SKIP THE DOWNSTREAM-MOST NODE
# CONFIRMED: THIS SEQUENTIAL LIST IS NEEDED
i = 0
for node in H.nodes():
    H.nodes[node]['id'] = i
    i += 1

upstream_segment_IDs = []
for node in H.nodes():
    _par_id = []
    parents = list(H.predecessors(node))
    for _node in parents:
        _par_id.append( H.nodes[_node]['id'] )
    upstream_segment_IDs.append(_par_id)   # [] if no parents

downstream_segment_IDs = []
for node in H.nodes():
    _id = []
    child = list(H.successors(node))
    for _node in child:
        _id.append( H.nodes[_node]['id'] )
    # Outlet has no downstream node in the subset
    #if child[0] == outlet_node:
    #    child = []
    downstream_segment_IDs.append( _id )

# Actually, I think we want downstream IDs to be the node itself
# This means that our segment names are set by their upstream ID
downstream_segment_IDs = []
for node in H.nodes():
    _id = []
    if node == outlet_node:
        pass
    else:
        _id.append( H.nodes[node]['id'] )
    downstream_segment_IDs.append( _id )
"""


# WE RATHER MUST LOOP OVER EDGES.
# NEED TO FIX ON PERROT EXAMPLE
# Name each edge by its upstream marker.
# Walk up from it to get the two upstream rivers\

# First, ID for each segment, in sequence.
# This is the index in Python for that segment.
i = 0
for edge in H.edges():
    H.edges[edge]['id'] = i
    i += 1

upstream_segment_IDs = []
downstream_segment_IDs = []
for edge in H.edges():
    parent, child = edge   # edge = (u, v)
    downstream_segment = list(H.out_edges(child))
    if len(downstream_segment) > 1:
        print("ERROR! ERROR! ERROR! ERROR! ERROR! ERROR!")
    elif len(downstream_segment) == 0:
        downstream_segment_IDs.append([])
    else:
        downstream_segment = downstream_segment[0]
        downstream_segment_IDs.append( [H.edges[downstream_segment]['id']] )
        # G, to get that downstream-most one
        #downstream_segment_IDs.append( [G.edges[downstream_segment]['id']] )
for edge in H.edges():
    parent, child = edge   # edge = (u, v)
    upstream_segments = list(G.in_edges(parent))
    if len(upstream_segments) == 0:
        upstream_segment_IDs.append([])
    else:
        _upseg_ids = []
        for _upseg in upstream_segments:
            _upseg_ids.append( H.edges[_upseg]['id'] )
        upstream_segment_IDs.append( _upseg_ids )    



###################################################
# VARIABLES: THOSE ON LINKS + FROM THE NODE ABOVE #
###################################################

# Use G to get the full 10
def network_variable_prepend_list(varname):
    new_values = []
    for edge in H.edges():
        parent, child = edge   # edge = (u, v)
        px = H.nodes[parent][varname]
        arr = H.edges[edge][varname]
        new_values.append( np.array(px + arr) )
        #print(px)
    return new_values

# Downstream distance (s in network, x in GRLP)
x = network_variable_prepend_list('s')
# Flip x direction
for i in range(len(x)):
    x[i] *= -1

# Elevation (z)
z = network_variable_prepend_list('z')

# Discharge
# Start with drainage area
A_in_list = network_variable_prepend_list('A')
# Then just modify it simply
Q = []
# Drainage areas might go negative (offmap flow)
# I THINK THIS WAS CAUSING INSTABILITY
# IT WAS! PERHAPS MORE THAN THE SPIKY INPUTS...!
for _A in A_in_list:
    Q.append(np.abs(_A)/1E7)

# Valley width: let's leave constant for now. Defined above (_B)
# Actually, a bit more realism: let's make small (Trout, Beaver) and large
# (mainstem Whitewater) valleys.
B = []
for _A in A_in_list:
    _Aa = np.abs(_A)
    if (_Aa < 1E8).all():
        B.append(_A*0 + 200)
    else:
        B.append(_A*0 + 900)
    
#######################
# BOUNDARY CONDITIONS #
#######################

# Base level: from mouth node
x_bl = (H.nodes.get(outlet_node)['s'][0] - 2) * -1 # x direction flipped; IS SET TO 0 WHEN THIS IS ALSO THE LOCATION OF THE CELL ABOVE
z_bl = H.nodes.get(outlet_node)['z'][0]

# Upstream boundary condition: 1.5% grade
S0 = []
for node in H.nodes():
    parents = list(G.predecessors(node))
    if len(parents) == 0:
        _u, _v = list(H.out_edges(node))[0]
        _s = H.edges[_u, _v]["s"]
        _z = H.edges[_u, _v]["z"]
        _S = np.abs( ( _z[0] - _z[-1] ) / ( _s[0] - _s[-1] ) )
        S0.append( _S )


##############################
# INSTANTIATE NETWORK OBJECT #
##############################

# Instantiate network object
net = grlp.Network()
self = net # For debugging, etc.

# Initialize network by passing x, z, etc. information to model
net.initialize(
                config_file = None,
                x_bl = x_bl,
                z_bl = z_bl,
                S0 = S0,
                Q_s_0 = None,
                upstream_segment_IDs = upstream_segment_IDs,
                downstream_segment_IDs = downstream_segment_IDs,
                x = x,
                z = z,
                Q = Q,
                B = B,
                overwrite=False
                )

# Should do this above
net.set_niter(4)
net.get_z_lengths()

# For testing
#net.evolve_threshold_width_river_network(nt=100, dt=dt*10)

#"""
# For plotting
# WHEN RUN FOR NT=10, GET BACKWARDS SLOPE ON TRIBUTARY
# THIS IS WHERE WE NEED TO ADD IN CLOSED BASINS AS ANOTHER SEGMENT TYPE
# Was becoming singular because I had flow R to L instead of L to R

# Start
for lp in net.list_of_LongProfile_objects:
    # If not downstream-most segment
    if len( lp.downstream_segment_IDs ) > 0:
        for _id in lp.downstream_segment_IDs:
            dsseg = net.list_of_LongProfile_objects[_id]
            _xjoin = [lp.x[-1], dsseg.x[0]]
            _zjoin = [lp.z[-1], dsseg.z[0]]
            plt.plot(_xjoin, _zjoin, '-', color='b', linewidth=4, alpha=.5)
    else:
        plt.plot(lp.x_ext[0][-2:], lp.z_ext[0][-2:], '-', color='b', linewidth=4, alpha=.5)
    plt.plot(lp.x, lp.z, '-', color='b', linewidth=4, alpha=1.)#, label=lp.)

# 10 years of evolution under constant bankfull conditions
for ts in range(40):
    net.evolve_threshold_width_river_network(nt=1, dt=dt*400)
    #lp = net.list_of_LongProfile_objects[2]
    for lp in net.list_of_LongProfile_objects:
        # If not downstream-most segment
        if len( lp.downstream_segment_IDs ) > 0:
            for _id in lp.downstream_segment_IDs:
                dsseg = net.list_of_LongProfile_objects[_id]
                _xjoin = [lp.x[-1], dsseg.x[0]]
                _zjoin = [lp.z[-1], dsseg.z[0]]
                plt.plot(_xjoin, _zjoin, '-', color='.5', linewidth=1, alpha=.5)
        else:
            plt.plot(lp.x_ext[0][-2:], lp.z_ext[0][-2:], 'b-', linewidth=1, alpha=.5)
        plt.plot(lp.x, lp.z, '-', color='.5', linewidth=1, alpha=1.)#, label=lp.)
for lp in net.list_of_LongProfile_objects:
    # If not downstream-most segment
    if len( lp.downstream_segment_IDs ) > 0:
        for _id in lp.downstream_segment_IDs:
            dsseg = net.list_of_LongProfile_objects[_id]
            _xjoin = [lp.x[-1], dsseg.x[0]]
            _zjoin = [lp.z[-1], dsseg.z[0]]
            plt.plot(_xjoin, _zjoin, '-', color='.5', linewidth=4, alpha=.5)
    else:
        plt.plot(lp.x_ext[0][-2:], lp.z_ext[0][-2:], '-', color='.5', linewidth=4, alpha=.5)
    plt.plot(lp.x, lp.z, '-', color='.5', linewidth=4, alpha=1.)#, label=lp.)


# NOW LET'S PERTURB IT

# 0.2 years, up to 2 years
base_level_rise_rate = 0.02 # m / yr
for ts in range(100):
    # Assuming 1-year time steps: BRITTLE
    net.update_z_ext_external_downstream(z_bl + base_level_rise_rate * (ts + 1))    
    net.evolve_threshold_width_river_network(nt=1, dt=dt)
    #lp = net.list_of_LongProfile_objects[2]
    if ts % 10 == 0:
        for lp in net.list_of_LongProfile_objects:
            # If not downstream-most segment
            if len( lp.downstream_segment_IDs ) > 0:
                for _id in lp.downstream_segment_IDs:
                    dsseg = net.list_of_LongProfile_objects[_id]
                    _xjoin = [lp.x[-1], dsseg.x[0]]
                    _zjoin = [lp.z[-1], dsseg.z[0]]
                    plt.plot(_xjoin, _zjoin, '-', color='r', linewidth=1, alpha=.5)
            else:
                plt.plot(lp.x_ext[0][-2:], lp.z_ext[0][-2:], '-', color='r', linewidth=1, alpha=.5)
            plt.plot(lp.x, lp.z, '-', color='r', linewidth=1, alpha=1.)#, label=lp.)
for lp in net.list_of_LongProfile_objects:
    # If not downstream-most segment
    if len( lp.downstream_segment_IDs ) > 0:
        for _id in lp.downstream_segment_IDs:
            dsseg = net.list_of_LongProfile_objects[_id]
            _xjoin = [lp.x[-1], dsseg.x[0]]
            _zjoin = [lp.z[-1], dsseg.z[0]]
            plt.plot(_xjoin, _zjoin, '-', color='r', linewidth=4, alpha=.5)
    else:
        plt.plot(lp.x_ext[0][-2:], lp.z_ext[0][-2:], '-', color='r', linewidth=4, alpha=.5)
    plt.plot(lp.x, lp.z, '-', color='r', linewidth=4, alpha=1.)#, label=lp.)

plt.show()

