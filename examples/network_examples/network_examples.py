from grlp import *

# --- Some generic properties
Q_max = 16.         # discharge at outlet
Qs_max = 0.001163   # sediment flux at outlet
dx = 1.e3           # x increment


# --- Build a simple network: main trunk with single-segment tributaries

# First instantiate class to generate lists describing topology
# Specify:
#   - total length of trunk
#   - length of trunk segments
#   - length of tributary segments [default same as trunk segments]
net_topo = Simple_Network(nx_total=101, nx_trunk_seg=20, nx_trib_seg=10)

# Next call function to set up grlp's Network object
# Specify:
#   - dx
#   - topology from net_topo
#   - Q & Qs
#   - optionally evolve for some time aiming for steady state
net = set_up_network_object(
    nx_list=net_topo.nxs,
    dx=dx, 
    upstream_segment_list=net_topo.upstream_segment_IDs, 
    downstream_segment_list=net_topo.downstream_segment_IDs,
    Q_max=Q_max, 
    Qs_max=Qs_max,
    evolve=False)

# Plot simple visualisation
# Doesn't work every time!
# OK for simple networks
__ = plot_network(net)


# ---- Build a topologically random network: algorithm from Shreve (1974, NRR).

# First instantiate class to generate lists describing topology
# Topology and link (segment) lengths assigned random based on specified limits
# Specify:
#	- network magnitude (i.e. number of sources)
#	- minimum link (segment) length
#	- maximum link (segment) length
net_topo = Shreve_Random_Network(
	magnitude=20,
	min_link_length=4,
	max_link_length=8)

# Proceed as before
net = set_up_network_object(
    nx_list=net_topo.nxs,
    dx=dx, 
    upstream_segment_list=net_topo.upstream_segment_IDs, 
    downstream_segment_list=net_topo.downstream_segment_IDs,
    Q_max=Q_max, 
    Qs_max=Qs_max,
    evolve=False)
__ = plot_network(net)