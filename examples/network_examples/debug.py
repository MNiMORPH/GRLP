from grlp import *

# ---- River properties
x0 = 10.e3
L = 100.e3
mean_Q = 10.
mean_Qs = 0.001
B = 98.1202038813591

# Get random network topology
net_topo = Shreve_Random_Network(magnitude=10)

# Get setup parameters
nx_list, dx, Q_in, Qs_in = get_simple_network_setup_params(
    net_topo.upstream_segment_IDs,
    net_topo.downstream_segment_IDs,
    L,
    mean_Q,
    mean_Qs)
downstream_segment_list = net_topo.downstream_segment_IDs
upstream_segment_list = net_topo.upstream_segment_IDs

# ---- Some parameters for use during set up
segments = []
sources = [i for i in range(len(nx_list)) if not upstream_segment_list[i]]

# ---- Loop over segments filling lists for network
x_ls = []
z_ls = []
Q_ls = []
B_ls = []
for i,nx in enumerate(nx_list):
    
    # Set up x domain
    down_IDs = downstream_IDs(downstream_segment_list, i)[1:]
    down_nx = sum([nx_list[j] for j in down_IDs])
    x0 = - down_nx - nx_list[i]
    x1 = x0 + nx_list[i]
    # x = np.arange( (x0-1), (x1+1), 1. ) * dx
    x = np.arange( x0, x1, 1. ) * dx
    x_ls.append(x)
    
    # set width
    B_ls.append(B)
    
    # Set initial z
    S0 = (Qs_in/(0.041*Q_in))**(6./7.)
    z = (x.max()-x)*S0
    if downstream_segment_list[i]:
        # if not mouth, reset downstream elevation to that of downstream segment
        z += z_ls[downstream_segment_list[i][0]][0] + dx*S0
    z_ls.append(z)
    
    # Set discharge
    if i in sources:
        # if segment is a source, set Q input values
        Q_ls.append(np.full(len(x), Q_in))
    else:
        # otherwise set based on number of upstream sources
        up_IDs = upstream_IDs(upstream_segment_list, i)
        num_sources = len([j for j in up_IDs if not upstream_segment_list[j]])
        Q_ls.append(np.full(len(x), Q_in*num_sources))
        
# ---- Update x coordinates to run from 0 at furthest upstream point, record max
x_min = min([min(x) for x in x_ls])
x_max = max([max(x) for x in x_ls])
for i,nx in enumerate(nx_list):
    x_ls[i] -= x_min

net = Network()
net.initialize(
    config_file = None,
    x_bl =x_max,
    z_bl = 0.,
    S0 = None,
    Q_s_0 = Qs_in,
    upstream_segment_IDs = upstream_segment_list,
    downstream_segment_IDs = downstream_segment_list,
    x = x_ls,
    z = z_ls,
    Q = Q_ls,
    B = B_ls,
    overwrite = False
    )
net.set_niter(3)
net.get_z_lengths()
net.evolve_threshold_width_river_network(nt=1000, dt=3.15e10)