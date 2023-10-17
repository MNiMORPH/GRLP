from grlp import *

# ---- River properties
x0 = 10.e3
L = 100.e3
mean_Q = 10.
mean_Qs = 0.001
B = 98.1202038813591

# ---- Basic lp object to get k_Qs
lp = LongProfile()
lp.basic_constants()
lp.bedload_lumped_constants()
lp.set_hydrologic_constants()

# ---- Arrays
dx = 1.e3
x = [np.arange(0., L+dx, dx)]
S0 = [(mean_Qs/(lp.k_Qs*mean_Q))**(6./7.)]
upstream_segment_IDs = [[]]
downstream_segment_IDs = [[]]
z = [(x[0].max()-x[0])*S0]
Q = [np.full(len(x),mean_Q)]
B = [np.full(len(x),B)]

# ---- Network object
net = Network()
net.initialize(
    config_file = None,
    x_bl =L+dx,
    z_bl = 0.,
    S0 = S0,
    upstream_segment_IDs = upstream_segment_IDs,
    downstream_segment_IDs = downstream_segment_IDs,
    x = x,
    z = z,
    Q = Q,
    B = B,
    overwrite = False
    )
net.set_niter(3)
net.get_z_lengths()
net.evolve_threshold_width_river_network(nt=100, dt=3.15e11)

plt.plot(net.list_of_LongProfile_objects[0].x, net.list_of_LongProfile_objects[0].z)
plt.show()

for seg in net.list_of_LongProfile_objects: seg.compute_Q_s()