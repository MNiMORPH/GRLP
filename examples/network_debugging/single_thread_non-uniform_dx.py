from grlp import *

# ---- River properties
x0 = 10.e3
L = 100.e3
mean_Q = 10.
mean_Qs = 0.001
B = 98.1202038813591
z_bl = -10.

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
z = [(x[0].max()-x[0])*S0 + z_bl]
Q = [np.full(len(x),mean_Q)]
B = [np.full(len(x),B)]

# ---- Network object
net = Network()
net.initialize(
    config_file = None,
    x_bl =L+(dx/2.), # <---- last bit different dx
    z_bl = z_bl,
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
net.evolve_threshold_width_river_network(nt=1000, dt=3.15e11)
for seg in net.list_of_LongProfile_objects: seg.compute_Q_s()

plt.plot(net.list_of_LongProfile_objects[0].x/1.e3, net.list_of_LongProfile_objects[0].z)
plt.xlabel("Downstream distance [km]")
plt.ylabel("Elevation [m]")
plt.title("Long profile")
plt.show()

plt.plot(net.list_of_LongProfile_objects[0].x_ext[0][-3:]/1.e3, net.list_of_LongProfile_objects[0].z_ext[0][-3:])
plt.xlabel("Downstream distance [km]")
plt.ylabel("Elevation [m]")
plt.title("Long profile: last three nodes")
plt.show()

plt.plot(net.list_of_LongProfile_objects[0].x/1.e3, net.list_of_LongProfile_objects[0].Q_s)
plt.xlabel("Downstream distance [km]")
plt.ylabel(r"$Q_{s,out}$ [m$^3$ s$^{-1}$]")
plt.title("Sediment discharge")
plt.show()