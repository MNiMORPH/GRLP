import numpy as np
from matplotlib import pyplot as plt
import importlib
plt.ion()
#plt.ioff()

import grlp
importlib.reload(grlp)

S0 = 0.015
P_xB = 0.2
z1 = 0+7.5 # Interesting -- my RHS always has to be 0 (not quite right)

dt = 3.15E7/2.

nseg = 3

segments = []
for i in range(nseg):
    segments.append(grlp.LongProfile())

Qlist = [10., 60., 70.]
upstream_segment_list = [[], [], [0,1]]
downstream_segment_list = [[2], [2], []]

i = 0
for lp in segments:
    lp.set_ID(i)
    lp.set_upstream_segment_IDs(upstream_segment_list[i])
    lp.set_downstream_segment_IDs(downstream_segment_list[i])
    lp.set_intermittency(1)
    lp.basic_constants()
    lp.bedload_lumped_constants()
    lp.set_hydrologic_constants()
    # Local or global x -- really just care about diffs for solver
    dx=500.
    nx=6
    x0=0
    _x = np.arange(x0, x0+dx*nx, dx)
    #_x[-3] += 1E-6
    x_ext = np.hstack((_x[0]-dx, _x, _x[-1]+dx))
    lp.set_x(x_ext=x_ext)
    lp.set_z(S0=-S0, z1=z1)
    lp.set_niter(10)
    #lp.set_z_bl(z1)
    lp.set_Q(Qlist[i])
    lp.set_B(100.)
    # START HERE
    if len(upstream_segment_list[i]) == 0:
        #Qs0 = 0.015 #lp.k_Qs * lp.Q[0] * (1*S0)**(7/6.)
        Qs0 = lp.k_Qs * lp.Q[0] * (1*S0)**(7/6.)
        lp.set_Qs_input_upstream(Qs0)
    i += 1
    lp.set_uplift_rate(0)
    

segments[0].z += segments[-1].z_ext[0] - 7.5 # BAND-AID
segments[0].z_ext += segments[-1].z_ext[0] - 7.5 # BAND-AID
segments[1].z += segments[-1].z_ext[0] - 7.5 # BAND-AID
segments[1].z_ext += segments[-1].z_ext[0] - 7.5 # BAND-AID

x0 = segments[0].x[-1]
segments[0].x += segments[-1].x_ext[0] - x0
segments[0].x_ext += segments[-1].x_ext[0] - x0
x0 = segments[1].x[-1]
segments[1].x += segments[-1].x_ext[0] - x0
segments[1].x_ext += segments[-1].x_ext[0] - x0

segments[-1].z_ext[-1] = -10

i = 0
for lp in segments:
    if len(downstream_segment_list[i]) == 0:
        lp.set_z_bl(lp.z_ext[-1])
    i += 1

net = grlp.Network(segments)
net.get_z_lengths()
net.set_niter(3)
net.build_ID_list()
#net.set_dQ()
segments[0].dQ[-1] = 0.#60/2.
segments[1].dQ[-1] = 0.#10/2.
segments[2].dQ[0] = 0.
net.evolve_threshold_width_river_network(nt=10, dt=dt)

plt.figure()
for lp in segments:
    plt.plot(lp.x_ext, lp.z_ext, '--', linewidth=4, alpha=.5)
    #plt.plot(lp.x, lp.z, '--', linewidth=4, alpha=.5, label=str(lp.Q[0]))
    plt.legend()

for lp in segments:
    lp.compute_Q_s()
    print(lp.Q_s)
    

net.list_of_LongProfile_objects[0].set_Qs_input_upstream(0.02)
net.evolve_threshold_width_river_network(nt=100, dt=dt)


plt.figure()
for lp in segments:
    plt.plot(lp.x_ext, lp.z_ext, '--', linewidth=4, alpha=.5)
    #plt.plot(lp.x, lp.z, '--', linewidth=4, alpha=.5, label=str(lp.Q[0]))
    plt.legend()
    
net.list_of_LongProfile_objects[-1].set_z_bl(-30)
net.list_of_LongProfile_objects[-1].set_bcr_Dirichlet()
net.evolve_threshold_width_river_network(nt=160, dt=dt/2.)

plt.figure()
for lp in segments:
    plt.plot(lp.x_ext, lp.z_ext, '--', linewidth=4, alpha=.5)
    #plt.plot(lp.x, lp.z, '--', linewidth=4, alpha=.5, label=str(lp.Q[0]))
    plt.legend()
    
