# Figure UpliftSubsidence in paper

import numpy as np
from matplotlib import pyplot as plt
plt.ion()

import grlp
reload(grlp)

S0 = 0.015
P_xB = 0.2
z1 = 0+7.5 # Interesting -- my RHS always has to be 0 (not quite right)

dt = 3.15E22

nseg = 1

segments = []
for i in range(nseg):
    segments.append(grlp.LongProfile())

# This should be built from network in the future
#Qlist = [10., 15., 25., 20., 45.]
#upstream_segment_list = [[], [], [0,1], [], [2,3]]
#downstream_segment_list = [[2], [2], [4], [4], []]

Qlist = [10., 10.]
upstream_segment_list = [[], [0]]
downstream_segment_list = [[1], []]

Qlist = [10.]
upstream_segment_list = [[]]
downstream_segment_list = [[]]


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
    _x[-3] += 1E-6
    x_ext = np.hstack((_x[0]-dx, _x, _x[-1]+dx))
    lp.set_x(x_ext=x_ext)
    lp.set_z(S0=-S0, z1=z1)
    lp.set_niter()
    #lp.set_z_bl(z1)
    lp.set_Q(Qlist[i])
    lp.set_B(100.)
    # START HERE
    if len(upstream_segment_list[i]) == 0:
        Qs0 = lp.k_Qs * lp.Q[0] * S0**(7/6.)
        lp.set_Qs_input_upstream(Qs0)
    i += 1
    lp.set_uplift_rate(0)
    


# z adjustment: should also be done with network code
# "zext" includes first/last cells of neighboring streams
# UPDATE WHOLE SET OF Z, NOT THIS
#segments[3].set_z_bl(segments[4].z[0])
#segments[2].set_z_bl(segments[4].z[0])
#segments[1].set_z_bl(segments[2].z[0])
#segments[0].set_z_bl(segments[2].z[0])
"""
segments[3].z += segments[4].z[0]
segments[3].z_ext += segments[4].z[0]
segments[2].z += segments[4].z[0]
segments[2].z_ext += segments[4].z[0]
segments[1].z += segments[2].z[0]
segments[1].z_ext += segments[2].z[0]
segments[0].z += segments[2].z[0]
segments[0].z_ext += segments[2].z[0]

#segments[4].x += segments[4].x_ext[0] - segments[4].x[-1]
#segments[4].x_ext += segments[4].x[0] - segments[4].x[-1]
x0 = segments[3].x_ext[-1]
segments[3].x += segments[4].x_ext[0] - x0 + dx
segments[3].x_ext += segments[4].x[0] - x0 + dx
x0 = segments[2].x_ext[-1]
segments[2].x += segments[4].x_ext[0] - x0 + dx
segments[2].x_ext += segments[4].x[0] - x0 + dx
x0 = segments[1].x_ext[-1]
segments[1].x += segments[2].x_ext[0] - x0 + 2*dx
segments[1].x_ext += segments[2].x[0] - x0 + 2*dx
x0 = segments[0].x_ext[-1]
segments[0].x += segments[2].x_ext[0] - x0 + 2*dx
segments[0].x_ext += segments[2].x[0] - x0 + 2*dx
"""

"""
segments[0].z += segments[1].z_ext[0]
segments[0].z_ext += segments[1].z_ext[0]

x0 = segments[0].x[-1]
segments[0].x += segments[1].x_ext[0] - x0
segments[0].x_ext += segments[1].x_ext[0] - x0
"""

i = 0
for lp in segments:
    if len(downstream_segment_list[i]) == 0:
        lp.set_z_bl(lp.z_ext[-1])
    i += 1

i = 0
for lp in segments:
    lp.build_LHS_coeff_C0(dt=dt)
    lp.build_matrices()
    i += 1

net = grlp.Network(segments)
#net.build_block_diagonal_matrix_core()
#net.add_block_diagonal_matrix_upstream_boundary_conditions()
#net.add_block_diagonal_matrix_downstream_boundary_conditions()
net.get_z_lengths()
net.set_niter()
net.evolve_threshold_width_river_network(nt=1, dt=dt)
#self = net



plt.figure()
for lp in segments:
    plt.plot(lp.x_ext, lp.z_ext, '--', linewidth=4, alpha=.5)



#lp.evolve_threshold_width_river(1, 1E22)
#lp.analytical_threshold_width()

"""
fig = plt.figure(figsize=(5,3))
ax1 = fig.add_subplot(1,1,1)
plt.xlabel('Downstream distance [km]', fontsize=14, fontweight='bold')
plt.ylabel('Elevation [m]', fontsize=14, fontweight='bold')
  
#lp.slope_area()
ax1.plot(lp.x/1000., lp.z, '-', color='.6', linewidth=6, label='Numerical')
ax1.plot(lp.x/1000., lp.zanalytical, '-', color='0', linewidth=2, label='Analytical')
plt.legend(loc='upper right')
plt.show()
plt.tight_layout()

#plt.savefig('Uplift_Subsidence.svg')
"""
