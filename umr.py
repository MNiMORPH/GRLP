import numpy as np
from matplotlib import pyplot as plt

import grlp
reload(grlp)

S0 = 1E-3
P_xB = 0.8

lp = grlp.LongProfile()
self = lp

lp.basic_constants()
lp.bedload_lumped_constants()
lp.set_hydrologic_constants()

lp.set_x(dx=1000, nx=1100, x0=3*110*1000)
lp.set_z(S0=-S0, z1=100)
lp.set_A(k_xA=1.)
lp.set_Q(q_R=0.01, A_R=1E7)
lp.set_B(k_xB=10./np.max(lp.x**P_xB), P_xB=P_xB)
#lp.set_B(k_xB=100., P_xB=0.)
lp.set_uplift_rate(0)
lp.set_niter()
Qs0 = lp.k_Qs * lp.Q[0] * (100*S0)**(7/6.)
lp.set_Qs_input_upstream(Qs0)
lp.set_bcr_Dirichlet(0)
#lp.set_uplift_rate(0.01/3.15E7)
lp.evolve_threshold_width_river(150, 1E9)
lp.analytical_threshold_width(P_xB=P_xB)
lp.compute_Q_s()


x0 = lp.x.copy()
z0 = lp.z.copy()

plt.ion()
plt.figure(figsize=(12,6))
plt.plot(lp.x/1000., lp.z, '0.6', linewidth=6)
plt.plot(lp.x/1000., lp.zanalytical, 'k', linewidth=2)
plt.xlabel('Downstream distance [km]', fontsize=26)
plt.ylabel('Elevation [m]', fontsize=26)
plt.tick_params(axis='both', which='major', labelsize=16)
plt.tight_layout()

# Glacial inputs
lp.Q += 1E3
lp.set_Qs_input_upstream(Qs0*2.)
for i in range(10):
    lp.evolve_threshold_width_river(1, 1E9)
    plt.plot(lp.x/1000., lp.z, '0.3', linewidth=2)
    plt.pause(.1)
# Now only excess water and a bit of sediment
#lp.Q += 5E3
for i in range(10):
    lp.set_Qs_input_upstream(Qs0*(2-(i+1)/10.))
    lp.evolve_threshold_width_river(1, 1E9)
    plt.plot(lp.x/1000., lp.z, '0.3', linewidth=2)
    plt.pause(.1)
for i in range(10):
    lp.evolve_threshold_width_river(1, 1E9)
    plt.plot(lp.x/1000., lp.z, '0.3', linewidth=2)
    plt.pause(.1)
# Now cut off the water
self.Q -= 1E3
lp.set_Qs_input_upstream(Qs0*1)
for i in range(10):
    lp.evolve_threshold_width_river(1, 1E9)
    plt.plot(lp.x/1000., lp.z, '0.3', linewidth=2)
    plt.pause(.1)

