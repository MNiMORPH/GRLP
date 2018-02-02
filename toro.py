import numpy as np
from matplotlib import pyplot as plt

import grlp
reload(grlp)

S0 = 1E-2
P_xB = 0.2
#z1 = 1000
z1 = 0

Qamp = 0.5
dt = 1E9
Qperiod = 20000.*3.15E7
nt = 2 * int(np.ceil(Qperiod/dt))
Bmax = 600.

lp = grlp.LongProfile()
self = lp

self.bcr = z1

lp.basic_constants()
lp.bedload_lumped_constants()
lp.set_hydrologic_constants()

#lp.set_uplift_rate(0.1)

lp.set_x(dx=1000, nx=60, x0=4E4)
lp.set_z(S0=-S0, z1=z1)
lp.set_A(k_xA=1.)
lp.set_Q(q_R=0.002, A_R=1E4)
lp.set_B(k_xB=Bmax/np.max(lp.x**P_xB), P_xB=P_xB)
#lp.set_B(k_xB=100., P_xB=0.)
lp.set_uplift_rate(0)
lp.set_niter()
lp.set_bcr_Dirichlet(z1)
#lp.set_uplift_rate(0.01/3.15E7)

Qs0 = lp.k_Qs * lp.Q[0] * (10*S0)**(7/6.)

x0 = lp.x.copy()
Q0 = lp.Q.copy()

lp.Q = Q0*.5
lp.set_Qs_input_upstream(Qs0)
lp.evolve_threshold_width_river(150, 1E11)
z_min_eq = lp.z.copy()

lp.Q = Q0 * 1.5
lp.set_Qs_input_upstream(Qs0)
lp.evolve_threshold_width_river(150, 1E11)
z_max_eq = lp.z.copy()

lp.Q = Q0
lp.set_Qs_input_upstream(Qs0)
lp.evolve_threshold_width_river(150, 1E11)
z0 = lp.z.copy()

lp.analytical_threshold_width(P_xB=P_xB)
lp.compute_Q_s()


plt.ion()
plt.figure(figsize=(12,6))
plt.plot(lp.x/1000., lp.z, '0.7', linewidth=6)
#plt.plot(lp.x/1000., lp.zanalytical, 'k', linewidth=2)
plt.plot(lp.x/1000., z_max_eq, '--', color='.5')
plt.plot(lp.x/1000., z_min_eq, '--', color='.5')
plt.xlabel('Downstream distance [km]', fontsize=26)
plt.ylabel('Elevation [m]', fontsize=26)
plt.tick_params(axis='both', which='major', labelsize=16)
plt.ylim((0, 400))
plt.tight_layout()

ts = np.arange(nt)
t = dt * ts
Qmult = np.sin(t / Qperiod * 2 * np.pi) * Qamp + 1

zmax = z0.copy()
zmin = z0.copy()

zall = []

for i in range(nt):
    lp.Q = Q0 * Qmult[i]
    lp.set_Qs_input_upstream(Qs0)
    lp.compute_Q_s() # slope too
    lp.evolve_threshold_width_river(1, dt)
    zmax = np.maximum(zmax, lp.z)
    zmin = np.minimum(zmin, lp.z)
    zall.append(lp.z.copy())
    if i%5 == 0:
        plt.cla()
        plt.xlabel('Downstream distance [km]', fontsize=26)
        plt.ylabel('Elevation [m]', fontsize=26)
        plt.tick_params(axis='both', which='major', labelsize=16)
        plt.plot(lp.x/1000., z0, '0.7', linewidth=6)
        plt.plot(lp.x/1000., lp.z, '0.3', linewidth=2)
        plt.ylim((0, 400))
        plt.pause(0.01)
plt.plot(lp.x/1000., zmax, '--', color='.0')
plt.plot(lp.x/1000., zmin, '--', color='.0')
plt.plot(lp.x/1000., z_max_eq, '--', color='.5')
plt.plot(lp.x/1000., z_min_eq, '--', color='.5')

# long profile with envelope of total change



"""
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
"""
