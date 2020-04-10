# Figure UpliftSubsidence in paper

import numpy as np
from matplotlib import pyplot as plt
plt.ion()

import grlp
reload(grlp)

S0 = 0.015
P_xB = 0.2
z1 = 0

Qamp = 0.5
dt = 3.15E7 * 1E2
nt = int(100)
Bmax = 250.

lp = grlp.LongProfile()
self = lp

self.bcr = z1

lp.basic_constants()
lp.bedload_lumped_constants()
lp.set_hydrologic_constants()

lp.set_x(dx=500, nx=180, x0=10E3)
lp.set_z(S0=-S0, z1=z1)
lp.set_A(k_xA=1.)
lp.set_Q(k_xQ=1.433776163432246e-05, P_xQ=7/4.*0.7)
lp.set_B(k_xB=Bmax/np.max(lp.x**P_xB), P_xB=P_xB)
lp.set_niter()
lp.set_z_bl(z1)
Qs0 = lp.k_Qs * lp.Q[0] * S0**(7/6.)
lp.set_Qs_input_upstream(Qs0)

fig = plt.figure(figsize=(6,3))
ax1 = fig.add_subplot(1,1,1)
plt.xlabel('Downstream distance [km]', fontsize=14, fontweight='bold')
plt.ylabel('Elevation [m]', fontsize=14, fontweight='bold')
plt.tight_layout()

# Starting case
U = 0.
lp.set_uplift_rate(U/3.15E7)
lp.evolve_threshold_width_river(10, 1E14)
ax1.plot(lp.x/1000., lp.z - lp.z[0] + 500, color='.5', linewidth=3)

# Transient
U = 1E-3
lp.set_uplift_rate(U/3.15E7)
for i in range(10):
    lp.evolve_threshold_width_river(1, 3E11)
    ax1.plot(lp.x/1000., lp.z - lp.z[0] + 500, color='.5', linewidth=1)

# New equilibrium
lp.evolve_threshold_width_river(1, 1E14)
ax1.plot(lp.x/1000., lp.z - lp.z[0] + 500, color='0', linewidth=3)

