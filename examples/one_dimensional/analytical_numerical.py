import numpy as np
from matplotlib import pyplot as plt
import importlib
import sys
#plt.ion()

import grlp

# Reload -- for testing purposes
if sys.version_info.major == 3:
    importlib.reload(grlp)
else:
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

lp.set_intermittency(1)

self.bcr = z1

lp.basic_constants()
lp.bedload_lumped_constants()
lp.set_hydrologic_constants()

lp.set_x(dx=500., nx=180, x0=10E3)
lp.set_z(S0=-S0, z1=z1)
lp.set_A(k_xA=1.)
lp.set_Q(k_xQ=1.433776163432246e-05, P_xQ=7/4.*0.7)
lp.set_B(k_xB=Bmax/np.max(lp.x**P_xB), P_xB=P_xB)
lp.set_niter(3)
lp.set_z_bl(z1)
Qs0 = lp.k_Qs * lp.Q[0] * S0**(7/6.)
lp.set_Qs_input_upstream(Qs0)

lp.set_uplift_rate(0)
lp.evolve_threshold_width_river(1, 1E15)
lp.analytical_threshold_width()

fig = plt.figure(figsize=(5,3))
ax1 = fig.add_subplot(1,1,1)
plt.xlabel('Downstream distance [km]', fontsize=14, fontweight='bold')
plt.ylabel('Elevation [m]', fontsize=14, fontweight='bold')

lp.slope_area()
ax1.plot(lp.x/1000., lp.z, '-', color='.6', linewidth=6, label='Numerical')
ax1.plot(lp.x/1000., lp.zanalytical, '-', color='0', linewidth=2, label='Analytical')
plt.legend(loc='upper right')
plt.tight_layout()
plt.show()

#plt.savefig('Uplift_Subsidence.svg')
