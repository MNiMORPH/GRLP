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

lp.set_intermittency(1)

self.bcr = z1

lp.basic_constants()
lp.bedload_lumped_constants()
lp.set_hydrologic_constants()

#x_ext = np.array([0, 400, 800, 1200, 1600, 1700, 1900, 1950, 3000, 3500, 5000, 5200]).astype(float)

#x_ext = np.linspace(0,500*9,10)
#x_ext = np.hstack((x_ext, np.max(x_ext)+500+1E-6))

dx=500.
nx=180
x0=10E3
x = np.arange(x0, x0+dx*nx, dx)
dx = dx
dx_isscalar = True
#x = np.arange(x0, x0+dx*nx, 50 + 900 * np.random.rand())
x_ext = np.hstack((x[0]-dx, x, x[-1]+dx))
#x_ext[-3] += 5000
#x_ext[-2] += 15000
#x_ext[-1] += 100
x_ext[-1] += 1E-6
#x_ext += np.linspace(-400, 400, len(x_ext))

lp.set_x(x_ext=x_ext)
lp.set_z(S0=-S0, z1=z1)
lp.set_A(k_xA=1.)
lp.set_Q(k_xQ=1.433776163432246e-05, P_xQ=7/4.*0.7)
lp.set_B(k_xB=Bmax/np.max(lp.x**P_xB), P_xB=P_xB)
lp.set_niter()
lp.set_z_bl(z1)
Qs0 = lp.k_Qs * lp.Q[0] * S0**(7/6.)
lp.set_Qs_input_upstream(Qs0)

lp.set_uplift_rate(0)
lp.evolve_threshold_width_river(1, 1E22)
lp.analytical_threshold_width()

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
