import numpy as np
from matplotlib import pyplot as plt
plt.ion()

import grlp
reload(grlp)

#S0 = 1E-2
S0 = 0.01#0.096316051787322637
P_xB = 0.2
z1 = 0

Qamp = 0.5
dt = 3.15E7 * 1E2
nt = int(100)
Bmax = 600.

lp = grlp.LongProfile()
self = lp

self.bcr = z1

lp.basic_constants()
lp.bedload_lumped_constants()
lp.set_hydrologic_constants()

#S0 = (lp.z[1] - lp.z[0])/lp.dx

lp.set_x(dx=1000, nx=75, x0=5E3)
lp.set_z(S0=-S0, z1=z1)
lp.set_A(k_xA=1.)
lp.set_Q(k_xQ=1.433776163432246e-05, P_xQ=7/4.*0.7)
#lp.set_B(k_xB=Bmax/np.max(lp.x**P_xB), P_xB=P_xB)
lp.set_B(k_xB=10, P_xB=0)
#lp.set_B(B=lp.Q*2)
#lp.set_B(B=10*np.ones(len(lp.Q)))
lp.set_uplift_rate(0.000/3.15E7)
lp.set_niter()
#lp.set_bcr_Dirichlet(z1)
lp.set_z_bl(z1)

Qs0 = lp.k_Qs * lp.Q[0] * S0**(7/6.)

x0 = lp.x.copy()
Q0 = lp.Q.copy()

lp.Q = Q0
lp.set_Qs_input_upstream(Qs0)
lp.evolve_threshold_width_river(200, 1E14)
z0 = lp.z.copy()

lp.analytical_threshold_width()
lp.compute_Q_s()


fig = plt.figure(figsize=(12,5))
ax1 = plt.subplot(111)
ax1.cla()
ax1.set_xlabel('Downstream distance [km]', fontsize=26)
ax1.set_ylabel('Elevation [m]', fontsize=26)
ax1.tick_params(axis='both', which='major', labelsize=16)
ax1.plot(lp.x/1000., lp.z, '0.7', linewidth=6)
ax1.plot(lp.x/1000., lp.zanalytical, '0.3', linewidth=2)
#ax1.set_ylim((0, 400))
#plt.savefig('LP_'+'%06d' %(t[i]/3.15E7)+'.png')
plt.tight_layout()
plt.show()

#print self.z[0] - self.z[-1]
