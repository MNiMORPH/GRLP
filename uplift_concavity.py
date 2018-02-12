import numpy as np
from matplotlib import pyplot as plt
plt.ion()

import grlp
reload(grlp)

#S0 = 1E-2
S0 = 0.01#0.096316051787322637
P_xB = 0.1
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
lp.set_B(k_xB=Bmax/np.max(lp.x**P_xB), P_xB=P_xB)
#lp.set_B(k_xB=10, P_xB=0)
lp.set_niter()
lp.set_z_bl(z1)
Qs0 = lp.k_Qs * lp.Q[0] * S0**(7/6.)
lp.set_Qs_input_upstream(Qs0)

plt.figure()
Uall = np.linspace(-0.0001,0.0001,51)
thetaall = []
R2all = []
for U in Uall:
    print U
    #plt.cla()
    lp.set_uplift_rate(U/3.15E7)
    lp.evolve_threshold_width_river(10, 1E14)
    lp.set_uplift_rate(0)
    lp.evolve_threshold_width_river(10, 1E14)
    self.slope_area()
    thetaall.append(self.theta)
    R2all.append(self.thetaR2)
    plt.plot(lp.x, lp.z - lp.z[-1], 'k-')
    plt.pause(0.01)

plt.figure()
plt.plot(Uall, thetaall, 'k-')
plt.show()

#print self.z[0] - self.z[-1]

