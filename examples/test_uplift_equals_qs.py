import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
plt.ion()

import grlp
reload(grlp)

S0 = 0.01
P_xB = 0.0
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

Q_s_out_1 = []
Q_s_out_2 = []
Q_s_out_3 = []

time_transient_kyr = np.arange(30, 601, 30)

z0 = lp.z.copy()
z_ext_0 = lp.z_ext.copy()

lp.intermittency = 1.

# Return to old values
lp.z = z0.copy()
lp.z_ext = z_ext_0.copy()
lp.set_z_bl(z1)
lp.set_Qs_input_upstream(0)

lp.set_niter(3)

# Calculate
U = 1E-3
lp.set_uplift_rate(U/3.15E7)
for i in range(100):
    lp.evolve_threshold_width_river(5, 1E4*3.15E7)
    lp.slope_area(verbose=False)
    lp.compute_Q_s()
    Q_s_out_2.append(lp.Q_s[-1] * lp.intermittency)
    print np.mean(lp.dz_dt/lp.U)




