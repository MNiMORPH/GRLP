import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
plt.ion()

import grlp
#reload(grlp)

S0 = 0.01
P_xB = 0.0
z1 = 0

Qamp = 0.5
dt = 3.15E7 * 1E3
nt = 50
Bmax = 250.

lp = grlp.LongProfile()
self = lp

self.bcr = z1

lp.basic_constants()
lp.bedload_lumped_constants()
lp.set_hydrologic_constants()

lp.set_x(dx=500, nx=180, x0=30E3)
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

"""
# Return to old values
lp.z = z0.copy()
lp.z_ext = z_ext_0.copy()
lp.set_z_bl(z1)
lp.set_Qs_input_upstream(0)
"""

lp.set_niter(3)

# Calculate
lambda_sternberg = 1E-4
U = 1E-9
lp.set_uplift_rate(U/3.15E7)
lp.set_Sternberg_gravel_loss(lambda_sternberg)

# Have to integrate this with Sternberg instead of leaving detached
uplift_Qs_divisor = (lp.A**.5 * lp.U)[0] / lp.Q_s_0
excess_Qs_in = lp.A**.5 * lp.U / uplift_Qs_divisor - lp.Q_s_0
excess_Qs_in_valley_equivalent = excess_Qs_in / ( (1-lp.lambda_p) * lp.B )
#lp.set_source_sink_distributed(excess_Qs_in_valley_equivalent/250000.)

# Actually, it should just be constant everywhere based on hillslope inputs
# And in fact, should relate to differential drainage area... but
# let's just say constant is a better approximation.
#lp.set_source_sink_distributed(lp.Q_s_0/5E4 / ( (1-lp.lambda_p) * lp.B) )

# More process based
#dA = np.diff(lp.A)
dA = (lp.A_ext[2:] - lp.A_ext[:-2])/2.
L_channel = dA**(4/7.) # Hack, pretending it works on an area-slice
# Production and decay
sed_input = 1E1*dA * U/3.15E7 * np.exp(-lambda_sternberg * L_channel)
#upstream_boundary_sed_input = 1E-1 * dA[0] * U/3.15E7 # approx hack need ext
#sed_input = np.hstack (( 0, sed_input ))
#lp.set_source_sink_distributed(lp.Q_s_0/lp.A** / ( (1-lp.lambda_p) * lp.B) )
lp.set_source_sink_distributed( sed_input / ( (1-lp.lambda_p) * lp.B) )



for i in range(nt):
    lp.evolve_threshold_width_river(5, dt)
    lp.slope_area(verbose=False)
    lp.compute_Q_s()
    Q_s_out_2.append(lp.Q_s[-1] * lp.intermittency)
    print( np.mean(lp.dz_dt/lp.U) )

plt.plot(lp.x, lp.z)
#plt.plot(lp.x, z02)


