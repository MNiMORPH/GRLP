import numpy as np
from matplotlib import pyplot as plt

import grlp
reload(grlp)

S0 = 1E-2
P_xB = 0.8

lp = grlp.LongProfile()
self = lp

lp.basic_constants()
lp.bedload_lumped_constants()
lp.set_hydrologic_constants()

lp.set_x(dx=1000, nx=50, x0=2000)
lp.set_z(S0=-S0)
lp.set_A(k_xA=1.)
lp.set_Q(q_R=0.01, A_R=1E5)
lp.set_B(k_xB=10./np.max(lp.x**P_xB), P_xB=P_xB)
#lp.set_B(k_xB=100., P_xB=0.)
lp.set_uplift_rate(0)
lp.set_niter()
Qs0 = lp.k_Qs * lp.Q[0] * (10*S0)**(7/6.)
lp.set_Qs_input_upstream(Qs0)
lp.set_bcr_Dirichlet(0)
#lp.set_uplift_rate(0.01/3.15E7)
lp.evolve_threshold_width_river(100, 1E8)
lp.compute_Q_s()

plt.ion()
plt.plot(self.x/1000., self.z)
