import numpy as np

import grlp
reload(grlp)

S0 = 1E-2

lp = grlp.LongProfile()
self = lp

lp.basic_constants()
lp.bedload_lumped_constants()
lp.set_hydrologic_constants()

lp.set_x(dx=1000, nx=20, x0=2000)
lp.set_z(S0=-S0)
lp.set_A(k_xA=1.)
lp.set_Q(q_R=0.001, A_R=1E6)
#lp.set_B(k_xB=10./np.max(lp.x**0.8), P_xB=0.8)
lp.set_B(k_xB=100., P_xB=0.)
lp.set_uplift_rate(0)
lp.set_niter()
lp.set_Qs_input_upstream(lp.k_Qs * lp.Q[0] * 10*S0**(7/6.))
lp.set_bcr_Dirichlet()
lp.evolve_threshold_width_river(500, 2E9)
lp.compute_Q_s()

plt.ion()
plt.plot(self.x/1000., self.z)
