#! /usr/bin/python3

# Figure UpliftSubsidence in paper

import numpy as np
from matplotlib import pyplot as plt
#plt.ion()

import grlp

Q = 700. # m3/s <-- mean annual @ Lago Argentino outlet, not bankfull
#Q*= 5

z0 = 180 # From modern surface near Lago Argentino outlet
z1 = 0
L = 250E3
S0 = (z0-z1)/L
B = 10000 # m: valley width, measured, approx uniform
U = 1E-4 # m/yr: from terrace age and height above modern floodplain

lp = grlp.LongProfile()
self = lp

self.bcr = z1

lp.basic_constants()
lp.bedload_lumped_constants()
lp.set_hydrologic_constants()

lp.set_x(dx=1000, nx=251, x0=0)
lp.set_z(S0=-S0, z1=z1)
lp.set_Q(Q)
lp.set_B(B)
lp.set_niter()
lp.set_z_bl(z1)

# Initial sed supply is the slope that we prescribe
Qs0 = lp.k_Qs * lp.Q[0] * S0**(7/6.)
lp.set_Qs_input_upstream(Qs0)

fig = plt.figure(figsize=(6,3))
ax1 = fig.add_subplot(1,1,1)
plt.xlabel('Downstream distance [km]', fontsize=14, fontweight='bold')
plt.ylabel('Elevation [m]', fontsize=14, fontweight='bold')
plt.tight_layout()

# Evolve for 1 year
lp.set_uplift_rate(U/3.15E7)
lp.evolve_threshold_width_river(1, 1*3.15E7)

# Let's try just comparing starting points
#"""
# Starting steady state
lp.set_uplift_rate(U/3.15E7)
lp.evolve_threshold_width_river(10, 1E14)
ax1.plot(lp.x/1000., lp.z, color='k', alpha=1, linewidth=3)

# Transient -- double sediment supply
ts = 1000 # years
nts = 50
lp.set_Qs_input_upstream(Qs0*4)
for i in range(nts):
    lp.evolve_threshold_width_river(1, ts * 3.15E7)
    #ax1.plot(lp.x/1000., lp.z, color='r', alpha=1-0.75/nts*i, linewidth=1)
    ax1.plot(lp.x/1000., lp.z, color='b', alpha=0.1, linewidth=1)
ax1.plot(lp.x/1000., lp.z, color='b', alpha=1, linewidth=3)

# Transient -- sediment supply to 0
lp.set_Qs_input_upstream(0)
for i in range(nts):
    lp.evolve_threshold_width_river(1, ts * 3.15E7)
    #ax1.plot(lp.x/1000., lp.z, color='r', alpha=1-0.75/nts*i, linewidth=1)
    ax1.plot(lp.x/1000., lp.z, color='r', alpha=0.1, linewidth=1)
ax1.plot(lp.x/1000., lp.z, color='r', alpha=1, linewidth=3)

plt.show()
#"""

