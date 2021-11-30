# run_grlp.py
# ADW, commented 29 November 2018

####################################
# Quick example to set up GRLP     #
# See Wickert and Schildgen, 2018  # 
# This example reproduces Figure 2 #
# in Wickert and Schildgen, 2019   # 
####################################

# Import numerical and plotting libraries
import numpy as np
from matplotlib import pyplot as plt
import importlib

# Import the GRLP module
# "reload(grlp)" is to refresh grlp if you are running this interactively
# while modifying the grlp library
import grlp
importlib.reload(grlp)

# Instantiate the long profile object
lp = grlp.LongProfile()
# Uncomment this so you can refer to lp as "self" for debugging in the 
# grlp module (e.g., while modifying it)
# self = lp

# Uncomment this if you want to enable interactive plotting
# plt.ion()

# S0 is the upstream-end slope that determines the sediment input to the 
# catchment
S0 = 1.5E-2
# Valley width: B = k_xB * x**P_xB (k_xB defined below)
P_xB = 0.2

# Intermittency: What fraction of the total time is the river experiencing a
# geomorphically-effective flood? This assumes a binary on--off state, common 
# for gravel-bed rivers with floodplains (see Blom et al., 2017)
lp.set_intermittency(0.8)

# Utility functions to create constants defined in the W&S 2018 paper
lp.basic_constants()
lp.bedload_lumped_constants()
lp.set_hydrologic_constants()

# Set up the x domain
lp.set_x(dx=1000, nx=90, x0=10000)

# Set up a starting set of channel-bed elevations (z) on a uniform slope (S0)
lp.set_z(S0=-S0)

# Set up transfer functions between drainage area (A), discharge (Q), and
# valley width (B). Use help(lp.set_A), etc., to see the other options 
# available for each of these functions.
lp.set_A(k_xA=1.)
lp.set_Q(k_xQ = 1.43e-5, P_xQ = 49/40)
# k_xB = 250./np.max(lp.x**P_xB)
lp.set_B(k_xB=25, P_xB=P_xB)

# Set the uplift rate [m/s]; positive upwards
lp.set_uplift_rate(0)

# Set up the number of iterations in semi-implicit solver; defaults to 3.
lp.set_niter()

# Input sediment discharge: this is set based on your defined S0, above.
# (this ficticious boundary-condition slope is the transport slope for the
#  amount of sediment being supplied)
Qs0 = lp.k_Qs * lp.Q[0] * (S0)**(7/6.)
lp.set_Qs_input_upstream(Qs0)

# Set the base level (the set_z function already assumes that z_bl starts at 0)
lp.set_z_bl(0.)

# Numerical and analytical solutions.
# Numerical: (number of time steps, length of time step [s])
lp.evolve_threshold_width_river(5, 1E12)
lp.analytical_threshold_width()
lp.compute_Q_s()

# Plot
plt.figure(figsize=(12,6))
plt.plot(lp.x/1000., lp.z, '0.6', linewidth=6, label='Numerical')
plt.plot(lp.x/1000., lp.zanalytical, 'k', linewidth=2, label='Analytical')
plt.xlabel('Downstream distance [km]', fontsize=26)
plt.ylabel('Elevation [m]', fontsize=26)
plt.tick_params(axis='both', which='major', labelsize=16)
plt.legend()
plt.tight_layout()
plt.show()
