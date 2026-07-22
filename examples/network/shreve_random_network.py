#! /usr/bin/env python3
"""
Generate a synthetic river network and run GRLP on it.

A Shreve (1974) random network is generated, given discharges and widths, and
evolved toward steady state with the GRLP threshold-width solver. This is the
quickest way to set up and run a GRLP network -- no DEM required -- so it is a
good first thing to try after installing GRLP.
"""

import random

import numpy as np
from matplotlib import pyplot as plt

import grlp

# Seed for a reproducible network; remove these two lines to draw a new random
# network each run.
random.seed(7)
np.random.seed(7)

# ---- Generate a Shreve random network -------------------------------------
# magnitude    : number of channel heads (network "size")
# max_length   : length of the longest source-to-outlet path [m]
# mean_discharge: sets segment discharges from drainage area
net, topo = grlp.generate_random_network(
    magnitude=8,
    max_length=2.0e4,
    mean_discharge=10.,
)
net.set_niter(3)
net.get_z_lengths()

# ---- Evolve toward steady state -------------------------------------------
net.evolve_threshold_width_river_network(nt=200, dt=3.15e11)

# ---- Long profiles --------------------------------------------------------
plt.figure()
for lp in net.list_of_LongProfile_objects:
    # connect each segment's end to its downstream neighbour (or the base-level
    # ghost at the outlet)
    if lp.downstream_segment_IDs:
        ds = net.list_of_LongProfile_objects[lp.downstream_segment_IDs[0]]
        plt.plot([lp.x[-1], ds.x[0]], [lp.z[-1], ds.z[0]], 'k-', linewidth=1)
    else:
        plt.plot([lp.x[-1], lp.x_ghost_downstream], [lp.z[-1], lp.z_bl],
                 'k-', linewidth=1)
    plt.plot(lp.x, lp.z, '-')
plt.xlabel('Downstream distance [m]')
plt.ylabel('Elevation [m]')
plt.title('Long profiles at steady state')

# ---- Network planform -----------------------------------------------------
plt.figure()
net.plot()   # draws the branching planform and shows both figures
