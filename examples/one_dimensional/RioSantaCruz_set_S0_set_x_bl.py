#! /usr/bin/env python3
"""
Río Santa Cruz: force the boundaries by *slope* (set_S0) and move base level
*along the continental shelf* (set_x_bl).

Two v3 boundary features on a real river (Patagonia, Argentina; parameters after
Ruby et al., 2026, AGU Advances, doi:10.1029/2025AV002035):

  * set_S0 -- force the UPSTREAM boundary with a slope measured from a DEM
    (here 180 m over 250 km = 7.2e-4), rather than back-deriving a sediment
    supply Q_s_0. set_S0 is the slope-first entry point; it applies the same
    boundary sediment supply that set_Qs_input_upstream would for that slope.

  * set_x_bl -- move the DOWNSTREAM base level along the continental-shelf
    gradient. When sea level falls by dz, the shoreline steps seaward by
    dz / shelf_slope, so the outlet's x AND z change together, coupled by the
    shelf slope. The Patagonian shelf is gentle (0.0009 m/m; Kennett, 1982), so
    a modest sea-level fall drives a large horizontal migration -- base level is
    controlled strongly by the shelf gradient. What the river feels is the
    downstream SLOPE: the mouth now grades at the shelf slope (0.0009), steeper
    than the river's own (7.2e-4), so a wave of incision propagates upstream
    (the regression phase in Ruby et al., 2026).

We force a sea-level fall to a lowstand and back (a regression): the mouth only
ever migrates seaward, which keeps the outlet downstream of the last node.
"""

import numpy as np
from matplotlib import pyplot as plt

import grlp

# ---- Río Santa Cruz parameters (field / DEM-derived) ----------------------
Q  = 700.          # m3/s, mean annual at the Lago Argentino outlet
z0 = 180.          # m, modern surface elevation at the head
z1 = 0.            # m, modern base level (sea level)
L  = 250e3         # m, length
B  = 10000.        # m, valley width (~uniform)
U  = 1e-4 / 3.15e7 # m/s, uplift from terrace ages
S0 = (z0 - z1) / L                  # 7.2e-4 -- the DEM-measured equilibrium slope
shelf_slope = 0.0009               # m/m, continental shelf gradient (Kennett, 1982)

lp = grlp.LongProfile()
lp.basic_constants()
lp.bedload_lumped_constants()
lp.set_hydrologic_constants()
lp.set_x(dx=1000, nx=251, x0=0)
lp.set_z(S0=-S0, z1=z1)
lp.set_Q(Q)
lp.set_B(B)
lp.set_niter(3)
lp.set_uplift_rate(U)
lp.set_z_bl(z1)

# ---- Upstream boundary by slope (set_S0) ----------------------------------
# The equilibrium slope is known from the DEM, so set the boundary from it
# directly. Confirm set_S0 matches the slope-to-supply route it replaces:
lp.set_S0(S0)
Qs0_from_S0 = lp.Q_s_0
Qs0_manual  = lp.k_Qs * lp.Q[0] * S0 ** (7 / 6.)
print("set_S0 boundary supply matches set_Qs_input_upstream: %s"
      % np.isclose(Qs0_from_S0, Qs0_manual))

# ---- Bring to steady state ------------------------------------------------
lp.evolve_threshold_width_river(10, 1e14)
z_eq = lp.z.copy()

# ---- Sea-level fall along the shelf (set_x_bl + set_z_bl) ------------------
amp     = 30.                    # m, sea-level amplitude (Ruby et al., 2026)
period  = 100e3 * 3.15e7         # s, 100 kyr
dt      = 2e3 * 3.15e7           # s, 2 kyr steps
nsteps  = 50                     # one 100-kyr regression cycle

x_bl0 = lp.x_ghost_downstream    # modern outlet position
for i in range(nsteps):
    t = (i + 1) * dt
    # sea level: fall to -amp at mid-cycle, recover to 0 (stays <= modern)
    dz = -amp / 2. * (1 - np.cos(2 * np.pi * t / period))
    lp.set_z_bl(z1 + dz)                        # base level drops
    lp.set_x_bl(x_bl0 - dz / shelf_slope)       # ...and steps seaward on the shelf
    lp.evolve_threshold_width_river(1, dt)

# ---- Profile response (the slope is the point; outlet position not plotted) -
plt.figure(figsize=(7, 3))
plt.plot(lp.x / 1000., z_eq, 'k-', lw=2, label='steady state (set_S0)')
plt.plot(lp.x / 1000., lp.z, 'b-', lw=2, label='after a shelf-coupled sea-level fall')
plt.xlabel('Downstream distance [km]', fontsize=12)
plt.ylabel('Elevation [m]', fontsize=12)
plt.legend()
plt.tight_layout()
plt.show()
