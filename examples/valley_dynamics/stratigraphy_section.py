# stratigraphy_section.py
#
# Valley-dynamics demo (3 of 3): record the valley-fill stratigraphy, and expose
# the overbank-deposition feedback through an internal river property rather than
# through the forcing.
#
# The forcing is held simple -- a steady base-level rise -- and the discharge and
# valley width are uniform.  The one thing that varies along the river is the
# channel's lateral mobility (migration rate zetadot), increasing downstream.
# Because the channel-deposit fraction f_ch keys on how fast the channel reworks
# the valley relative to how fast it aggrades, the slow-migrating upstream reach
# preserves overbank-rich deposits (low f_ch) while the fast-migrating downstream
# reach reworks the whole valley into channel-rich deposits (f_ch -> 1).  So a
# single internal gradient, under constant forcing, prints a downstream facies
# gradient into the valley fill.
#
# The recorder (grlp.stratigraphy.StratRecorder, via set_stratigraphic_recording)
# stores a per-node (z, f_ch) polyline and compresses it on the fly with
# vertical-metric Douglas-Peucker, so the stored size is set by the f_ch
# tolerance, not the time step.

import numpy as np
from matplotlib import pyplot as plt

import grlp

YEAR = 3.15e7
RISE_RATE = 1.0e-3 / YEAR   # steady base-level rise [m/s] -- constant forcing
GRAIN_D = 0.05
NT, DT = 120, 1.0e11
TOL = 0.02
# Lateral migration rate: slow upstream -> fast downstream (the internal gradient).
ZETADOT_UP, ZETADOT_DOWN = 5.0e-11, 1.0e-7


def build_steady_profile():
    S0 = 1.5e-2
    lp = grlp.LongProfile()
    lp.set_intermittency(0.8)
    lp.basic_constants(); lp.bedload_lumped_constants(); lp.set_hydrologic_constants()
    lp.set_x(dx=1000.0, nx=90, x0=10000.0)
    lp.set_z(S0=-S0)
    lp.set_A(k_xA=1.0)
    lp.set_Q(Q=10.0)                       # uniform discharge along the domain
    lp.set_B(B=200.0)
    lp.set_uplift_rate(0.0); lp.set_niter(3)
    lp.set_Qs_input_upstream(lp.k_Qs * lp.Q[0] * S0 ** (7 / 6.0))
    lp.set_z_bl(0.0)
    lp.evolve_threshold_width_river(nt=60, dt=1.0e13)
    return lp


lp = build_steady_profile()
lp.D = GRAIN_D
zetadot = np.logspace(np.log10(ZETADOT_UP), np.log10(ZETADOT_DOWN), len(lp.x))
lp.set_lateral_migration_rate(zetadot)     # per-node (spatially variable) mobility
lp.set_valley_dynamics(partition_by_aggradation=True)
lp.set_stratigraphic_recording(tol=TOL, record_thickness=True)

# Constant forcing: steady base-level rise throughout.
for _ in range(NT):
    lp.set_z_bl(lp.z_bl + RISE_RATE * DT)
    lp.evolve_threshold_width_river(nt=1, dt=DT)

raw = sum(len(zc) for zc in lp.strat.z)
simp = sum(len(zc) for zc, _ in lp.strat.columns())
print('raw vertices %d -> compressed %d (%.1fx)' % (raw, simp, raw / simp))

fig, (ax_sec, ax_map) = plt.subplots(1, 2, figsize=(13, 6),
                                     gridspec_kw={'width_ratios': [2.4, 1]})
ax, pc = lp.strat.plot_section(lp.x / 1e3, ax=ax_sec)
fig.colorbar(pc, ax=ax_sec, label='channel-deposit fraction f_ch')
ax_sec.set_xlabel('Downstream distance [km]')
ax_sec.set_ylabel('Elevation [m]')
ax_sec.set_title('Constant forcing, uniform Q: a downstream migration-rate gradient\n'
                 'prints a facies gradient (overbank-rich -> channel-rich)')

# The internal control: migration rate and the f_ch it produces, vs distance.
f_ch = np.broadcast_to(np.asarray(lp.f_ch, dtype=float), lp.x.shape)
ax_map.plot(lp.x / 1e3, f_ch, 'C0-', lw=2.5, label='f_ch (facies)')
ax_map.set_xlabel('Downstream distance [km]')
ax_map.set_ylabel('channel-deposit fraction f_ch', color='C0')
ax_map.tick_params(axis='y', labelcolor='C0')
ax_map.set_ylim(0, 1.05)
ax2 = ax_map.twinx()
ax2.semilogy(lp.x / 1e3, zetadot, 'C3--', lw=2, label='migration rate')
ax2.set_ylabel('migration rate zetadot [m/s]', color='C3')
ax2.tick_params(axis='y', labelcolor='C3')
ax_map.set_title('The internal control')

fig.tight_layout()
fig.savefig('stratigraphy_section.png', dpi=150)
print('Wrote stratigraphy_section.png')
plt.show()
