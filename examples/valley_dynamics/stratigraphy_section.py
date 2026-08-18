# stratigraphy_section.py
#
# Valley-dynamics demo (3 of 3): record the valley-fill stratigraphy.
#
# Everything is uniform in space and constant in time -- uniform discharge,
# valley width, and lateral migration rate, and a simple base-level history
# (steady rise, then held).  The one thing that changes is the aggradation rate:
# it is high while base level rises and decays once base level is held.  Because
# the channel-deposit fraction f_ch keys on how fast the channel reworks the
# valley relative to how fast it aggrades, f_ch sits in the mid-range (a
# channel/overbank mix) while base level rises and climbs toward 1 (channel-rich)
# as aggradation wanes.  So the valley fill records a vertical facies transition,
# uniform along the valley.
#
# The migration rate and aggradation rate are chosen so f_ch sits in a clearly
# visible mid-range rather than pinned near its confined floor (low migration) or
# saturated at 1 (high migration).
#
# The recorder (grlp.stratigraphy.StratRecorder, via set_stratigraphic_recording)
# stores a per-node (z, f_ch) polyline and compresses it on the fly with
# vertical-metric Douglas-Peucker, so the stored size is set by the f_ch
# tolerance, not the time step.

import numpy as np
from matplotlib import pyplot as plt

import grlp

YEAR = 3.15e7
# Values chosen (uniform in space, constant in time) so the f_ch dynamics are
# clearly visible: the migration rate sets where f_ch sits during aggradation --
# too low and it is pinned near its confined floor, too high and it saturates at
# 1.  ZETADOT = 1.5e-8 with this aggradation rate puts it in the visible mid-range.
RISE_RATE = 2.0e-3 / YEAR   # steady base-level rise [m/s]
ZETADOT = 1.5e-8           # uniform migration rate
GRAIN_D = 0.05
N_FORCE, N_RELAX, DT = 40, 90, 1.0e11
TOL = 0.02


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
lp.set_lateral_migration_rate(ZETADOT)     # uniform along the domain
lp.set_valley_dynamics(partition_by_aggradation=True)
lp.set_stratigraphic_recording(tol=TOL, record_thickness=True)

# Base level rises steadily, then is held.
for _ in range(N_FORCE):
    lp.set_z_bl(lp.z_bl + RISE_RATE * DT)
    lp.evolve_threshold_width_river(nt=1, dt=DT)
for _ in range(N_RELAX):
    lp.evolve_threshold_width_river(nt=1, dt=DT)

raw = sum(len(zc) for zc in lp.strat.z)
simp = sum(len(zc) for zc, _ in lp.strat.columns())
_, _, vol_err = lp.strat.volume_closure()
print('raw vertices %d -> compressed %d (%.1fx)' % (raw, simp, raw / simp))
print('volume closure max error: %.3f m' % vol_err)

fig, (ax_sec, ax_col) = plt.subplots(1, 2, figsize=(13, 6),
                                     gridspec_kw={'width_ratios': [2.4, 1]})
ax, pc = lp.strat.plot_section(lp.x / 1e3, ax=ax_sec)
fig.colorbar(pc, ax=ax_sec, label='channel-deposit fraction f_ch')
ax_sec.set_xlabel('Downstream distance [km]')
ax_sec.set_ylabel('Elevation [m]')
ax_sec.set_title('Uniform Q, valley width, and migration rate; simple base-level history\n'
                 'f_ch rises upward: mixed channel/overbank (aggrading) -> channel-rich (aggradation wanes)')

# one column: f_ch changing with elevation, and the kept Douglas-Peucker vertices
i0 = 45
zc_raw, fc_raw = np.asarray(lp.strat.z[i0]), np.asarray(lp.strat.f[i0])
zc, fc = lp.strat.columns()[i0]
ax_col.plot(fc_raw, zc_raw, '-', color='0.7', lw=1, label='raw (%d)' % len(zc_raw))
ax_col.plot(fc, zc, 'o-', color='C3', lw=1.5, ms=4, label='DP kept (%d)' % len(zc))
ax_col.set_xlabel('f_ch')
ax_col.set_ylabel('Elevation [m]')
ax_col.set_xlim(0, 1.05)
ax_col.set_title('Column at x = %.0f km' % (lp.x[i0] / 1e3))
ax_col.legend(fontsize=9)

fig.tight_layout()
fig.savefig('stratigraphy_section.png', dpi=150)
print('Wrote stratigraphy_section.png')
plt.show()
