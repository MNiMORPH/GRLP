# stratigraphy_section.py
#
# Valley-dynamics demo (3 of 3): record the valley-fill stratigraphy, and make
# the overbank-deposition feedback visible as cyclic facies.
#
# With uniform discharge the steady profile is a straight line and there is no
# downstream gradient, so the overbank feedback has to be revealed in *time*.
# The channel-deposit fraction f_ch keys on the aggradation rate (fast rise ->
# gravel confined to the channel belt, low f_ch, overbank-rich; slow rise ->
# channel reworks the whole valley, f_ch -> 1, channel-rich).  So a cyclic base
# level -- alternating fast and slow rise -- lays down alternating overbank-rich
# and channel-rich bands: the feedback recorded as stratigraphy, uniform along
# the valley but banded in elevation.
#
# The recorder (grlp.stratigraphy.StratRecorder, via set_stratigraphic_recording)
# stores a per-node (z, f_ch) polyline and compresses it on the fly with
# vertical-metric Douglas-Peucker, so the stored size is set by the f_ch
# tolerance, not the time step.
#
# NB the accurate transient f_ch needs the in-Picard valley coupling (the default
# once valley dynamics are active) or a finer dt; the between-step coupling's
# one-step f_ch lag inflates transient features.

import numpy as np
from matplotlib import pyplot as plt

import grlp

YEAR = 3.15e7
R0 = 1.0e-3 / YEAR         # reference base-level rise rate [m/s]
FAST, SLOW = 2.0 * R0, 0.15 * R0  # rise rate in the fast / slow half-cycles
N_FAST, N_SLOW = 15, 55    # steps in each half-cycle (slow longer -> visible band)
N_CYCLE = 4                # number of fast/slow cycles
ZETADOT = 1.0e-9
GRAIN_D = 0.05
DT = 1.0e11
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
lp.set_lateral_migration_rate(ZETADOT)
lp.set_valley_dynamics(partition_by_aggradation=True)
lp.set_stratigraphic_recording(tol=TOL, record_thickness=True)

# Cyclic base level: alternating fast and slow rise.
for _ in range(N_CYCLE):
    for _ in range(N_FAST):
        lp.set_z_bl(lp.z_bl + FAST * DT)
        lp.evolve_threshold_width_river(nt=1, dt=DT)
    for _ in range(N_SLOW):
        lp.set_z_bl(lp.z_bl + SLOW * DT)
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
ax_sec.set_title('Cyclic base-level rise, constant Q: overbank feedback as banded facies\n'
                 '(dark = overbank-rich / fast; light = channel-rich / slow; %.1fx compression)'
                 % (raw / simp))

# one column: f_ch cycling with elevation, and the kept Douglas-Peucker vertices
i0 = 45
zc_raw, fc_raw = np.asarray(lp.strat.z[i0]), np.asarray(lp.strat.f[i0])
zc, fc = lp.strat.columns()[i0]
ax_col.plot(fc_raw, zc_raw, '-', color='0.7', lw=1, label='raw (%d)' % len(zc_raw))
ax_col.plot(fc, zc, 'o-', color='C3', lw=1.5, ms=3, label='DP kept (%d)' % len(zc))
ax_col.set_xlabel('f_ch')
ax_col.set_ylabel('Elevation [m]')
ax_col.set_title('Column at x = %.0f km' % (lp.x[i0] / 1e3))
ax_col.legend(fontsize=9)

fig.tight_layout()
fig.savefig('stratigraphy_section.png', dpi=150)
print('Wrote stratigraphy_section.png')
plt.show()
