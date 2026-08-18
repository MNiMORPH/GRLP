# stratigraphy_section.py
#
# Valley-dynamics demo (3 of 3): record the valley-fill stratigraphy.
#
# During the aggradation-with-overbank-deposition experiment, record the
# channel-deposit fraction f_ch as a function of elevation and downstream
# distance -- a synthetic valley-fill section.  The recorder (grlp.stratigraphy
# .StratRecorder, attached via set_stratigraphic_recording) stores a per-node
# (z, f_ch) polyline and compresses it on the fly with vertical-metric
# Douglas-Peucker line simplification, so the stored size is set by the f_ch
# tolerance, not the time step.
#
# The section reads as an upward transition from overbank-dominated (low f_ch,
# during forcing) to channel-dominated (f_ch -> 1, as base level holds and
# aggradation ceases).
#
# NB For quantitative transient stratigraphy, use the in-Picard valley coupling
# (set_valley_coupling('in_picard')) or a finer time step: with the default
# between-step coupling the one-step f_ch lag, amplified by the overbank
# feedback, inflates transient features (a spurious sub-steady f_ch minimum) at
# coarse dt.  The gross section is robust either way.

import numpy as np
from matplotlib import pyplot as plt

import grlp

YEAR = 3.15e7
RISE_RATE = 1.0e-3 / YEAR
ZETADOT = 1.0e-9
GRAIN_D = 0.05
N_FORCE, N_RELAX, NT_SNAP, DT = 6, 6, 10, 1.0e11
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

# Phase 1: base-level rise (aggradation).  Phase 2: base level held (relaxation).
for _ in range(N_FORCE):
    for _ in range(NT_SNAP):
        lp.set_z_bl(lp.z_bl + RISE_RATE * DT)
        lp.evolve_threshold_width_river(nt=1, dt=DT)
for _ in range(N_RELAX):
    lp.evolve_threshold_width_river(nt=NT_SNAP, dt=DT)

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
ax_sec.set_title('Valley-fill section f_ch(x, z)\n(%.1fx compression at tol = %.2f)'
                 % (raw / simp, TOL))

# one column, showing the kept Douglas-Peucker vertices
i0 = 20
zc_raw, fc_raw = np.asarray(lp.strat.z[i0]), np.asarray(lp.strat.f[i0])
zc, fc = lp.strat.columns()[i0]
ax_col.plot(fc_raw, zc_raw, '.-', color='0.7', lw=1, ms=3, label='raw (%d)' % len(zc_raw))
ax_col.plot(fc, zc, 'o-', color='C3', lw=2, ms=5, label='DP kept (%d)' % len(zc))
ax_col.set_xlabel('f_ch')
ax_col.set_ylabel('Elevation [m]')
ax_col.set_title('Column at x = %.0f km' % (lp.x[i0] / 1e3))
ax_col.legend(fontsize=9)

fig.tight_layout()
fig.savefig('stratigraphy_section.png', dpi=150)
print('Wrote stratigraphy_section.png')
plt.show()
