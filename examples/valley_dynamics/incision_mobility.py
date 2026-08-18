# incision_mobility.py
#
# Base-level fall (incision) with a fixed valley width vs an adaptive one, across
# a range of channel lateral mobilities.
#
# When the bed incises, the adaptive valley narrows toward the channel width b;
# lateral migration (the channel mobility) widens it back toward its prescribed
# width B_max.  So the valley width during and after an incision episode is a
# race between narrowing and widening, set by the migration rate:
#
#   * low mobility  -> the valley narrows sharply and stays narrow, and
#   * high mobility -> widening keeps pace, so the valley barely narrows.
#
# The left panel shows the incision itself (the long profile falling to its new
# equilibrium); the right panel shows the valley width at the mid-point through
# time for the fixed case and for the adaptive case at each mobility.

import numpy as np
from matplotlib import pyplot as plt

import grlp

YEAR = 3.15e7
FALL_RATE = 2.0e-3 / YEAR
GRAIN_D = 0.05
N_FALL, N_HOLD, DT = 40, 120, 1.0e11
ZETADOTS = [3.0e-10, 1.0e-9, 3.0e-9, 1.0e-8, 3.0e-8]   # range of channel mobilities
MID = 45
# NB the mobility trend below (low -> narrow, high -> wide) is robust, but the
# specific valley widths are dt-sensitive: the abrupt B -> b narrowing competes
# with the rate-limited widening, so coarse dt under-narrows, and a valley driven
# very close to the channel width becomes numerically stiff near equilibrium.
# We stay in the stable regime here; use a finer dt for quantitative widths.


def build_steady_profile():
    S0 = 1.5e-2
    lp = grlp.LongProfile()
    lp.set_intermittency(0.8)
    lp.basic_constants(); lp.bedload_lumped_constants(); lp.set_hydrologic_constants()
    lp.set_x(dx=1000.0, nx=90, x0=10000.0)
    lp.set_z(S0=-S0)
    lp.set_A(k_xA=1.0)
    lp.set_Q(Q=10.0)
    lp.set_B(B=200.0)
    lp.set_uplift_rate(0.0); lp.set_niter(3)
    lp.set_Qs_input_upstream(lp.k_Qs * lp.Q[0] * S0 ** (7 / 6.0))
    lp.set_z_bl(0.0)
    lp.evolve_threshold_width_river(nt=60, dt=1.0e13)
    return lp


def run_case(zetadot=None, n_snap=0):
    """Base-level fall then hold.  zetadot=None -> fixed valley width.  Returns
    the mid-node valley-width series, and (if n_snap) evenly spaced z snapshots."""
    lp = build_steady_profile()
    lp.D = GRAIN_D
    if zetadot is not None:
        lp.set_lateral_migration_rate(zetadot)
        lp.set_valley_dynamics(narrow_by_incision=True, widen_by_migration=True)
    n_total = N_FALL + N_HOLD
    snap_every = max(1, n_total // n_snap) if n_snap else 0
    B_mid = [np.broadcast_to(lp.B, lp.x.shape)[MID]]
    snaps = [lp.z.copy()] if n_snap else []
    for k in range(n_total):
        if k < N_FALL:
            lp.set_z_bl(lp.z_bl - FALL_RATE * DT)
        lp.evolve_threshold_width_river(nt=1, dt=DT)
        B_mid.append(np.broadcast_to(lp.B, lp.x.shape)[MID])
        if n_snap and (k + 1) % snap_every == 0:
            snaps.append(lp.z.copy())
    return np.array(B_mid), snaps


t = np.arange(N_FALL + N_HOLD + 1) * DT / YEAR / 1e3   # kyr
fixed_B, fixed_snaps = run_case(zetadot=None, n_snap=7)

fig, (ax_z, ax_B) = plt.subplots(1, 2, figsize=(13, 6))

# Left: the incision (fixed-width run), long profile falling to equilibrium.
x = build_steady_profile().x / 1e3
shades = plt.cm.viridis(np.linspace(0.15, 0.9, len(fixed_snaps)))
for zs, c in zip(fixed_snaps, shades):
    ax_z.plot(x, zs, color=c, lw=1.8)
ax_z.set_xlabel('Downstream distance [km]')
ax_z.set_ylabel('Bed elevation [m]')
ax_z.set_title('Base-level fall: incision to a new equilibrium\n(colour: light = early, dark = late)')

# Right: valley width at the mid-point through time -- fixed vs adaptive, swept
# over migration rate.
ax_B.plot(t, fixed_B, 'k--', lw=2, label='fixed width')
cols = plt.cm.plasma(np.linspace(0.1, 0.85, len(ZETADOTS)))
for zd, c in zip(ZETADOTS, cols):
    B_mid, _ = run_case(zetadot=zd)
    ax_B.plot(t, B_mid, color=c, lw=2, label='adaptive, %.0e m/s' % zd)
ax_B.axvline(N_FALL * DT / YEAR / 1e3, color='0.7', ls=':', lw=1)
ax_B.text(N_FALL * DT / YEAR / 1e3, 205, ' base level held', fontsize=8, color='0.4')
ax_B.set_xlabel('Time [kyr]')
ax_B.set_ylabel('Valley width at mid-point [m]')
ax_B.set_title('Valley width: fixed vs adaptive, across channel mobilities\n'
               '(narrow while incising, re-widen by migration when it stops)')
ax_B.legend(fontsize=9, loc='center right', title='migration rate')

fig.tight_layout()
fig.savefig('incision_mobility.png', dpi=150)
print('Wrote incision_mobility.png')
plt.show()
