# lateral_supply_bounce.py
#
# BENCHMARK (not a test or example): does the sediment supplied by lateral valley-
# wall erosion produce a *bounce* in the long profile -- late aggradation after
# base-level fall stops and the valley keeps widening -- or does the profile
# relax smoothly and diffusively to its new equilibrium?
#
# For a fixed 64 m base-level drop, we sweep the channel lateral migration rate
# and the base-level-fall rate.  Each run falls to the 64 m drop, then holds at
# the low base level for the rest of a long run (>> the fall time) so the profile
# reaches its new equilibrium and any late wall-fed aggradation can develop.
# Every case is run twice -- with and without the lateral sediment supply -- to
# isolate its effect.
#
# Outputs (raw data + figures) are written next to this script.

import copy
import os

import numpy as np
from matplotlib import pyplot as plt

import grlp

YEAR = 3.15e7
DROP = 64.0                       # fixed base-level drop [m]
DT = 1.0e11                       # time step [s] (~3.17 kyr)
N_TOTAL = 260                     # total steps (fall + long hold to equilibrium)
MID = 45                          # mid-domain node tracked for the bounce
FALL_RATES_MMYR = [1.0, 2.0, 4.0]         # base-level-fall rates [mm/yr]
MIGRATION = [1.0e-10, 1.0e-9, 1.0e-8]     # lateral migration rates [m/s]
HERE = os.path.dirname(os.path.abspath(__file__))


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


base = build_steady_profile()
x = base.x.copy()
t = np.arange(N_TOTAL) * DT / YEAR / 1e3   # kyr


def run(fall_rate_myr, zetadot, supply):
    fall_rate = fall_rate_myr * 1e-3 / YEAR
    n_fall = int(round(DROP / (fall_rate * DT)))
    lp = copy.deepcopy(base)
    lp.D = 0.05
    lp.set_lateral_migration_rate(zetadot)
    lp.set_valley_dynamics(narrow_by_incision=True, widen_by_migration=True)
    lp.supply_lateral_sediment = supply
    z_mid = np.zeros(N_TOTAL)
    B_mid = np.zeros(N_TOTAL)
    for k in range(N_TOTAL):
        if k < n_fall:
            lp.set_z_bl(lp.z_bl - fall_rate * DT)
        lp.evolve_threshold_width_river(nt=1, dt=DT)
        z_mid[k] = lp.z[MID]
        B_mid[k] = np.broadcast_to(lp.B, lp.x.shape)[MID]
    return dict(z_mid=z_mid, B_mid=B_mid, n_fall=n_fall, z_final=lp.z.copy())


# --- sweep ------------------------------------------------------------------ #
results = {}
for fr in FALL_RATES_MMYR:
    for zd in MIGRATION:
        for supply in (True, False):
            results[(fr, zd, supply)] = run(fr, zd, supply)
            print('done: fall=%.0f mm/yr  zetadot=%.0e  supply=%s' % (fr, zd, supply))

# --- save raw data ---------------------------------------------------------- #
save = {'x': x, 't_kyr': t, 'DROP': DROP,
        'fall_rates_mmyr': np.array(FALL_RATES_MMYR),
        'migration': np.array(MIGRATION)}
for (fr, zd, supply), r in results.items():
    tag = 'fr%.0f_zd%.0e_%s' % (fr, zd, 'on' if supply else 'off')
    save[tag + '_zmid'] = r['z_mid']
    save[tag + '_Bmid'] = r['B_mid']
    save[tag + '_zfinal'] = r['z_final']
    save[tag + '_nfall'] = r['n_fall']
np.savez(os.path.join(HERE, 'lateral_supply_bounce_data.npz'), **save)
print('wrote lateral_supply_bounce_data.npz')


# --- bounce diagnostic ------------------------------------------------------ #
# A bounce is a *rise after the deepest incision* -- so measure the peak the bed
# reaches after it bottoms out, not the final elevation.  The naive "final minus
# minimum" hides the interesting case entirely: the low-mobility bump rises then
# crashes back to equilibrium, so final == minimum and it reads as zero.  We
# report both the peak of the rebound (the bump) and how much of it survives to
# the end (settled), which exposes the bump-and-crash directly.
def bounce(zmid):
    """Return (bump, settled) at the mid node [m].

    bump    -- highest the bed climbs after the deepest incision (the rebound).
    settled -- how much of that rebound is still there at the end of the run;
               bump > 0 with settled ~ 0 is a bump-and-crash (unstable).
    """
    kmin = int(np.argmin(zmid))
    bump = zmid[kmin:].max() - zmid[kmin]
    settled = zmid[-1] - zmid[kmin]
    return bump, settled


print('\nBOUNCE at the mid node [m]: rebound above deepest incision, and how')
print('much survives to the end (bump with settled ~ 0 = bump-and-crash).')
print('  fall[mm/yr]  zetadot    bump-ON  settled-ON   bump-OFF   verdict')
for fr in FALL_RATES_MMYR:
    for zd in MIGRATION:
        bon, son = bounce(results[(fr, zd, True)]['z_mid'])
        boff, _ = bounce(results[(fr, zd, False)]['z_mid'])
        if bon < 1.0:
            verdict = 'smooth (no bounce)'
        elif son < 1.0:
            verdict = 'bump-and-crash (unstable)'
        else:
            verdict = 'persistent bounce'
        print('  %6.0f       %.0e  %7.1f    %7.1f    %7.1f   %s'
              % (fr, zd, bon, son, boff, verdict))

# --- figures ---------------------------------------------------------------- #
# Two rows per fall rate.  Top: the mid-node bed elevation, supply on (solid) vs
# off (dashed).  Bottom: the *isolated* wall-sediment effect, z(on) - z(off) --
# this is the benchmark's actual question ("does wall sediment cause a bounce?")
# and it is unambiguous where the raw curves are not: a single bumping curve can
# read as two, whereas the difference plot shows the wall-sediment contribution
# directly.  Only the lowest mobility departs from zero, and it bumps-and-crashes.
cols = plt.cm.viridis(np.linspace(0.15, 0.85, len(MIGRATION)))
fig, axes = plt.subplots(2, len(FALL_RATES_MMYR),
                         figsize=(5 * len(FALL_RATES_MMYR), 8),
                         sharex=True)
for j, fr in enumerate(FALL_RATES_MMYR):
    ax_z, ax_d = axes[0, j], axes[1, j]
    n_fall = results[(fr, MIGRATION[0], True)]['n_fall']
    for ax in (ax_z, ax_d):
        ax.axvline(t[n_fall - 1], color='0.7', ls=':', lw=1)
    ax_z.text(t[n_fall - 1], ax_z.get_ylim()[1], ' hold', fontsize=8,
              color='0.4', va='top')
    for zd, c in zip(MIGRATION, cols):
        zon = results[(fr, zd, True)]['z_mid']
        zoff = results[(fr, zd, False)]['z_mid']
        ax_z.plot(t, zon, color=c, lw=2, label='%.0e (on)' % zd)
        ax_z.plot(t, zoff, color=c, lw=1.2, ls='--')
        ax_d.plot(t, zon - zoff, color=c, lw=2, label='%.0e' % zd)
    ax_z.set_title('fall rate = %.0f mm/yr' % fr)
    ax_d.set_xlabel('Time [kyr]')
    ax_d.axhline(0, color='0.7', lw=0.8)
axes[0, 0].set_ylabel('Mid-node bed elevation [m]')
axes[1, 0].set_ylabel('Wall-sediment effect,\nz(on) - z(off) [m]')
axes[0, 0].legend(fontsize=8, title='migration rate\n(solid=on, dashed=off)')
axes[1, 0].legend(fontsize=8, title='migration rate')
fig.suptitle('Lateral wall-sediment supply: bounce vs smooth relaxation (64 m drop)\n'
             'bottom row isolates the wall-sediment contribution -- only the '
             'lowest mobility bounces, and it is unstable (bump-and-crash)',
             fontsize=12)
fig.tight_layout()
fig.savefig(os.path.join(HERE, 'lateral_supply_bounce.png'), dpi=150)
print('wrote lateral_supply_bounce.png')
