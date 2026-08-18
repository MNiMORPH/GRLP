# aggradation_deposition.py
#
# Valley-dynamics demo (1 of 2): aggradation with and without suspended-load
# (overbank) deposition.
#
# A gravel long profile is spun up to steady state and then forced to aggrade
# by a rising base level (following Fergus McNab's width.py: a base-level ramp
# is the natural aggradation driver -- the upstream sediment supply is fixed and
# backfills from the downstream boundary).  We compare two rivers:
#
#   * WITHOUT overbank deposition -- the channel gravel fills the whole valley
#     width B (the channel-deposit fraction f_ch = 1), and
#   * WITH overbank deposition    -- the laterally migrating channel reworks
#     only part of the valley, so its gravel is confined to a channel belt of
#     width f_ch * B < B while the rest of the valley receives suspended-load
#     fines.  f_ch is set by the Wickert et al. (2013) deposit partition, closed
#     with the Turowski crossing time, and evolves with the aggradation rate.
#
# With the sediment supply fixed by the upstream boundary, confining the gravel
# to a narrower belt makes the same supply fill a smaller width, so the river
# WITH overbank deposition aggrades its channel bed *faster*.  This mirrors the
# incision demo, where narrowing concentrates the gravel and deepens incision:
# in both, confining the gravel speeds the bed's response to base level.  The
# lower panel shows f_ch, the mechanism behind the difference.
#
# NB the lateral migration rate here is finite but small: it sets f_ch through
# the deposit partition without widening the valley (no channel-belt coefficient
# k0 is set, so the widening rule is inactive).  A migration rate of exactly
# zero would collapse f_ch to b/B ~ 0.01, an extreme end-member -- not a useful
# illustration.

import numpy as np
from matplotlib import pyplot as plt

import grlp

YEAR = 3.15e7


def build_steady_profile():
    """A canonical gravel long profile evolved to steady state (no forcing).

    Mirrors examples/run_grlp.py (a reproduction of Wickert & Schildgen, 2019,
    Fig. 2): a 90-node domain with power-law discharge and valley width.  The
    incision demo uses the identical spin-up, so the two start from the same
    river.
    """
    S0 = 1.5e-2
    lp = grlp.LongProfile()
    lp.set_intermittency(0.8)
    lp.basic_constants()
    lp.bedload_lumped_constants()
    lp.set_hydrologic_constants()
    lp.set_x(dx=1000.0, nx=90, x0=10000.0)
    lp.set_z(S0=-S0)
    lp.set_A(k_xA=1.0)
    lp.set_Q(k_xQ=1.43e-5, P_xQ=49 / 40.0)
    lp.set_B(k_xB=25.0, P_xB=0.2)          # this B is now the *valley* width
    lp.set_uplift_rate(0.0)
    lp.set_niter(3)
    lp.set_Qs_input_upstream(lp.k_Qs * lp.Q[0] * S0 ** (7 / 6.0))
    lp.set_z_bl(0.0)
    lp.evolve_threshold_width_river(nt=60, dt=1.0e13)
    return lp


# Forcing and sampling (mirror of the incision demo: same rate and cadence,
# opposite sign of base-level change).
RISE_RATE = 1.0e-3 / YEAR   # base-level rise rate [m/s] -> aggradation
ZETADOT = 1.0e-9           # lateral migration rate [m/s] -> sets f_ch, no widening
GRAIN_D = 0.05             # grain size [m], needed to resolve channel width & depth
N_SNAP = 6                 # number of long-profile snapshots
NT_SNAP = 10               # solver steps between snapshots
DT = 1.0e11                # time step [s]


def run(with_deposition):
    """Evolve an aggrading profile, returning (x, list-of-z-snapshots, f_ch)."""
    lp = build_steady_profile()
    lp.D = GRAIN_D
    if with_deposition:
        lp.set_lateral_migration_rate(ZETADOT)
        lp.set_valley_dynamics(partition_by_aggradation=True)
    snapshots = [lp.z.copy()]
    # Step one at a time so base level rises smoothly rather than in jumps.
    for _ in range(N_SNAP):
        for _ in range(NT_SNAP):
            lp.set_z_bl(lp.z_bl + RISE_RATE * DT)
            lp.evolve_threshold_width_river(nt=1, dt=DT)
        snapshots.append(lp.z.copy())
    # f_ch defaults to the scalar 1.; broadcast to the profile for plotting.
    f_ch = np.broadcast_to(np.asarray(lp.f_ch, dtype=float), lp.x.shape).copy()
    return lp.x, snapshots, f_ch


x, z_off, fch_off = run(with_deposition=False)
_, z_on, fch_on = run(with_deposition=True)

times_kyr = np.arange(N_SNAP + 1) * NT_SNAP * DT / YEAR / 1.0e3

fig, (ax_prof, ax_fch) = plt.subplots(
    2, 1, figsize=(9, 8), height_ratios=[2, 1], sharex=True)

# Long profiles: colour ramps light (early) -> dark (late).
shades = plt.cm.viridis(np.linspace(0.15, 0.9, N_SNAP + 1))
for zs, sh in zip(z_off, shades):
    ax_prof.plot(x / 1e3, zs, color=sh, lw=1.5, ls='--')
for zs, sh in zip(z_on, shades):
    ax_prof.plot(x / 1e3, zs, color=sh, lw=2.5)
ax_prof.plot([], [], 'k--', lw=1.5, label='no overbank deposition (f_ch = 1)')
ax_prof.plot([], [], 'k-', lw=2.5, label='with overbank deposition (f_ch < 1)')
ax_prof.set_ylabel('Bed elevation [m]', fontsize=13)
ax_prof.set_title(
    'Aggradation: overbank deposition confines the gravel and speeds bed rise\n'
    '(colour: light = early, dark = late; %.0f kyr total)' % times_kyr[-1],
    fontsize=12)
ax_prof.legend(fontsize=11, loc='upper right')

# f_ch: the mechanism.  f_ch = 1 without deposition; < 1 with it.
ax_fch.plot(x / 1e3, fch_off, 'k--', lw=1.5,
            label='no overbank deposition')
ax_fch.plot(x / 1e3, fch_on, 'k-', lw=2.5,
            label='with overbank deposition')
ax_fch.set_ylim(0, 1.05)
ax_fch.set_xlabel('Downstream distance [km]', fontsize=13)
ax_fch.set_ylabel('Channel-deposit\nfraction f_ch [-]', fontsize=13)
ax_fch.legend(fontsize=11, loc='lower right')

fig.tight_layout()
fig.savefig('aggradation_deposition.png', dpi=150)
print('Wrote aggradation_deposition.png')
print('bed rise at x0 -- no deposition: %.1f m ; with deposition: %.1f m'
      % (z_off[-1][0] - z_off[0][0], z_on[-1][0] - z_on[0][0]))
plt.show()
