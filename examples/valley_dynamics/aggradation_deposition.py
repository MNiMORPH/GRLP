# aggradation_deposition.py
#
# Valley-dynamics demo (1 of 2): aggradation with and without suspended-load
# (overbank) deposition.
#
# A gravel long profile is spun up to steady state and then run in two phases:
#
#   Phase 1 (forcing) -- base level rises at a constant rate, driving
#     aggradation (following Fergus McNab's width.py: a base-level ramp is the
#     natural aggradation driver -- the upstream sediment supply is fixed and
#     backfills from the downstream boundary).
#   Phase 2 (relaxation) -- base level is held fixed and the profile relaxes the
#     rest of the way toward its new equilibrium.
#
# In each phase we compare two rivers:
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
# WITH overbank deposition aggrades its channel bed *faster* during forcing.
# Once base level stops (Phase 2), aggradation decays toward equilibrium and the
# deposit partition switches off (f_ch -> 1), so the two rivers converge.  The
# lower panel shows f_ch at the end of forcing, when the feedback is strongest.
#
# NB the lateral migration rate here is finite but small: it sets f_ch through
# the deposit partition without widening the valley (widen_by_migration is left
# off, so the widening rule is inactive).  A migration rate of exactly zero would
# collapse f_ch to b/B ~ 0.01, an extreme end-member -- not a useful illustration.
#
# For simplicity in these tests the initial valley width is uniform along the
# whole domain.

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
    lp.set_B(B=200.0)                      # uniform valley width B (= B_max)
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
N_FORCE = 6                # long-profile snapshots during forcing
N_RELAX = 6                # long-profile snapshots during relaxation
NT_SNAP = 10               # solver steps between snapshots
DT = 1.0e11                # time step [s]


def run(with_deposition):
    """Evolve through forcing then relaxation.

    Returns (x, forcing-snapshots, relaxation-snapshots, f_ch-at-forcing-end).
    The forcing list includes the initial steady profile.
    """
    lp = build_steady_profile()
    lp.D = GRAIN_D
    if with_deposition:
        lp.set_lateral_migration_rate(ZETADOT)
        lp.set_valley_dynamics(partition_by_aggradation=True)

    def f_ch_now():
        # f_ch defaults to the scalar 1.; broadcast to the profile for plotting.
        return np.broadcast_to(np.asarray(lp.f_ch, dtype=float), lp.x.shape).copy()

    # Phase 1: base level rises, one step at a time so it climbs smoothly.
    forcing = [lp.z.copy()]
    for _ in range(N_FORCE):
        for _ in range(NT_SNAP):
            lp.set_z_bl(lp.z_bl + RISE_RATE * DT)
            lp.evolve_threshold_width_river(nt=1, dt=DT)
        forcing.append(lp.z.copy())
    f_ch_force = f_ch_now()          # deposit fraction at the end of forcing

    # Phase 2: base level held fixed; the profile relaxes toward equilibrium.
    relaxation = []
    for _ in range(N_RELAX):
        lp.evolve_threshold_width_river(nt=NT_SNAP, dt=DT)
        relaxation.append(lp.z.copy())
    f_ch_final = f_ch_now()          # deposit fraction at the end of the experiment

    return lp.x, forcing, relaxation, f_ch_force, f_ch_final


x, force_off, relax_off, fch_off_force, fch_off_final = run(with_deposition=False)
_, force_on, relax_on, fch_on_force, fch_on_final = run(with_deposition=True)

t_force_kyr = N_FORCE * NT_SNAP * DT / YEAR / 1.0e3
t_total_kyr = (N_FORCE + N_RELAX) * NT_SNAP * DT / YEAR / 1.0e3

fig, (ax_prof, ax_fch) = plt.subplots(
    2, 1, figsize=(9, 8), height_ratios=[2, 1], sharex=True)

# Long profiles.  Forcing phase in a cool ramp, relaxation phase in a warm ramp;
# within each, light -> dark is early -> late.  Solid = with deposition, dashed
# = without.
cool = plt.cm.viridis(np.linspace(0.15, 0.9, N_FORCE + 1))
warm = plt.cm.autumn_r(np.linspace(0.35, 0.95, N_RELAX))
for zs, c in zip(force_off, cool):
    ax_prof.plot(x / 1e3, zs, color=c, lw=1.5, ls='--')
for zs, c in zip(force_on, cool):
    ax_prof.plot(x / 1e3, zs, color=c, lw=2.5)
for zs, c in zip(relax_off, warm):
    ax_prof.plot(x / 1e3, zs, color=c, lw=1.5, ls='--')
for zs, c in zip(relax_on, warm):
    ax_prof.plot(x / 1e3, zs, color=c, lw=2.5)
ax_prof.plot([], [], 'k-', lw=2.5, label='with overbank deposition (f_ch < 1)')
ax_prof.plot([], [], 'k--', lw=1.5, label='no overbank deposition (f_ch = 1)')
ax_prof.plot([], [], color=cool[-1], lw=6, label='base level rising (0-%.0f kyr)' % t_force_kyr)
ax_prof.plot([], [], color=warm[-1], lw=6,
             label='base level held, relaxing (%.0f-%.0f kyr)' % (t_force_kyr, t_total_kyr))
ax_prof.set_ylabel('Bed elevation [m]', fontsize=13)
ax_prof.set_title(
    'Aggradation: overbank deposition confines the gravel and speeds bed rise\n'
    'then base level holds and the profile relaxes to equilibrium',
    fontsize=12)
ax_prof.legend(fontsize=10, loc='upper right')

# f_ch at two times, colour-matched to the profile phases: end of forcing (cool)
# vs end of experiment (warm).  Without deposition f_ch = 1 at both times; with
# deposition it is < 1 at the end of forcing and relaxes back toward 1 as the
# aggradation ceases.
ax_fch.plot(x / 1e3, fch_off_force, color='0.6', ls='--', lw=1.5,
            label='no deposition (f_ch = 1, both times)')
ax_fch.plot(x / 1e3, fch_on_force, color=cool[-1], lw=2.5,
            label='with deposition, end of forcing')
ax_fch.plot(x / 1e3, fch_on_final, color=warm[-1], lw=2.5,
            label='with deposition, end of experiment')
ax_fch.set_ylim(0, 1.05)
ax_fch.set_xlabel('Downstream distance [km]', fontsize=13)
ax_fch.set_ylabel('Channel-deposit\nfraction f_ch [-]', fontsize=12)
ax_fch.legend(fontsize=10, loc='lower right')

fig.tight_layout()
fig.savefig('aggradation_deposition.png', dpi=150)
print('Wrote aggradation_deposition.png')
print('bed rise at x0 -- no deposition: end of forcing %.1f m, final %.1f m'
      % (force_off[-1][0] - force_off[0][0], relax_off[-1][0] - force_off[0][0]))
print('bed rise at x0 -- with deposition: end of forcing %.1f m, final %.1f m'
      % (force_on[-1][0] - force_on[0][0], relax_on[-1][0] - force_on[0][0]))
plt.show()
