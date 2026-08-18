# incision_narrowing.py
#
# Valley-dynamics demo (2 of 2): incision with and without valley narrowing.
#
# This is the mirror image of aggradation_deposition.py.  A gravel long profile
# is spun up to steady state and then run in two phases:
#
#   Phase 1 (forcing) -- base level falls at a constant rate, driving incision
#     (the same base-level knob as the aggradation demo, opposite sign).
#   Phase 2 (relaxation) -- base level is held fixed and the profile relaxes the
#     rest of the way toward its new, lower equilibrium.
#
# In each phase we compare two rivers:
#
#   * WITHOUT narrowing -- the valley keeps its full width B as the bed incises,
#     and
#   * WITH narrowing    -- where the bed incises the channel abandons its
#     floodplain and, with vertical walls, the valley collapses to the channel
#     width b while the walls grow by the incision depth (Turowski-style vertical
#     entrenchment).
#
# The narrowed valley concentrates the gravel storage into a much smaller width
# (b << B), so the channel bed responds faster to the falling base level and
# incises more during forcing.
#
# Lateral migration is now switched on (the "channel lateral mobility"): it
# widens the valley back toward its prescribed width B_max (Turowski et al.,
# 2025), competing with the incision-driven narrowing.  At this migration rate it
# is gentle -- it keeps the valley from fully collapsing during forcing and lets
# it re-widen a little once incision ceases in Phase 2.  The lower panel shows the
# valley width B at both the end of forcing and the end of the experiment, so the
# re-widening is visible.
#
# For simplicity in these tests the initial valley width is uniform along the
# whole domain, so B_max is a single value everywhere.

import numpy as np
from matplotlib import pyplot as plt

import grlp

YEAR = 3.15e7


def build_steady_profile():
    """A canonical gravel long profile evolved to steady state (no forcing).

    Identical spin-up to aggradation_deposition.py (Wickert & Schildgen, 2019,
    Fig. 2 geometry), so the two demos start from the same river.
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


# Forcing and sampling (mirror of the aggradation demo: same rate and cadence,
# opposite sign of base-level change).
FALL_RATE = 1.0e-3 / YEAR   # base-level fall rate [m/s] -> incision
ZETADOT = 1.0e-9           # lateral migration rate [m/s] -> valley widening
GRAIN_D = 0.05             # grain size [m], needed to resolve channel width & depth
N_FORCE = 6                # long-profile snapshots during forcing
N_RELAX = 6                # long-profile snapshots during relaxation
NT_SNAP = 10               # solver steps between snapshots
DT = 1.0e11                # time step [s]


def run(with_narrowing):
    """Evolve through forcing then relaxation.

    Returns (x, forcing-snapshots, relaxation-snapshots, B-at-forcing-end,
    B-at-experiment-end).  The forcing list includes the initial steady profile.
    """
    lp = build_steady_profile()
    lp.D = GRAIN_D
    # Lateral migration (the channel mobility) widens the valley back toward
    # B_max in both cases; the toggle is the incision-narrowing feedback.
    lp.set_lateral_migration_rate(ZETADOT)
    lp.set_valley_dynamics(widen_by_migration=True)
    if with_narrowing:
        lp.set_valley_dynamics(narrow_by_incision=True)

    def B_now():
        return np.broadcast_to(np.asarray(lp.B, dtype=float), lp.x.shape).copy()

    # Phase 1: base level falls, one step at a time so it drops smoothly.
    forcing = [lp.z.copy()]
    for _ in range(N_FORCE):
        for _ in range(NT_SNAP):
            lp.set_z_bl(lp.z_bl - FALL_RATE * DT)
            lp.evolve_threshold_width_river(nt=1, dt=DT)
        forcing.append(lp.z.copy())
    B_force = B_now()                # valley width at the end of forcing

    # Phase 2: base level held fixed; the profile relaxes toward equilibrium.
    relaxation = []
    for _ in range(N_RELAX):
        lp.evolve_threshold_width_river(nt=NT_SNAP, dt=DT)
        relaxation.append(lp.z.copy())
    B_final = B_now()                # valley width at the end of the experiment

    return lp.x, forcing, relaxation, B_force, B_final


x, force_wide, relax_wide, B_wide_force, B_wide_final = run(with_narrowing=False)
_, force_narrow, relax_narrow, B_narrow_force, B_narrow_final = run(with_narrowing=True)

t_force_kyr = N_FORCE * NT_SNAP * DT / YEAR / 1.0e3
t_total_kyr = (N_FORCE + N_RELAX) * NT_SNAP * DT / YEAR / 1.0e3

fig, (ax_prof, ax_B) = plt.subplots(
    2, 1, figsize=(9, 8), height_ratios=[2, 1], sharex=True)

# Long profiles.  Forcing phase in a cool ramp, relaxation phase in a warm ramp;
# within each, light -> dark is early -> late.  Solid = with narrowing, dashed
# = without.
cool = plt.cm.viridis(np.linspace(0.15, 0.9, N_FORCE + 1))
warm = plt.cm.autumn_r(np.linspace(0.35, 0.95, N_RELAX))
for zs, c in zip(force_wide, cool):
    ax_prof.plot(x / 1e3, zs, color=c, lw=1.5, ls='--')
for zs, c in zip(force_narrow, cool):
    ax_prof.plot(x / 1e3, zs, color=c, lw=2.5)
for zs, c in zip(relax_wide, warm):
    ax_prof.plot(x / 1e3, zs, color=c, lw=1.5, ls='--')
for zs, c in zip(relax_narrow, warm):
    ax_prof.plot(x / 1e3, zs, color=c, lw=2.5)
ax_prof.plot([], [], 'k-', lw=2.5, label='with narrowing (valley -> channel width b)')
ax_prof.plot([], [], 'k--', lw=1.5, label='no narrowing (valley stays width B)')
ax_prof.plot([], [], color=cool[-1], lw=6, label='base level falling (0-%.0f kyr)' % t_force_kyr)
ax_prof.plot([], [], color=warm[-1], lw=6,
             label='base level held, relaxing (%.0f-%.0f kyr)' % (t_force_kyr, t_total_kyr))
ax_prof.set_ylabel('Bed elevation [m]', fontsize=13)
ax_prof.set_title(
    'Incision: narrowing concentrates the gravel and deepens incision\n'
    'then base level holds and the profile relaxes to equilibrium',
    fontsize=12)
ax_prof.legend(fontsize=10, loc='upper right')

# Valley width at two times, colour-matched to the profile phases: end of
# forcing (cool) vs end of experiment (warm).  Without narrowing B keeps its full
# width (migration does not widen a valley already wider than W0); with narrowing
# it collapses during incision, then re-widens by migration once incision ceases.
ax_B.semilogy(x / 1e3, B_wide_force, color='0.6', ls='--', lw=1.5,
              label='no narrowing (both times)')
ax_B.semilogy(x / 1e3, B_narrow_force, color=cool[-1], lw=2.5,
              label='with narrowing, end of forcing')
ax_B.semilogy(x / 1e3, B_narrow_final, color=warm[-1], lw=2.5,
              label='with narrowing, end of experiment')
ax_B.set_xlabel('Downstream distance [km]', fontsize=13)
ax_B.set_ylabel('Valley width B [m]', fontsize=12)
ax_B.legend(fontsize=10, loc='center right')

fig.tight_layout()
fig.savefig('incision_narrowing.png', dpi=150)
print('Wrote incision_narrowing.png')
print('bed drop at x0 -- no narrowing: end of forcing %.1f m, final %.1f m'
      % (force_wide[0][0] - force_wide[-1][0], force_wide[0][0] - relax_wide[-1][0]))
print('bed drop at x0 -- with narrowing: end of forcing %.1f m, final %.1f m'
      % (force_narrow[0][0] - force_narrow[-1][0], force_narrow[0][0] - relax_narrow[-1][0]))
plt.show()
