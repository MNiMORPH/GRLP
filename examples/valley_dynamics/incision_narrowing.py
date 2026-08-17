# incision_narrowing.py
#
# Valley-dynamics demo (2 of 2): incision with and without valley narrowing.
#
# This is the mirror image of aggradation_deposition.py.  A gravel long profile
# is spun up to steady state and then forced to incise by a falling base level
# (in GRLP a falling base level and a positive uplift rate are exact mirror
# images: they create relief at the same rate, with the same profile response).
# We compare two rivers:
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
# incises more.  The lower panel shows the valley width B, the mechanism behind
# the difference.
#
# The lateral migration rate is zero here (no widening), isolating the narrowing
# process.

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
    lp.set_B(k_xB=25.0, P_xB=0.2)          # this B is now the *valley* width
    lp.set_uplift_rate(0.0)
    lp.set_niter(3)
    lp.set_Qs_input_upstream(lp.k_Qs * lp.Q[0] * S0 ** (7 / 6.0))
    lp.set_z_bl(0.0)
    lp.evolve_threshold_width_river(nt=60, dt=1.0e13)
    return lp


# Forcing and sampling (mirror of the aggradation demo: same rate and cadence).
FALL_RATE = 1.0e-3 / YEAR   # base-level fall rate [m/s] (= mirror of uplift)
GRAIN_D = 0.05             # grain size [m], needed to resolve channel width & depth
N_SNAP = 6                 # number of long-profile snapshots
NT_SNAP = 10               # solver steps between snapshots
DT = 1.0e11                # time step [s]


def run(with_narrowing):
    """Evolve an incising profile, returning (x, z-snapshots, B, base level)."""
    lp = build_steady_profile()
    lp.D = GRAIN_D
    if with_narrowing:
        lp.set_lateral_migration_rate(0.0)     # no widening
        lp.set_valley_dynamics(narrow_by_incision=True)
    snapshots = [lp.z.copy()]
    z_bl = [lp.z_bl]
    # Step one at a time so base level falls smoothly rather than in jumps.
    for _ in range(N_SNAP):
        for _ in range(NT_SNAP):
            lp.set_z_bl(lp.z_bl - FALL_RATE * DT)
            lp.evolve_threshold_width_river(nt=1, dt=DT)
        snapshots.append(lp.z.copy())
        z_bl.append(lp.z_bl)
    B = np.broadcast_to(np.asarray(lp.B, dtype=float), lp.x.shape).copy()
    return lp.x, snapshots, B, z_bl


x, z_wide, B_wide, zbl = run(with_narrowing=False)
_, z_narrow, B_narrow, _ = run(with_narrowing=True)

times_kyr = np.arange(N_SNAP + 1) * NT_SNAP * DT / YEAR / 1.0e3

fig, (ax_prof, ax_B) = plt.subplots(
    2, 1, figsize=(9, 8), height_ratios=[2, 1], sharex=True)

# Long profiles: colour ramps light (early) -> dark (late).
shades = plt.cm.viridis(np.linspace(0.15, 0.9, N_SNAP + 1))
for zs, sh in zip(z_wide, shades):
    ax_prof.plot(x / 1e3, zs, color=sh, lw=1.5, ls='--')
for zs, sh in zip(z_narrow, shades):
    ax_prof.plot(x / 1e3, zs, color=sh, lw=2.5)
ax_prof.plot([], [], 'k--', lw=1.5, label='no narrowing (valley stays width B)')
ax_prof.plot([], [], 'k-', lw=2.5, label='with narrowing (valley -> channel width b)')
ax_prof.set_ylabel('Bed elevation [m]', fontsize=13)
ax_prof.set_title(
    'Incision: narrowing concentrates the gravel and deepens incision\n'
    '(colour: light = early, dark = late; %.0f kyr total)' % times_kyr[-1],
    fontsize=12)
ax_prof.legend(fontsize=11, loc='upper right')

# Valley width: the mechanism.  Full B without narrowing; collapses to b with it.
ax_B.semilogy(x / 1e3, B_wide, 'k--', lw=1.5, label='no narrowing')
ax_B.semilogy(x / 1e3, B_narrow, 'k-', lw=2.5, label='with narrowing')
ax_B.set_xlabel('Downstream distance [km]', fontsize=13)
ax_B.set_ylabel('Valley width B [m]', fontsize=13)
ax_B.legend(fontsize=11, loc='center right')

fig.tight_layout()
fig.savefig('incision_narrowing.png', dpi=150)
print('Wrote incision_narrowing.png')
print('bed drop at x0 -- no narrowing: %.1f m ; with narrowing: %.1f m'
      % (z_wide[0][0] - z_wide[-1][0], z_narrow[0][0] - z_narrow[-1][0]))
plt.show()
