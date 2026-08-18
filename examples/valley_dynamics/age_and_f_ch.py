# age_and_f_ch.py
#
# Record two deposit attributes for the same valley fill -- the channel-deposit
# fraction f_ch and the deposit age (the model time at which each layer was laid
# down) -- and show them side by side, for two runs:
#
#   1. overbank deposition ON  (f_ch < 1, the partition feedback active), and
#   2. the "standard" run       (overbank off, so f_ch = 1 everywhere).
#
# Both are driven by a period of base-level rise, then base level is held until
# the profile relaxes to equilibrium.  The f_ch panel shows the facies; the age
# panel is the chronostratigraphy -- and because aggradation is fast during the
# rise and slows to nothing at equilibrium, most of the elapsed time is recorded
# in a thin interval near the top (a condensed section).

import numpy as np
from matplotlib import pyplot as plt

import grlp

YEAR = 3.15e7
RISE_RATE = 2.0e-3 / YEAR
ZETADOT = 1.5e-8
GRAIN_D = 0.05
N_RISE, N_HOLD, DT = 30, 150, 1.0e11   # rise, then hold to equilibrium
TOL = 0.02


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


def run_case(overbank):
    lp = build_steady_profile()
    lp.D = GRAIN_D
    if overbank:
        lp.set_lateral_migration_rate(ZETADOT)
        lp.set_valley_dynamics(partition_by_aggradation=True)
    lp.set_stratigraphic_recording(tol=TOL, record_age=True, record_thickness=True)
    for _ in range(N_RISE):
        lp.set_z_bl(lp.z_bl + RISE_RATE * DT)
        lp.evolve_threshold_width_river(nt=1, dt=DT)
    for _ in range(N_HOLD):
        lp.evolve_threshold_width_river(nt=1, dt=DT)
    return lp


def make_figure(lp, label, fname):
    t_final = lp.t
    fig, (ax_f, ax_a) = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

    # f_ch section
    _, pcf = lp.strat.plot_section(lp.x / 1e3, attribute='f_ch', ax=ax_f)
    fig.colorbar(pcf, ax=ax_f, label='channel-deposit fraction f_ch')
    ax_f.set_xlabel('Downstream distance [km]')
    ax_f.set_ylabel('Elevation [m]')
    ax_f.set_title('f_ch (facies)')

    # age section: age before the end of the run, in kyr (younger at top)
    Xa, Za, Va = lp.strat.section(lp.x / 1e3, attribute='age')
    age_kyr = (t_final - Va) / YEAR / 1e3
    pca = ax_a.pcolormesh(Xa, Za, age_kyr, shading='gouraud', cmap='plasma')
    fig.colorbar(pca, ax=ax_a, label='deposit age [kyr before end]')
    ax_a.set_xlabel('Downstream distance [km]')
    ax_a.set_title('deposit age (chronostratigraphy)')

    fig.suptitle(label, fontsize=13)
    fig.tight_layout()
    fig.savefig(fname, dpi=150)
    print('Wrote %s' % fname)


lp_ov = run_case(overbank=True)
make_figure(lp_ov, 'Overbank deposition ON (f_ch < 1)', 'age_and_f_ch_overbank.png')

lp_std = run_case(overbank=False)
make_figure(lp_std, 'Standard run (overbank off, f_ch = 1)', 'age_and_f_ch_standard.png')

plt.show()
