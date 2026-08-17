"""
Valley dynamics: channel-deposit fraction f_ch coupling (#19).

The channel-deposit fraction ``f_ch`` reduces the effective width the channel
gravel occupies: the Exner storage *and* transport act over ``f_ch * B`` (the
solver divides the flux conductance by ``f_ch * B``, not ``B``).  The key
regression property is that ``f_ch = 1`` must reproduce the standard
constant-width run *exactly* -- a smaller ``f_ch`` fills a narrower width and so
must aggrade faster.

These tests exercise the full valley-dynamics path (``update_valley`` computes
``f_ch`` from the deposit closure each step) rather than poking ``f_ch`` by hand:

  * in the high-migration limit the closure gives ``f_ch -> 1`` (the exponential
    underflows), so with no valley widening the profile must match a standard
    run to machine precision; and
  * a small migration rate gives ``f_ch < 1``, which must change the profile --
    guarding against the coupling being silently reverted.
"""

import numpy as np

import grlp
from conftest import make_long_profile


_YEAR = 3.15e7
_AGGRADATION = 1.0e-4 / _YEAR   # U > 0 aggrades the bed in GRLP's sign convention


def _spun_up_aggrading():
    """A steady single-segment profile with an aggradation forcing applied."""
    lp = make_long_profile()
    lp.evolve_threshold_width_river(nt=40, dt=1.0e12)
    lp.set_uplift_rate(_AGGRADATION)
    return lp


def test_f_ch_unity_reproduces_standard_run():
    """f_ch = 1 (high-migration limit, no widening) matches the standard run."""
    # Standard run: no valley dynamics, so f_ch stays at its default of 1.
    lp_std = _spun_up_aggrading()
    lp_std.evolve_threshold_width_river(nt=10, dt=1.0e11)
    z_std = lp_std.z.copy()

    # Same run with the deposit-partition path active in the high-migration
    # limit: the closure drives f_ch -> 1 exactly (the exponential underflows),
    # and with no channel-belt coefficient set there is no widening, so the
    # effective width f_ch * B equals B and the profile must match to machine
    # precision.
    lp_val = _spun_up_aggrading()
    lp_val.D = 0.05
    lp_val.set_lateral_migration_rate(1.0e6)   # huge -> f_ch -> 1
    lp_val.set_valley_dynamics(partition_by_aggradation=True)
    lp_val.evolve_threshold_width_river(nt=10, dt=1.0e11)

    assert np.allclose(lp_val.f_ch, 1.0), \
        "the high-migration limit should give f_ch == 1 everywhere"
    assert np.allclose(z_std, lp_val.z, rtol=0.0, atol=1.0e-9), \
        "f_ch = 1 must reproduce the standard constant-width run exactly"


def test_f_ch_below_unity_changes_the_profile():
    """f_ch < 1 (a finite migration rate) actually couples: the profile differs.

    Guards against the f_ch -> effective-width coupling being reverted or
    orphaned (folding f_ch into storage alone cancels in the volume-first
    row-scaling and has no effect)."""
    lp_std = _spun_up_aggrading()
    lp_std.evolve_threshold_width_river(nt=10, dt=1.0e11)
    z_std = lp_std.z.copy()

    lp_val = _spun_up_aggrading()
    lp_val.D = 0.05
    lp_val.set_lateral_migration_rate(3.0e-10)   # small -> f_ch < 1
    lp_val.set_valley_dynamics(partition_by_aggradation=True)
    lp_val.evolve_threshold_width_river(nt=10, dt=1.0e11)

    assert lp_val.f_ch.min() < 0.9, \
        "a small migration rate should give f_ch meaningfully below 1"
    assert not np.allclose(z_std, lp_val.z), \
        "f_ch < 1 must change the profile (the gravel fills a narrower width)"


def _evolve_Y_network(f_ch=None):
    """A confluence (Y) network under uplift, optionally with a uniform f_ch.

    Two tributaries (0, 1) join into a trunk (2); the trunk's first node is the
    confluence.  Uplift makes every cell -- interiors *and* the junction cell --
    a sediment source over its land area, so the confluence assembly is exercised
    (its ``land_area`` is built from the tributaries' terminal ``f_ch * B`` and
    the confluence node's own ``f_ch * B``).  Returns the three segment ``z``.
    """
    nseg = 30
    dx = 1000.0
    xt = dx * np.arange(1, nseg + 1.0)
    xtr = dx * np.arange(nseg + 1, 2 * nseg + 1.0)
    U = 5.0e-4 / _YEAR                              # ~0.5 mm/yr uplift -> aggrades
    net = grlp.Network()
    net.initialize(
        x_bl=dx * (2 * nseg + 1), z_bl=0.0, S0=[1.5e-2, 1.5e-2], Q_s_0=None,
        upstream_segment_IDs=[[], [], [0, 1]],
        downstream_segment_IDs=[[2], [2], []],
        x=[xt.copy(), xt.copy(), xtr.copy()], z=[np.zeros(nseg)] * 3,
        Q=[10.0 * np.ones(nseg), 10.0 * np.ones(nseg), 20.0 * np.ones(nseg)],
        B=[100.0 * np.ones(nseg)] * 3,
    )
    for lp in net.list_of_LongProfile_objects:
        lp.set_intermittency(1.0)
        lp.set_uplift_rate(U)
        if f_ch is not None:
            lp.f_ch = f_ch * np.ones(len(lp.x))
    net.set_niter(4)
    net.get_z_lengths()
    net.evolve_threshold_width_river_network(nt=6000, dt=1.0e12)
    return [lp.z.copy() for lp in net.list_of_LongProfile_objects]


def test_f_ch_unity_reproduces_standard_network_run():
    """Explicit f_ch = 1 everywhere reproduces the standard confluence run.

    Guards the confluence f_ch coupling from corrupting the junction cell at
    unity: the tributary-terminal and confluence-node land areas fold in f_ch,
    which must be a no-op when f_ch = 1."""
    z_std = _evolve_Y_network(f_ch=None)
    z_one = _evolve_Y_network(f_ch=1.0)
    for i, (a, b) in enumerate(zip(z_std, z_one)):
        assert np.array_equal(a, b), \
            "f_ch = 1 must reproduce the standard network run exactly (seg %d)" % i


def test_f_ch_below_unity_couples_at_the_confluence():
    """f_ch < 1 changes the network profile, the confluence node included.

    The confluence cell's sediment source is generated over its land area; with
    f_ch < 1 that area shrinks to f_ch * B, so the junction elevation must
    differ.  Guards the confluence coupling (solver ``land_area`` / ``A_confluence``
    terms) from being reverted or left on plain ``B``."""
    z_std = _evolve_Y_network(f_ch=None)
    z_val = _evolve_Y_network(f_ch=0.4)
    # every segment shifts, and in particular the confluence node (trunk[0])
    assert not np.allclose(z_std[2], z_val[2]), \
        "f_ch < 1 must change the trunk profile"
    assert abs(z_std[2][0] - z_val[2][0]) > 1.0, \
        "f_ch < 1 must move the confluence node elevation (couples at the junction)"
