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


def _spun_up(D=0.05):
    """A steady single-segment profile with grain size set (for b / h)."""
    lp = make_long_profile()
    lp.evolve_threshold_width_river(nt=40, dt=1.0e12)
    lp.D = D
    return lp


def _evolve_aggrading(lp, nt=10, dt=1.0e11, rate=1.0e-3 / _YEAR):
    """Drive aggradation by raising base level, one step at a time.

    Base-level rise is the aggradation driver (a positive uplift rate would make
    the channel *incise* into the rising valley -- see the uplift tests below),
    so it is what actually exercises the overbank deposit partition."""
    for _ in range(nt):
        lp.set_z_bl(lp.z_bl + rate * dt)
        lp.evolve_threshold_width_river(nt=1, dt=dt)


def test_f_ch_unity_reproduces_standard_run():
    """f_ch = 1 (high-migration limit, no widening) matches the standard run."""
    # Standard run: no valley dynamics, so f_ch stays at its default of 1.
    lp_std = _spun_up()
    _evolve_aggrading(lp_std)
    z_std = lp_std.z.copy()

    # Same run with the deposit-partition path active in the high-migration
    # limit: the closure drives f_ch -> 1 exactly (the exponential underflows),
    # and with no channel-belt coefficient set there is no widening, so the
    # effective width f_ch * B equals B and the profile must match to machine
    # precision.
    lp_val = _spun_up()
    lp_val.set_lateral_migration_rate(1.0e6)   # huge -> f_ch -> 1
    lp_val.set_valley_dynamics(partition_by_aggradation=True)
    _evolve_aggrading(lp_val)

    assert np.allclose(lp_val.f_ch, 1.0), \
        "the high-migration limit should give f_ch == 1 everywhere"
    assert np.allclose(z_std, lp_val.z, rtol=0.0, atol=1.0e-9), \
        "f_ch = 1 must reproduce the standard constant-width run exactly"


def test_f_ch_below_unity_changes_the_profile():
    """f_ch < 1 (a finite migration rate) actually couples: the profile differs.

    Guards against the f_ch -> effective-width coupling being reverted or
    orphaned (folding f_ch into storage alone cancels in the volume-first
    row-scaling and has no effect)."""
    lp_std = _spun_up()
    _evolve_aggrading(lp_std)
    z_std = lp_std.z.copy()

    lp_val = _spun_up()
    lp_val.set_lateral_migration_rate(3.0e-10)   # small -> f_ch < 1
    lp_val.set_valley_dynamics(partition_by_aggradation=True)
    _evolve_aggrading(lp_val)

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


def test_uplift_drives_incision_not_overbank_deposition():
    """Under uplift the channel incises into the rising valley, so the deposit
    partition must NOT fire: f_ch stays 1.

    The absolute bed rises under uplift (uplift adds rock faster than the river
    erodes it), so a gate on the raw dz/dt would spuriously switch on overbank
    deposition.  The rule keys off dz/dt - U (the bed's motion relative to the
    uplifting valley floor), which is negative here -- the river incises."""
    lp = _spun_up()
    lp.set_uplift_rate(1.0e-3 / _YEAR)          # uplift -> incision into rising valley
    lp.set_lateral_migration_rate(1.0e-9)
    lp.set_valley_dynamics(partition_by_aggradation=True)
    lp.evolve_threshold_width_river(nt=20, dt=1.0e11)
    assert np.allclose(lp.f_ch, 1.0), \
        "uplift must not trigger overbank deposition (the river incises)"


def test_uplift_drives_valley_narrowing():
    """Under uplift the channel entrenches into the rising valley, so narrowing
    fires and the valley collapses toward the channel width b -- even though the
    absolute bed elevation is rising."""
    lp = _spun_up()
    lp.set_uplift_rate(1.0e-3 / _YEAR)
    lp.set_lateral_migration_rate(0.0)
    lp.set_valley_dynamics(narrow_by_incision=True)
    B0 = lp.B.copy()
    lp.evolve_threshold_width_river(nt=20, dt=1.0e11)
    assert lp.B.max() < B0.min(), \
        "uplift must entrench the channel: the valley narrows toward b"


def test_spatially_variable_uplift_splits_narrowing_and_deposition():
    """Spatially variable U: an uplifting reach incises (narrows, no overbank)
    while a subsiding reach aggrades relative to its dropping floor (overbank
    deposition, f_ch < 1).  Exercises the per-node dz/dt - U gate."""
    lp = _spun_up()
    n = len(lp.x)
    U = 1.0e-3 / _YEAR
    lp.set_uplift_rate(np.where(np.arange(n) < n // 2, +U, -U))
    lp.set_lateral_migration_rate(1.0e-9)
    lp.set_valley_dynamics(narrow_by_incision=True, partition_by_aggradation=True)
    B0 = lp.B.copy()
    lp.evolve_threshold_width_river(nt=20, dt=1.0e11)
    up = slice(0, n // 2)
    down = slice(n // 2, n)
    assert np.allclose(lp.f_ch[up], 1.0), \
        "the uplifting reach must not deposit overbank (it incises)"
    assert lp.B[up].max() < B0[up].min(), \
        "the uplifting reach must narrow toward b"
    assert lp.f_ch[down].min() < 0.9, \
        "the subsiding reach must deposit overbank (f_ch < 1)"


def test_subsidence_drives_aggradation_not_narrowing():
    """Conjugate of the uplift case: subsidence under a constant base level makes
    the channel aggrade relative to its dropping floor, so overbank deposition
    fires (f_ch < 1) and the valley does NOT narrow.  Confirms dz/dt - U has the
    right sign for U < 0 too."""
    lp = _spun_up()
    lp.set_uplift_rate(-1.0e-3 / _YEAR)          # subsidence -> aggradation
    lp.set_lateral_migration_rate(1.0e-9)
    lp.set_valley_dynamics(narrow_by_incision=True, partition_by_aggradation=True)
    B0 = lp.B.copy()
    lp.evolve_threshold_width_river(nt=20, dt=1.0e11)
    assert lp.f_ch.min() < 0.9, \
        "subsidence must trigger overbank deposition (f_ch < 1)"
    assert np.allclose(lp.B, B0), \
        "subsidence must not narrow the valley (the channel aggrades, not incises)"
