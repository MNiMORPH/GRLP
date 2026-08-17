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
