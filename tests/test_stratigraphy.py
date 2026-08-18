"""
Valley-fill stratigraphic recorder (#19).

The StratRecorder stores per-node (z, f_ch) polylines compressed with
vertical-metric Douglas-Peucker.  The properties that matter:

  * the compression bounds the f_ch reconstruction error by the tolerance, and
    the stored vertex count is set by the f_ch shape, not the time step;
  * incision scours the record (unconformities);
  * with record_thickness the compressed section conserves deposit volume; and
  * recording is a passive diagnostic -- a run with it on is byte-identical to
    one without.
"""

import numpy as np

import grlp
from grlp.stratigraphy import StratRecorder, douglas_peucker_vertical
from conftest import make_long_profile


_YEAR = 3.15e7


def test_dp_reduces_a_linear_ramp_to_its_endpoints():
    z = np.linspace(0.0, 10.0, 60)
    f = 0.3 + 0.05 * z
    assert douglas_peucker_vertical(z, f, 0.02) == [0, len(z) - 1]


def test_dp_holds_a_sharp_inflection():
    z = np.linspace(0.0, 10.0, 61)
    f = np.full_like(z, 0.3)
    f[30] = 0.9                                   # a spike well above tolerance
    idx = douglas_peucker_vertical(z, f, 0.02)
    assert 30 in idx, "the sharp inflection must be kept"


def test_dp_reconstruction_error_within_tolerance():
    rng = np.random.default_rng(0)
    z = np.sort(rng.uniform(0, 100, 400))
    f = 0.5 + 0.4 * np.sin(z / 10.0)
    tol = 0.02
    idx = douglas_peucker_vertical(z, f, tol)
    f_rec = np.interp(z, z[idx], f[idx])
    assert np.max(np.abs(f_rec - f)) <= tol + 1e-12


def test_stored_vertex_count_is_dt_independent():
    """Refining the sampling (more raw points on the same f_ch(z) curve) leaves
    the compressed vertex count essentially unchanged -- storage is set by the
    composition tolerance, not the step count."""
    def compressed_count(nsamp):
        z = np.linspace(0.0, 100.0, nsamp)
        f = 0.2 + 0.6 * np.sin(z / 30.0) ** 2
        rec = StratRecorder(1, tol=0.02, max_raw=None)
        for zi, fi in zip(z, f):
            rec.record(np.array([zi]), fi)
        return len(rec.columns()[0][0])
    counts = [compressed_count(n) for n in (200, 800, 3200)]
    assert max(counts) - min(counts) <= 2, \
        "compressed count should not grow with sampling density: %r" % counts


def test_incision_scours_the_record():
    """Aggrade, then incise below part of the column: the scoured samples are
    removed (an unconformity) and the top of the record sits at the new bed."""
    rec = StratRecorder(1, tol=0.0, max_raw=None)
    for z in np.arange(1.0, 21.0):                # aggrade 1 -> 20
        rec.record(np.array([z]), 0.5)
    rec.record(np.array([12.0]), 0.5)             # incise back to 12
    zc, _ = rec.columns()[0]
    assert zc[-1] == 12.0, "record top must follow the incised bed"
    assert zc.max() <= 12.0, "everything above the new bed must be scoured"


def test_volume_closure_within_tolerance():
    """The compressed section's integral(f_ch dz) matches the exact accumulated
    sum(f_ch dz) to within ~tol * thickness."""
    rec = StratRecorder(1, tol=0.02, max_raw=None, record_thickness=True)
    z = 0.0
    for k in range(300):
        z += 1.0
        rec.record(np.array([z]), 0.2 + 0.6 * k / 300.0)
    compressed, exact, maxdiff = rec.volume_closure()
    assert maxdiff <= 0.02 * z, \
        "volume closure %.3f exceeds tol*thickness %.3f" % (maxdiff, 0.02 * z)


def _aggrading_profile(record=False):
    lp = make_long_profile()
    lp.evolve_threshold_width_river(nt=40, dt=1.0e12)
    lp.D = 0.05
    lp.set_lateral_migration_rate(1.0e-9)
    lp.set_valley_dynamics(partition_by_aggradation=True)
    if record:
        lp.set_stratigraphic_recording(tol=0.02)
    return lp


def test_recording_does_not_perturb_the_solve():
    """Recording is a passive diagnostic: the bed evolution is byte-identical
    with it on or off."""
    lp_off = _aggrading_profile(record=False)
    lp_on = _aggrading_profile(record=True)
    for _ in range(10):
        for lp in (lp_off, lp_on):
            lp.set_z_bl(lp.z_bl + 1.0e-3 / _YEAR * 1.0e11)
            lp.evolve_threshold_width_river(nt=1, dt=1.0e11)
    assert np.array_equal(lp_off.z, lp_on.z), \
        "stratigraphic recording must not change the solve"


def _uniform_demo_profile(nx=60):
    """The valley-demo build (uniform Q and B, straight steady profile), the
    setting the stratigraphy examples use."""
    S0 = 1.5e-2
    lp = grlp.LongProfile()
    lp.set_intermittency(0.8)
    lp.basic_constants(); lp.bedload_lumped_constants(); lp.set_hydrologic_constants()
    lp.set_x(dx=1000.0, nx=nx, x0=10000.0)
    lp.set_z(S0=-S0)
    lp.set_A(k_xA=1.0)
    lp.set_Q(Q=10.0)
    lp.set_B(B=200.0)
    lp.set_uplift_rate(0.0)
    lp.set_niter(3)
    lp.set_Qs_input_upstream(lp.k_Qs * lp.Q[0] * S0 ** (7 / 6.0))
    lp.set_z_bl(0.0)
    lp.evolve_threshold_width_river(nt=40, dt=1.0e13)
    return lp


def test_stratigraphy_golden_demo_column():
    """Golden master for the tuned uniform stratigraphy demo (examples/
    valley_dynamics/stratigraphy_section.py): the mid-node f_ch column reproduces
    the recorded mixed-body -> channel-rich-cap profile.  Guards the whole
    valley-dynamics + recorder path (deposit-partition closure, in-Picard default,
    Douglas-Peucker compression) against numerical drift."""
    lp = _uniform_demo_profile(nx=60)
    lp.D = 0.05
    lp.set_lateral_migration_rate(1.5e-8)
    lp.set_valley_dynamics(partition_by_aggradation=True)
    lp.set_stratigraphic_recording(tol=0.02, record_thickness=True)
    rate, dt = 2.0e-3 / _YEAR, 1.0e11
    for _ in range(30):                          # base-level rise
        lp.set_z_bl(lp.z_bl + rate * dt)
        lp.evolve_threshold_width_river(nt=1, dt=dt)
    for _ in range(40):                          # held (relaxation)
        lp.evolve_threshold_width_river(nt=1, dt=dt)

    zc, fc = lp.strat.columns()[30]
    zn = (zc - zc[0]) / (zc[-1] - zc[0])
    deciles = np.array([np.interp(q, zn, fc) for q in (0.1, 0.3, 0.5, 0.7, 0.9)])
    golden = np.array([0.4815, 0.4101, 0.4092, 0.4082, 0.5368])
    np.testing.assert_allclose(
        deciles, golden, atol=2.0e-3,
        err_msg="stratigraphy f_ch column drifted from the golden master")
    assert abs(zc[-1] - 640.0919) < 0.05, "column top elevation drifted"
    _, _, maxd = lp.strat.volume_closure()
    assert maxd < 0.02 * (zc[-1] - zc[0]), "volume closure exceeds tol*thickness"


def test_recording_builds_a_section():
    """A short aggradation run produces a reconstructable f_ch(x, z) section with
    f_ch < 1 (overbank deposition recorded)."""
    lp = _aggrading_profile(record=True)
    for _ in range(20):
        lp.set_z_bl(lp.z_bl + 1.0e-3 / _YEAR * 1.0e11)
        lp.evolve_threshold_width_river(nt=1, dt=1.0e11)
    _, _, F = lp.strat.section(lp.x)
    assert np.nanmin(F) < 0.9, "the section should record overbank deposition"
    assert np.nanmax(F) <= 1.0 + 1e-9
