"""
Characterization (golden-master) tests.

These re-run the configurations in ``characterization_configs`` and assert that
every recorded output array still matches the stored reference to tight
tolerance. On the same machine the model is bit-reproducible; the tolerance here
allows only the small floating-point differences a different platform/BLAS might
introduce, so any genuine change in behavior fails loudly.

There is **one golden set per time-integration scheme**: the original first-order
backward Euler (``characterization_reference.npz``) and the current default
second-order BDF2 iterated to convergence (``characterization_reference_bdf2.npz``).
The steady-state arrays are essentially identical between the two schemes; the
transient snapshots differ by the second-order accuracy gain (BDF2 lands several
times closer to a fine-step reference). Keeping both pins each scheme so a change
to either fails loudly.

Their job is to guard forthcoming work that makes the valley width ``B`` a
dynamic function of the other variables: when that code is configured to
reproduce a prescribed ``B``, these tests confirm the numerical results are
unchanged from what the prescribed-``B`` model produces today. The reference
deliberately includes transient snapshots and ``B``-loaded quantities
(diffusivity, uplift profiles, arbitrary ``B(x)``, networks with per-segment
widths), because the no-uplift equilibrium profile alone is insensitive to ``B``.

To intentionally update the references after a knowing change, run
``python tests/generate_characterization_reference.py`` and review the diff.
"""

import os

import numpy as np
import pytest

from characterization_configs import run_all

# scheme -> reference file
REFERENCE_FILES = {
    "euler": "characterization_reference.npz",
    "bdf2": "characterization_reference_bdf2.npz",
}

# Tight enough to catch any real behavioral change (which would be orders of
# magnitude larger), loose enough to tolerate cross-platform float noise.
RTOL = 1e-9
ATOL = 1e-12


def _path(scheme):
    return os.path.join(os.path.dirname(__file__), REFERENCE_FILES[scheme])


def _reference_keys(scheme):
    """Keys to parametrize over; empty if that reference has not been generated."""
    path = _path(scheme)
    if not os.path.exists(path):
        return []
    with np.load(path) as npz:
        return sorted(npz.files)


@pytest.fixture(scope="module")
def references():
    refs = {}
    for scheme in REFERENCE_FILES:
        path = _path(scheme)
        if os.path.exists(path):
            with np.load(path) as npz:
                refs[scheme] = {key: npz[key] for key in npz.files}
    return refs


@pytest.fixture(scope="module")
def currents():
    return {scheme: run_all(scheme=scheme) for scheme in REFERENCE_FILES}


@pytest.mark.parametrize("scheme", sorted(REFERENCE_FILES))
def test_reference_and_current_have_same_keys(references, currents, scheme):
    # A drift in the set of configs/arrays should be an explicit, reviewed
    # reference regeneration, not a silent mismatch.
    if scheme not in references:
        pytest.skip("%s reference missing; regenerate with "
                    "python tests/generate_characterization_reference.py"
                    % scheme)
    assert set(currents[scheme]) == set(references[scheme])


@pytest.mark.parametrize(
    "scheme,key",
    [(s, k) for s in sorted(REFERENCE_FILES) for k in _reference_keys(s)])
def test_matches_reference(references, currents, scheme, key):
    ref = references[scheme][key]
    cur = currents[scheme][key]
    assert cur.shape == ref.shape, (
        "shape changed for %s/%s: %s vs %s" % (scheme, key, cur.shape, ref.shape)
    )
    np.testing.assert_allclose(
        cur, ref, rtol=RTOL, atol=ATOL,
        err_msg="characterization mismatch for %s/%s" % (scheme, key),
    )
