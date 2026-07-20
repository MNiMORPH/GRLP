"""
Characterization (golden-master) tests.

These re-run the configurations in ``characterization_configs`` and assert that
every recorded output array still matches the reference in
``characterization_reference.npz`` to tight tolerance. On the same machine the
model is bit-reproducible; the tolerance here allows only the small
floating-point differences a different platform/BLAS might introduce, so any
genuine change in behavior fails loudly.

Their job is to guard forthcoming work that makes the valley width ``B`` a
dynamic function of the other variables: when that code is configured to
reproduce a prescribed ``B``, these tests confirm the numerical results are
unchanged from what the prescribed-``B`` model produces today. The reference
deliberately includes transient snapshots and ``B``-loaded quantities
(diffusivity, uplift profiles, arbitrary ``B(x)``, networks with per-segment
widths), because the no-uplift equilibrium profile alone is insensitive to ``B``.

To intentionally update the reference after a knowing change, run
``python tests/generate_characterization_reference.py`` and review the diff.
"""

import os

import numpy as np
import pytest

from characterization_configs import run_all

REFERENCE_PATH = os.path.join(os.path.dirname(__file__),
                              "characterization_reference.npz")

# Tight enough to catch any real behavioral change (which would be orders of
# magnitude larger), loose enough to tolerate cross-platform float noise.
RTOL = 1e-9
ATOL = 1e-12


def _reference_keys():
    """Keys to parametrize over; empty if the reference has not been generated."""
    if not os.path.exists(REFERENCE_PATH):
        return []
    with np.load(REFERENCE_PATH) as npz:
        return sorted(npz.files)


@pytest.fixture(scope="module")
def reference():
    if not os.path.exists(REFERENCE_PATH):
        pytest.skip(
            "characterization_reference.npz missing; generate it with "
            "python tests/generate_characterization_reference.py"
        )
    with np.load(REFERENCE_PATH) as npz:
        return {key: npz[key] for key in npz.files}


@pytest.fixture(scope="module")
def current():
    return run_all()


def test_reference_and_current_have_same_keys(reference, current):
    # A drift in the set of configs/arrays should be an explicit, reviewed
    # reference regeneration, not a silent mismatch.
    assert set(current) == set(reference)


@pytest.mark.parametrize("key", _reference_keys())
def test_matches_reference(reference, current, key):
    ref = reference[key]
    cur = current[key]
    assert cur.shape == ref.shape, (
        "shape changed for %s: %s vs %s" % (key, cur.shape, ref.shape)
    )
    np.testing.assert_allclose(
        cur, ref, rtol=RTOL, atol=ATOL,
        err_msg="characterization mismatch for %s" % key,
    )
