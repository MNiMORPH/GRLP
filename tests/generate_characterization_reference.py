"""
Regenerate the characterization reference recording.

Run this ONLY when you have deliberately, knowingly changed the model's
numerical output and want the golden-master tests to adopt the new values:

    python tests/generate_characterization_reference.py

It overwrites tests/characterization_reference.npz with the current outputs of
every configuration in characterization_configs.run_all(). Review the resulting
diff (git diff --stat, and spot-check magnitudes) before committing: an
unexpected change here is exactly the regression these tests exist to catch.
"""

import os

import numpy as np

from characterization_configs import run_all

# One golden set per time-integration scheme: the original first-order backward
# Euler, and the current default second-order BDF2 (iterated to convergence).
# The steady-state arrays are essentially identical between the two; the
# transient snapshots differ by the second-order accuracy gain.
REFERENCE_PATHS = {
    "euler": os.path.join(os.path.dirname(__file__),
                          "characterization_reference.npz"),
    "bdf2": os.path.join(os.path.dirname(__file__),
                         "characterization_reference_bdf2.npz"),
}


def main():
    for scheme, path in REFERENCE_PATHS.items():
        data = run_all(scheme=scheme)
        np.savez_compressed(path, **data)
        print("Wrote %d %s reference arrays to %s"
              % (len(data), scheme, path))


if __name__ == "__main__":
    main()
