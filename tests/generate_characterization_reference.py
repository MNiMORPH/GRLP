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

REFERENCE_PATH = os.path.join(os.path.dirname(__file__),
                              "characterization_reference.npz")


def main():
    data = run_all()
    np.savez_compressed(REFERENCE_PATH, **data)
    print("Wrote %d reference arrays to %s" % (len(data), REFERENCE_PATH))


if __name__ == "__main__":
    main()
