#!/usr/bin/env bash
#
# Build the GRLP wheel and serve the standalone interactive demo in JupyterLite.
# GRLP runs entirely in your browser (Pyodide) -- no server, no local kernel.
#
# One-time toolchain install:
#     pip install jupyterlite-core jupyterlite-pyodide-kernel
#
set -euo pipefail
cd "$(dirname "$0")"

# 1. Build the GRLP wheel from the repository and put it in the wheelhouse.
mkdir -p pypi
python3 -m pip wheel .. --no-deps -w pypi

# 2. Point piplite (JupyterLite's installer) at the freshly built wheel, so
#    `%pip install grlp` in the notebook resolves it. Version-agnostic.
WHEEL="$(ls pypi/grlp-*.whl | head -1)"
cat > jupyter_lite_config.json <<EOF
{
  "PipliteAddon": {
    "piplite_urls": ["${WHEEL}"]
  }
}
EOF
echo "Using wheel: ${WHEEL}"

# 3. Build and serve JupyterLite (Pyodide downloads on first run). Open
#    interactive_single_segment.ipynb in the browser tab that appears and run it.
jupyter lite serve --contents interactive_single_segment.ipynb
