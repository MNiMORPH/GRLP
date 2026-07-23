#!/usr/bin/env bash
#
# Build the GRLP wheel and serve the standalone interactive demos in JupyterLite.
# GRLP runs entirely in your browser (Pyodide) -- no server, no local kernel.
#
# One-time toolchain install (the *frontends* jupyterlab_widgets and ipympl must
# be in this build env so `jupyter lite build` bundles the widget + matplotlib
# managers into the app; jupyter-server is needed to add the notebooks):
#
#     pip install jupyterlite-core jupyterlite-pyodide-kernel \
#                 jupyterlab_widgets ipympl jupyter-server
#
set -euo pipefail
cd "$(dirname "$0")"

# Guard: the widget/plot frontends must be present, or the sliders and the live
# canvas silently fail to render. (The notebooks' own `%pip install` only covers
# the in-browser Python kernel, not these frontends.)
for pkg in jupyterlab_widgets ipympl; do
  if ! python3 -c "import ${pkg}" 2>/dev/null; then
    echo "ERROR: ${pkg} is not installed in this build environment." >&2
    echo "Install the frontends so JupyterLite bundles them:" >&2
    echo "    pip install jupyterlab_widgets ipympl" >&2
    exit 1
  fi
done

# Clean any previous build state (a stale config or output can collide).
rm -rf _output .jupyterlite.doit.db jupyter_lite_config.json pypi

# Build the GRLP wheel into the pypi/ wheelhouse. jupyterlite-pyodide-kernel
# auto-discovers and indexes wheels placed in pypi/, so `%pip install grlp` in
# the notebooks resolves this wheel in the browser. (Do NOT also list it in a
# jupyter_lite_config.json piplite_urls -- that double-targets the same file.)
mkdir -p pypi
python3 -m pip wheel .. --no-deps -w pypi
echo "Built wheel: $(ls pypi/grlp-*.whl)"

# Build and serve JupyterLite (Pyodide downloads on first run). Two notebooks:
#   interactive_single_segment.ipynb       -- equilibrium response (sliders)
#   interactive_single_segment_live.ipynb  -- live transient (runs continuously)
jupyter lite serve \
    --contents interactive_single_segment.ipynb \
    --contents interactive_single_segment_live.ipynb
