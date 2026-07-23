# Standalone interactive GRLP demos (experimental)

Browser-based, slider-driven demos of a single gravel-bed river responding to its
**water input**, **sediment input**, and **base level** — GRLP running entirely
in your browser via Pyodide / JupyterLite. Two notebooks:

- **`interactive_single_segment.ipynb`** — *equilibrium* response: each slider
  change recomputes the steady-state profile (Lane's balance).
- **`interactive_single_segment_live.ipynb`** — *live transient*: the model runs
  continuously; drag the sliders and watch the profile aggrade/incise toward the
  new boundary conditions in real time.

Deliberately kept **out of the ReadTheDocs site** for now: standalone, until
tested enough.

## Run it

One-time toolchain install (the **frontends** — `jupyterlab_widgets` for the
sliders, `ipympl` for the live canvas — must be in this env so `jupyter lite
build` bundles them; `jupyter-server` lets it add the notebooks):

```sh
pip install jupyterlite-core jupyterlite-pyodide-kernel \
            jupyterlab_widgets ipympl jupyter-server
```

These frontends are separate from the `ipywidgets`/`ipympl` the notebooks install
into the in-browser Python kernel — you need both sides. `run_demo.sh` checks for
the frontends and errors early if they are missing.

Then:

```sh
./run_demo.sh
```

This builds the GRLP wheel, points JupyterLite's installer at it, and serves the
site (Pyodide downloads on first run). In the browser tab that opens, open either
notebook, run all cells, and use the sliders.

## What to check (test checklist)

- [ ] `%pip install grlp networkx ipywidgets [ipympl]` resolves (GRLP from the
      bundled wheel).
- [ ] **Equilibrium:** the three sliders (Water Q, Sediment Qs, Base level)
      redraw the profile; Lane's balance reads correctly (more sediment →
      steeper; more water → gentler; lower base level → profile drops).
- [ ] **Live:** the profile evolves continuously; dragging sediment up makes it
      aggrade/steepen over ~seconds; dropping base level sends an incision wave
      upstream; Running/Reset behave.
- [ ] Frame rate is acceptable in-browser (native step is ~ms).

## Files

- `interactive_single_segment.ipynb` — equilibrium demo (sliders + GRLP compute).
- `interactive_single_segment_live.ipynb` — live transient demo (async run loop).
- `run_demo.sh` — build the wheel + launch JupyterLite.
- `pypi/`, `_output/`, `.jupyterlite.doit.db` — build artifacts (git-ignored).
  `pypi/` is the wheelhouse: `jupyterlite-pyodide-kernel` auto-discovers and
  indexes the `grlp` wheel there, so `%pip install grlp` resolves it.

## Status / notes

- **Verified natively:** the compute logic (correct Lane's-balance response,
  ~9 ms/update) and the pure-Python `grlp` wheel. **Not yet verified:** the
  in-browser (WASM) execution and slider responsiveness — that is what this
  standalone package is for.
- If `%pip install grlp` does not resolve, confirm the freshly built wheel is in
  `pypi/` (the auto-discovered wheelhouse) and that the JupyterLite build indexed
  it (look for `piplite:copy:whl` in the build log).
- Once tested and trusted, the route to embed this in RTD is `jupyterlite-sphinx`
  (a `{jupyterlite}`/`{voici}` directive). We are intentionally not doing that
  yet.
