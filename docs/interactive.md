# Interactive demo

Try GRLP live in your browser — no installation required. The model runs entirely
on your machine via [Pyodide](https://pyodide.org) (Python compiled to
WebAssembly); nothing is sent to a server.

This demo evolves a **single-segment** gravel river through time. Press **▶ Run**,
then drag the sliders *while the simulation runs* and watch the long profile
respond:

- **Water discharge $Q$** — more water lowers the equilibrium slope (gentler
  profile).
- **Bed-load sediment input $Q_s$** — more sediment steepens it.
- **Base level** — raising or lowering the downstream boundary aggrades or
  incises the profile.

These are the three controls of Lane's balance, played out dynamically. Use
**Set to equilibrium** to jump to the steady state for the current settings.

```{note}
The demo takes ~10–30 s to start the first time while your browser downloads the
Python runtime, then runs smoothly. Best viewed on a full-width screen.
```

<iframe src="_static/interactive/grlp_panel.html"
        width="100%" height="760" style="border: none;"
        title="GRLP interactive demo"></iframe>

The same demo is available as a Jupyter notebook (for classroom or notebook use)
in the project's `interactive_demo/` directory.
