"""
Interactive GRLP demo — a single gravel river adjusting in real time.

This is the source for the live, in-browser demo embedded in the documentation.
It is compiled to a standalone WebAssembly app (no server, no install) with::

    panel convert grlp_panel.py --to pyodide-worker --out <dir> \\
        --requirements <grlp wheel> numpy scipy networkx

GRLP then runs entirely in the browser via Pyodide. Press play and drag the
sliders while it runs to change the boundary conditions and watch the long
profile respond transiently — Lane's balance, played out in time.

A Jupyter-notebook version of the same demo lives alongside this file
(interactive_single_segment*.ipynb) for classroom / notebook use.
"""
import numpy as np
import panel as pn
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource

import grlp

pn.extension()

YEAR = 31556926.        # seconds per year
DT = 2 * YEAR           # 2 years advanced per animation frame (small step -> smooth)


def make_equilibrium(Qw, Qs, zbl):
    """A single segment started at steady state for the given inputs."""
    lp = grlp.LongProfile()
    lp.basic_constants()
    lp.bedload_lumped_constants()
    lp.set_hydrologic_constants()
    lp.set_x(dx=1000., nx=60, x0=1000.)
    lp.set_z(S0=-1e-2, z1=zbl)
    lp.set_Q(Qw)
    lp.set_B(100.)
    lp.set_niter(3)
    lp.set_uplift_rate(0.)
    lp.set_z_bl(zbl)
    lp.set_Qs_input_upstream(Qs)
    lp.evolve_threshold_width_river(nt=10, dt=1e13)
    return lp


Q0, QS0, ZBL0 = 100., 0.02, 0.
# `sim` (NOT `state`, which collides with panel's global pn.state) holds the
# evolving model between animation frames.
sim = {"lp": make_equilibrium(Q0, QS0, ZBL0), "t": 0.}
_lp = sim["lp"]


def _bl_xy(zbl):
    return {"x": [_lp.x.min() / 1000., _lp.x.max() / 1000.], "z": [zbl, zbl]}


profile = ColumnDataSource(data={"x": _lp.x / 1000., "z": _lp.z})
baselevel = ColumnDataSource(data=_bl_xy(ZBL0))

# The plot fills whatever column it is embedded in -- the documentation page, a
# course page, a projected slide -- while holding its shape. `scale_width`
# scales height with width, so ASPECT_RATIO fixes the vertical exaggeration:
# a reader on an ultrawide monitor and one on a laptop see the same river at
# the same apparent steepness. `stretch_width` would not do this; it pins the
# height, so a wider window silently flattens the profile.
#
# Growth is deliberately unbounded: the figure takes whatever width it is
# given. Holding the ratio means the height follows, so on a very wide screen
# it becomes tall -- about 1340 px in a 2400 px container. That is the price of
# a constant vertical exaggeration, and it is the right way round for a figure
# people read slopes off. Set `max_width` on the figure to cap it.
ASPECT_RATIO = 680. / 380.    # the proportions this figure was designed at

fig = figure(height=380, sizing_mode="scale_width",
             aspect_ratio=ASPECT_RATIO,
             title="t = 0.0 kyr",
             x_axis_label="Downstream distance [km]",
             y_axis_label="Elevation [m]")
fig.line("x", "z", source=profile, line_width=3)
fig.line("x", "z", source=baselevel, line_width=1, line_dash="dashed",
         color="gray", legend_label="base level")
fig.y_range.start, fig.y_range.end = -120, 1300
fig.legend.location = "top_right"

Qw = pn.widgets.FloatSlider(name="Water discharge  Q  [m³/s]",
                            start=20, end=600, step=20, value=Q0)
Qs = pn.widgets.FloatSlider(name="Bed-load sediment input  Qₛ  [m³/s]",
                            start=0.005, end=0.06, step=0.005, value=QS0,
                            format="0.000")
zbl = pn.widgets.FloatSlider(name="Base level  [m]",
                             start=-100, end=100, step=5, value=ZBL0)

run = pn.widgets.Toggle(name="▶ Run", value=False)
reset = pn.widgets.Button(name="Set to equilibrium", button_type="primary")


def step(event=None):
    """Advance one frame, reading the sliders as live boundary conditions."""
    lp = sim["lp"]
    lp.set_Q(Qw.value)
    lp.set_Qs_input_upstream(Qs.value)
    lp.set_z_bl(zbl.value)
    lp.evolve_threshold_width_river(nt=1, dt=DT)
    sim["t"] += DT
    profile.data = {"x": lp.x / 1000., "z": lp.z}
    baselevel.data = _bl_xy(zbl.value)
    fig.title.text = "t = %.1f kyr" % (sim["t"] / (1000. * YEAR))


def do_reset(event=None):
    """Restart from equilibrium for the current slider settings."""
    sim["lp"] = make_equilibrium(Qw.value, Qs.value, zbl.value)
    sim["t"] = 0.
    lp = sim["lp"]
    profile.data = {"x": lp.x / 1000., "z": lp.z}
    baselevel.data = _bl_xy(zbl.value)
    fig.title.text = "t = 0.0 kyr"


# A single play/pause button drives the animation via a periodic callback (the
# timer only runs while the toggle is on).
_ticker = pn.state.add_periodic_callback(step, period=33, start=False)  # ~30 fps


def toggle_run(event):
    if event.new:
        _ticker.start()
        run.name = "⏸ Pause"
    else:
        _ticker.stop()
        run.name = "▶ Run"


run.param.watch(toggle_run, "value")
reset.on_click(do_reset)

pn.Column(
    pn.pane.Markdown(
        "### Watch a gravel river adjust\n"
        "Press **▶** to run, then drag the sliders while it plays: more "
        "**sediment** aggrades and steepens the profile; more **water** lowers "
        "its slope; dropping **base level** sends an incision wave upstream. "
        "**Set to equilibrium** jumps to the steady state for the current "
        "settings."),
    pn.Row(run, reset),
    Qw, Qs, zbl,
    fig,
    sizing_mode="stretch_width",
).servable(title="GRLP interactive demo")
