"""
Deterministic model configurations captured by the characterization
(golden-master) tests.

This module is the single source of truth shared by
``generate_characterization_reference.py`` (which records the current outputs to
``characterization_reference.npz``) and ``test_characterization.py`` (which
re-runs these configurations and compares against that recording).

Purpose: pin the *current* numerical behavior of the model for a range of
prescribed valley-width (``B``) configurations -- uniform, power-law, arbitrary
array, with and without uplift and intermittency, plus grain-size-dependent
channel width/depth and networks -- so that forthcoming work making ``B`` a
dynamic function of the other variables can be checked for an exact regression
to the prescribed-``B`` behavior. Because the no-uplift equilibrium profile is
independent of ``B``, transient snapshots and ``B``-loaded quantities
(diffusivity, the uplift profile, arbitrary ``B(x)``) are deliberately included.

Each ``run_*`` function returns a dict of named float arrays. ``run_all``
flattens these into ``"<config>__<array>"`` keys.
"""

import numpy as np

import grlp
from network_helpers import NETWORK_TOPOLOGIES, run_topology_arrays


# Evolution schedule for single-segment configs: cumulative (nt, dt) legs.
# The first two legs leave the model mid-transient (B acts on the *rate*); the
# last drives it to steady state.
_SINGLE_SCHEDULE = [
    ("z_t0", 2, 5e11),
    ("z_t1", 6, 5e11),
    ("z_eq", 400, 1e13),
]


# --------------------------------------------------------------------------- #
# Single-segment builders
# --------------------------------------------------------------------------- #

def _build_single(S0=1.5e-2, k_xQ=1.43e-5, P_xQ=7 / 4.0 * 0.7, k_xA=1.0,
                  intermittency=0.8, U_mm_yr=0.0, nx=90, dx=1000.0, x0=10000.0,
                  B_spec=("power", 25.0, 0.2), D=None):
    """
    Build (but do not evolve) a single-segment LongProfile.

    B_spec selects the width field:
      ("power", k_xB, P_xB)   -> B = k_xB * x**P_xB
      ("uniform", value)      -> B = value (scalar, broadcast)
      ("array", callable)     -> B = callable(x)  (arbitrary B(x))
    """
    lp = grlp.LongProfile()
    lp.set_intermittency(intermittency)
    lp.basic_constants()
    lp.bedload_lumped_constants()
    lp.set_hydrologic_constants()
    lp.set_x(dx=dx, nx=nx, x0=x0)
    lp.set_z(S0=-S0, z1=0.0)
    lp.set_A(k_xA=k_xA)
    lp.set_Q(k_xQ=k_xQ, P_xQ=P_xQ)
    kind = B_spec[0]
    if kind == "power":
        lp.set_B(k_xB=B_spec[1], P_xB=B_spec[2])
    elif kind == "uniform":
        lp.set_B(B=B_spec[1])
    elif kind == "array":
        lp.set_B(B=B_spec[1](lp.x))
    else:
        raise ValueError("unknown B_spec kind: " + kind)
    lp.set_uplift_rate(U_mm_yr * 1e-3 / 3.15e7)
    lp.set_niter(3)
    lp.set_z_bl(0.0)
    Qs0 = lp.k_Qs * lp.Q[0] * S0 ** (7 / 6.0)
    lp.set_Qs_input_upstream(Qs0)
    if D is not None:
        lp.D = D
    return lp


def run_single_segment(**kwargs):
    """Evolve a single-segment config through the schedule and snapshot it."""
    lp = _build_single(**kwargs)
    out = {}
    for name, nt, dt in _SINGLE_SCHEDULE:
        lp.evolve_threshold_width_river(nt=nt, dt=dt)
        out[name] = lp.z.copy()
    # Quantities evaluated at the final (steady) state.
    lp.compute_Q_s()
    out["Q_s"] = lp.Q_s.copy()
    out["S"] = lp.S.copy()
    lp.compute_diffusivity()
    out["diffusivity"] = lp.diffusivity.copy()
    # Grain-size-dependent channel geometry, when a grain size is set.
    if getattr(lp, "D", None) is not None:
        lp.compute_channel_width()
        lp.compute_flow_depth()
        out["channel_width"] = lp.b.copy()
        out["flow_depth"] = lp.h.copy()
    return out


# --------------------------------------------------------------------------- #
# Network builders
# --------------------------------------------------------------------------- #

_NET_NT = 2000
_NET_DT = 1e9
_NSEG = 4
_DX = 2000.0


def _seg_x(start_mult):
    return _DX * np.arange(start_mult, start_mult + _NSEG, dtype=float)


def run_confluence(B_lists=None, S0=0.015, Qh=5.0):
    """
    Symmetric Y network (two heads -> trunk). B_lists, if given, is a list of
    three width arrays (one per segment); otherwise uniform B=100.
    """
    x = [_seg_x(1), _seg_x(1), _seg_x(_NSEG + 1)]
    if B_lists is None:
        B_lists = [100.0 * np.ones(_NSEG) for _ in range(3)]
    net = grlp.Network()
    net.initialize(
        x_bl=_DX * (2 * _NSEG + 1),
        z_bl=0.0,
        S0=[S0, S0],
        Q_s_0=None,
        upstream_segment_IDs=[[], [], [0, 1]],
        downstream_segment_IDs=[[2], [2], []],
        x=x,
        z=[np.zeros(_NSEG) for _ in range(3)],
        Q=[Qh * np.ones(_NSEG), Qh * np.ones(_NSEG), 2 * Qh * np.ones(_NSEG)],
        B=B_lists,
    )
    net.set_niter(3)
    net.get_z_lengths()
    net.evolve_threshold_width_river_network(nt=_NET_NT, dt=_NET_DT)
    out = {}
    out["z_all"] = np.hstack([lp.z for lp in net.list_of_LongProfile_objects])
    for lp in net.list_of_LongProfile_objects:
        S = np.abs(np.diff(lp.z) / np.diff(lp.x))
        Q_mid = (lp.Q[:-1] + lp.Q[1:]) / 2.0
        out["Qs_seg%d" % lp.ID] = lp.k_Qs * Q_mid * S ** (7 / 6.0)
    return out


def run_chain(S0=0.015, Qh=5.0):
    """Two segments in series, uniform discharge."""
    x = [_seg_x(1), _seg_x(_NSEG + 1)]
    net = grlp.Network()
    net.initialize(
        x_bl=_DX * (2 * _NSEG + 1),
        z_bl=0.0,
        S0=[S0],
        Q_s_0=None,
        upstream_segment_IDs=[[], [0]],
        downstream_segment_IDs=[[1], []],
        x=x,
        z=[np.zeros(_NSEG), np.zeros(_NSEG)],
        Q=[Qh * np.ones(_NSEG), Qh * np.ones(_NSEG)],
        B=[100.0 * np.ones(_NSEG), 100.0 * np.ones(_NSEG)],
    )
    net.set_niter(3)
    net.get_z_lengths()
    net.evolve_threshold_width_river_network(nt=_NET_NT, dt=_NET_DT)
    return {"z_all": np.hstack([lp.z for lp in net.list_of_LongProfile_objects])}


# --------------------------------------------------------------------------- #
# The catalog of configurations
# --------------------------------------------------------------------------- #

def _wavy_B(x):
    # A non-power-law, strictly positive B(x): a mean widening trend with a
    # superimposed oscillation. Representative of an arbitrary dynamic B array.
    xr = (x - x[0]) / (x[-1] - x[0])
    return 40.0 + 30.0 * xr + 10.0 * np.sin(4.0 * np.pi * xr)


# Single-segment configs: name -> kwargs for run_single_segment.
SINGLE_CONFIGS = {
    "uniform_B":            dict(B_spec=("uniform", 100.0)),
    "powerlaw_B":           dict(B_spec=("power", 25.0, 0.2)),
    "powerlaw_B_steep":     dict(B_spec=("power", 5.0, 0.5)),
    "arbitrary_B":          dict(B_spec=("array", _wavy_B)),
    "uplift_uniform_B":     dict(B_spec=("uniform", 100.0), U_mm_yr=10.0,
                                 intermittency=1.0),
    "uplift_powerlaw_B":    dict(B_spec=("power", 25.0, 0.2), U_mm_yr=10.0,
                                 intermittency=1.0),
    "uplift_arbitrary_B":   dict(B_spec=("array", _wavy_B), U_mm_yr=10.0,
                                 intermittency=1.0),
    "intermittent_B":       dict(B_spec=("power", 25.0, 0.2), intermittency=0.6),
    "grainsize_uniform_B":  dict(B_spec=("uniform", 100.0), D=0.05),
    "grainsize_arbitrary_B": dict(B_spec=("array", _wavy_B), D=0.05),
    "fine_grid_B":          dict(B_spec=("power", 25.0, 0.2), nx=180, dx=500.0),
}


def _wavy_seg(x):
    xr = (x - x[0]) / (x[-1] - x[0] + 1.0)
    return 80.0 + 40.0 * xr


# Network configs: name -> (callable, kwargs).
NETWORK_CONFIGS = {
    "confluence_uniform_B": (run_confluence, {}),
    "chain_uniform_B":      (run_chain, {}),
    "confluence_varying_B": (
        run_confluence,
        dict(B_lists=[_wavy_seg(_seg_x(1)),
                      _wavy_seg(_seg_x(1)) * 1.1,
                      _wavy_seg(_seg_x(_NSEG + 1)) + 20.0]),
    ),
}


def run_all():
    """Run every configuration; return a flat {'<config>__<array>': arr} dict."""
    data = {}
    for name, kwargs in SINGLE_CONFIGS.items():
        for key, arr in run_single_segment(**kwargs).items():
            data["%s__%s" % (name, key)] = arr
    for name, (fn, kwargs) in NETWORK_CONFIGS.items():
        for key, arr in fn(**kwargs).items():
            data["%s__%s" % (name, key)] = arr
    # Golden master for the shared correctness topologies (asymmetric, unequal
    # dx, multi-level, balanced tree). Safe to pin because Group A confirms the
    # current solver is correct on them.
    for name, spec in NETWORK_TOPOLOGIES.items():
        for key, arr in run_topology_arrays(spec).items():
            data["topology_%s__%s" % (name, key)] = arr
    return data
