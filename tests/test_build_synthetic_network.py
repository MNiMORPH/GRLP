"""
Smoke tests for the shipped synthetic-network generator (build_synthetic_network).

This module is exported from the package (grlp.generate_random_network,
grlp.Shreve_Random_Network) and is a natural "install and build-and-run" entry
point, but it had no test and silently broke on v3 when the de-pad removed the
padded arrays. These guard the generate -> compute -> evolve -> plot path.
"""
import random

import numpy as np
import pytest

import grlp


def _make_network(seed=12345, magnitude=6):
    random.seed(seed)
    np.random.seed(seed)
    net, topo = grlp.generate_random_network(
        magnitude=magnitude, max_length=1.5e4, mean_discharge=10.)
    return net, topo


def test_generate_random_network_runs_and_evolves():
    net, topo = _make_network()
    segs = net.list_of_LongProfile_objects
    assert len(segs) > 1                                  # a branching network
    assert isinstance(net, grlp.Network)

    # compute_Q_s over the network sets S and Q_s on every segment
    net.compute_Q_s()
    assert np.all(np.isfinite(segs[0].Q_s))

    # it evolves on the v3 walking solver
    z0_before = segs[0].z.copy()
    net.set_niter(3)
    net.get_z_lengths()
    net.evolve_threshold_width_river_network(nt=10, dt=3.15e11)
    assert np.all(np.isfinite(segs[0].z))
    assert not np.array_equal(segs[0].z, z0_before)       # something happened


def test_network_plot_returns_planform():
    net, topo = _make_network()
    planform = net.plot(show=False)
    assert len(planform) == len(net.list_of_LongProfile_objects)
    # each segment has matching-length x/y planform coordinates
    for i, seg in enumerate(net.list_of_LongProfile_objects):
        assert len(planform[i]["x"]) == len(planform[i]["y"])
    assert net.max_topological_length >= 1


def test_find_hack_parameters():
    net, topo = _make_network(magnitude=8)
    out = net.find_hack_parameters()
    assert set(("k", "p", "d", "Q")).issubset(out.keys())
    assert np.isfinite(out["k"]) and np.isfinite(out["p"])
    # non-dimensional variant (fixed net.sources -> channel-head IDs)
    nd = net.find_hack_parameters_non_dim()
    assert set(("k_src", "p_src", "k_seg", "p_seg")).issubset(nd.keys())
    assert np.isfinite(nd["k_src"]) and np.isfinite(nd["p_src"])
