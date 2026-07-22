"""
Smoke tests for the DEM-network input cleaner (preprocessing/network_preprocessing.py).

That module is a preprocessing step, not part of the installed grlp package, so we
add its directory to the path explicitly.
"""
import os
import sys
import copy

import numpy as np
import networkx as nx
import pytest

sys.path.insert(
    0, os.path.join(os.path.dirname(__file__), "..", "preprocessing")
)
import network_preprocessing as npp  # noqa: E402


def _tiny_network():
    """Two heads (1, 2) draining to an outlet (0); each edge carries x/y/z/A/s
    arrays sampled parent->child with a deliberate elevation spike to clean."""
    G = nx.DiGraph()
    G.add_node(0, z=[0.0],  s=[0.0])
    G.add_node(1, z=[20.0], s=[400.0])
    G.add_node(2, z=[25.0], s=[500.0])

    def edge(z_hi, s_hi, n=40):
        s = np.linspace(s_hi, 0.0, n)
        z = np.linspace(z_hi, 0.0, n).astype(float)
        z[n // 2] += 8.0                      # spike
        x = np.linspace(0, 1, n)
        y = np.linspace(0, 1, n)
        A = np.linspace(1, 5, n)
        return dict(x=x.tolist(), y=y.tolist(), z=z.tolist(),
                    A=A.tolist(), s=s.tolist())

    G.add_edge(1, 0, **edge(20.0, 400.0))
    G.add_edge(2, 0, **edge(25.0, 500.0))
    return G


def test_preprocess_network_monotone_and_coarsen():
    G = _tiny_network()
    n_before = len(G.edges[1, 0]["z"])

    H = npp.preprocess_network(G, node_attr="z", edge_attr="z",
                               smooth_window=7, keep_fraction=0.25)

    zc = np.asarray(H.edges[1, 0]["z"])
    # coarsened to ~keep_fraction of the samples
    assert len(zc) == max(1, round(0.25 * n_before))
    # elevation is monotone non-increasing downstream (parent higher than child)
    assert np.all(np.diff(zc) <= 1e-9)
    # the input graph is not modified
    assert len(G.edges[1, 0]["z"]) == n_before


def test_clean_edge_values_despike_smooth_in_place():
    G = _tiny_network()
    n = len(G.edges[1, 0]["z"])
    spike_i = n // 2
    spike_before = G.edges[1, 0]["z"][spike_i]   # local +8 bump above the trend

    npp.clean_edge_values(G, "z", "z", despike_window=3, despike_k=2.0,
                          alpha=0.5, n_iters=50)

    z = np.asarray(G.edges[1, 0]["z"])
    assert len(z) == n                          # length preserved
    assert z[spike_i] < spike_before - 4.0      # the local spike is knocked down
    assert z[0] == pytest.approx(20.0)          # endpoints clamped to node values
    assert z[-1] == pytest.approx(0.0)


def test_monotone_envelope_edge_clamps_a_spike():
    y = [10.0, 9.0, 20.0, 7.0, 6.0]           # spike at index 2
    out = npp.monotone_envelope_edge(y, parent_val=10.0, child_val=6.0)
    assert np.all(np.diff(out) <= 1e-9)        # non-increasing
    assert out[0] == pytest.approx(10.0)


def test_json_round_trip(tmp_path):
    G = npp.preprocess_network(_tiny_network(), node_attr="z", edge_attr="z",
                               smooth_window=7, keep_fraction=0.25)
    path = tmp_path / "net.json"
    npp.save_network_json(G, str(path))
    H = npp.load_network_json(str(path))
    assert H.number_of_edges() == G.number_of_edges()
    assert H.number_of_nodes() == G.number_of_nodes()
