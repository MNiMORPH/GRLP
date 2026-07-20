"""
The NetworkX topology graph must faithfully encode the network topology.

Phase B of the Network-primary refactor introduces ``Network.build_graph``: a
directed graph whose edges are river segments (carrying ``segment_id`` and the
``LongProfile``) and whose nodes are junctions. These tests assert the graph
agrees with the existing per-segment upstream/downstream ID lists, so that
topology queries can be migrated onto the graph without changing behavior.
"""

import networkx as nx
import pytest

from network_helpers import NETWORK_TOPOLOGIES, build_network


@pytest.fixture(scope="module", params=list(NETWORK_TOPOLOGIES),
                ids=list(NETWORK_TOPOLOGIES))
def net(request):
    spec = NETWORK_TOPOLOGIES[request.param]
    return build_network(spec["x"], spec["Q"], spec["up"], spec["down"],
                        spec["x_bl"], evolve=False)


def _edge_for(G, seg_id):
    for u, v, d in G.edges(data=True):
        if d["segment_id"] == seg_id:
            return u, v
    raise KeyError(seg_id)


def _graph_downstream_ids(G, seg_id):
    """Segment IDs from seg_id to the outlet, following out-edges."""
    _, node = _edge_for(G, seg_id)
    ids = [seg_id]
    while G.out_degree(node) > 0:
        (_, nxt, data), = list(G.out_edges(node, data=True))
        ids.append(data["segment_id"])
        node = nxt
    return ids


def _graph_upstream_ids(G, seg_id):
    """Segment IDs at or upstream of seg_id (edges within the ancestor set)."""
    u, _ = _edge_for(G, seg_id)
    upstream_nodes = nx.ancestors(G, u) | {u}
    return {seg_id} | {
        d["segment_id"]
        for a, b, d in G.edges(data=True)
        if a in upstream_nodes and b in upstream_nodes
    }


def test_graph_is_a_dag_with_single_outlet(net):
    G = net.graph
    assert nx.is_directed_acyclic_graph(G)
    outlets = [n for n in G.nodes if G.out_degree(n) == 0]
    sources = [n for n in G.nodes if G.in_degree(n) == 0]
    assert len(outlets) == 1
    assert len(sources) == len(net.list_of_channel_head_segment_IDs)


def test_one_edge_per_segment_carrying_the_object(net):
    G = net.graph
    segs = net.list_of_LongProfile_objects
    assert G.number_of_edges() == len(segs)
    seen = {}
    for _, _, d in G.edges(data=True):
        seen[d["segment_id"]] = d["segment"]
    assert set(seen) == {lp.ID for lp in segs}
    for lp in segs:
        assert seen[lp.ID] is lp


def test_graph_downstream_matches_id_lists(net):
    G = net.graph
    for lp in net.list_of_LongProfile_objects:
        assert _graph_downstream_ids(G, lp.ID) == net.find_downstream_IDs(lp.ID)


def test_graph_upstream_matches_id_lists(net):
    G = net.graph
    for lp in net.list_of_LongProfile_objects:
        assert _graph_upstream_ids(G, lp.ID) == set(net.find_upstream_IDs(lp.ID))


def test_confluence_nodes_have_expected_degree(net):
    # A confluence node (named for the outflow segment c) has in-degree equal to
    # c's number of tributaries and out-degree 1.
    G = net.graph
    for lp in net.list_of_LongProfile_objects:
        if lp.upstream_segment_IDs:
            node = ("jcn", lp.ID)
            assert G.in_degree(node) == len(lp.upstream_segment_IDs)
            assert G.out_degree(node) == 1
