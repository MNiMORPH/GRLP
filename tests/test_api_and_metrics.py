"""
Previously untested public API: response-timescale helpers, network Strahler /
Horton metrics, and the synthetic-network builder.
"""

import random

import numpy as np
import pytest

import grlp
from network_helpers import build_network
from test_spectral import make_uniform_profile


# --------------------------------------------------------------------------- #
# Single-segment response-timescale helpers
# --------------------------------------------------------------------------- #

def test_equilibration_time_is_L2_over_diffusivity():
    lp = make_uniform_profile()
    lp.compute_equilibration_time()
    assert lp.equilibration_time == pytest.approx(
        lp.L ** 2 / lp.diffusivity.mean(), rel=1e-12
    )


def test_wavenumber_formula():
    lp = make_uniform_profile()
    lp.compute_equilibration_time()  # sets self.L
    for n in range(4):
        assert lp.compute_wavenumber(n) == pytest.approx(
            (2 * n + 1) * np.pi / 2.0 / lp.L, rel=1e-12
        )


def test_series_coefficient_positive_and_decreasing_in_n():
    lp = make_uniform_profile()
    lp.compute_diffusivity()
    period = lp.equilibration_time
    coeffs = [lp.compute_series_coefficient(n, period) for n in range(6)]
    assert all(c > 0 for c in coeffs)
    assert all(coeffs[i + 1] < coeffs[i] for i in range(len(coeffs) - 1))


def test_e_folding_time_formula():
    lp = make_uniform_profile()
    tau = lp.compute_e_folding_time(0)
    expected = 1.0 / lp.diffusivity.mean() / lp.compute_wavenumber(0) ** 2
    assert tau == pytest.approx(expected, rel=1e-12)


# --------------------------------------------------------------------------- #
# Network Strahler / Horton metrics on a balanced binary tree
# --------------------------------------------------------------------------- #

D = 2000.0


def _balanced_tree():
    """4 heads -> 2 mids (order 2) -> trunk (order 3); Q doubles each order."""
    def h():
        return D * np.arange(1, 5, dtype=float)

    x = [h(), h(), h(), h(),
         D * np.arange(5, 9, dtype=float),
         D * np.arange(5, 9, dtype=float),
         D * np.arange(9, 13, dtype=float)]
    Q = [5 * np.ones(4)] * 4 + [10 * np.ones(4), 10 * np.ones(4), 20 * np.ones(4)]
    up = [[], [], [], [], [0, 1], [2, 3], [4, 5]]
    down = [[4], [4], [5], [5], [6], [6], []]
    net = build_network(x, Q, up, down, D * 13, S0=0.015)
    net.compute_Q_s()
    return net


@pytest.fixture(scope="module")
def balanced_net():
    return _balanced_tree()


def test_strahler_orders_balanced_tree(balanced_net):
    balanced_net.compute_strahler_orders()
    assert list(balanced_net.segment_orders) == [1, 1, 1, 1, 2, 2, 3]


def test_horton_bifurcation_ratio(balanced_net):
    balanced_net.compute_network_properties()
    # Orders 1,2,3 have 4,2,1 streams -> bifurcation ratio 2.
    assert balanced_net.bifurcation_ratio == pytest.approx(2.0, rel=1e-9)


def test_horton_length_ratio_equal_segments(balanced_net):
    balanced_net.compute_network_properties()
    # All segments are the same length -> length ratio 1.
    assert balanced_net.length_ratio == pytest.approx(1.0, rel=1e-9)


def test_horton_discharge_ratio(balanced_net):
    balanced_net.compute_network_properties()
    # Discharge doubles each order (5 -> 10 -> 20) -> ratio 2.
    assert balanced_net.discharge_ratio == pytest.approx(2.0, rel=1e-9)


# --------------------------------------------------------------------------- #
# Synthetic-network builder (Shreve random topology)
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("magnitude", [4, 8, 16])
def test_shreve_link_counts(magnitude):
    random.seed(0)
    net = grlp.Shreve_Random_Network(magnitude=magnitude)
    up, down = net.upstream_segment_IDs, net.downstream_segment_IDs
    # A Shreve network of magnitude m has m sources, m-1 internal junctions,
    # 2m-1 links total, and a single outlet.
    assert len(up) == 2 * magnitude - 1
    assert sum(len(u) == 0 for u in up) == magnitude
    assert sum(len(d) == 0 for d in down) == 1


@pytest.mark.parametrize("magnitude", [4, 8, 16])
def test_shreve_is_valid_tree(magnitude):
    random.seed(0)
    net = grlp.Shreve_Random_Network(magnitude=magnitude)
    # Every segment drains to at most one downstream segment.
    assert all(len(d) <= 1 for d in net.downstream_segment_IDs)
    # Internal junctions are binary (exactly two tributaries).
    for u in net.upstream_segment_IDs:
        assert len(u) in (0, 2)


def test_shreve_deterministic_with_seed():
    random.seed(42)
    a = grlp.Shreve_Random_Network(magnitude=8)
    random.seed(42)
    b = grlp.Shreve_Random_Network(magnitude=8)
    assert a.upstream_segment_IDs == b.upstream_segment_IDs
    assert a.downstream_segment_IDs == b.downstream_segment_IDs
