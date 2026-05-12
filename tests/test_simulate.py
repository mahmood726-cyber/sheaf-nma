import numpy as np
from sheafnma.simulate import (
    xoshiro128ss, normal_random,
    generate_planted_network,
)
from sheafnma.core import gii


def test_xoshiro_is_deterministic():
    r1 = xoshiro128ss(42)
    r2 = xoshiro128ss(42)
    vals1 = [r1() for _ in range(10)]
    vals2 = [r2() for _ in range(10)]
    assert vals1 == vals2


def test_xoshiro_seeds_differ():
    r1 = xoshiro128ss(42)
    r2 = xoshiro128ss(43)
    assert [r1() for _ in range(5)] != [r2() for _ in range(5)]


def test_normal_random_within_bounds_99pct():
    r = xoshiro128ss(7)
    vals = [normal_random(r) for _ in range(1000)]
    assert -5 < min(vals) and max(vals) < 5
    assert abs(np.mean(vals)) < 0.2


def test_planted_zero_tau_gives_consistent_network():
    net = generate_planted_network(
        n_treatments=4, n_studies_per_edge=3, tau_inc=0.0, seed=1,
    )
    assert gii(net) < 0.5  # noise floor


def test_planted_high_tau_gives_inconsistent_network():
    net = generate_planted_network(
        n_treatments=4, n_studies_per_edge=3, tau_inc=1.0, seed=1,
    )
    assert gii(net) > 2.0


def test_planted_returns_known_inconsistent_edge():
    net = generate_planted_network(
        n_treatments=4, n_studies_per_edge=3, tau_inc=1.0, seed=1,
    )
    assert "planted_edge_idx" in net
    assert 0 <= net["planted_edge_idx"] < len(net["edges"])
