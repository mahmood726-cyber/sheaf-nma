import numpy as np
from sheafnma.simulate import (
    xoshiro128ss, normal_random,
    generate_planted_network,
    splitmix32,
)
from sheafnma.core import gii


def _ref_splitmix32(seed):
    """Reference splitmix32 copied verbatim from generate_figures.py lines 27-38."""
    a = seed & 0xFFFFFFFF
    def next_val():
        nonlocal a
        a = (a + 0x9E3779B9) & 0xFFFFFFFF
        t = a ^ (a >> 16)
        t = (t * 0x21F0AAAD) & 0xFFFFFFFF
        t = t ^ (t >> 15)
        t = (t * 0x735A2D97) & 0xFFFFFFFF
        return (t ^ (t >> 15)) & 0xFFFFFFFF
    return next_val


def _ref_xoshiro128ss(seed):
    """Reference xoshiro128ss copied verbatim from generate_figures.py lines 41-62."""
    sm = _ref_splitmix32(seed)
    s = [sm(), sm(), sm(), sm()]

    def next_val():
        result = ((s[1] * 5) & 0xFFFFFFFF)
        result = (((result << 7) | (result >> 25)) & 0xFFFFFFFF)
        result = (result * 9) & 0xFFFFFFFF
        t = (s[1] << 9) & 0xFFFFFFFF

        s[2] ^= s[0]
        s[3] ^= s[1]
        s[1] ^= s[2]
        s[0] ^= s[3]

        s[2] ^= t
        s[3] = ((s[3] << 11) | (s[3] >> 21)) & 0xFFFFFFFF

        return result / 4294967296.0

    return next_val


def test_xoshiro_matches_generate_figures_byte_identical():
    """xoshiro128ss(42) in simulate.py produces the SAME stream as
    generate_figures.py's xoshiro128ss(42).

    This is the byte-identicality guarantee for the v0.1 -> v0.2
    refactor: Task 5 reduces generate_figures.py to a thin wrapper
    around sheafnma.simulate, and the Phase-1 figures must not shift.
    """
    r_new = xoshiro128ss(42)
    r_old = _ref_xoshiro128ss(42)
    new_vals = [r_new() for _ in range(20)]
    old_vals = [r_old() for _ in range(20)]
    assert new_vals == old_vals, (
        f"xoshiro divergence at first "
        f"{next(i for i, (a, b) in enumerate(zip(new_vals, old_vals)) if a != b)}"
        f" samples"
    )


def test_splitmix32_matches_generate_figures_byte_identical():
    """splitmix32 in simulate.py produces byte-identical output to generate_figures.py."""
    sm_new = splitmix32(99)
    sm_old = _ref_splitmix32(99)
    new_vals = [sm_new() for _ in range(20)]
    old_vals = [sm_old() for _ in range(20)]
    assert new_vals == old_vals, (
        f"splitmix32 divergence at index "
        f"{next(i for i, (a, b) in enumerate(zip(new_vals, old_vals)) if a != b)}"
    )


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
