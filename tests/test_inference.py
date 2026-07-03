import pytest

from sheafnma.core import gii
from sheafnma.inference import gii_permutation_test


def _triangle(effect_ac=2.0):
    return {
        "nodes": ["A", "B", "C"],
        "edges": [
            {"treat1": "A", "treat2": "B", "effect": 1.0, "se": 1.0},
            {"treat1": "B", "treat2": "C", "effect": 1.0, "se": 1.0},
            {"treat1": "A", "treat2": "C", "effect": effect_ac, "se": 1.0},
        ],
    }


def test_permutation_returns_expected_keys():
    res = gii_permutation_test(_triangle(), n_perm=50, seed=1)
    for k in ("gii_observed", "p_value", "n_perm", "null_mean", "null_q95"):
        assert k in res
    assert res["n_perm"] == 50


def test_observed_gii_matches_core():
    net = _triangle(effect_ac=5.0)
    res = gii_permutation_test(net, n_perm=20, seed=1)
    assert res["gii_observed"] == pytest.approx(gii(net))


def test_permutation_is_deterministic_for_fixed_seed():
    a = gii_permutation_test(_triangle(5.0), n_perm=100, seed=7)
    b = gii_permutation_test(_triangle(5.0), n_perm=100, seed=7)
    assert a == b


def test_p_value_in_open_unit_interval():
    res = gii_permutation_test(_triangle(5.0), n_perm=200, seed=3)
    # add-one correction => strictly in (0, 1]
    assert 0.0 < res["p_value"] <= 1.0


def test_permutation_validates_disconnected():
    disc = {
        "nodes": ["A", "B", "C", "D"],
        "edges": [
            {"treat1": "A", "treat2": "B", "effect": 1.0, "se": 1.0},
            {"treat1": "C", "treat2": "D", "effect": 2.0, "se": 1.0},
        ],
    }
    with pytest.raises(ValueError, match="disconnected"):
        gii_permutation_test(disc, n_perm=10)


def test_permutation_does_not_mutate_input():
    net = _triangle(5.0)
    effects_before = [e["effect"] for e in net["edges"]]
    gii_permutation_test(net, n_perm=30, seed=2)
    effects_after = [e["effect"] for e in net["edges"]]
    assert effects_before == effects_after


def test_bad_n_perm_raises():
    with pytest.raises(ValueError):
        gii_permutation_test(_triangle(), n_perm=0)
