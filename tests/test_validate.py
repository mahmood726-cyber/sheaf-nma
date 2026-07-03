import pytest

from sheafnma.validate import (
    NetworkValidationError,
    connected_components,
    is_connected,
    validate_network,
)


def _triangle():
    return {
        "nodes": ["A", "B", "C"],
        "edges": [
            {"treat1": "A", "treat2": "B", "effect": 1.0, "se": 1.0},
            {"treat1": "B", "treat2": "C", "effect": 1.0, "se": 1.0},
            {"treat1": "A", "treat2": "C", "effect": 2.0, "se": 1.0},
        ],
    }


def _disconnected():
    return {
        "nodes": ["A", "B", "C", "D"],
        "edges": [
            {"treat1": "A", "treat2": "B", "effect": 1.0, "se": 1.0},
            {"treat1": "C", "treat2": "D", "effect": 2.0, "se": 1.0},
        ],
    }


def test_connected_triangle_is_one_component():
    comps = connected_components(_triangle())
    assert comps == [["A", "B", "C"]]
    assert is_connected(_triangle()) is True


def test_disconnected_two_components():
    comps = connected_components(_disconnected())
    assert comps == [["A", "B"], ["C", "D"]]
    assert is_connected(_disconnected()) is False


def test_validate_passes_connected_and_returns_same_object():
    net = _triangle()
    out = validate_network(net)
    assert out is net  # returned unchanged for inline use


def test_validate_rejects_disconnected():
    with pytest.raises(NetworkValidationError, match="disconnected"):
        validate_network(_disconnected())


def test_validate_can_skip_connectivity():
    # structural checks still run, but disconnected is allowed
    out = validate_network(_disconnected(), require_connected=False)
    assert out is not None


def test_validate_rejects_self_loop():
    net = _triangle()
    net["edges"][0]["treat2"] = "A"  # A vs A
    with pytest.raises(NetworkValidationError, match="self-loop"):
        validate_network(net)


def test_validate_rejects_nonpositive_se():
    net = _triangle()
    net["edges"][0]["se"] = 0.0
    with pytest.raises(NetworkValidationError, match="non-positive|se"):
        validate_network(net)


def test_validate_rejects_missing_key():
    net = _triangle()
    del net["edges"][1]["effect"]
    with pytest.raises(NetworkValidationError, match="missing required key"):
        validate_network(net)


def test_validate_rejects_empty_edges():
    with pytest.raises(NetworkValidationError, match="zero edges"):
        validate_network({"nodes": [], "edges": []})


def test_validate_rejects_non_finite_effect():
    net = _triangle()
    net["edges"][0]["effect"] = float("inf")
    with pytest.raises(NetworkValidationError, match="non-finite"):
        validate_network(net)


def test_isolated_declared_node_counts_as_component():
    net = _triangle()
    net["nodes"].append("Z")  # declared but never in any edge
    comps = connected_components(net)
    assert ["Z"] in comps
    assert is_connected(net) is False
