"""Unit tests for the heterogeneity-robust null calibration (sheafnma.calibrate)."""
import numpy as np
import pytest

from sheafnma.calibrate import calibrated_inconsistency, studentized_residuals


def _consistent_triangle():
    # A-B=1, B-C=1, A-C=2 — consistent by construction.
    return {
        "nodes": ["A", "B", "C"],
        "edges": [
            {"treat1": "A", "treat2": "B", "effect": 1.0, "se": 0.1},
            {"treat1": "B", "treat2": "C", "effect": 1.0, "se": 0.1},
            {"treat1": "A", "treat2": "C", "effect": 2.0, "se": 0.1},
        ],
    }


def _square_network():
    # 4 nodes, fully connected (6 edges), node truth t_i = -0.2*i, all loops
    # consistent: effect = truth[t2] - truth[t1].
    truth = {"A": 0.0, "B": -0.2, "C": -0.4, "D": -0.6}
    pairs = [("A", "B"), ("A", "C"), ("A", "D"),
             ("B", "C"), ("B", "D"), ("C", "D")]
    edges = [{"treat1": a, "treat2": b, "effect": truth[b] - truth[a], "se": 0.1}
             for a, b in pairs]
    return {"nodes": ["A", "B", "C", "D"], "edges": edges}


def test_consistent_network_not_flagged():
    res = calibrated_inconsistency(_square_network())
    assert res["inconsistent"] is False
    assert res["p_value"] > 0.05


def test_returns_required_fields():
    res = calibrated_inconsistency(_square_network())
    for k in ("p_value", "inconsistent", "phi", "df", "max_abs_t",
              "argmax_edge", "edge_t", "edge_p"):
        assert k in res
    assert res["df"] == 6 - (4 - 1)  # m - (n-1) = 3


def _hexagon_network():
    """6 fully-connected nodes (15 edges), all loops consistent by construction.

    A bigger network than the triangle/square gives the test enough residual
    degrees of freedom (df = 15 - 5 = 10) for the calibrated test to have real
    power, mirroring the harness's n_treat=6 regime.
    """
    import string
    nodes = list(string.ascii_uppercase[:6])
    truth = {nd: -0.2 * i for i, nd in enumerate(nodes)}
    pairs = [(nodes[i], nodes[j])
             for i in range(6) for j in range(i + 1, 6)]
    edges = [{"treat1": a, "treat2": b, "effect": truth[b] - truth[a], "se": 0.1}
             for a, b in pairs]
    return {"nodes": nodes, "edges": edges}


def test_inconsistent_edge_is_flagged_and_localized():
    net = _hexagon_network()
    # break the A-C edge by a large known amount
    ac_idx = next(i for i, e in enumerate(net["edges"])
                  if (e["treat1"], e["treat2"]) == ("A", "C"))
    net["edges"][ac_idx]["effect"] += 1.0
    res = calibrated_inconsistency(net)
    assert res["inconsistent"] is True
    assert res["argmax_edge"] == ac_idx


def test_phi_tracks_overdispersion():
    """phi ~ 1 when SEs are adequate; > 1 when heterogeneity inflates scatter."""
    rng = np.random.default_rng(0)
    truth = {"A": 0.0, "B": -0.2, "C": -0.4, "D": -0.6, "E": -0.8}
    pairs = [(a, b) for i, a in enumerate("ABCDE") for b in "ABCDE"[i + 1:]]

    def build(noise_sd):
        edges = []
        for a, b in pairs:
            eff = truth[b] - truth[a] + rng.normal(0, noise_sd)
            edges.append({"treat1": a, "treat2": b, "effect": eff, "se": 0.1})
        return {"nodes": list("ABCDE"), "edges": edges}

    phi_lo = np.mean([studentized_residuals(build(0.02))["phi"] for _ in range(50)])
    phi_hi = np.mean([studentized_residuals(build(0.30))["phi"] for _ in range(50)])
    assert phi_hi > phi_lo
    assert phi_hi > 1.5  # heavy heterogeneity clearly inflates dispersion


def test_saturated_network_returns_no_detection():
    """A triangle has df = 3 - 2 = 1; still testable. A single edge (df<1)
    must not raise and must report not-inconsistent."""
    net = {"nodes": ["A", "B"],
           "edges": [{"treat1": "A", "treat2": "B", "effect": 1.0, "se": 0.1}]}
    res = calibrated_inconsistency(net)
    assert res["inconsistent"] is False
    assert res["p_value"] == 1.0
