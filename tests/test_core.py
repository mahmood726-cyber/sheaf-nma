import numpy as np
import pytest
from sheafnma.core import build_coboundary, sheaf_laplacian, gii, edge_residuals


def _triangle_network():
    # 3 treatments A/B/C, 3 edges with unit precision (SE=1)
    return {
        "nodes": ["A", "B", "C"],
        "edges": [
            {"treat1": "A", "treat2": "B", "effect": 1.0, "se": 1.0},
            {"treat1": "B", "treat2": "C", "effect": 1.0, "se": 1.0},
            {"treat1": "A", "treat2": "C", "effect": 2.0, "se": 1.0},  # consistent
        ],
    }


def test_coboundary_shape_and_signs():
    net = _triangle_network()
    F = build_coboundary(net)
    assert F.shape == (3, 3)
    # Edge 0 (A->B): F[0,A]=-1, F[0,B]=+1, F[0,C]=0
    assert F[0, 0] == -1.0
    assert F[0, 1] == 1.0
    assert F[0, 2] == 0.0


def test_coboundary_precision_weighting():
    net = _triangle_network()
    net["edges"][0]["se"] = 0.5  # higher precision (2x weight)
    F = build_coboundary(net)
    # weight = 1/SE = 2.0
    assert F[0, 0] == -2.0
    assert F[0, 1] == 2.0


def test_sheaf_laplacian_is_psd():
    net = _triangle_network()
    F = build_coboundary(net)
    L = sheaf_laplacian(F)
    assert L.shape == (3, 3)
    np.testing.assert_allclose(L, L.T, atol=1e-12)  # symmetric
    eigvals = np.linalg.eigvalsh(L)
    assert eigvals.min() > -1e-9  # PSD


def test_gii_is_zero_on_consistent_network():
    net = _triangle_network()  # effects A-B=1, B-C=1, A-C=2 — consistent
    g = gii(net)
    assert g < 1e-9


def test_gii_is_positive_on_inconsistent_network():
    net = _triangle_network()
    net["edges"][2]["effect"] = 5.0  # break A-C closure (should be 2)
    g = gii(net)
    assert g > 0.5


def test_edge_residuals_hand_verified_triangle():
    # Consistent triangle: residuals should be ~0 everywhere
    net = _triangle_network()
    r = edge_residuals(net)
    assert len(r) == 3
    np.testing.assert_allclose(r, np.zeros(3), atol=1e-9)
