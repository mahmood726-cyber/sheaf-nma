import numpy as np
import pytest
from sheafnma.core import build_coboundary, sheaf_laplacian, gii, edge_residuals, _solve_node_estimates


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
    np.testing.assert_allclose(g, 3.0, atol=1e-9)


def test_edge_residuals_hand_verified_triangle():
    # Consistent triangle: residuals should be ~0 everywhere
    net = _triangle_network()
    r = edge_residuals(net)
    assert len(r) == 3
    np.testing.assert_allclose(r, np.zeros(3), atol=1e-9)


def test_gii_scale_invariance_under_uniform_se_doubling():
    """If every (effect, SE) is doubled together, the precision-weighted
    data d = effect/se is unchanged. F is halved (entries are 1/se), so
    L is quartered and b is halved, giving x_red doubled and F@x
    unchanged. GII = sum((d - F@x)²) must therefore equal the GII of the
    SE=1 case.

    This discriminates effect/se from effect*se: with effect*se the data
    vector would scale by 4 instead of being invariant, breaking the
    invariance.

    Reference: the SE=1 inconsistent triangle of
    test_gii_is_positive_on_inconsistent_network has effects (1,1,5) and
    GII=3.0. We double everything to (2,2,10) at SE=2 and assert GII==3.0.
    """
    net = {
        "nodes": ["A", "B", "C"],
        "edges": [
            {"treat1": "A", "treat2": "B", "effect": 2.0, "se": 2.0},
            {"treat1": "B", "treat2": "C", "effect": 2.0, "se": 2.0},
            {"treat1": "A", "treat2": "C", "effect": 10.0, "se": 2.0},
        ],
    }
    g = gii(net)
    np.testing.assert_allclose(g, 3.0, atol=1e-9)


def test_solve_raises_value_error_on_disconnected_network():
    """Disconnected network (two components: A-B and C-D, no bridge)
    should produce a ValueError with a clear message, not a raw
    LinAlgError from numpy."""
    net = {
        "nodes": ["A", "B", "C", "D"],
        "edges": [
            {"treat1": "A", "treat2": "B", "effect": 1.0, "se": 1.0},
            {"treat1": "C", "treat2": "D", "effect": 2.0, "se": 1.0},
        ],
    }
    with pytest.raises(ValueError, match="disconnected"):
        gii(net)
