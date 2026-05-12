"""Cellular-sheaf math for NMA inconsistency localization.

Lifted from generate_figures.py::sheaf_analysis and decomposed.
No I/O, no figure generation, no RNG — pure NumPy.
"""
from __future__ import annotations

import numpy as np


def build_coboundary(network: dict) -> np.ndarray:
    """Coboundary map F: edges x nodes with precision-weighted signed entries.

    For edge e = (treat1, treat2) with SE_e:
      F[e, treat1] = -1/SE_e
      F[e, treat2] = +1/SE_e
    """
    nodes = network["nodes"]
    edges = network["edges"]
    n = len(nodes)
    m = len(edges)
    idx = {nd: i for i, nd in enumerate(nodes)}
    F = np.zeros((m, n))
    for e_idx, edge in enumerate(edges):
        w = 1.0 / edge["se"]
        F[e_idx, idx[edge["treat1"]]] = -w
        F[e_idx, idx[edge["treat2"]]] = +w
    return F


def sheaf_laplacian(F: np.ndarray) -> np.ndarray:
    """L = F^T F. Symmetric PSD."""
    return F.T @ F


def _solve_node_estimates(F: np.ndarray, d: np.ndarray) -> np.ndarray:
    """WLS solve with reference node = 0 (fix node 0)."""
    n = F.shape[1]
    L = sheaf_laplacian(F)
    Ftd = F.T @ d
    if n == 1:
        return np.array([0.0])
    L_red = L[1:, 1:]
    b_red = Ftd[1:]
    try:
        x_red = np.linalg.solve(L_red, b_red)
    except np.linalg.LinAlgError as exc:
        raise ValueError(
            "Network graph is disconnected or degenerate: reduced "
            "sheaf Laplacian is singular. SheafNMA v0.2 requires a "
            "connected network with no zero/negative SE values."
        ) from exc
    x = np.zeros(n)
    x[1:] = x_red
    return x


def _precision_weighted_data(network: dict) -> np.ndarray:
    return np.array([edge["effect"] / edge["se"] for edge in network["edges"]])


def edge_residuals(network: dict) -> np.ndarray:
    """Per-edge residuals r_e = d_e - (Fx)_e where x is the WLS node-estimate."""
    F = build_coboundary(network)
    d = _precision_weighted_data(network)
    x = _solve_node_estimates(F, d)
    return d - F @ x


def gii(network: dict) -> float:
    """Global Inconsistency Index = sum of squared edge residuals."""
    r = edge_residuals(network)
    return float(np.sum(r ** 2))
