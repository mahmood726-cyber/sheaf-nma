"""Heterogeneity-robust null calibration for SheafNMA inconsistency detection.

WHY THIS EXISTS
---------------
The raw Global Inconsistency Index ``gii`` (``core.gii``) and the fixed-effect
loop-closure comparators (``comparators.bucher_loop_closure`` /
``design_by_treatment_chi2``) judge inconsistency against a null that assumes
the supplied edge SEs are the *whole* story. They are not: real NMA contrasts
carry between-study heterogeneity (tau^2) that an inverse-variance pooled SE
does not capture. With unmodeled heterogeneity the precision-weighted residuals
are over-dispersed, so the raw GII and the FE z-statistics are inflated and the
detector OVER-DETECTS — it flags consistent-but-heterogeneous networks as
inconsistent far above the nominal alpha (measured type-I 0.16-0.91 at
noise_sd 0.15-0.30 vs target 0.05; the HTML's absolute GII>=2 "RED" cutoff fired
on 99.7% of consistent heterogeneous networks).

THE FIX (same family as the meta-entropy null-calibration fix)
--------------------------------------------------------------
Judge the statistic against its *proper* null instead of an absolute threshold:

1. Fit the consistent model by WLS (reference node fixed); take the
   precision-weighted residuals r = d - F x.
2. Studentize each residual by its leverage: the hat matrix H projects onto the
   fitted (consistent) subspace, so under the null
   r_e ~ N(0, phi * (1 - h_ee)) and rr_e = r_e / sqrt(1 - h_ee) is homoscedastic.
3. Estimate the *multiplicative over-dispersion* phi robustly from the studentized
   residuals. phi absorbs the unmodeled tau^2: phi ~ 1 when SEs are adequate, and
   phi > 1 exactly when heterogeneity inflates the scatter. The estimator is
   iterated with outlier trimming so a single genuinely-inconsistent edge cannot
   inflate phi and mask itself (which would destroy power).
4. Compare the maximum standardized residual against its calibrated null
   (standard normal, Bonferroni over the m edges). The returned ``p_value`` is the
   calibrated inconsistency p-value; ``edge_p`` localizes it.

This makes the test heterogeneity-INVARIANT: type-I stays at nominal alpha
whether noise_sd is 0.05, 0.15 or 0.30, while power under a genuine planted
inconsistency is preserved and monotone.

The module is ADDITIVE: ``core.gii`` and the comparators are unchanged so all
existing exact-value tests still hold. Callers that want a yes/no decision under
heterogeneity should use :func:`calibrated_inconsistency` here.
"""
from __future__ import annotations

import numpy as np
from scipy import stats as sp_stats

from .core import build_coboundary, _solve_node_estimates, _precision_weighted_data

# E[z^2 | z^2 < 9] for z ~ N(0,1). Used to de-bias the trimmed mean of squared
# studentized residuals back to phi (the trim at z^2 < 9*phi removes ~0.27% of a
# clean normal's mass, so the surviving mean must be rescaled up by this factor).
# Precomputed via scipy.stats.chi2.expect(lambda x: x, lb=0, ub=9, args=(1,)) / chi2.cdf(9,1).
_TRUNC_MEAN_9 = 0.9733369246698526
# Median of a chi^2_1 deviate (= median of z^2). median(rr^2)/this is a robust
# first-pass scale estimate for phi.
_CHI2_1_MEDIAN = 0.4549364231195736


def _hat_leverages(F: np.ndarray) -> np.ndarray:
    """Diagonal of the hat matrix H = F_free (F_free^T F_free)^-1 F_free^T.

    F_free drops the reference-node column (node 0 is fixed at 0 in the WLS
    solve), matching ``core._solve_node_estimates``. h_ee in [0, 1) is the
    leverage of edge e on its own fitted value.
    """
    F_free = F[:, 1:]
    G = F_free.T @ F_free
    # Solve rather than invert; G is small (n-1 square) and PD for connected nets.
    H = F_free @ np.linalg.solve(G, F_free.T)
    h = np.clip(np.diag(H), 0.0, 1.0 - 1e-9)
    return h


def _robust_phi(rr2: np.ndarray) -> float:
    """Robust multiplicative over-dispersion phi from squared studentized resids.

    rr2 = rr_e^2 where rr_e ~ N(0, phi) under the consistent null. We want a phi
    estimate that (a) tracks heterogeneity-induced over-dispersion and (b) is NOT
    inflated by one or two genuinely inconsistent edges (which would mask them).

    Start from the median (highly outlier-resistant), then iterate a trimmed mean:
    keep rr2 < 9*phi (~3 sigma in z^2 units) and rescale by the truncated-normal
    correction so phi remains unbiased for a clean network.
    """
    rr2 = np.asarray(rr2, dtype=float)
    phi = max(float(np.median(rr2)) / _CHI2_1_MEDIAN, 1e-12)
    for _ in range(5):
        keep = rr2 < 9.0 * phi
        if int(keep.sum()) < 2:
            break
        new = float(np.mean(rr2[keep])) / _TRUNC_MEAN_9
        new = max(new, 1e-12)
        if abs(new - phi) <= 1e-9 * phi:
            phi = new
            break
        phi = new
    return phi


def studentized_residuals(network: dict) -> dict:
    """Heterogeneity-robust studentized edge residuals under the consistent model.

    Returns a dict with:
      ``residuals``  precision-weighted residuals r_e = d_e - (Fx)_e
      ``leverages``  hat-matrix diagonal h_ee
      ``phi``        robust multiplicative over-dispersion (absorbs unmodeled tau^2)
      ``t``          standardized residuals rr_e / sqrt(phi) (~N(0,1) under null)
      ``df``         residual degrees of freedom m - (n - 1)
    """
    F = build_coboundary(network)
    d = _precision_weighted_data(network)
    x = _solve_node_estimates(F, d)
    r = d - F @ x
    m, n = F.shape
    df = m - (n - 1)
    h = _hat_leverages(F)
    rr = r / np.sqrt(1.0 - h)            # homoscedastic studentized residual
    phi = _robust_phi(rr ** 2)
    t = rr / np.sqrt(phi)
    return {
        "residuals": r,
        "leverages": h,
        "phi": float(phi),
        "t": t,
        "df": int(df),
    }


def calibrated_inconsistency(network: dict, alpha: float = 0.05) -> dict:
    """Calibrated network-level inconsistency test, robust to between-study tau^2.

    Statistic: the maximum standardized residual max_e |t_e|, judged against a
    standard-normal null with a Bonferroni correction over the m edges. Because
    the over-dispersion phi is estimated from the data, the test calibrates to the
    nominal alpha whether or not heterogeneity is present.

    Returns dict:
      ``p_value``        Bonferroni network p-value (min over edges, capped at 1).
      ``inconsistent``   bool: p_value < alpha.
      ``phi``            estimated multiplicative over-dispersion.
      ``df``             residual degrees of freedom.
      ``max_abs_t``      max |t_e| (the localizing statistic).
      ``argmax_edge``    index of the most-inconsistent edge (localization).
      ``edge_t``         per-edge standardized residuals t_e.
      ``edge_p``         per-edge two-sided Bonferroni p-values.
    """
    s = studentized_residuals(network)
    t = s["t"]
    m = len(t)
    if s["df"] < 1 or m == 0:
        # Saturated / degenerate network: no residual room to test inconsistency.
        return {
            "p_value": 1.0,
            "inconsistent": False,
            "phi": s["phi"],
            "df": s["df"],
            "max_abs_t": float(np.max(np.abs(t))) if m else 0.0,
            "argmax_edge": int(np.argmax(np.abs(t))) if m else None,
            "edge_t": t,
            "edge_p": np.ones(m),
        }
    edge_p = np.clip(2.0 * sp_stats.norm.sf(np.abs(t)) * m, 0.0, 1.0)
    p_value = float(edge_p.min())
    return {
        "p_value": p_value,
        "inconsistent": bool(p_value < alpha),
        "phi": s["phi"],
        "df": s["df"],
        "max_abs_t": float(np.max(np.abs(t))),
        "argmax_edge": int(np.argmax(np.abs(t))),
        "edge_t": t,
        "edge_p": edge_p,
    }
