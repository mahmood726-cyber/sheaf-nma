"""Truth-recovery assertions for SheafNMA, vs known-truth networks.

Run: python -m pytest truth-recovery/test_truth_recovery.py -q
 or: python truth-recovery/test_truth_recovery.py

These encode what was MEASURED against the repo's own engine:
  (A) On consistent, low-noise networks the calibrated tests do not false-alarm.
  (B) Power rises with injected inconsistency (POWER is real).
  (C) BUG: the fixed-effect loop-closure tests AND the un-calibrated GII
      over-detect under within-edge heterogeneity (tau2 ignored) -> type-I
      far above 0.05 on networks that are CONSISTENT by construction.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[0]))
sys.path.insert(0, str(HERE))

from sheafnma.core import gii  # noqa: E402
from sheafnma.comparators import design_by_treatment_chi2, bucher_loop_closure  # noqa: E402
from dgp_known_truth import make_network  # noqa: E402


def _dbt_rate(noise, delta, N=600, seed=3):
    rng = np.random.default_rng(seed)
    k = 0
    for _ in range(N):
        net = make_network(6, 5, rng, inject_delta=delta, noise_sd=noise)
        if design_by_treatment_chi2(net)["p_value"] < 0.05:
            k += 1
    return k / N


def test_no_false_alarm_low_noise_consistent():
    """(A) Consistent + low within-edge noise -> DBT does not false-alarm."""
    rate = _dbt_rate(noise=0.05, delta=0.0)
    assert rate <= 0.05, f"DBT type-I {rate} should be <=0.05 at low noise"


def test_power_rises_with_injected_inconsistency():
    """(B) POWER: reject rate strictly increases as the known delta grows."""
    r02 = _dbt_rate(0.05, 0.2)
    r05 = _dbt_rate(0.05, 0.5)
    r10 = _dbt_rate(0.05, 1.0)
    assert r02 < r05 < r10, f"power not monotone: {r02},{r05},{r10}"
    assert r10 >= 0.95, f"power at delta=1.0 should be high, got {r10}"


def test_BUG_dbt_overdetects_under_heterogeneity():
    """(C) BUG: with unmodeled within-edge heterogeneity the fixed-effect DBT
    test over-rejects on CONSISTENT networks (type-I >> 0.05)."""
    rate = _dbt_rate(noise=0.15, delta=0.0)
    assert rate > 0.10, (
        f"expected over-detection from FE/no-tau2 SE; measured type-I={rate}")


def test_BUG_gii_uncalibrated_redflags_consistent_networks():
    """(C) BUG: the GII has no null reference; its hardcoded RED cutoff (>=2)
    fires on CONSISTENT-but-heterogeneous networks."""
    rng = np.random.default_rng(5)
    red = 0
    N = 500
    for _ in range(N):
        net = make_network(6, 5, rng, inject_delta=0.0, noise_sd=0.15)
        if gii(net) >= 2.0:
            red += 1
    frac = red / N
    assert frac > 0.5, (
        f"expected GII to over-flag consistent heterogeneous nets; got {frac}")


def test_planted_edge_is_largest_residual_when_clearly_inconsistent():
    """Sanity: when a large delta is injected, the planted edge carries the
    largest absolute residual (localization works in the easy regime)."""
    from sheafnma.core import edge_residuals
    rng = np.random.default_rng(9)
    hits = 0
    N = 200
    for _ in range(N):
        net = make_network(6, 5, rng, inject_delta=1.0, noise_sd=0.05)
        r = np.abs(edge_residuals(net))
        planted = None
        for i, e in enumerate(net["edges"]):
            if tuple(sorted((e["treat1"], e["treat2"]))) == net["inconsistent_pair"]:
                planted = i
                break
        if int(np.argmax(r)) == planted:
            hits += 1
    assert hits / N >= 0.85, f"localization weak: {hits/N}"


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    passed = 0
    for fn in fns:
        fn()
        print(f"PASS {fn.__name__}")
        passed += 1
    print(f"\n{passed}/{len(fns)} assertions passed")
