"""Truth-recovery harness for SheafNMA's inconsistency detector.

Wires the repo's OWN engine (sheafnma.core.edge_residuals, sheafnma.comparators
.design_by_treatment_chi2 / bucher_loop_closure) against KNOWN-TRUTH networks
produced by dgp_known_truth.make_network.

Measures, vs the known truth:
  1. TYPE-I of the calibrated DBT-chi2 test on CONSISTENT networks (target ~0.05).
  2. TYPE-I of the Bucher any-loop test on CONSISTENT networks (target ~0.05).
  3. POWER of both tests as injected inconsistency delta rises.
  4. Residual-threshold detector type-I/power, using a threshold FROZEN on a
     held-out consistent+inconsistent training grid (the repo's Youden design),
     so we can compare an absolute-threshold detector against the calibrated test.
  5. OVER-DETECTION probe of the HTML dashboard's RELATIVE normalization
     (score/max*100): under a fully consistent network, what fraction of
     networks show an edge at >=50% "inconsistent"? (relative scale always
     pushes the worst edge toward 100%).
"""
from __future__ import annotations

import sys, json
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from sheafnma.core import edge_residuals  # noqa: E402
from sheafnma.comparators import design_by_treatment_chi2, bucher_loop_closure  # noqa: E402
from dgp_known_truth import make_network  # noqa: E402


def dbt_flag(net) -> bool:
    return design_by_treatment_chi2(net)["p_value"] < 0.05


def bucher_flag(net) -> bool:
    return any(b["p_value"] < 0.05 for b in bucher_loop_closure(net))


def abs_residuals(net) -> np.ndarray:
    return np.abs(edge_residuals(net))


def relative_norm_max(net) -> float:
    """HTML dashboard's normalized %: score/maxScore*100, the MAX over edges.

    Always 100 by construction (max edge maps to 100), so instead we report the
    SECOND-highest normalized score -> how 'red' a non-worst edge looks even when
    the network is fully consistent. This is the over-detection surface."""
    r = abs_residuals(net)
    sc = r ** 2
    mx = sc.max()
    if mx <= 1e-15:
        return 0.0
    norm = sc / mx * 100.0
    norm_sorted = np.sort(norm)[::-1]
    return float(norm_sorted[1]) if len(norm_sorted) > 1 else float(norm_sorted[0])


def freeze_residual_threshold(rng, n_treat=6, n_studies=5, delta=0.5, reps=300):
    """Youden-J frozen threshold on a held-out training grid (repo's design):
    consistent + single-edge-inconsistent networks. Returns absolute residual t."""
    res_list, planted_list = [], []
    for _ in range(reps):
        # half consistent, half inconsistent training reps
        d = delta if rng.random() < 0.5 else 0.0
        net = make_network(n_treat, n_studies, rng, inject_delta=d)
        r = abs_residuals(net)
        planted = None
        if net["inconsistent_pair"] is not None:
            for i, e in enumerate(net["edges"]):
                if tuple(sorted((e["treat1"], e["treat2"]))) == net["inconsistent_pair"]:
                    planted = i
                    break
        res_list.append(r)
        planted_list.append(planted)
    flat = np.concatenate(res_list)
    cands = np.linspace(0.0, float(flat.max()), 60)
    best_j, best_t = -1.0, 0.0
    for t in cands:
        tp = fp = tn = fn = 0
        for r, planted in zip(res_list, planted_list):
            for i, v in enumerate(r):
                flagged = v > t
                is_p = (planted is not None and i == planted)
                if flagged and is_p: tp += 1
                elif flagged and not is_p: fp += 1
                elif not flagged and not is_p: tn += 1
                else: fn += 1
        sens = tp / (tp + fn) if (tp + fn) else 0.0
        spec = tn / (tn + fp) if (tn + fp) else 0.0
        j = sens + spec - 1.0
        if j > best_j:
            best_j, best_t = j, float(t)
    return best_t


def measure(reps=600, n_treat=6, n_studies=5, seed=20260615):
    rng = np.random.default_rng(seed)

    thr = freeze_residual_threshold(np.random.default_rng(seed + 1),
                                    n_treat, n_studies)

    deltas = [0.0, 0.2, 0.5, 1.0]
    out = {"threshold": thr, "n_treat": n_treat, "n_studies": n_studies,
           "reps": reps, "cells": []}

    for delta in deltas:
        dbt_flags = bucher_flags = 0
        sheaf_planted_flag = sheaf_clean_fp = 0
        planted_count = clean_edge_count = 0
        overdet50 = overdet90 = 0
        for _ in range(reps):
            net = make_network(n_treat, n_studies, rng, inject_delta=delta)
            if dbt_flag(net):
                dbt_flags += 1
            if bucher_flag(net):
                bucher_flags += 1
            r = abs_residuals(net)
            planted_idx = None
            if net["inconsistent_pair"] is not None:
                for i, e in enumerate(net["edges"]):
                    if tuple(sorted((e["treat1"], e["treat2"]))) == net["inconsistent_pair"]:
                        planted_idx = i
                        break
            for i, v in enumerate(r):
                flagged = v > thr
                if i == planted_idx:
                    planted_count += 1
                    if flagged:
                        sheaf_planted_flag += 1
                else:
                    clean_edge_count += 1
                    if flagged:
                        sheaf_clean_fp += 1
            # over-detection: relative-normalized 2nd edge
            if delta == 0.0:
                rn = relative_norm_max(net)
                if rn >= 50.0:
                    overdet50 += 1
                if rn >= 90.0:
                    overdet90 += 1

        cell = {
            "delta": delta,
            "dbt_reject_rate": dbt_flags / reps,
            "bucher_reject_rate": bucher_flags / reps,
            "sheaf_planted_sensitivity": (sheaf_planted_flag / planted_count) if planted_count else None,
            "sheaf_clean_edge_fp_rate": (sheaf_clean_fp / clean_edge_count) if clean_edge_count else None,
        }
        if delta == 0.0:
            cell["html_relnorm_overdetect_ge50pct"] = overdet50 / reps
            cell["html_relnorm_overdetect_ge90pct"] = overdet90 / reps
        out["cells"].append(cell)
    return out


if __name__ == "__main__":
    res = measure()
    print(json.dumps(res, indent=2))
    outp = Path(__file__).resolve().parent / "results.json"
    outp.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(f"\nWrote {outp}")
