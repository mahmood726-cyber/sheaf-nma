"""Power-simulation grid for SheafNMA vs DBT-χ² vs Bucher.

Grid (defaults): n_treatments ∈ {4,6,8} × n_studies_per_edge ∈ {3,5,10}
× τ_inc ∈ {0,0.2,0.5,1.0}. 1000 reps per cell.

Threshold for SheafNMA per-edge flagging:
  - Training cell: (n_treat=6, n_studies=5, τ_inc=0.5)
  - Frozen at Youden-maximising threshold over candidate set
    np.linspace(0, max_residual, 50)
  - Evaluation: remaining 35 cells, threshold reused.

Outputs:
  analysis/results/power_results.csv   per-(cell,method) sensitivity/specificity
  analysis/results/threshold.json      the frozen threshold + Youden-J trace

Determinism: each (cell, rep) uses seed = seed_base*1_000_003 +
cell_idx*10_007 + rep. test_reproducibility.py enforces byte-identical
re-runs.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from sheafnma.core import edge_residuals
from sheafnma.simulate import generate_planted_network
from sheafnma.comparators import bucher_loop_closure, design_by_treatment_chi2


def _cell_id(n_treat: int, n_studies: int, tau: float) -> str:
    return f"t{n_treat}_s{n_studies}_tau{tau:.2f}"


def _run_one_rep(n_treat: int, n_studies: int, tau: float, seed: int) -> dict:
    net = generate_planted_network(
        n_treatments=n_treat,
        n_studies_per_edge=n_studies,
        tau_inc=tau,
        seed=seed,
    )
    residuals = np.abs(edge_residuals(net))
    planted = net["planted_edge_idx"]
    dbt = design_by_treatment_chi2(net)
    bucher = bucher_loop_closure(net)
    return {
        "residuals": residuals,
        "planted_idx": planted,
        "n_edges": len(net["edges"]),
        "dbt_p": dbt["p_value"],
        "bucher_any_flag": int(any(b["p_value"] < 0.05 for b in bucher)),
    }


def _seed(seed_base: int, cell_idx: int, rep: int) -> int:
    return (seed_base * 1_000_003 + cell_idx * 10_007 + rep) % (2**31 - 1)


def _select_threshold(
    reps: int,
    residuals_all: list[np.ndarray],
    planted_all: list[int],
) -> tuple[float, list[dict]]:
    """Youden on the training cell: J = sensitivity + specificity - 1."""
    flat_r = np.concatenate(residuals_all) if residuals_all else np.array([1.0])
    max_val = float(flat_r.max()) if flat_r.size else 1.0
    candidates = np.linspace(0.0, max_val, 50)
    trace = []
    best_j = -1.0
    best_t = 0.0
    for t in candidates:
        tp = fp = tn = fn = 0
        for r, planted in zip(residuals_all, planted_all):
            for i, v in enumerate(r):
                flagged = bool(v > t)
                is_planted = (i == planted)
                if flagged and is_planted:
                    tp += 1
                elif flagged and not is_planted:
                    fp += 1
                elif not flagged and not is_planted:
                    tn += 1
                elif not flagged and is_planted:
                    fn += 1
        sens = tp / (tp + fn) if (tp + fn) else 0.0
        spec = tn / (tn + fp) if (tn + fp) else 0.0
        j = sens + spec - 1.0
        trace.append({
            "threshold": float(t),
            "sensitivity": float(sens),
            "specificity": float(spec),
            "youden_j": float(j),
        })
        if j > best_j:
            best_j = j
            best_t = float(t)
    return best_t, trace


def main() -> int:
    ap = argparse.ArgumentParser(
        description="36-cell power simulation: SheafNMA vs DBT-χ² vs Bucher"
    )
    ap.add_argument("--treatments", default="4,6,8",
                    help="Comma-separated list of n_treatments values")
    ap.add_argument("--studies", default="3,5,10",
                    help="Comma-separated list of n_studies_per_edge values")
    ap.add_argument("--taus", default="0,0.2,0.5,1.0",
                    help="Comma-separated list of tau_inc values")
    ap.add_argument("--reps", type=int, default=1000,
                    help="Number of Monte-Carlo replications per cell")
    ap.add_argument("--seed", type=int, default=42,
                    help="Base random seed (deterministic)")
    ap.add_argument("--out", default="analysis/results/power_results.csv",
                    help="Output CSV path")
    ap.add_argument("--threshold-out", default="analysis/results/threshold.json",
                    help="Output JSON path for frozen threshold + Youden trace")
    args = ap.parse_args()

    treats = [int(x) for x in args.treatments.split(",")]
    studies_list = [int(x) for x in args.studies.split(",")]
    taus = [float(x) for x in args.taus.split(",")]

    # Build ordered cell list — order is deterministic
    cells: list[tuple[int, int, int, float]] = []
    cell_idx = 0
    for nt in treats:
        for ns in studies_list:
            for tau in taus:
                cells.append((cell_idx, nt, ns, tau))
                cell_idx += 1

    total_cells = len(cells)
    print(f"Grid: {total_cells} cells × {args.reps} reps")

    # Training cell: (n_treat=6, n_studies=5, τ_inc=0.5) — held-out threshold selection
    training_signature = (6, 5, 0.5)
    training = next(
        (c for c in cells if (c[1], c[2], c[3]) == training_signature),
        cells[total_cells // 2],
    )
    print(
        f"Training cell: idx={training[0]}  n_treat={training[1]}  "
        f"n_studies={training[2]}  tau={training[3]}"
    )

    # --- Phase 1: Run training cell reps to select threshold ---
    train_residuals: list[np.ndarray] = []
    train_planted: list[int] = []
    for rep in range(args.reps):
        r = _run_one_rep(
            training[1], training[2], training[3],
            _seed(args.seed, training[0], rep),
        )
        train_residuals.append(r["residuals"])
        train_planted.append(r["planted_idx"])

    threshold, j_trace = _select_threshold(args.reps, train_residuals, train_planted)
    print(f"Selected sheaf threshold: {threshold:.4f}")

    # --- Phase 2: Evaluate ALL cells with frozen threshold ---
    rows: list[dict] = []
    for c_idx, nt, ns, tau in cells:
        sheaf_tp = sheaf_fp = sheaf_tn = sheaf_fn = 0
        dbt_pos = dbt_neg = 0
        bucher_pos = bucher_neg = 0

        for rep in range(args.reps):
            r = _run_one_rep(nt, ns, tau, _seed(args.seed, c_idx, rep))
            for i, v in enumerate(r["residuals"]):
                flagged = bool(v > threshold)
                is_planted = (i == r["planted_idx"])
                if flagged and is_planted:
                    sheaf_tp += 1
                elif flagged and not is_planted:
                    sheaf_fp += 1
                elif not flagged and not is_planted:
                    sheaf_tn += 1
                elif not flagged and is_planted:
                    sheaf_fn += 1

            if r["dbt_p"] < 0.05:
                dbt_pos += 1
            else:
                dbt_neg += 1

            if r["bucher_any_flag"]:
                bucher_pos += 1
            else:
                bucher_neg += 1

        sheaf_sens = sheaf_tp / (sheaf_tp + sheaf_fn) if (sheaf_tp + sheaf_fn) else 0.0
        sheaf_spec = sheaf_tn / (sheaf_tn + sheaf_fp) if (sheaf_tn + sheaf_fp) else 0.0
        dbt_rate = dbt_pos / (dbt_pos + dbt_neg) if (dbt_pos + dbt_neg) else 0.0
        bucher_rate = bucher_pos / (bucher_pos + bucher_neg) if (bucher_pos + bucher_neg) else 0.0

        rows.append({
            "cell_id": _cell_id(nt, ns, tau),
            "n_treatments": nt,
            "n_studies_per_edge": ns,
            "tau_inc": tau,
            "is_training_cell": int((nt, ns, tau) == training_signature),
            "sheaf_edge_sensitivity": round(sheaf_sens, 6),
            "sheaf_edge_specificity": round(sheaf_spec, 6),
            "dbt_reject_rate": round(dbt_rate, 6),
            "bucher_any_reject_rate": round(bucher_rate, 6),
            "sheaf_threshold": round(threshold, 6),
        })
        print(
            f"  cell {c_idx:2d}/{total_cells-1}  {_cell_id(nt, ns, tau)}"
            f"  sheaf_sens={sheaf_sens:.3f}  sheaf_spec={sheaf_spec:.3f}"
            f"  dbt_rate={dbt_rate:.3f}  bucher_rate={bucher_rate:.3f}"
        )

    # --- Write CSV ---
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {len(rows)} cell rows to {out_path}")

    # --- Write threshold JSON ---
    thresh_path = Path(args.threshold_out)
    thresh_path.parent.mkdir(parents=True, exist_ok=True)
    thresh_path.write_text(
        json.dumps(
            {
                "threshold": round(threshold, 6),
                "training_cell": {
                    "n_treatments": training[1],
                    "n_studies": training[2],
                    "tau_inc": training[3],
                },
                "youden_trace": j_trace,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    print(f"Wrote threshold metadata to {thresh_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
