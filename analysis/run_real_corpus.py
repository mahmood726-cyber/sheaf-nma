"""Apply SheafNMA + DBT-χ² + Bucher + netsplit to every NMA in corpus/data/.

Output: analysis/results/corpus_results.csv with one row per (NMA, edge)
where SheafNMA's per-edge residual and netmeta::netsplit's per-edge p-value
can both be compared.

Usage:
    python -m analysis.run_real_corpus
    python -m analysis.run_real_corpus --threshold 0.45   # custom sheaf threshold
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from sheafnma.io import load_csv
from sheafnma.core import build_coboundary, edge_residuals, gii
from sheafnma.comparators import bucher_loop_closure, design_by_treatment_chi2
from sheafnma.r_bridge import run_netsplit


def analyze_one(nma_id: str, network: dict, sheaf_threshold: float) -> list[dict]:
    G = gii(network)
    residuals = edge_residuals(network).tolist()
    dbt = design_by_treatment_chi2(network)
    bucher = bucher_loop_closure(network)
    try:
        ns = run_netsplit(network)
        ns_lookup = {
            tuple(sorted([e["treat1"], e["treat2"]])): e["p_value"]
            for e in ns["edges"]
        }
    except Exception as exc:  # R missing, netmeta error, etc.
        ns_lookup = {}
        print(f"  WARN netsplit unavailable for {nma_id}: {exc}")

    rows = []
    for i, edge in enumerate(network["edges"]):
        key = tuple(sorted([edge["treat1"], edge["treat2"]]))
        ns_p = ns_lookup.get(key)
        rows.append({
            "nma_id": nma_id,
            "edge_idx": i,
            "treat1": edge["treat1"],
            "treat2": edge["treat2"],
            "effect": edge["effect"],
            "se": edge["se"],
            "sheaf_residual": residuals[i],
            "sheaf_flagged": int(abs(residuals[i]) > sheaf_threshold),
            "nma_gii": G,
            "dbt_chi2": dbt["chi2"],
            "dbt_df": dbt["df"],
            "dbt_p": dbt["p_value"],
            "n_loops_flagged_bucher": int(sum(1 for b in bucher if b["p_value"] < 0.05)),
            "netsplit_p": ns_p if ns_p is not None else "",
            "netsplit_flagged": int(ns_p is not None and ns_p < 0.05),
        })
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--corpus-dir", default="corpus/data")
    ap.add_argument("--out", default="analysis/results/corpus_results.csv")
    ap.add_argument("--threshold", type=float, default=0.30,
                    help="sheaf-residual threshold for flagging (placeholder; "
                         "frozen value comes from power_simulation.py)")
    args = ap.parse_args()

    corpus = Path(args.corpus_dir)
    manifest = json.loads((corpus / "_manifest.json").read_text(encoding="utf-8"))
    all_rows: list[dict] = []
    for nma_id, _meta in manifest["datasets"].items():
        csv_path = corpus / f"{nma_id}.csv"
        if not csv_path.exists():
            print(f"SKIP {nma_id}: csv not found")
            continue
        print(f"Analyzing {nma_id} ...")
        net = load_csv(csv_path)
        all_rows.extend(analyze_one(nma_id, net, args.threshold))

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if not all_rows:
        print("No rows produced.")
        return 1
    cols = list(all_rows[0].keys())
    with out_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(all_rows)
    print(f"Wrote {len(all_rows)} rows to {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
