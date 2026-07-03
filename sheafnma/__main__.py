"""Command-line entry point: ``python -m sheafnma``.

Runs a SheafNMA analysis on a contrast-format CSV (or a bundled worked
example) with pure NumPy/SciPy — no R required. Prints the Global
Inconsistency Index, the design-by-treatment chi-squared, the top per-edge
residuals, and (optionally) a permutation p-value for the GII.

Examples
--------
    python -m sheafnma --example                 # bundled worked example
    python -m sheafnma path/to/contrasts.csv      # your own CSV
    python -m sheafnma --example --permute 2000    # add a GII permutation test
    python -m sheafnma path/to/net.csv --json      # machine-readable output
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import __version__
from .comparators import bucher_loop_closure, design_by_treatment_chi2
from .core import edge_residuals, gii
from .inference import gii_permutation_test
from .io import load_csv
from .validate import NetworkValidationError, connected_components, validate_network

# Youden-frozen per-edge residual threshold from the v0.2 corpus study
# (README "v0.2 headline"). Edges with |residual| >= this are flagged.
DEFAULT_EDGE_THRESHOLD = 2.6912


def _example_network() -> dict:
    """A small, self-contained inconsistent 4-treatment network.

    A-B-C-D fully-connected. The A-D contrast is planted inconsistent
    (walking A->B->C->D gives 1+1+1 = 3, but the direct A-D says 6), so the
    example demonstrably flags the A-D edge. Numbers are illustrative only.
    """
    return {
        "nodes": ["A", "B", "C", "D"],
        "edges": [
            {"study": "S1", "treat1": "A", "treat2": "B", "effect": 1.0, "se": 0.30},
            {"study": "S2", "treat1": "B", "treat2": "C", "effect": 1.0, "se": 0.30},
            {"study": "S3", "treat1": "C", "treat2": "D", "effect": 1.0, "se": 0.30},
            {"study": "S4", "treat1": "A", "treat2": "C", "effect": 2.0, "se": 0.30},
            {"study": "S5", "treat1": "B", "treat2": "D", "effect": 2.0, "se": 0.30},
            {"study": "S6", "treat1": "A", "treat2": "D", "effect": 6.0, "se": 0.30},
        ],
    }


def analyze(network: dict, threshold: float, permute: int, seed: int) -> dict:
    validate_network(network)
    resid = edge_residuals(network).tolist()
    g = gii(network)
    dbt = design_by_treatment_chi2(network)
    bucher = bucher_loop_closure(network)

    edges_out = []
    for e, r in zip(network["edges"], resid):
        edges_out.append({
            "treat1": e["treat1"],
            "treat2": e["treat2"],
            "effect": e["effect"],
            "se": e["se"],
            "residual": r,
            "flagged": abs(r) >= threshold,
        })

    result = {
        "n_treatments": len(network["nodes"]),
        "n_edges": len(network["edges"]),
        "gii": g,
        "design_by_treatment": dbt,
        "n_bucher_loops": len(bucher),
        "edge_threshold": threshold,
        "edges": edges_out,
        "n_flagged": sum(1 for e in edges_out if e["flagged"]),
    }
    if permute > 0:
        result["gii_permutation"] = gii_permutation_test(
            network, n_perm=permute, seed=seed
        )
    return result


def _print_human(res: dict) -> None:
    print(f"SheafNMA v{__version__}")
    print(f"  treatments={res['n_treatments']}  edges={res['n_edges']}  "
          f"loops={res['n_bucher_loops']}")
    print(f"  GII (Global Inconsistency Index) = {res['gii']:.6g}")
    dbt = res["design_by_treatment"]
    print(f"  design-by-treatment chi2={dbt['chi2']:.4g} "
          f"df={dbt['df']} p={dbt['p_value']:.4g}")
    if "gii_permutation" in res:
        pt = res["gii_permutation"]
        print(f"  GII permutation p={pt['p_value']:.4g} "
              f"(n_perm={pt['n_perm']}, null_q95={pt['null_q95']:.4g})")
    print(f"  per-edge residuals (threshold |r| >= {res['edge_threshold']:.4g}):")
    for e in sorted(res["edges"], key=lambda x: -abs(x["residual"])):
        mark = "  <-- FLAGGED" if e["flagged"] else ""
        print(f"    {e['treat1']:>4} vs {e['treat2']:<4}  "
              f"r={e['residual']:+8.4f}{mark}")
    print(f"  {res['n_flagged']} edge(s) flagged as localized inconsistency.")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        prog="python -m sheafnma",
        description="Cellular-sheaf inconsistency localization for NMA.",
    )
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("csv", nargs="?", help="contrast CSV: study,treat1,treat2,effect,se")
    src.add_argument("--example", action="store_true",
                     help="run the bundled worked example (no CSV needed)")
    ap.add_argument("--threshold", type=float, default=DEFAULT_EDGE_THRESHOLD,
                    help=f"per-edge residual flag threshold (default {DEFAULT_EDGE_THRESHOLD})")
    ap.add_argument("--permute", type=int, default=0, metavar="N",
                    help="run an N-permutation GII test (0 = skip, default)")
    ap.add_argument("--seed", type=int, default=42,
                    help="permutation-test RNG seed (default 42)")
    ap.add_argument("--json", action="store_true", help="emit JSON instead of text")
    ap.add_argument("--version", action="version",
                    version=f"sheafnma {__version__}")
    args = ap.parse_args(argv)

    if args.example:
        network = _example_network()
    else:
        path = Path(args.csv)
        if not path.exists():
            print(f"error: no such file: {path}", file=sys.stderr)
            return 2
        network = load_csv(path)

    try:
        res = analyze(network, args.threshold, args.permute, args.seed)
    except NetworkValidationError as exc:
        print(f"error: invalid network: {exc}", file=sys.stderr)
        comps = connected_components(network)
        if len(comps) > 1:
            print(f"  ({len(comps)} disconnected components found)", file=sys.stderr)
        return 2

    if args.json:
        print(json.dumps(res, indent=2))
    else:
        _print_human(res)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
