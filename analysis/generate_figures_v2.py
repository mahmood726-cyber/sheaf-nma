"""Paper figures from analysis/results/{power_results,corpus_results}.csv.

Three figures:
  fig_power_grid.png       — heatmap: sheaf edge sensitivity vs DBT reject-rate
                             across the 36-cell grid (one panel per τ).
  fig_corpus_agreement.png — 2x2 contingency: sheaf-flagged vs netsplit-flagged
                             edges, per-corpus and pooled.
  fig_gii_forest.png       — per-NMA GII values (sorted), with the n_edges_flagged
                             annotated on each bar.

All figures use the same color palette as the v0.1 dashboard for visual continuity.
"""
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def _load_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8") as f:
        return list(csv.DictReader(f))


def fig_power_grid(power_rows: list[dict], outpath: Path) -> None:
    taus = sorted({float(r["tau_inc"]) for r in power_rows})
    fig, axes = plt.subplots(1, len(taus), figsize=(4 * len(taus), 4), sharey=True)
    if len(taus) == 1:
        axes = [axes]
    for ax, tau in zip(axes, taus):
        cell = [r for r in power_rows if float(r["tau_inc"]) == tau]
        nt = sorted({int(r["n_treatments"]) for r in cell})
        ns = sorted({int(r["n_studies_per_edge"]) for r in cell})
        grid = np.zeros((len(nt), len(ns)))
        for r in cell:
            i = nt.index(int(r["n_treatments"]))
            j = ns.index(int(r["n_studies_per_edge"]))
            grid[i, j] = float(r["sheaf_edge_sensitivity"])
        im = ax.imshow(grid, vmin=0, vmax=1, cmap="viridis", aspect="auto")
        ax.set_xticks(range(len(ns)), [str(s) for s in ns])
        ax.set_yticks(range(len(nt)), [str(t) for t in nt])
        ax.set_xlabel("studies/edge")
        ax.set_ylabel("treatments")
        ax.set_title(f"τ_inc = {tau}")
    fig.colorbar(im, ax=axes, fraction=0.03)
    fig.suptitle(
        "SheafNMA per-edge sensitivity across (n_treatments × n_studies) × τ",
        fontsize=11,
    )
    outpath.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)


def fig_corpus_agreement(corpus_rows: list[dict], outpath: Path) -> None:
    # 2x2 sheaf_flagged x netsplit_flagged
    counts = {(0, 0): 0, (0, 1): 0, (1, 0): 0, (1, 1): 0}
    for r in corpus_rows:
        if r["netsplit_p"] in ("", None):
            continue
        a = int(r["sheaf_flagged"])
        b = int(r["netsplit_flagged"])
        counts[(a, b)] += 1
    fig, ax = plt.subplots(figsize=(4.5, 4))
    grid = np.array([[counts[(0, 0)], counts[(0, 1)]],
                     [counts[(1, 0)], counts[(1, 1)]]])
    im = ax.imshow(grid, cmap="Blues")
    for i in range(2):
        for j in range(2):
            ax.text(j, i, str(grid[i, j]), ha="center", va="center", fontsize=14)
    ax.set_xticks([0, 1], ["netsplit not flag", "netsplit flagged"])
    ax.set_yticks([0, 1], ["sheaf not flag", "sheaf flagged"])
    ax.set_title("Per-edge agreement: SheafNMA vs node-splitting")
    fig.colorbar(im, ax=ax, fraction=0.04)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)


def fig_gii_forest(corpus_rows: list[dict], outpath: Path) -> None:
    by_nma: dict[str, dict] = defaultdict(lambda: {"gii": 0.0, "flagged": 0})
    for r in corpus_rows:
        by_nma[r["nma_id"]]["gii"] = float(r["nma_gii"])
        by_nma[r["nma_id"]]["flagged"] += int(r["sheaf_flagged"])
    items = sorted(by_nma.items(), key=lambda kv: kv[1]["gii"])
    labels = [k for k, _ in items]
    giis = [v["gii"] for _, v in items]
    flagged = [v["flagged"] for _, v in items]

    fig, ax = plt.subplots(figsize=(6, max(4, 0.3 * len(labels))))
    ax.barh(range(len(labels)), giis, color="#2c7fb8")
    for i, f in enumerate(flagged):
        if f > 0:
            ax.text(giis[i], i, f"  {f} edge(s) flagged", va="center", fontsize=8)
    ax.set_yticks(range(len(labels)), labels)
    ax.set_xlabel("GII (Σ residual²)")
    ax.set_title("Per-NMA Global Inconsistency Index")
    outpath.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--power", default="analysis/results/power_results.csv")
    ap.add_argument("--corpus", default="analysis/results/corpus_results.csv")
    ap.add_argument("--outdir", default="paper/figures")
    args = ap.parse_args()

    out = Path(args.outdir)
    power_rows = _load_csv(Path(args.power))
    corpus_rows = _load_csv(Path(args.corpus))
    fig_power_grid(power_rows, out / "fig_power_grid.png")
    fig_corpus_agreement(corpus_rows, out / "fig_corpus_agreement.png")
    fig_gii_forest(corpus_rows, out / "fig_gii_forest.png")
    print(f"Wrote 3 figures to {out}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
