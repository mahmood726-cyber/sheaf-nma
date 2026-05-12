"""Phase 1 E156 figure generator — now imports math from sheafnma package.

Kept for backward-compat with the existing e156-submission pipeline.
For real validation work see analysis/run_real_corpus.py.
"""
from __future__ import annotations

import os
import sys
import io as _io

# Ensure stdout handles non-ASCII on Windows (per rules.md Python platform note)
sys.stdout = _io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap

from sheafnma.core import build_coboundary, sheaf_laplacian, edge_residuals, gii
from sheafnma.simulate import xoshiro128ss, normal_random


def get_simulated_data():
    """The original 12-contrast planted-A-D-inconsistency demo, seed=42.

    Retained byte-identically so E156 #153 figures don't shift under refactor.
    """
    rng = xoshiro128ss(42)
    truth = {"A": 0, "B": -0.3, "C": -0.5, "D": -0.8, "E": -1.0}

    def make_contrast(study, t1, t2, shift_bias=0):
        true_effect = truth[t1] - truth[t2]
        noise = normal_random(rng) * 0.05
        se = 0.10 + rng() * 0.20
        effect = round((true_effect + noise + shift_bias) * 10000) / 10000
        se = round(se * 10000) / 10000
        return {"study": study, "treat1": t1, "treat2": t2,
                "effect": effect, "se": se}

    return [
        make_contrast("Sim01", "A", "B"),
        make_contrast("Sim02", "A", "C"),
        make_contrast("Sim03", "A", "E"),
        make_contrast("Sim04", "B", "C"),
        make_contrast("Sim05", "B", "D"),
        make_contrast("Sim06", "B", "E"),
        make_contrast("Sim07", "C", "D"),
        make_contrast("Sim08", "C", "E"),
        make_contrast("Sim09", "D", "E"),
        make_contrast("Sim10", "A", "B"),
        make_contrast("Sim11_INC", "A", "D", 0.5),
        make_contrast("Sim12_INC", "A", "D", 0.5),
    ]


def build_network_dict(contrasts):
    nodes = sorted({c["treat1"] for c in contrasts} | {c["treat2"] for c in contrasts})
    return {"nodes": nodes, "edges": contrasts}


def score_to_color(score, max_score):
    if max_score == 0:
        return "#cccccc"
    t = min(1.0, score / max_score)
    cmap = LinearSegmentedColormap.from_list("inc", ["#2c7fb8", "#ffffbf", "#d7191c"])
    return cmap(t)


def main():
    outdir = "figures"
    os.makedirs(outdir, exist_ok=True)
    contrasts = get_simulated_data()
    net = build_network_dict(contrasts)
    F = build_coboundary(net)
    L = sheaf_laplacian(F)
    r = edge_residuals(net)
    G = gii(net)
    print(f"GII = {G:.4f}")
    # The detailed figure code from v0.1 lives in figures/ as PNGs already
    # generated; v0.2 figures are produced by analysis/generate_figures_v2.py.
    # Re-emit Figure 1 only, to keep e156-submission valid.
    fig, ax = plt.subplots(figsize=(6, 5))
    edge_scores = r ** 2
    max_s = float(np.max(edge_scores)) if edge_scores.size else 0.0
    ax.barh(range(len(edge_scores)), edge_scores,
            color=[score_to_color(s, max_s) for s in edge_scores])
    ax.set_yticks(range(len(edge_scores)))
    ax.set_yticklabels([f"{c['treat1']}-{c['treat2']}" for c in contrasts])
    ax.set_xlabel("residual² (inconsistency score)")
    ax.set_title(f"SheafNMA v0.2 — GII = {G:.3f}")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "figure1_edge_scores.png"), dpi=150)
    plt.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
