"""Known-truth NMA network generator for SheafNMA truth-recovery.

INDEPENDENT of the repo's own `sheafnma.simulate.generate_planted_network`.
We construct networks directly from a known node-truth vector so that, BY
CONSTRUCTION, every fundamental loop obeys the consistency identity

    d_AC = d_AB + d_BC   (here: effect_e = truth[t2] - truth[t1])

A CONSISTENT network has every edge's mean equal to that identity (+ sampling
noise). An INCONSISTENT network has a single known delta added to one design's
pooled effect, breaking exactly that one loop.

Uses numpy's own Generator (PCG64) for reproducibility, deliberately NOT the
repo's xoshiro, so the test is an out-of-distribution check on the detector.
"""
from __future__ import annotations

import string
import numpy as np


def make_network(
    n_treatments: int,
    n_studies_per_edge: int,
    rng: np.random.Generator,
    inject_delta: float = 0.0,
    inject_pair: tuple[str, str] | None = None,
    noise_sd: float = 0.05,
    se_lo: float = 0.10,
    se_hi: float = 0.20,
) -> dict:
    """Fully-connected NMA from known node truth.

    truth[node_i] = -0.2 * i  (so all loops are consistent by construction).
    If inject_delta != 0, add it to the pooled effect of one chosen pair,
    breaking that loop by a KNOWN amount.
    Returns dict with nodes, edges, truth, and inconsistent_pair (or None).
    """
    letters = string.ascii_uppercase
    nodes = [letters[i] for i in range(n_treatments)]
    truth = {nd: -0.2 * i for i, nd in enumerate(nodes)}

    pairs = [(nodes[i], nodes[j])
             for i in range(n_treatments) for j in range(i + 1, n_treatments)]

    if inject_delta != 0.0:
        if inject_pair is None:
            inject_pair = pairs[rng.integers(len(pairs))]
        inject_pair = tuple(sorted(inject_pair))

    edges = []
    inconsistent_pair = None
    for (a, b) in pairs:
        true_eff = truth[b] - truth[a]  # effect_e = truth[t2]-truth[t1]
        effects, ses = [], []
        for _ in range(n_studies_per_edge):
            se = se_lo + rng.random() * (se_hi - se_lo)
            noise = rng.normal(0.0, noise_sd)
            effects.append(true_eff + noise)
            ses.append(se)
        w = np.array([1.0 / s ** 2 for s in ses])
        pooled_eff = float((np.array(effects) * w).sum() / w.sum())
        pooled_se = float((1.0 / w.sum()) ** 0.5)
        if inject_delta != 0.0 and tuple(sorted((a, b))) == inject_pair:
            pooled_eff += inject_delta
            inconsistent_pair = inject_pair
        edges.append({
            "study": f"E_{a}{b}",
            "treat1": a, "treat2": b,
            "effect": pooled_eff, "se": pooled_se,
        })

    return {
        "nodes": nodes,
        "edges": edges,
        "truth": truth,
        "inconsistent_pair": inconsistent_pair,
    }
