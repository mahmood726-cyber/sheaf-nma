"""Deterministic RNG (xoshiro128**) + planted-inconsistency NMA generator.

The xoshiro implementation is lifted unchanged from
generate_figures.py for byte-identical determinism across the v0.1 -> v0.2
transition (the E156 #153 figures use seed=42; we don't want them shifting).
"""
from __future__ import annotations

import string
from typing import Callable


def splitmix32(seed: int) -> Callable[[], int]:
    state = [seed & 0xFFFFFFFF]
    def next32() -> int:
        state[0] = (state[0] + 0x9E3779B9) & 0xFFFFFFFF
        z = state[0]
        z = ((z ^ (z >> 16)) * 0x85EBCA6B) & 0xFFFFFFFF
        z = ((z ^ (z >> 13)) * 0xC2B2AE35) & 0xFFFFFFFF
        z = (z ^ (z >> 16)) & 0xFFFFFFFF
        return z
    return next32


def xoshiro128ss(seed: int) -> Callable[[], float]:
    sm = splitmix32(seed)
    s = [sm(), sm(), sm(), sm()]

    def rotl(x: int, k: int) -> int:
        return ((x << k) | (x >> (32 - k))) & 0xFFFFFFFF

    def nxt() -> float:
        result = (rotl((s[1] * 5) & 0xFFFFFFFF, 7) * 9) & 0xFFFFFFFF
        t = (s[1] << 9) & 0xFFFFFFFF
        s[2] ^= s[0]
        s[3] ^= s[1]
        s[1] ^= s[2]
        s[0] ^= s[3]
        s[2] ^= t
        s[3] = rotl(s[3], 11)
        return result / 4294967296.0  # uniform [0,1)
    return nxt


def normal_random(rng: Callable[[], float]) -> float:
    """Box-Muller standard normal."""
    import math
    u1 = max(rng(), 1e-12)
    u2 = rng()
    return math.sqrt(-2.0 * math.log(u1)) * math.cos(2.0 * math.pi * u2)


def generate_planted_network(
    n_treatments: int,
    n_studies_per_edge: int,
    tau_inc: float,
    seed: int,
) -> dict:
    """Random fully-connected NMA with planted inconsistency on one random edge.

    Truth: treatment k has true effect = -0.2*k vs treatment 0.
    Each edge gets n_studies_per_edge studies, IV-pooled.
    One randomly-chosen edge gets tau_inc added to its pooled effect.
    """
    rng = xoshiro128ss(seed)
    letters = string.ascii_uppercase
    nodes = [letters[i] for i in range(n_treatments)]
    truth = {nd: -0.2 * i for i, nd in enumerate(nodes)}

    # all pairs
    pairs = []
    for i in range(n_treatments):
        for j in range(i + 1, n_treatments):
            pairs.append((nodes[i], nodes[j]))

    # planted edge index
    planted_idx = int(rng() * len(pairs))

    edges = []
    for e_idx, (a, b) in enumerate(pairs):
        # n_studies_per_edge studies, IV-pooled
        effects = []
        ses = []
        for k in range(n_studies_per_edge):
            true_eff = truth[a] - truth[b]
            noise = normal_random(rng) * 0.05
            se = 0.10 + rng() * 0.10
            effects.append(true_eff + noise)
            ses.append(se)
        # inverse-variance pool
        weights = [1.0 / s ** 2 for s in ses]
        w_sum = sum(weights)
        pooled_effect = sum(e * w for e, w in zip(effects, weights)) / w_sum
        pooled_se = (1.0 / w_sum) ** 0.5
        if e_idx == planted_idx:
            pooled_effect += tau_inc
        edges.append({
            "study": f"E{e_idx:02d}",
            "treat1": a,
            "treat2": b,
            "effect": pooled_effect,
            "se": pooled_se,
        })

    return {"nodes": nodes, "edges": edges, "planted_edge_idx": planted_idx}
