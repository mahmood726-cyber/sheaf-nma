"""Canonical NMA inconsistency tests for comparator use against SheafNMA.

- design_by_treatment_chi2: global Wald χ² (Higgins et al. 2012 RSM §3)
- bucher_loop_closure:      per-loop direct-vs-indirect (Bucher 1997)
- enumerate_independent_loops: cycle basis of the contrast graph
"""
from __future__ import annotations

from collections import defaultdict
from typing import Iterable

import numpy as np
from scipy import stats as sp_stats


def _pool_edges_by_pair(network: dict) -> dict[tuple[str, str], dict]:
    """IV-pool multi-study edges by (sorted) treatment pair."""
    buckets: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for e in network["edges"]:
        key = tuple(sorted([e["treat1"], e["treat2"]]))
        # orient: effect is treat1→treat2; flip sign if we sorted-swap
        sign = 1.0 if (e["treat1"], e["treat2"]) == key else -1.0
        buckets[key].append({"effect": sign * e["effect"], "se": e["se"]})
    pooled = {}
    for key, rows in buckets.items():
        weights = np.array([1.0 / r["se"] ** 2 for r in rows])
        effects = np.array([r["effect"] for r in rows])
        w_sum = weights.sum()
        pooled[key] = {
            "effect": float((effects * weights).sum() / w_sum),
            "se": float((1.0 / w_sum) ** 0.5),
        }
    return pooled


def enumerate_independent_loops(network: dict) -> list[dict]:
    """Cycle basis via spanning-tree fundamental cycles."""
    nodes = list(network["nodes"])
    adj: dict[str, set[str]] = {nd: set() for nd in nodes}
    pooled = _pool_edges_by_pair(network)
    for (a, b) in pooled.keys():
        adj[a].add(b)
        adj[b].add(a)

    parent: dict[str, str | None] = {}
    visited: set[str] = set()
    tree_edges: set[tuple[str, str]] = set()
    if nodes:
        root = nodes[0]
        stack = [root]
        parent[root] = None
        while stack:
            u = stack.pop()
            if u in visited:
                continue
            visited.add(u)
            for v in adj[u]:
                if v not in visited and v not in parent:
                    parent[v] = u
                    tree_edges.add(tuple(sorted([u, v])))
                    stack.append(v)

    def path_to_root(x: str) -> list[str]:
        path = [x]
        while parent.get(x) is not None:
            x = parent[x]
            path.append(x)
        return path

    loops = []
    for (a, b) in pooled.keys():
        if (a, b) in tree_edges:
            continue
        # non-tree edge closes a fundamental cycle
        pa = path_to_root(a)
        pb = path_to_root(b)
        # find LCA
        set_pb = set(pb)
        lca = next(n for n in pa if n in set_pb)
        cycle = pa[: pa.index(lca) + 1] + list(reversed(pb[: pb.index(lca)]))
        loops.append({"nodes": cycle, "closing_edge": (a, b)})
    return loops


def bucher_loop_closure(network: dict) -> list[dict]:
    """For each independent loop, compute direct-minus-indirect difference, SE, z, p."""
    pooled = _pool_edges_by_pair(network)
    results = []
    for loop in enumerate_independent_loops(network):
        cycle = loop["nodes"]
        # Walk the cycle: sum signed pooled edge effects around the loop.
        signed_sum = 0.0
        var_sum = 0.0
        n = len(cycle)
        for i in range(n):
            a, b = cycle[i], cycle[(i + 1) % n]
            key = tuple(sorted([a, b]))
            edge = pooled[key]
            sign = 1.0 if (a, b) == key else -1.0
            signed_sum += sign * edge["effect"]
            var_sum += edge["se"] ** 2
        se = float(var_sum ** 0.5)
        z = signed_sum / se if se > 0 else 0.0
        p = float(2.0 * (1.0 - sp_stats.norm.cdf(abs(z))))
        results.append({
            "loop_nodes": cycle,
            "diff": float(signed_sum),
            "se": se,
            "z": float(z),
            "p_value": p,
        })
    return results


def design_by_treatment_chi2(network: dict) -> dict:
    """Higgins 2012 §3.2 global Wald χ² for design-by-treatment interaction.

    df = (number of designs) − (number of treatments) + 1.
    Test statistic: sum of squared standardized Bucher loop closures, summed
    over the fundamental cycle basis. This is the loop-based realisation of
    the DBT test (equivalent to the full Wald-χ² form for single-D stalks).
    """
    loops = bucher_loop_closure(network)
    chi2 = float(sum((r["z"]) ** 2 for r in loops))
    designs = {e.get("design") or tuple(sorted([e["treat1"], e["treat2"]]))
               for e in network["edges"]}
    n_designs = len(designs)
    n_treatments = len(network["nodes"])
    df = max(1, n_designs - n_treatments + 1)
    p = float(1.0 - sp_stats.chi2.cdf(chi2, df=df))
    return {"chi2": chi2, "df": df, "p_value": p, "n_designs": n_designs}
