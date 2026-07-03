"""Input validation for NMA networks (public entry-point guard).

NMA good practice — and the SheafNMA methodology — requires checking network
connectivity *before* any inconsistency statistic is interpreted: a
disconnected network has no shared reference frame across components, so its
sheaf Laplacian is singular and `gii` / `edge_residuals` are undefined.

This module exposes those checks as first-class, testable functions so callers
can fail closed with a clear message *before* the linear-algebra layer, rather
than surfacing a raw ``LinAlgError`` (or a solve-time ``ValueError``) deep in
``core``. Nothing here alters any scientific output — it only rejects inputs
that have no well-defined answer.
"""
from __future__ import annotations

from collections import defaultdict


class NetworkValidationError(ValueError):
    """Raised when a network cannot support a well-defined sheaf analysis."""


REQUIRED_EDGE_KEYS = ("treat1", "treat2", "effect", "se")


def _iter_edges(network: dict):
    if not isinstance(network, dict):
        raise NetworkValidationError(
            f"network must be a dict with 'nodes' and 'edges'; got {type(network).__name__}"
        )
    edges = network.get("edges")
    if edges is None:
        raise NetworkValidationError("network has no 'edges' key")
    if len(edges) == 0:
        raise NetworkValidationError("network has zero edges; nothing to analyse")
    return edges


def connected_components(network: dict) -> list[list[str]]:
    """Return the connected components of the treatment graph, sorted.

    Each component is a sorted list of node labels; the outer list is sorted by
    the smallest label in each component for determinism. Multi-study edges on
    the same pair collapse to a single graph edge.
    """
    edges = _iter_edges(network)
    adj: dict[str, set[str]] = defaultdict(set)
    seen_nodes: set[str] = set()
    for e in edges:
        a, b = e["treat1"], e["treat2"]
        seen_nodes.add(a)
        seen_nodes.add(b)
        if a != b:
            adj[a].add(b)
            adj[b].add(a)
    # include any isolated declared nodes too (edges never reference them, but
    # they are part of the treatment set)
    for nd in network.get("nodes", []):
        seen_nodes.add(nd)

    components: list[list[str]] = []
    unvisited = set(seen_nodes)
    while unvisited:
        start = min(unvisited)
        stack = [start]
        comp: set[str] = set()
        while stack:
            u = stack.pop()
            if u in comp:
                continue
            comp.add(u)
            unvisited.discard(u)
            for v in adj[u]:
                if v not in comp:
                    stack.append(v)
        components.append(sorted(comp))
    return sorted(components, key=lambda c: c[0] if c else "")


def is_connected(network: dict) -> bool:
    """True iff every treatment sits in a single connected component."""
    return len(connected_components(network)) <= 1


def validate_network(network: dict, *, require_connected: bool = True) -> dict:
    """Validate a network for sheaf analysis; raise on the first fatal problem.

    Checks (in order):
      1. structural: dict with a non-empty 'edges' list, each edge carrying the
         required keys with numeric, finite ``effect`` and strictly-positive,
         finite ``se``;
      2. no self-loops (``treat1 == treat2``);
      3. connectivity (only when ``require_connected``): a disconnected network
         has a singular sheaf Laplacian and no well-defined GII.

    Returns the (unmodified) ``network`` so it can be used inline::

        gii(validate_network(net))

    Raises ``NetworkValidationError`` (a subclass of ``ValueError``) with an
    actionable message. Does **not** mutate the input or any numeric value.
    """
    edges = _iter_edges(network)

    for i, e in enumerate(edges):
        for key in REQUIRED_EDGE_KEYS:
            if key not in e:
                raise NetworkValidationError(
                    f"edge {i} is missing required key '{key}'; has {sorted(e)}"
                )
        a, b = e["treat1"], e["treat2"]
        if a == b:
            raise NetworkValidationError(
                f"edge {i} is a self-loop ({a!r}=={b!r}); NMA contrasts must "
                f"compare two distinct treatments"
            )
        try:
            eff = float(e["effect"])
            se = float(e["se"])
        except (TypeError, ValueError) as exc:
            raise NetworkValidationError(
                f"edge {i} has non-numeric effect/se: "
                f"effect={e['effect']!r}, se={e['se']!r}"
            ) from exc
        if not _isfinite(eff):
            raise NetworkValidationError(f"edge {i} has non-finite effect {eff!r}")
        if not _isfinite(se) or se <= 0.0:
            raise NetworkValidationError(
                f"edge {i} has non-positive or non-finite se {se!r}; "
                f"standard errors must be > 0"
            )

    if require_connected:
        comps = connected_components(network)
        if len(comps) > 1:
            summary = "; ".join(
                "{" + ", ".join(c) + "}" for c in comps
            )
            raise NetworkValidationError(
                f"network is disconnected into {len(comps)} components: {summary}. "
                f"SheafNMA requires a single connected network - split it or add "
                f"a bridging comparison before analysis."
            )
    return network


def _isfinite(x: float) -> bool:
    return x == x and x not in (float("inf"), float("-inf"))
