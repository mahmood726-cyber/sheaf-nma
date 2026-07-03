"""Reference distribution for the Global Inconsistency Index (GII).

The v0.2 methods paper (§4.4 Limitations) noted: *"the current GII does not
have a formal reference distribution, so p-values for the global test are not
yet available. A chi-squared approximation ... is a natural candidate but
requires validation."* This module supplies the paper's suggested calibration
route — a **permutation test** — as a concrete, deterministic function.

Under the null of global consistency, every closed loop in the network should
close to within sampling error, i.e. the *signs* attached to each contrast are
exchangeable around each loop. We approximate that null by repeatedly flipping
the orientation (sign of the effect) of a random subset of edges and recomputing
the GII: sign flips preserve each edge's magnitude and SE (so the marginal
precision structure is untouched) while destroying any coherent loop closure,
generating a null spread of GII values. The observed GII is compared against
that permutation distribution.

This is a *calibration aid / sensitivity tool*, not a replacement for the
design-by-treatment test or node-splitting — it gives the global sheaf statistic
an empirical p-value on the same algebraic object. No planted numbers or fixture
outputs are changed; the observed GII is computed by the unmodified ``core.gii``.
"""
from __future__ import annotations

from copy import deepcopy

from .core import gii
from .simulate import xoshiro128ss
from .validate import validate_network


def gii_permutation_test(
    network: dict,
    n_perm: int = 1000,
    seed: int = 42,
    *,
    validate: bool = True,
) -> dict:
    """Permutation p-value for the observed GII.

    Parameters
    ----------
    network : dict
        Standard sheafnma network (``nodes`` + ``edges``).
    n_perm : int
        Number of sign-flip permutations (default 1000).
    seed : int
        Seed for the deterministic xoshiro128** RNG (default 42), so repeated
        calls with the same arguments return byte-identical results.
    validate : bool
        Run :func:`sheafnma.validate.validate_network` first (default True).

    Returns
    -------
    dict with keys:
        ``gii_observed``  – GII of the input network (from ``core.gii``);
        ``p_value``       – ``(1 + #{perm >= observed}) / (n_perm + 1)``
                            (add-one correction, never exactly 0);
        ``n_perm``        – permutations actually run;
        ``null_mean``     – mean of the permutation GII distribution;
        ``null_q95``      – 95th percentile of the null distribution.
    """
    if validate:
        validate_network(network)
    if n_perm < 1:
        raise ValueError("n_perm must be >= 1")

    observed = gii(network)
    rng = xoshiro128ss(seed)
    m = len(network["edges"])

    null_gii: list[float] = []
    ge = 0
    for _ in range(n_perm):
        perm = deepcopy(network)
        for e in perm["edges"]:
            if rng() < 0.5:
                e["effect"] = -e["effect"]
        g = gii(perm)
        null_gii.append(g)
        if g >= observed:
            ge += 1

    p_value = (1 + ge) / (n_perm + 1)
    null_sorted = sorted(null_gii)
    q95 = _percentile(null_sorted, 0.95)
    null_mean = sum(null_gii) / len(null_gii)
    return {
        "gii_observed": float(observed),
        "p_value": float(p_value),
        "n_perm": int(n_perm),
        "null_mean": float(null_mean),
        "null_q95": float(q95),
    }


def _percentile(sorted_vals: list[float], q: float) -> float:
    """Linear-interpolation percentile on an already-sorted list."""
    if not sorted_vals:
        return float("nan")
    if len(sorted_vals) == 1:
        return sorted_vals[0]
    pos = q * (len(sorted_vals) - 1)
    lo = int(pos)
    hi = min(lo + 1, len(sorted_vals) - 1)
    frac = pos - lo
    return sorted_vals[lo] * (1 - frac) + sorted_vals[hi] * frac
