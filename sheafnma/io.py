"""Standard CSV schema for NMA contrasts: study, treat1, treat2, effect, se [, design]."""
from __future__ import annotations

import csv
from pathlib import Path
from typing import Union

REQUIRED_COLUMNS = ("study", "treat1", "treat2", "effect", "se")
OPTIONAL_COLUMNS = ("design",)


def load_csv(path: Union[str, Path]) -> dict:
    """Load a contrast-format NMA CSV. Fail closed on missing columns / empty input."""
    path = Path(path)
    with path.open(encoding="utf-8") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames or []
        for col in REQUIRED_COLUMNS:
            if col not in fieldnames:
                raise KeyError(
                    f"missing required column '{col}' in {path}; "
                    f"found: {fieldnames}"
                )
        rows = list(reader)
    if not rows:
        raise ValueError(f"no contrast rows in {path}")

    edges = []
    nodes: set[str] = set()
    for row in rows:
        edge = {
            "study": row["study"],
            "treat1": row["treat1"],
            "treat2": row["treat2"],
            "effect": float(row["effect"]),
            "se": float(row["se"]),
        }
        if "design" in fieldnames and row.get("design"):
            edge["design"] = row["design"]
        edges.append(edge)
        nodes.add(row["treat1"])
        nodes.add(row["treat2"])
    return {"nodes": sorted(nodes), "edges": edges}


def save_csv(network: dict, path: Union[str, Path]) -> None:
    """Inverse of load_csv. Writes utf-8, LF line-endings."""
    path = Path(path)
    cols = list(REQUIRED_COLUMNS)
    if any("design" in e for e in network["edges"]):
        cols.append("design")
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for e in network["edges"]:
            w.writerow({c: e.get(c, "") for c in cols})
