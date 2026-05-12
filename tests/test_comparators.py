import math
import pytest
from sheafnma.comparators import (
    design_by_treatment_chi2,
    bucher_loop_closure,
    enumerate_independent_loops,
)


def _consistent_triangle():
    return {
        "nodes": ["A", "B", "C"],
        "edges": [
            {"study": "S1", "treat1": "A", "treat2": "B",
             "effect": 1.0, "se": 0.1, "design": "AB"},
            {"study": "S2", "treat1": "B", "treat2": "C",
             "effect": 1.0, "se": 0.1, "design": "BC"},
            {"study": "S3", "treat1": "A", "treat2": "C",
             "effect": 2.0, "se": 0.1, "design": "AC"},
        ],
    }


def _inconsistent_triangle():
    n = _consistent_triangle()
    n["edges"][2]["effect"] = 5.0  # break the loop
    return n


def test_enumerate_loops_triangle():
    loops = enumerate_independent_loops(_consistent_triangle())
    assert len(loops) == 1
    loop = loops[0]
    assert set(loop["nodes"]) == {"A", "B", "C"}


def test_bucher_consistent_triangle_p_near_one():
    n = _consistent_triangle()
    results = bucher_loop_closure(n)
    assert len(results) == 1
    assert results[0]["p_value"] > 0.5


def test_bucher_inconsistent_triangle_p_small():
    n = _inconsistent_triangle()
    results = bucher_loop_closure(n)
    assert results[0]["p_value"] < 0.001


def test_dbt_chi2_consistent_triangle_p_high():
    n = _consistent_triangle()
    result = design_by_treatment_chi2(n)
    assert "chi2" in result
    assert "df" in result
    assert "p_value" in result
    assert result["p_value"] > 0.10


def test_dbt_chi2_inconsistent_triangle_p_low():
    n = _inconsistent_triangle()
    result = design_by_treatment_chi2(n)
    assert result["p_value"] < 0.05


def test_dbt_chi2_df_formula():
    """df = number of designs − number of treatments + 1 (Higgins 2012 §3.2)."""
    n = _consistent_triangle()  # 3 designs, 3 treatments → df = 1
    result = design_by_treatment_chi2(n)
    assert result["df"] == 1
