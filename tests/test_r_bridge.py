import json
import shutil
import subprocess
import pytest

from sheafnma.r_bridge import (
    rscript_path,
    list_netmeta_datasets,
    load_netmeta_dataset,
    run_netsplit,
)

pytestmark = pytest.mark.skipif(
    not shutil.which("Rscript")
    and not shutil.which(r"C:\Program Files\R\R-4.5.2\bin\Rscript.exe"),
    reason="Rscript not on PATH",
)


def test_rscript_path_resolves():
    p = rscript_path()
    assert p is not None
    assert "Rscript" in str(p)


def test_list_netmeta_datasets_returns_non_empty():
    names = list_netmeta_datasets()
    assert isinstance(names, list)
    assert len(names) >= 5
    assert any("Senn" in n or "smoking" in n.lower() or "Dong" in n for n in names)


def test_load_senn_diabetes_has_required_columns():
    net = load_netmeta_dataset("Senn2013")
    assert "nodes" in net
    assert "edges" in net
    assert len(net["edges"]) > 0
    for e in net["edges"]:
        assert "treat1" in e
        assert "treat2" in e
        assert "effect" in e
        assert "se" in e


def test_run_netsplit_returns_per_edge_p_values():
    net = load_netmeta_dataset("Senn2013")
    result = run_netsplit(net)
    assert "edges" in result
    assert len(result["edges"]) > 0
    for edge in result["edges"]:
        assert "treat1" in edge
        assert "treat2" in edge
        assert "p_value" in edge
        assert 0.0 <= edge["p_value"] <= 1.0
