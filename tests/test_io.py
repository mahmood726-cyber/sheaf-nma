import io
import pytest
from sheafnma.io import load_csv, save_csv, REQUIRED_COLUMNS


CSV_OK = """study,treat1,treat2,effect,se
Sim01,A,B,0.30,0.10
Sim02,A,C,0.50,0.12
"""


def test_load_csv_basic(tmp_path):
    p = tmp_path / "n.csv"
    p.write_text(CSV_OK, encoding="utf-8")
    net = load_csv(p)
    assert sorted(net["nodes"]) == ["A", "B", "C"]
    assert len(net["edges"]) == 2
    assert net["edges"][0]["effect"] == 0.30


def test_load_csv_missing_column_fails_closed(tmp_path):
    bad = "study,treat1,treat2,effect\nSim01,A,B,0.30\n"  # no `se`
    p = tmp_path / "bad.csv"
    p.write_text(bad, encoding="utf-8")
    with pytest.raises(KeyError, match="se"):
        load_csv(p)


def test_load_csv_empty_fails_closed(tmp_path):
    p = tmp_path / "empty.csv"
    p.write_text("study,treat1,treat2,effect,se\n", encoding="utf-8")
    with pytest.raises(ValueError, match="no contrast rows"):
        load_csv(p)


def test_save_csv_roundtrip(tmp_path):
    net = {
        "nodes": ["A", "B"],
        "edges": [{"study": "S1", "treat1": "A", "treat2": "B",
                   "effect": 0.5, "se": 0.1}],
    }
    p = tmp_path / "out.csv"
    save_csv(net, p)
    net2 = load_csv(p)
    assert net2["edges"][0]["effect"] == 0.5
    assert net2["edges"][0]["se"] == 0.1


def test_required_columns_constant():
    assert "study" in REQUIRED_COLUMNS
    assert "treat1" in REQUIRED_COLUMNS
    assert "treat2" in REQUIRED_COLUMNS
    assert "effect" in REQUIRED_COLUMNS
    assert "se" in REQUIRED_COLUMNS
