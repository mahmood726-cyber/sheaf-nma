import json
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def _run(*args):
    return subprocess.run(
        [sys.executable, "-m", "sheafnma", *args],
        cwd=REPO, capture_output=True, text=True,
    )


def test_cli_example_runs_and_flags_planted_edge():
    p = _run("--example")
    assert p.returncode == 0, p.stderr
    assert "GII" in p.stdout
    assert "FLAGGED" in p.stdout
    # planted inconsistency is on the A-D contrast
    assert "A vs D" in p.stdout


def test_cli_json_is_parseable():
    p = _run("--example", "--json")
    assert p.returncode == 0, p.stderr
    data = json.loads(p.stdout)
    assert data["n_treatments"] == 4
    assert data["n_edges"] == 6
    assert data["n_flagged"] == 1
    assert data["gii"] > 0


def test_cli_permute_adds_permutation_block():
    p = _run("--example", "--json", "--permute", "50", "--seed", "5")
    assert p.returncode == 0, p.stderr
    data = json.loads(p.stdout)
    assert "gii_permutation" in data
    assert data["gii_permutation"]["n_perm"] == 50
    assert 0.0 < data["gii_permutation"]["p_value"] <= 1.0


def test_cli_rejects_disconnected_csv(tmp_path):
    csv = tmp_path / "disc.csv"
    csv.write_text(
        "study,treat1,treat2,effect,se\nS1,A,B,1.0,0.3\nS2,C,D,2.0,0.3\n",
        encoding="utf-8",
    )
    p = _run(str(csv))
    assert p.returncode == 2
    assert "disconnected" in p.stderr


def test_cli_missing_file_returns_2():
    p = _run("does_not_exist_12345.csv")
    assert p.returncode == 2
    assert "no such file" in p.stderr


def test_cli_roundtrips_a_written_csv(tmp_path):
    csv = tmp_path / "net.csv"
    csv.write_text(
        "study,treat1,treat2,effect,se\n"
        "S1,A,B,1.0,1.0\nS2,B,C,1.0,1.0\nS3,A,C,2.0,1.0\n",
        encoding="utf-8",
    )
    p = _run(str(csv), "--json")
    assert p.returncode == 0, p.stderr
    data = json.loads(p.stdout)
    # consistent triangle => GII ~ 0
    assert data["gii"] < 1e-6
