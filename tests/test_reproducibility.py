"""Two consecutive small-grid power_simulation runs must be byte-identical."""
import subprocess
import sys
from pathlib import Path


def test_power_simulation_is_deterministic(tmp_path):
    repo = Path(__file__).resolve().parents[1]
    out1 = tmp_path / "r1.csv"
    out2 = tmp_path / "r2.csv"
    cmd_base = [
        sys.executable, "-m", "analysis.power_simulation",
        "--reps", "5",
        "--treatments", "4",
        "--studies", "3",
        "--taus", "0,0.5",
        "--seed", "42",
    ]
    subprocess.run(cmd_base + ["--out", str(out1), "--threshold-out", str(tmp_path / "t1.json")], cwd=repo, check=True)
    subprocess.run(cmd_base + ["--out", str(out2), "--threshold-out", str(tmp_path / "t2.json")], cwd=repo, check=True)
    assert out1.read_bytes() == out2.read_bytes()
