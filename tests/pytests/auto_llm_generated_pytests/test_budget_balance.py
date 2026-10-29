"""Check budget failures, completeness, and per-level tolerances in small NetCDF files."""

import importlib.util
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest
from netCDF4 import Dataset


CHECKER_PATH = Path(__file__).resolve().parents[3] / "postprocessing/check_budgets_balance/checkBudget.py"
spec = importlib.util.spec_from_file_location("budget_checker", CHECKER_PATH)
checker = importlib.util.module_from_spec(spec)
spec.loader.exec_module(checker)


def test_each_level_uses_its_own_tolerance(capsys):
    """A large earlier value must not hide a later imbalance."""
    left = np.array([0.0, 1e8, 2.0])
    right = np.array([0.0, 1e8, 0.0])
    checker.numFails = 0
    assert not checker.dispError(left, right, (left - right)[1:], 1e-5, 2, 1,
                                 "vm", "m s^{-2}", True)
    assert "vm fails at t= 2 and z= 3" in capsys.readouterr().out


def test_budget_components_and_completeness(tmp_path, monkeypatch, capsys):
    """Find component and completeness faults in the selected column and timestep."""
    path = tmp_path / "case_stats.nc"
    with Dataset(path, "w") as f:
        for name, size in (("time", 3), ("zt", 3), ("col", 2)):
            f.createDimension(name, size)
        time = f.createVariable("time", "f8", ("time",))
        time.units = "seconds"
        time[:] = [0, 300, 600]
        bt = f.createVariable("rtm_bt", "f8", ("time", "zt", "col"))
        bt.units = "kg kg^{-1} s^{-1}"
        bt[:] = 0
        component = f.createVariable("rtm_ta", "f8", ("time", "zt", "col"))
        component.long_name = "rtm budget: transport"
        component[:] = 0
        component[1, 2, 0] = 0.001
        component[:, :, 1] = 1.0
        state = f.createVariable("rtm", "f8", ("time", "zt", "col"))
        state[:] = 0
        state[1:, 1, 0] = 0.03

    monkeypatch.setattr(checker, "FILEPATH", str(tmp_path) + "/")
    checker.numFails = 0
    assert checker.checkNetcdfBudgets(path.name, 1)
    assert not checker.checkNetcdfBudgets(path.name, 2)
    output = capsys.readouterr().out
    assert "rtm completeness test fails at t= 2 and z= 2" in output
    assert "rtm fails at t= 2 and z= 3" in output


def test_nonfinite_values_do_not_pass():
    """A non-finite budget value must not silently pass the vector comparison."""
    left = np.array([0.0, np.nan])
    right = np.zeros(2)
    with pytest.raises(ValueError, match="Non-finite value"):
        checker.dispError(left, right, (left - right)[1:], 1e-5, 2, 1,
                          "vm", "m s^{-2}", True)


def test_budget_checker_detects_corruption(tmp_path):
    """The real CLI accepts a balanced budget and rejects two independent faults."""
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    checker_dir = tmp_path / "postprocessing" / "check_budgets_balance"
    checker_dir.mkdir(parents=True)
    path = output_dir / "example_stats.nc"

    with Dataset(path, "w") as f:
        for name, size in (("time", 3), ("zt", 3), ("col", 1)):
            f.createDimension(name, size)
        time = f.createVariable("time", "f8", ("time",))
        time.units = "seconds"
        time[:] = [0, 300, 600]
        bt = f.createVariable("rtm_bt", "f8", ("time", "zt", "col"))
        bt.units = "kg kg^{-1} s^{-1}"
        bt[:] = np.zeros((3, 3, 1))
        for name, value in (("rtm_ta", 0.001), ("rtm_f", -0.001)):
            term = f.createVariable(name, "f8", ("time", "zt", "col"))
            term.long_name = "rtm budget: test component"
            term[:] = np.full((3, 3, 1), value)
        state = f.createVariable("rtm", "f8", ("time", "zt", "col"))
        state[:] = np.zeros((3, 3, 1))

    def run_checker():
        """Run the checker on the isolated output and capture its status and messages."""
        return subprocess.run(
            [sys.executable, str(CHECKER_PATH), "all", "0"], cwd=checker_dir,
            text=True, capture_output=True, check=False,
        )

    control = run_checker()
    assert control.returncode == 0, control.stdout + control.stderr
    assert "Budgets successfully balance!" in control.stdout

    with Dataset(path, "r+") as f:
        f["rtm_ta"][1, 2, 0] = 0.002
    imbalance = run_checker()
    assert imbalance.returncode == 1, imbalance.stdout + imbalance.stderr
    assert "rtm fails at t= 2 and z= 3" in imbalance.stdout
    assert "Budgets fail to balance" in imbalance.stdout

    with Dataset(path, "r+") as f:
        f["rtm_ta"][1, 2, 0] = 0.001
        f["rtm"][1:, 2, 0] = 0.03
    incomplete = run_checker()
    assert incomplete.returncode == 1, incomplete.stdout + incomplete.stderr
    assert "rtm completeness test fails at t= 2 and z= 3" in incomplete.stdout
    assert "Budgets fail to balance" in incomplete.stdout
