"""Structured bindiff results and the runners that consume them."""

import importlib.util
import json
import shutil
import subprocess
import sys
from pathlib import Path

import netCDF4
import numpy as np
import pytest

from tests.run_jax_vs_fortran_cases import _first_failing_timestep
from tests.run_python_vs_fortran_cases import _average_earliest_timestep


BINDIFF = Path(__file__).resolve().parents[2] / "run_scripts" / "run_bindiff_all.py"


def load_bindiff():
    spec = importlib.util.spec_from_file_location("bindiff_report_test", BINDIFF)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_case(path, *, right):
    with netCDF4.Dataset(path, "w") as nc:
        nc.createDimension("time", 2)
        nc.createDimension("col", 2 if right else 1)
        nc.createDimension("empty_dim", 0)
        nc.createVariable("time", "f8", ("time",))[:] = [60, 120]
        for name, values in {
            "all_zero": [0, 0],
            "identical": [2, 2],
            "within_tolerance": [1 + 1e-8, 1] if right else [1, 1],
            "different": [1, 1.01] if right else [1, 1],
        }.items():
            nc.createVariable(name, "f8", ("time",))[:] = values
        nc.createVariable("shape", "f8", ("time", "col"))[:] = np.ones((2, 2 if right else 1))
        nc.createVariable("text", "S1", ("time",))[:] = np.array([b"a", b"b"])
        nc.createVariable("empty", "f8", ("empty_dim",))
        nc.createVariable("only_right" if right else "only_left", "f8", ("time",))[:] = [1, 1]


def test_json_reports_paths_logs_and_exclusive_variable_categories(tmp_path):
    left, right = tmp_path / "left", tmp_path / "right"
    left.mkdir()
    right.mkdir()
    write_case(left / "case_stats.nc", right=False)
    write_case(right / "case_stats.nc", right=True)
    shutil.copyfile(left / "case_stats.nc", left / "case_zt.nc")
    shutil.copyfile(left / "case_stats.nc", right / "case_zt.nc")

    bindiff = load_bindiff()
    bindiff.outFilePath = str(tmp_path / "logs")
    Path(bindiff.outFilePath).mkdir()
    report_path = tmp_path / "result.json"
    _, differs, _, passed, failed = bindiff.find_diffs_in_all_files(
        str(left), str(right), "replace", 0, 1e-6, 1e-6, False,
        result_json=str(report_path), strict=True,
    )

    assert (differs, passed, failed) == (True, [], ["case"])
    report = json.loads(report_path.read_text())
    assert report["inputs"] == [str(left.resolve()), str(right.resolve())]
    assert report["result_json"] == str(report_path.resolve())
    assert report["status"] == "diff"
    assert report["average_earliest_timestep"] == 1.0
    case = report["cases"]["case"]
    assert case["status"] == "diff"
    assert case["log"] == str((tmp_path / "logs" / "case_diff.log").resolve())
    assert Path(case["log"]).is_file()
    files = {item["name"]: item for item in case["files"]}
    stats = files["case_stats.nc"]
    assert stats["inputs"] == [str((left / "case_stats.nc").resolve()), str((right / "case_stats.nc").resolve())]
    assert stats["first_failing_prefix"] == {"record": 1, "saved_records": 2, "output_time": 120.0}
    assert stats["variables"] == {
        "all_zero": ["all_zero"],
        "identical": ["identical", "time"],
        "within_tolerance": ["within_tolerance"],
        "different": ["different"],
        "shape_mismatch": ["shape"],
        "non_numeric": ["text"],
        "empty": ["empty"],
        "only_left": ["only_left"],
        "only_right": ["only_right"],
    }
    categorized = [name for names in stats["variables"].values() for name in names]
    assert len(categorized) == len(set(categorized))
    assert files["case_zt.nc"]["comparison"] == "byte_identical"
    assert files["case_zt.nc"]["variables"] == {}


def test_json_records_no_comparable_files_error(tmp_path):
    left, right = tmp_path / "left", tmp_path / "right"
    left.mkdir()
    right.mkdir()
    report_path = tmp_path / "result.json"
    with pytest.raises(SystemExit) as exit_info:
        load_bindiff().find_diffs_in_all_files(
            str(left), str(right), None, 0, 1e-7, 1e-7, False,
            result_json=str(report_path),
        )
    assert exit_info.value.code == 2
    report = json.loads(report_path.read_text())
    assert report["status"] == "error"
    assert report["issues"] == ["no_comparable_files"]
    assert report["cases"] == {}


def test_json_records_time_length_mismatch(tmp_path):
    left, right = tmp_path / "left", tmp_path / "right"
    left.mkdir()
    right.mkdir()
    for folder, ntime in ((left, 2), (right, 3)):
        with netCDF4.Dataset(folder / "case_stats.nc", "w") as nc:
            nc.createDimension("time", ntime)
            nc.createVariable("signal", "f8", ("time",))[:] = np.ones(ntime)
    report_path = tmp_path / "result.json"
    _, differs, _, _, failed = load_bindiff().find_diffs_in_all_files(
        str(left), str(right), None, 0, 1e-7, 1e-7, False,
        result_json=str(report_path),
    )
    assert differs and failed == ["case"]
    file = json.loads(report_path.read_text())["cases"]["case"]["files"][0]
    assert file["status"] == "diff"
    assert file["issues"] == ["time_length_mismatch"]
    assert file["first_failing_prefix"] is None


def test_json_records_no_comparable_stats_error(tmp_path):
    left, right = tmp_path / "left", tmp_path / "right"
    left.mkdir()
    right.mkdir()
    for folder, value in ((left, 1.0), (right, 2.0)):
        with netCDF4.Dataset(folder / "case_stats.nc", "w") as nc:
            nc.createDimension("time", 1)
            nc.createVariable("time", "f8", ("time",))[:] = [value]
    report_path = tmp_path / "result.json"
    _, differs, _, _, failed = load_bindiff().find_diffs_in_all_files(
        str(left), str(right), None, 0, 1e-7, 1e-7, False,
        result_json=str(report_path),
    )
    assert differs and failed == ["case"]
    file = json.loads(report_path.read_text())["cases"]["case"]["files"][0]
    assert file["status"] == "diff"
    assert file["issues"] == ["no_comparable_stats"]


@pytest.mark.parametrize("problem", ["missing_input_directory", "same_input_directory"])
def test_cli_writes_json_for_invalid_inputs(tmp_path, problem):
    left = tmp_path / "left"
    left.mkdir()
    right = tmp_path / "missing" if problem == "missing_input_directory" else left
    report_path = tmp_path / "result.json"
    run = subprocess.run(
        [sys.executable, str(BINDIFF), "-result_json", str(report_path), str(left), str(right)],
        capture_output=True, text=True,
    )
    assert run.returncode == 2
    report = json.loads(report_path.read_text())
    assert report["status"] == "error"
    assert report["issues"] == [problem]


def test_first_failing_prefix_is_not_first_pointwise_difference(tmp_path):
    left, right = tmp_path / "left", tmp_path / "right"
    left.mkdir()
    right.mkdir()
    for directory, values in (
        (left, [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]),
        (right, [[2e-6, 0.0], [3e-6, 3e-6], [3e-6, 3e-6]]),
    ):
        with netCDF4.Dataset(directory / "case_stats.nc", "w") as nc:
            nc.createDimension("time", 3)
            nc.createDimension("column", 2)
            nc.createDimension("bounds", 2)
            nc.createVariable("time", "f8", ("time",))[:] = [90.0, 240.0, 450.0]
            nc.createVariable("time_bnds", "f8", ("time", "bounds"))[:, :] = np.array(
                [[0.0, 180.0], [180.0, 300.0], [300.0, 600.0]]
            )
            nc.createVariable("signal", "f8", ("time", "column"))[:, :] = np.array(values)

    report = tmp_path / "result.json"
    run = subprocess.run(
        [sys.executable, str(BINDIFF), "-verbose", "0", "-case", "case", "-threshold", "1e-6",
         "-percent_threshold", "1e-7", "-result_json", str(report), str(left), str(right)],
        capture_output=True, text=True,
    )

    assert run.returncode == 1, run.stdout + run.stderr
    result = json.loads(report.read_text())
    case = result["cases"]["case"]
    assert result["inputs"] == [str(left.resolve()), str(right.resolve())]
    assert result["status"] == case["status"] == "diff"
    assert case["files"][0]["first_failing_prefix"] == {
        "record": 1, "saved_records": 3, "output_time": 300.0
    }
    assert case["files"][0]["variables"]["different"] == ["signal"]
    assert _first_failing_timestep(report, "case", 10, (0.0, 60.0)) == 5
    assert _average_earliest_timestep(report, "case") == 0.0
    assert _average_earliest_timestep(report) == 0.0
