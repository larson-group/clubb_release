"""Restart CLI forwarding and all-column bit-for-bit comparisons."""

import numpy as np
import pytest
from netCDF4 import Dataset

from tests import run_restart_test as runner


def _write_output(path, values, times=(300.0,)):
    with Dataset(path, "w") as ds:
        ds.createDimension("time", len(times))
        ds.createVariable('time', 'f8', ('time',))[:] = times
        ds.createDimension("zt", values.shape[0])
        ds.createDimension("col", values.shape[1])
        ds.createVariable("thlm", "f8", ("time", "zt", "col"))[:] = values


def test_final_comparison_checks_every_column(tmp_path, monkeypatch):
    output = tmp_path / "output"
    restart = tmp_path / "restart"
    output.mkdir()
    restart.mkdir()
    monkeypatch.setattr(runner, "OUTPUT_DIR", output)
    monkeypatch.setattr(runner, "RESTART_DIR", restart)
    values = np.array([[1.0, 2.0], [3.0, 4.0]])
    _write_output(restart / "bomex_stats.nc", values)
    _write_output(output / "bomex_stats.nc", values)
    assert runner._compare_final_timestep("bomex", "thlm", all_columns=True)
    values[0, 1] += 1.0
    _write_output(output / "bomex_stats.nc", values)
    assert not runner._compare_final_timestep("bomex", "thlm", all_columns=True)
    assert runner._compare_final_timestep("bomex", "thlm")


@pytest.mark.parametrize("argv,option", [
    (["bomex", "-jax"], "-jax"),
    (["-jax", "bomex"], "-jax"),
    (["bomex", "-jax=gpu"], "-jax=gpu"),
    (["bomex"], None),
])
def test_cli_selects_both_runners_and_uses_effective_duration(tmp_path, monkeypatch, argv, option):
    output = tmp_path / "output"
    restart = tmp_path / "restart"
    model = tmp_path / "input" / "case_setups" / "bomex_model.in"
    model.parent.mkdir(parents=True)
    model.write_text("time_initial=0.\ntime_final=21600.\n")
    monkeypatch.setattr(runner, "CLUBB_ROOT", tmp_path)
    monkeypatch.setattr(runner, "OUTPUT_DIR", output)
    monkeypatch.setattr(runner, "RESTART_DIR", restart)
    monkeypatch.setattr(runner.sys, "argv", ["run_restart_test.py", *argv, "-max_iters", "10"])
    calls = []

    def run_scm(case, options, override=None):
        calls.append((case, options[:], override))
        output.mkdir(exist_ok=True)
        (output / "bomex.in").write_text("time_initial=0.\ntime_final=600.\n")
        _write_output(output / "bomex_stats.nc", np.ones((2, 2)))
        return 0, ""

    monkeypatch.setattr(runner, "_run_scm", run_scm)
    assert runner.main() == 0
    assert len(calls) == 2
    for _, options, _ in calls:
        assert options == ([option] if option else []) + ["-max_iters", "10"]
    assert calls[0][2] == "l_restart=.false."
    assert "time_restart = 300.0" in calls[1][2]
    assert not restart.exists()


@pytest.mark.parametrize('times,expected', [
    ([60., 120., 180., 240., 300.], 120.),
    ([120., 240.], 120.),
    ([60., 120., 180., 240., 300., 360.], 180.),
])
def test_restart_uses_saved_interior_midpoint(tmp_path, monkeypatch, times, expected):
    monkeypatch.setattr(runner, 'OUTPUT_DIR', tmp_path)
    _write_output(tmp_path / 'bomex_stats.nc', np.ones((2, 1)), times)
    assert runner._saved_restart_time('bomex', 0., times[-1]) == expected


def test_split_restart_time_is_present_on_every_grid(tmp_path, monkeypatch):
    monkeypatch.setattr(runner, 'OUTPUT_DIR', tmp_path)
    for grid, times in [('zt', [60., 120., 180., 240.]),
                        ('zm', [120., 240.]), ('sfc', [120., 240.])]:
        _write_output(tmp_path / f'bomex_{grid}.nc', np.ones((2, 1)), times)
    assert runner._saved_restart_time('bomex', 0., 240.) == 120.


def test_restart_without_interior_output_has_clear_error(tmp_path, monkeypatch):
    monkeypatch.setattr(runner, 'OUTPUT_DIR', tmp_path)
    _write_output(tmp_path / 'bomex_stats.nc', np.ones((2, 1)), [300.])
    with pytest.raises(ValueError, match='strictly inside'):
        runner._saved_restart_time('bomex', 0., 300.)
