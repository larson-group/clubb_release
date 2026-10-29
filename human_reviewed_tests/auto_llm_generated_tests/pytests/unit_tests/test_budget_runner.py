"""Check the budget runner's worker forwarding and failure propagation."""

import importlib.util
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import pytest


SCRIPT = Path(__file__).resolve().parents[4] / "human_reviewed_tests/check_budget_balance.py"
spec = importlib.util.spec_from_file_location("budget_runner", SCRIPT)
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


@pytest.mark.parametrize("worker_args, expected", [([], None), (["-workers", "3"], "3")])
def test_workers_and_cases_forwarded(tmp_path, monkeypatch, worker_args, expected):
    """Forward the selected worker count and existing cases before checking budgets."""
    scripts = tmp_path / "run_scripts"
    scripts.mkdir()
    (scripts / "run_scm_all.py").touch()
    budget = tmp_path / "postprocessing/check_budgets_balance"
    budget.mkdir(parents=True)
    calls = []

    def fake_run(cmd, **kwargs):
        calls.append((cmd, kwargs))
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(sys, "argv", [str(SCRIPT), str(tmp_path), *worker_args])
    monkeypatch.setattr(runner.subprocess, "run", fake_run)
    assert runner.main() == 0
    command, options = calls[0]
    assert command[1] == str(scripts / "run_scm_all.py")
    if expected is None:
        assert "-workers" not in command
    else:
        assert command[command.index("-workers") + 1] == expected
    assert command[command.index("-cases") + 1].split(",") == runner.CASES
    assert command[command.index("-override") + 1] == "l_lmm_stepping=.false.,stats_tsamp=300.0"
    assert options["check"] is True
    assert calls[1][0][1:] == ["checkBudget.py", "all", "0"]
    assert calls[1][1]["cwd"] == str(budget)


def test_failed_cases_do_not_run_budget_checker(tmp_path, monkeypatch):
    """A failing concurrent case run stops the budget checker."""
    scripts = tmp_path / "run_scripts"
    scripts.mkdir()
    (scripts / "run_scm_all.py").touch()
    (tmp_path / "postprocessing/check_budgets_balance").mkdir(parents=True)
    calls = []

    def fake_run(cmd, **kwargs):
        calls.append(cmd)
        raise subprocess.CalledProcessError(1, cmd)

    monkeypatch.setattr(sys, "argv", [str(SCRIPT), str(tmp_path)])
    monkeypatch.setattr(runner.subprocess, "run", fake_run)
    with pytest.raises(subprocess.CalledProcessError):
        runner.main()
    assert len(calls) == 1
