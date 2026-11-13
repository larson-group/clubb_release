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
    calls = []
    checks = []

    def fake_run(cmd, **kwargs):
        calls.append((cmd, kwargs))
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(sys, "argv", [str(SCRIPT), str(tmp_path), *worker_args])
    monkeypatch.setattr(runner.subprocess, "run", fake_run)
    monkeypatch.setattr(runner, "check_budgets", lambda *args: checks.append(args) or 0)
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
    assert command[command.index("-output_dir") + 1] == str(tmp_path / "output")
    assert checks == [(tmp_path / "output", 0)]


def test_failed_cases_do_not_run_budget_checker(tmp_path, monkeypatch):
    """A failing concurrent case run stops the budget checker."""
    scripts = tmp_path / "run_scripts"
    scripts.mkdir()
    (scripts / "run_scm_all.py").touch()
    calls = []

    def fake_run(cmd, **kwargs):
        calls.append(cmd)
        raise subprocess.CalledProcessError(1, cmd)

    monkeypatch.setattr(sys, "argv", [str(SCRIPT), str(tmp_path)])
    monkeypatch.setattr(runner.subprocess, "run", fake_run)
    monkeypatch.setattr(runner, "check_budgets", lambda *args: pytest.fail("Checked failed case output"))
    with pytest.raises(subprocess.CalledProcessError):
        runner.main()
    assert len(calls) == 1


def test_budget_failure_propagates(tmp_path, monkeypatch):
    """Return the in-process checker failure after a successful SCM run."""
    scripts = tmp_path / "run_scripts"
    scripts.mkdir()
    (scripts / "run_scm_all.py").touch()
    monkeypatch.setattr(sys, "argv", [str(SCRIPT), str(tmp_path)])
    monkeypatch.setattr(runner.subprocess, "run", lambda *args, **kwargs: SimpleNamespace(returncode=0))
    monkeypatch.setattr(runner, "check_budgets", lambda *args: 1)
    assert runner.main() == 1


def test_check_only_selects_output_and_iteration(tmp_path, monkeypatch):
    """Check existing output without requiring a source tree or launching cases."""
    checks = []
    monkeypatch.setattr(sys, "argv", [str(SCRIPT), "-check_only", "-output_dir", str(tmp_path),
                                    "-iteration", "2"])
    monkeypatch.setattr(runner.subprocess, "run", lambda *args, **kwargs: pytest.fail("Ran SCM cases"))
    monkeypatch.setattr(runner, "check_budgets", lambda *args: checks.append(args) or 0)
    assert runner.main() == 0
    assert checks == [(tmp_path, 2)]
