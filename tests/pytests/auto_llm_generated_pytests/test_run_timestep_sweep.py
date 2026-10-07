"""Exercise timestep selection, runner forwarding, and failure reporting."""

import ast
import re
import subprocess
import sys
import threading

import pytest

from tests import run_timestep_sweep as harness


def test_native_defaults_preserve_the_existing_sweep(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["run_timestep_sweep.py"])
    args = harness.parse_args()
    assert args.cases == list(harness.DEFAULT_CASES)
    assert args.timesteps == [600, 1200, 1800, 2400, 3000]
    assert not args.jax
    assert args.max_iters is None
    assert args.out_root is None
    assert args.workers == 4


@pytest.mark.parametrize("jax_option", [
    "-jax", "-jax=cpu", "-jax=gpu,xla_prealloc",
    "-jax=gpu,device=GPU-aaaaaaaa-bbbb-cccc-dddd-eeeeeeeeeeee,prealloc_gpu_mem=false",
])
def test_selected_timesteps_forward_to_jax(monkeypatch, tmp_path, jax_option):
    monkeypatch.setattr(sys, "argv", [
        "run_timestep_sweep.py", jax_option, "-cases", "bomex,atex",
        "-timesteps", "60,300", "-max_iters", "8", "-workers", "1",
        "-output_root", str(tmp_path),
    ])
    calls = []

    def run(cmd, **kwargs):
        calls.append((cmd, kwargs))
        return subprocess.CompletedProcess(cmd, 0, stdout="completed\n")

    monkeypatch.setattr(harness.subprocess, "run", run)
    assert harness.main() == 0
    assert len(calls) == 4
    for (case, dt), (cmd, kwargs) in zip(
        [("bomex", 60), ("bomex", 300), ("atex", 60), ("atex", 300)], calls,
    ):
        assert cmd[-1] == case
        assert cmd.count(jax_option) == 1
        assert cmd[cmd.index("-dt_main") + 1] == str(dt)
        assert cmd[cmd.index("-dt_rad") + 1] == str(dt)
        assert cmd[cmd.index("-max_iters") + 1] == "8"
        assert cmd[cmd.index("-tout") + 1] == "0"
        assert cmd[cmd.index("-output_dir") + 1] == str(tmp_path / case / f"dt{dt}")
        assert kwargs["cwd"] == harness.CLUBB_ROOT


@pytest.mark.parametrize('backend', [[], ['-jax']])
def test_failure_reports_stability_limit_without_failing_sweep(monkeypatch, capsys, backend):
    monkeypatch.setattr(sys, "argv", [
        "run_timestep_sweep.py", *backend, "-cases", "bomex,atex",
        "-timesteps", "60,120,300", "-workers", "1",
    ])
    calls = []

    def run(cmd, **_kwargs):
        assert '-output_dir' not in cmd
        calls.append((cmd[-1], cmd[cmd.index("-dt_main") + 1]))
        returncode = 1 if calls[-1] == ("bomex", "120") else 0
        return subprocess.CompletedProcess(cmd, returncode, stdout="model failure details\n")

    monkeypatch.setattr(harness.subprocess, "run", run)
    assert harness.main() == 0
    assert calls == [
        ("bomex", "60"), ("bomex", "120"),
        ("atex", "60"), ("atex", "120"), ("atex", "300"),
    ]
    output = capsys.readouterr().out
    assert "model failure details" in output
    assert "FAIL @ dt = 120" in output


@pytest.mark.parametrize("options", [
    ["-cases", ""], ["-timesteps", ""], ["-timesteps", "60,0"],
    ["-timesteps", "60,abc"], ["-max_iters", "0"], ["-jax", "-jax=cpu"],
    ["-workers", "0"], ["-cases", "bomex,bomex"],
])
def test_invalid_selection_is_rejected(monkeypatch, options):
    monkeypatch.setattr(sys, "argv", ["run_timestep_sweep.py", *options])
    with pytest.raises(SystemExit) as exc:
        harness.parse_args()
    assert exc.value.code == 2


@pytest.mark.parametrize("value", ["timestep-audit", "output/timestep-audit"])
def test_output_labels_use_the_shared_output_owner(monkeypatch, value):
    monkeypatch.setattr(sys, "argv", ["run_timestep_sweep.py", "-output_root", value])
    assert harness.parse_args().out_root == harness.CLUBB_ROOT / "output" / "timestep-audit"


@pytest.mark.parametrize("options", [["-output_root", "../outside"], ["-time", "60"]])
def test_invalid_paths_and_abbreviated_options_are_rejected(monkeypatch, options):
    monkeypatch.setattr(sys, "argv", ["run_timestep_sweep.py", *options])
    with pytest.raises(SystemExit) as exc:
        harness.parse_args()
    assert exc.value.code == 2


def test_jenkins_jax_sweep_uses_shared_cases_and_native_run_defaults():
    """Keep CI selection aligned without importing comparison/JAX dependencies."""
    comparison = ast.parse((harness.CLUBB_ROOT / "tests" / "run_jax_vs_fortran_cases.py").read_text())
    assignment = next(
        node for node in comparison.body
        if isinstance(node, ast.Assign)
        and any(isinstance(target, ast.Name) and target.id == "DEFAULT_CASES" for target in node.targets)
    )
    comparison_names = {ast.literal_eval(call.args[0]) for call in assignment.value.elts}
    expected = [case for case in harness.DEFAULT_CASES if case in comparison_names]
    pipeline = (harness.CLUBB_ROOT / "jenkins_tests" / "clubb_timestep" / "Jenkinsfile").read_text()
    block = pipeline.split("stage('JAX timestep sweep')", 1)[1].split("stage('Fortran timestep sweep')", 1)[0]
    selected = re.search(r'-cases "(.*?)"', block, re.S).group(1)
    assert [case.strip() for case in selected.split(",")] == expected
    assert "-timesteps" not in block and "-max_iters" not in block
    assert "run_timestep_sweep.py -jax=cpu" in block
    assert "-workers 4" in block


def test_parallel_cases_keep_timestep_order_and_reports_together(monkeypatch, capsys):
    """All four jobs overlap, but only one timestep per case runs at once."""
    monkeypatch.setattr(sys, "argv", [
        "run_timestep_sweep.py", "-cases", "bomex,atex,arm,fire",
        "-timesteps", "600,1200", "-workers", "4",
    ])
    barrier = threading.Barrier(4)
    calls = {case: [] for case in ["bomex", "atex", "arm", "fire"]}

    def run(cmd, **kwargs):
        case, dt = cmd[-1], int(cmd[cmd.index("-dt_main") + 1])
        calls[case].append(dt)
        if dt == 600:
            barrier.wait(timeout=5)
        return subprocess.CompletedProcess(cmd, 0, stdout="done\n")

    monkeypatch.setattr(harness.subprocess, "run", run)
    assert harness.main() == 0
    assert all(values == [600, 1200] for values in calls.values())
    output = capsys.readouterr().out
    for case in calls:
        assert (
            f"---------------- Running {case} ----------------\n"
            "--- PASS @ dt = 600\n--- PASS @ dt = 1200\n"
        ) in output


def test_expensive_cases_are_scheduled_first(monkeypatch):
    monkeypatch.setattr(sys, "argv", [
        "run_timestep_sweep.py", "-cases", "bomex,lba,rico,atex", "-workers", "1",
    ])
    cases = []
    monkeypatch.setattr(harness, "run_case", lambda case, args, stop: cases.append(case) or case)
    assert harness.main() == 0
    assert cases == ["rico", "lba", "bomex", "atex"]


def test_stop_request_prevents_another_model_launch(monkeypatch):
    monkeypatch.setattr(sys, "argv", [
        "run_timestep_sweep.py", "-cases", "bomex", "-timesteps", "600,1200",
    ])
    args = harness.parse_args()
    stop = threading.Event()
    calls = []

    def run(cmd, **kwargs):
        calls.append(cmd)
        stop.set()  # An interruption arrives while the current model finishes.
        return subprocess.CompletedProcess(cmd, 0, stdout="done\n")

    monkeypatch.setattr(harness.subprocess, "run", run)
    report = harness.run_case("bomex", args, stop)
    assert len(calls) == 1
    assert "PASS @ dt = 600" in report and "1200" not in report
