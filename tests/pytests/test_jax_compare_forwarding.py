"""The JAX/Fortran harness forwards one model option set to both drivers."""

import sys

import pytest

from tests import run_jax_vs_fortran_cases as harness


def test_single_dash_options_forward_without_a_separator(monkeypatch):
    monkeypatch.setattr(sys, "argv", [
        "run_jax_vs_fortran_cases.py", "-cases", "rico", "-jobs", "2",
        "-max_iters", "3", "-dt_main", "60", "-stats", "input/stats/standard_stats.in",
    ])
    args = harness.parse_args()
    assert args.cases == ["rico"]
    assert args.jobs == 2
    assert args.run_scm_args == [
        "-max_iters", "3", "-dt_main", "60", "-stats", "input/stats/standard_stats.in"
    ]


@pytest.mark.parametrize("option", ["-out_dir", "-jax", "-override"])
def test_comparison_owned_model_options_are_rejected(monkeypatch, option):
    monkeypatch.setattr(sys, "argv", ["run_jax_vs_fortran_cases.py", "-cases", "rico", option, "value"])
    with pytest.raises(SystemExit) as exc:
        harness.parse_args()
    assert exc.value.code == 2


def test_both_runs_receive_the_same_forwarded_options(tmp_path):
    commands = []

    class Supervisor:
        def run_and_log(self, cmd, _cwd, log_path):
            commands.append(cmd)
            log_path.parent.mkdir(parents=True, exist_ok=True)
            if "-jax" in cmd:
                log_path.write_text("Completed 3 timesteps\n")
            elif cmd[1].endswith("run_scm.py"):
                log_path.write_text("iteration: 3 / 3 -- time = 180\n")
            return 0

    task = harness.TaskCtx(
        config=harness.CaseConfig("rico", max_iters=3, dt_main=60),
        flag_data=harness.FlagData("default", None),
        repo_root=tmp_path,
        run_scm_args=["-max_iters", "3", "-dt_main", "60", "-debug", "0"],
        bindiff_threshold=1e-7,
        bindiff_percent_threshold=1e-7,
        run_output_root=tmp_path / "results",
        supervisor=Supervisor(),
    )
    result = harness._run_case_w_flags(task)
    jax, fortran, bindiff = commands
    assert result.status == "match"
    assert "-strict" in bindiff
    assert jax.count("-max_iters") == fortran.count("-max_iters") == 1
    assert jax.count("-dt_main") == fortran.count("-dt_main") == 1
    start = jax.index("-multicol") + 2
    assert task.run_scm_args == jax[start:start + len(task.run_scm_args)]
    assert task.run_scm_args == fortran[start:start + len(task.run_scm_args)]
    assert "-jax" in jax and "-jax" not in fortran
    assert jax[-1] == fortran[-1] == "rico"
