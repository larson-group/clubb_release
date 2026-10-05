"""The JAX/Fortran harness forwards one model option set to both drivers."""

import sys

import pytest

from tests import run_jax_vs_fortran_cases as harness
from tests.pytests.jax_comparison_inputs import run_comparison


def test_single_dash_options_forward_without_a_separator(monkeypatch):
    monkeypatch.setattr(sys, "argv", [
        "run_jax_vs_fortran_cases.py", "-cases", "rico", "-workers", "2",
        "-max_iters", "3", "-dt_main", "60", "-stats", "input/stats/standard_stats.in",
    ])
    args = harness.parse_args()
    assert args.cases == ["rico"]
    assert args.jobs == 2
    assert args.run_scm_args == [
        "-max_iters", "3", "-dt_main", "60", "-stats", "input/stats/standard_stats.in"
    ]


@pytest.mark.parametrize("option", ["-output_dir", "-jax", "-override"])
def test_comparison_owned_model_options_are_rejected(monkeypatch, option):
    monkeypatch.setattr(sys, "argv", ["run_jax_vs_fortran_cases.py", "-cases", "rico", option, "value"])
    with pytest.raises(SystemExit) as exc:
        harness.parse_args()
    assert exc.value.code == 2


def test_both_runs_receive_the_same_forwarded_options(tmp_path):
    config = harness.CaseConfig("rico", max_iters=3, dt_main=60)
    result, commands, task = run_comparison(tmp_path, config)
    jax, fortran, bindiff = commands
    assert result.status == "match"
    assert "-strict" in bindiff
    assert bindiff[bindiff.index("-percent_threshold") + 1] == str(config.percent_threshold)
    assert jax.count("-max_iters") == fortran.count("-max_iters") == 1
    assert jax.count("-dt_main") == fortran.count("-dt_main") == 1
    start = jax.index("-multicol") + 2
    assert task.run_scm_args == jax[start:start + len(task.run_scm_args)]
    assert task.run_scm_args == fortran[start:start + len(task.run_scm_args)]
    assert "-jax" in jax and "-jax" not in fortran
    assert jax[-1] == fortran[-1] == config.source_case


def test_comparison_rejects_empty_explicit_case_list():
    with pytest.raises(SystemExit) as error:
        harness.parse_args(["-cases", ","])
    assert error.value.code == 2
