"""Shared CLI defaults, model forwarding, and comparison command behavior."""

import argparse
from pathlib import Path

import pytest

from tuner import system_defaults
from utilities.create_case_namelist import add_namelist_arguments, parse_forwarded_args
from run_scripts import run_clubb_w_varying_flags as flags
import importlib.util
import sys

_spec = importlib.util.spec_from_file_location(
    "arg_python_comparison", Path(__file__).resolve().parents[2] / "tests/run_python_vs_fortran_cases.py"
)
comparison = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = comparison
_spec.loader.exec_module(comparison)


@pytest.mark.parametrize("cpus, expected", [(1, 1), (2, 1), (3, 1), (8, 4), (32, 16)])
def test_worker_default_is_half_available_logical_cpus(monkeypatch, cpus, expected):
    monkeypatch.setattr(system_defaults.os, "sched_getaffinity", lambda _: set(range(cpus)), raising=False)
    assert system_defaults.default_max_workers() == expected


def test_worker_default_has_a_portable_fallback(monkeypatch):
    def unavailable(_):
        raise OSError("affinity unavailable")
    monkeypatch.setattr(system_defaults.os, "sched_getaffinity", unavailable, raising=False)
    monkeypatch.setattr(system_defaults.os, "cpu_count", lambda: 7)
    assert system_defaults.default_max_workers() == 3


def test_model_values_are_not_consumed_as_the_optional_case():
    parser = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
    parser.add_argument("case", nargs="?")
    parser.add_argument("-workers", type=int)
    args, forwarded = parse_forwarded_args(parser, [
        "-max_iters", "5", "-stats_tstart", "-1", "-stats_tend", "120",
        "-override", '{"all":{"thlm_sponge_damp_settings%l_sponge_damping":false}}',
        "-workers", "2", "bomex", "-override", "time_final=60",
    ])
    assert args.case == "bomex" and args.workers == 2
    assert forwarded == ["-max_iters", "5", "-stats_tstart", "-1", "-stats_tend", "120",
                         "-override", '{"all":{"thlm_sponge_damp_settings%l_sponge_damping":false}}',
                         "-override", "time_final=60"]
    model = argparse.ArgumentParser(add_help=False)
    add_namelist_arguments(model)
    parsed = model.parse_args(forwarded)
    assert parsed.stats_tstart == -1 and len(parsed.override) == 2


def test_varying_flags_forwards_model_settings_and_retains_override_precedence(monkeypatch, tmp_path):
    monkeypatch.setattr(flags.sys, "argv", ["flags", "-max_iters", "5", "-tout", "0",
                                          "-override", "debug_level=1", "bomex"])
    args = flags.get_cli_args()
    tasks = flags.build_tasks(str(tmp_path), {"default": {"debug_level": 0}}, ["bomex"], args)
    command = tasks[0]["cmd"]
    assert args.case_name == "bomex"
    assert command[-1] == "bomex"
    assert command[command.index("-max_iters") + 1] == "5"
    assert command[command.index("-tout") + 1] == "0"
    assert command[-3:-1] == ["-override", "debug_level=1"]


def test_python_comparison_forwards_settings_to_both_backends(monkeypatch, tmp_path):
    args = comparison.parse_args(["-cases", "bomex,rico", "arm", "-workers", "2", "-debug", "0"])
    assert args.cases == ["bomex", "rico", "arm"]
    assert args.run_scm_args == ("-debug", "0")
    commands = []
    monkeypatch.setattr(comparison, "_run_and_log", lambda cmd, *rest: commands.append(cmd) or 0)
    comparison._run_case("bomex", tmp_path, Path("stats.in"), 5, 1e-7,
                         tmp_path / "python", tmp_path / "fortran", tmp_path,
                         args.run_scm_args)
    for command in commands[:2]:
        assert command[command.index("-debug") + 1] == "0"
        assert command[command.index("-max_iters") + 1] == "5"
    assert "-python" in commands[0] and "-python" not in commands[1]


@pytest.mark.parametrize("option", ["-python", "-multicol=8", "-output_dir=elsewhere"])
def test_python_comparison_keeps_control_of_backend_and_output(option):
    with pytest.raises(SystemExit):
        comparison.parse_args([option])


def test_branch_comparison_translates_only_option_tokens_for_old_clones(tmp_path):
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "arg_branch_compare", Path(__file__).resolve().parents[2] / "tests/run_bindiff_w_flags.py"
    )
    branch = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(branch)
    scripts = ["run_scripts/run_clubb_w_varying_flags.py", "run_scripts/run_scm.py"]
    for file, declarations in zip(scripts, [
        ["--flag-config-file", "--priority-cases", "-nproc"],
        ["-params", "-out_dir", "-override"],
    ]):
        path = tmp_path / file
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("import argparse\np = argparse.ArgumentParser()\n" +
                        "\n".join(f"p.add_argument({name!r})" for name in declarations))
    arguments = ["-flag_config_file", "flags.json", "-priority_cases", "-workers=2",
                 "-params_file", "parameters.in", "-output_dir=/tmp/output",
                 "-override", '{"name":"-workers"}']
    assert branch.clone_arguments(tmp_path, arguments, scripts) == [
        "--flag-config-file", "flags.json", "--priority-cases", "-nproc=2",
        "-params", "parameters.in", "-out_dir=/tmp/output", "-override", '{"name":"-workers"}',
    ]
    current = Path(__file__).resolve().parents[2]
    assert branch.clone_arguments(current, arguments, [*scripts, "utilities/create_case_namelist.py"]) == arguments


def test_python_comparison_rejects_empty_explicit_case_list():
    with pytest.raises(SystemExit) as error:
        comparison.parse_args(["-cases", ","])
    assert error.value.code == 2
