"""SCM forwards opaque JAX values; only the wrapper validates their meaning."""
from utilities.output_paths import REPO_ROOT as _REPO_ROOT

import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

from clubb_jax import run_jax
from clubb_jax.backends import metal, rocm
from run_scripts import run_scm


@pytest.mark.parametrize("value", [None, "", "gpu,xla_prealloc", "FUTURE,option=some value"])
@pytest.mark.parametrize("case_first", [False, True])
def test_scm_main_passes_opaque_options_to_case_launch(value, case_first, tmp_path, monkeypatch):
    token = "-jax" if value is None else f"-jax={value}"
    args = ["arm", token] if case_first else [token, "arm"]
    monkeypatch.setattr(sys, "argv", ["run_scm.py", *args, "-stats", "none"])
    monkeypatch.setattr(run_scm, "resolve_output_dir", lambda _: tmp_path)
    namelist = str(tmp_path / "arm.in")
    monkeypatch.setattr(run_scm, "create_case_namelist", lambda *_: namelist)
    captured = []

    def launch(command, cwd, case, input_file, *args, **kwargs):
        captured.append(command)
        assert case == "arm" and input_file == namelist
        return 0

    monkeypatch.setattr(run_scm, "run_case", launch)
    assert run_scm.main() == 0
    assert captured[0][1:] == ([] if value is None else [f"-options={value}"])


@pytest.mark.parametrize("value", [None, "", "cpu", "gpu,xla_prealloc",
    "FUTURE,option=some value", "gpu,unknown", "literal;$NOT_EXPANDED", "gpu\n,typo"])
def test_runner_forwards_value_verbatim_without_interpreting_it(value, tmp_path, monkeypatch):
    # No JAX driver/runtime is needed by the SCM runner, only the wrapper entrypoint.
    wrapper = tmp_path / "clubb_jax" / "run_jax.py"
    wrapper.parent.mkdir()
    wrapper.touch()
    monkeypatch.setattr(run_scm, "CLUBB_ROOT", str(tmp_path))
    command, _, _ = run_scm.choose_run_command(SimpleNamespace(
        exe=None, python=False, jax=True, jax_options=value, gdb=False,
    ))
    assert command == [str(wrapper)] + ([] if value is None else [f"-options={value}"])


@pytest.mark.parametrize("value", ["", "tpu", "gpu,", "gpu,typo", "cpu,xla_prealloc",
    "gpu,xla_prealloc,xla_prealloc", "gpu,,xla_prealloc", "gpu\n,typo"])
def test_bad_options_are_rejected_by_wrapper_before_setup(value, tmp_path):
    import os

    wrapper = _REPO_ROOT / "clubb_jax" / "run_jax.py"
    env = os.environ | {"CLUBB_JAX_VENV": str(tmp_path / "venv"),
                        "CLUBB_JAX_TOOLS_DIR": str(tmp_path / "tools")}
    result = subprocess.run([str(wrapper), f"-options={value}", "-init_env"],
                            env=env, capture_output=True, text=True, timeout=10)
    assert result.returncode != 0
    assert "ERROR:" in result.stderr
    assert not (tmp_path / "venv").exists()
    assert not (tmp_path / "tools").exists()


@pytest.mark.parametrize("options", [
    ["-options=cpu", "-profile=gpu"],
    ["-options=cpu", "-options=gpu"],
    ["-options=gpu,xla_prealloc", "-xla_prealloc"],
])
def test_wrapper_rejects_conflicting_or_duplicate_settings(options):
    wrapper = _REPO_ROOT / "clubb_jax" / "run_jax.py"
    result = subprocess.run([str(wrapper), *options, "-info=json"],
                            capture_output=True, text=True, timeout=10)
    assert result.returncode != 0
    assert "ERROR:" in result.stderr


def test_wrapper_rejects_preallocation_on_cpu():
    wrapper = _REPO_ROOT / "clubb_jax" / "run_jax.py"
    result = subprocess.run([str(wrapper), "-profile=cpu", "-xla_prealloc", "-info=json"],
                            capture_output=True, text=True, timeout=10)
    assert result.returncode != 0
    assert "CUDA-only" in result.stderr


@pytest.mark.parametrize(
    ("system", "expected"),
    [("Darwin", "metal"), ("Linux", "cuda13")],
)
def test_gpu_profile_resolves_to_the_host_native_backend(system, expected, monkeypatch):
    monkeypatch.setattr(metal.platform, "system", lambda: system)
    monkeypatch.setattr(rocm, "has_gpu", lambda: False)
    values, driver_args = run_jax.parse_launcher_args(["-profile=gpu", "arm.in"])

    accelerator, profile = run_jax.resolve_accelerator(values)

    assert (accelerator, profile) == (expected, "gpu")
    assert driver_args == ["arm.in"]


@pytest.mark.parametrize("host_version", [(3, 9), (3, 10)])
@pytest.mark.parametrize("accelerator,required_jax", [("cpu", "0.11.0"), ("metal", "0.4.34")])
def test_discovery_without_target_venv_uses_old_host_python(
    tmp_path, monkeypatch, host_version, accelerator, required_jax
):
    monkeypatch.setattr(run_jax.sys, "version_info", host_version)
    monkeypatch.setattr(run_jax.sys, "executable", "/host/python")
    monkeypatch.setattr(run_jax, "_find_python", lambda _: None)
    venv = tmp_path / "missing-venv"
    assert run_jax._inspection_python(venv, accelerator) == (
        Path("/host/python"), "3.12", required_jax
    )
    assert not venv.exists()


def test_discovery_plans_packages_for_target_but_executes_on_host(tmp_path, monkeypatch):
    target = tmp_path / "python3.11"
    monkeypatch.setattr(run_jax.sys, "executable", "/host/python")
    monkeypatch.setattr(run_jax, "_find_python", lambda _: target)
    monkeypatch.setattr(run_jax, "_python_version", lambda _: (3, 11))
    assert run_jax._inspection_python(tmp_path / "missing", "cpu") == (
        Path("/host/python"), "3.11", "0.10.0"
    )


def test_unsuitable_setup_python_does_not_block_discovery(tmp_path, monkeypatch):
    def reject_override(_):
        raise run_jax.LauncherError("Python override is too old for JAX")

    monkeypatch.setattr(run_jax, "_find_python", reject_override)
    python, planned, jax = run_jax._inspection_python(tmp_path / "missing", "metal")
    assert python == Path(sys.executable)
    assert (planned, jax) == ("3.12", "0.4.34")


@pytest.mark.parametrize("options", [
    ["-jax", "-jax=gpu"], ["-jax=cpu", "-jax"],
    ["-jax", "-python"], ["-jax=cpu", "-exe", "model"],
])
def test_scm_rejects_duplicate_or_conflicting_implementations(options, monkeypatch):
    monkeypatch.setattr(sys, "argv", ["run_scm.py", *options, "arm"])
    with pytest.raises(SystemExit) as exc:
        run_scm.main()
    assert exc.value.code == 2


def test_comparison_help_does_not_prepare_environment(tmp_path):
    import os

    root = _REPO_ROOT
    env = os.environ | {"CLUBB_JAX_VENV": str(tmp_path / "venv"),
                        "CLUBB_JAX_TOOLS_DIR": str(tmp_path / "tools")}
    result = subprocess.run([sys.executable, str(root / "tests/run_jax_vs_fortran_cases.py"), "-h"],
                            env=env, capture_output=True, text=True, timeout=10)
    assert result.returncode == 0, result.stderr
    assert "-cases" in result.stdout
    assert not (tmp_path / "venv").exists()
    assert not (tmp_path / "tools").exists()


@pytest.mark.parametrize("accelerator", ["cpu", "cuda13", "rocm", "metal"])
@pytest.mark.parametrize("custom_venv", [None, "custom-runtime"])
def test_environment_setup_uses_launcher_paths_and_restarts_once(
    accelerator, custom_venv, tmp_path, monkeypatch
):
    marker = "_CLUBB_JAX_ENVIRONMENT_PYTHON"
    monkeypatch.delenv(marker, raising=False)
    monkeypatch.setenv("CLUBB_JAX_ACCELERATOR", accelerator)
    monkeypatch.delenv("CLUBB_JAX_VENV", raising=False)
    if custom_venv:
        monkeypatch.setenv("CLUBB_JAX_VENV", custom_venv)
    monkeypatch.setattr(run_jax, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(sys, "argv", ["compare.py", "-cases", "arm"])
    monkeypatch.setattr(sys, "executable", "/host/python")
    setup = []
    launched = []

    def initialize(command):
        setup.append(command)
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(run_jax.subprocess, "run", initialize)
    monkeypatch.setattr(run_jax.os, "execve", lambda *args: launched.append(args))
    run_jax.ensure_environment()
    assert setup == [[str(run_jax.SCRIPT_DIR / "run_jax.py"), "-init_env"]]
    default = {"cpu": ".venv-jax", "cuda13": ".venv-jax-cuda13", "rocm": ".venv-jax-rocm", "metal": ".venv-jax-metal"}
    python = str(tmp_path / (custom_venv or default[accelerator]) / "bin/python")
    assert launched[0][:2] == (python, [python, *sys.argv])
    assert launched[0][2][marker] == python
    monkeypatch.setattr(sys, "executable", python)
    run_jax.ensure_environment()
    assert len(setup) == len(launched) == 1


def test_environment_setup_failure_stops_before_restart(monkeypatch):
    monkeypatch.delenv("_CLUBB_JAX_ENVIRONMENT_PYTHON", raising=False)
    monkeypatch.setattr(run_jax.subprocess, "run", lambda *args: SimpleNamespace(returncode=7))
    monkeypatch.setattr(run_jax.os, "execve", lambda *args: pytest.fail("must not restart"))
    with pytest.raises(SystemExit) as exc:
        run_jax.ensure_environment()
    assert exc.value.code == 7


@pytest.mark.parametrize("module", [
    "clubb_jax.src.clubb_standalone",
    "tuner.tune_clubb",
    "clubb_jax.src.clubb_standalone_loss",
])
def test_managed_launcher_dispatches_selected_entry_module(module, tmp_path, monkeypatch):
    python = tmp_path / "bin" / "python"
    calls = []
    monkeypatch.setattr(run_jax, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(run_jax, "_prepare_environment", lambda *args: python)
    monkeypatch.setattr(run_jax, "_print_runtime_summary", lambda *args: None)
    monkeypatch.setattr(run_jax.os, "chdir", lambda *args: None)
    monkeypatch.setattr(run_jax.os, "execvpe", lambda *args: calls.append(args))
    assert run_jax.main(["-options=cpu", f"-module={module}", "run.in"]) == 0
    executable, command, env = calls[0]
    assert executable == str(python)
    assert command == [str(python), "-m", module, "run.in"]
    assert env["_CLUBB_JAX_ENVIRONMENT_PYTHON"] == str(python)


@pytest.mark.parametrize("options", [
    ["-module=unrecognized"],
    ["-module=tuner.tune_clubb", "-module=clubb_jax.src.clubb_standalone_loss"],
])
def test_managed_launcher_rejects_bad_entry_modules(options):
    with pytest.raises(run_jax.LauncherError, match="module"):
        run_jax.parse_launcher_args(options)
