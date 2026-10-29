"""Backend policy stays shared between launcher setup and read-only inspection."""

import hashlib
import importlib
import json
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

from clubb_jax import backends, run_jax, runtime_info
from clubb_jax.backends.common import GpuSelection


@pytest.mark.parametrize(
    "accelerator,requirements,venv,platform",
    [
        ("cpu", "requirements.txt", ".venv-jax", "cpu"),
        ("cuda13", "requirements-cuda13.txt", ".venv-jax-cuda13", "cuda,cpu"),
        ("rocm", "requirements-rocm.txt", ".venv-jax-rocm", "rocm,cpu"),
        ("metal", "requirements-metal.txt", ".venv-jax-metal", "METAL,cpu"),
    ],
)
def test_backend_runtime_paths(accelerator, requirements, venv, platform):
    paths = run_jax.runtime_paths(accelerator)
    assert paths == (
        run_jax.SCRIPT_DIR / requirements,
        run_jax.REPO_ROOT / venv,
        platform,
    )


@pytest.mark.parametrize("accelerator", backends.BACKENDS)
@pytest.mark.parametrize("inspect_only", [False, True])
def test_launcher_dispatches_backend_hooks(
    accelerator, inspect_only, tmp_path, monkeypatch
):
    backend = backends.get_backend(accelerator)
    events = []
    python = tmp_path / "bin" / "python"
    monkeypatch.delenv("CLUBB_JAX_VENV", raising=False)
    monkeypatch.delenv("CLUBB_JAX_TOOLS_DIR", raising=False)
    monkeypatch.setattr(run_jax, "REPO_ROOT", tmp_path)

    def configure(env, tools_dir, preallocate):
        assert env["JAX_PLATFORMS"] == backend.PLATFORM
        assert env["CLUBB_JAX_PROFILE"] == backend.PROFILE
        assert tools_dir == tmp_path / ".clubb-jax-tools"
        assert preallocate is False
        env["TEST_BACKEND"] = accelerator
        events.append("configure")

    def inspect(*args, **kwargs):
        assert args[1] == accelerator
        assert args[-1]["TEST_BACKEND"] == accelerator
        events.append("inspect")
        return 0

    def prepare_run(env):
        assert env["TEST_BACKEND"] == accelerator
        events.append("prepare_run")

    def prepare_environment(*args):
        assert args[0] == accelerator
        assert args[-1]["TEST_BACKEND"] == accelerator
        events.append("prepare_environment")
        return python

    def verify(executable, selected, env):
        assert executable == python
        assert selected == accelerator
        assert env["TEST_BACKEND"] == accelerator
        events.append("verify")

    monkeypatch.setattr(backend, "configure_environment", configure, raising=False)
    monkeypatch.setattr(backend, "prepare_run", prepare_run, raising=False)
    monkeypatch.setattr(run_jax, "_run_inspection", inspect)
    monkeypatch.setattr(run_jax, "_prepare_environment", prepare_environment)
    monkeypatch.setattr(run_jax, "_print_runtime_summary", lambda *args: None)
    monkeypatch.setattr(run_jax, "_verify_backend", verify)
    option = "-info=json" if inspect_only else "-init_env"
    assert run_jax.main([f"-accelerator={accelerator}", option]) == 0
    if inspect_only:
        assert events == ["configure", "inspect"]
    else:
        assert events == [
            "configure",
            *(["inspect"] if backend.PROFILE == "gpu" else []),
            "prepare_run",
            "prepare_environment",
            "verify",
        ]


@pytest.mark.parametrize("accelerator", backends.BACKENDS)
def test_setup_and_inspection_share_expected_packages(
    accelerator, tmp_path, monkeypatch
):
    backend = backends.get_backend(accelerator)
    expected = {"test-backend-package": "1.2.3"}
    monkeypatch.setattr(backend, "expected_packages", lambda _: expected)
    monkeypatch.setattr(
        backend, "inspect_devices",
        lambda: ([], GpuSelection((), None, order_known=False), True, ""),
    )
    requirements = tmp_path / "requirements.txt"
    requirements.write_text("test-backend-package==1.2.3\n")
    (tmp_path / ".clubb-jax-requirements.sha256").write_text(
        hashlib.sha256(requirements.read_bytes()).hexdigest()
    )
    packages = dict(expected)
    monkeypatch.setattr(
        runtime_info, "_installed_runtime", lambda _: {"packages": packages}
    )
    info = runtime_info.inspect_runtime(
        backend.PROFILE, accelerator, requirements, tmp_path, "0.11.0"
    )
    assert info["status"] == "ready"
    packages["test-backend-package"] = "wrong-version"
    info = runtime_info.inspect_runtime(
        backend.PROFILE, accelerator, requirements, tmp_path, "0.11.0"
    )
    assert info["status"] == "setup_required"
    commands = []

    def run(command, **kwargs):
        commands.append(command)
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr(run_jax.subprocess, "run", run)
    assert run_jax._packages_are_ready(Path("python"), "0.11.0", accelerator)
    assert json.loads(commands[0][-1]) == expected


@pytest.mark.parametrize("minor", [12, 13, 14])
def test_rocm_install_arguments_select_only_matching_python_wheels(
    minor, tmp_path, monkeypatch
):
    wheelhouse = tmp_path / "wheels"
    wheelhouse.mkdir()
    pjrt = wheelhouse / "jax_rocm7_pjrt-0.11.0-py3-none-any.whl"
    pjrt.touch()
    plugins = {}
    for supported in (12, 13, 14):
        plugins[supported] = (
            wheelhouse / f"jax_rocm7_plugin-0.11.0-cp3{supported}-manylinux.whl"
        )
        plugins[supported].touch()
    monkeypatch.setattr(backends.rocm, "prepare_wheelhouse", lambda _: wheelhouse)
    assert set(backends.rocm.install_arguments((3, minor), tmp_path)) == {
        str(pjrt), str(plugins[minor]),
    }


@pytest.mark.parametrize("accelerator", backends.BACKENDS)
def test_environment_setup_calls_install_hook_only_when_needed(
    accelerator, tmp_path, monkeypatch
):
    backend = backends.get_backend(accelerator)
    requirements = tmp_path / "requirements.txt"
    requirements.write_text("jax==0.11.0\n")
    venv = tmp_path / "venv"
    python = venv / "bin" / "python"
    python.parent.mkdir(parents=True)
    python.touch()
    tools_dir = tmp_path / "tools"
    hooks = []
    commands = []

    def install_arguments(version, directory):
        hooks.append((version, directory))
        return ["backend-wheel.whl"]

    def run(command, **kwargs):
        commands.append(command)
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr(backend, "install_arguments", install_arguments, raising=False)
    monkeypatch.setattr(run_jax, "_python_is_compatible", lambda *args: True)
    monkeypatch.setattr(run_jax, "_python_version", lambda _: (3, 12))
    monkeypatch.setattr(run_jax, "_packages_are_ready", lambda *args: True)
    monkeypatch.setattr(run_jax, "_ensure_uv", lambda *args: Path("uv"))
    monkeypatch.setattr(run_jax.subprocess, "run", run)
    for _ in range(2):
        assert run_jax._prepare_environment(
            accelerator, requirements, venv, tools_dir, {}
        ) == python
    assert hooks == [((3, 12), tools_dir)]
    assert commands == [[
        "uv", "pip", "install", "--python", str(python),
        "-r", str(requirements), "backend-wheel.whl",
    ]]


def test_importing_backends_has_no_runtime_or_setup_side_effects():
    script = """
import subprocess
import sys
import urllib.request
def unexpected(*args, **kwargs):
    raise AssertionError('import must not inspect hardware or download packages')
subprocess.run = unexpected
subprocess.check_output = unexpected
urllib.request.urlopen = unexpected
from clubb_jax import backends, run_jax, runtime_info
assert set(backends.BACKENDS) == {'cpu', 'cuda13', 'rocm', 'metal'}
assert not any(name == 'jax' or name.startswith('jax.') for name in sys.modules)
"""
    subprocess.run([sys.executable, "-c", script], check=True, cwd=run_jax.REPO_ROOT)


@pytest.mark.parametrize("accelerator", ["cuda13", "rocm", "metal"])
def test_host_callback_support_does_not_allow_cpu_model_fallback(accelerator, monkeypatch):
    backend = backends.get_backend(accelerator)
    assert backend.PLATFORM.split(",")[-1] == "cpu"
    assert backend.PLATFORM.split(",")[0] != "cpu"
    device = SimpleNamespace(platform="cpu", id=0, device_kind="CPU")
    fake_jax = SimpleNamespace(
        default_backend=lambda: "cpu", devices=lambda *args: [device],
        __version__="test",
    )
    monkeypatch.setitem(sys.modules, "jax", fake_jax)
    monkeypatch.setitem(sys.modules, "jaxlib", SimpleNamespace(__version__="test"))
    monkeypatch.setattr(sys, "argv", ["verify", backend.EXPECTED_BACKEND])
    original_import = importlib.import_module

    def import_component(name, *args):
        if name.startswith("jax_rocm7_plugin."):
            return SimpleNamespace()
        return original_import(name, *args)

    monkeypatch.setattr(importlib, "import_module", import_component)
    monkeypatch.setattr(
        run_jax.subprocess, "run", lambda command, **kwargs: exec(command[2], {})
    )
    with pytest.raises(AssertionError, match="JAX initialized cpu"):
        run_jax._verify_backend("python", accelerator, {})
