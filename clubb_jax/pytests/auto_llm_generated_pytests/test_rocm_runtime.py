import hashlib
import importlib
import subprocess
import sys
import zipfile
from types import SimpleNamespace

import pytest

from clubb_jax import run_jax, runtime_info
from clubb_jax.backends import rocm, metal, native_gpu_accelerator


def test_gpu_profile_detects_amd_without_nvidia(monkeypatch):
    monkeypatch.setattr(metal.platform, "system", lambda: "Linux")
    monkeypatch.setattr(run_jax.shutil, "which", lambda _: None)
    monkeypatch.setattr(rocm, "has_gpu", lambda: True)
    values, _ = run_jax.parse_launcher_args(["-options=gpu"])
    assert run_jax.resolve_accelerator(values) == ("rocm", "gpu")


def test_nvidia_remains_the_default_on_a_mixed_gpu_host(monkeypatch):
    monkeypatch.setattr(metal.platform, "system", lambda: "Linux")
    monkeypatch.setattr(run_jax.shutil, "which", lambda _: "/usr/bin/nvidia-smi")
    monkeypatch.setattr(rocm, "has_gpu", lambda: True)
    assert native_gpu_accelerator() == "cuda13"


def test_explicit_rocm_uses_an_isolated_environment():
    values, _ = run_jax.parse_launcher_args(["-accelerator=rocm"])
    assert run_jax.resolve_accelerator(values) == ("rocm", "gpu")
    requirements, venv, platform = run_jax.runtime_paths("rocm")
    assert requirements.name == "requirements-rocm.txt"
    assert venv.name == ".venv-jax-rocm"
    assert platform == "rocm,cpu"


@pytest.mark.parametrize("version,compatible", [((3, 11), False), ((3, 12), True),
                                                   ((3, 14), True), ((3, 15), False)])
def test_rocm_python_versions_match_the_pinned_wheels(version, compatible, monkeypatch):
    monkeypatch.setattr(run_jax, "_python_version", lambda _: version)
    assert run_jax._python_is_compatible("python", "rocm") is compatible


def test_rocm_query_reads_gpu_agents_and_ignores_the_cpu(tmp_path, monkeypatch):
    metadata = tmp_path / ".info"
    metadata.mkdir()
    (metadata / "version").write_text("7.2.4\n")
    monkeypatch.setenv("ROCM_PATH", str(tmp_path))
    monkeypatch.setattr(rocm.platform, "system", lambda: "Linux")
    monkeypatch.setattr(runtime_info.os, "access", lambda *_: True)
    monkeypatch.setattr(runtime_info.subprocess, "run", lambda *args, **kwargs:
                        subprocess.CompletedProcess(args, 0, """
Agent 1
  Name: AMD CPU
  Device Type: CPU
Agent 2
  Name: gfx1151
  Marketing Name: AMD Radeon 8060S Graphics
  Device Type: GPU
  Name: amdgcn-amd-amdhsa--gfx1151
""", ""))
    gpus, error = rocm.query_gpus()
    assert error == ""
    assert len(gpus) == 1
    assert gpus[0]["name"] == "AMD Radeon 8060S Graphics"
    assert gpus[0]["gfx_target"] == "gfx1151"
    assert gpus[0]["driver_version"] == "7.2.4"


def test_rocm_query_rejects_incompatible_runtime_before_setup(tmp_path, monkeypatch):
    metadata = tmp_path / ".info"
    metadata.mkdir()
    (metadata / "version").write_text("7.14.1\n")
    monkeypatch.setenv("ROCM_PATH", str(tmp_path))
    monkeypatch.setattr(rocm.platform, "system", lambda: "Linux")
    gpus, error = rocm.query_gpus()
    assert not gpus
    assert "require ROCm 7.2.x" in error


def test_rocm_inspection_requires_both_plugin_packages(tmp_path, monkeypatch):
    gpu = {"index": "0", "name": "Radeon", "gfx_target": "gfx1151", "driver_version": "7.2.4"}
    monkeypatch.setattr(rocm, "query_gpus", lambda: ([gpu], ""))
    requirements = tmp_path / "requirements.txt"
    requirements.write_text("jax==0.11.0\n")
    venv = tmp_path / "venv"
    venv.mkdir()
    (venv / ".clubb-jax-requirements.sha256").write_text(
        hashlib.sha256(requirements.read_bytes()).hexdigest()
    )
    packages = {"jax": "0.11.0", "jaxlib": "0.11.0", "jax-rocm7-plugin": "0.11.0"}
    monkeypatch.setattr(runtime_info, "_installed_runtime", lambda _: {"packages": packages})
    info = runtime_info.inspect_runtime("gpu", "rocm", requirements, venv, "0.11.0")
    assert info["status"] == "setup_required"
    assert info["selectable"] is True
    assert "AMD GPU: Radeon (gfx1151)" in runtime_info.format_human(info)
    packages["jax-rocm7-pjrt"] = "0.11.0"
    assert runtime_info.inspect_runtime("gpu", "rocm", requirements, venv, "0.11.0")["status"] == "ready"


def test_rocm_inspection_reports_hidden_host_devices(tmp_path, monkeypatch):
    monkeypatch.setattr(rocm, "query_gpus", lambda: ([], "No /dev/kfd access"))
    requirements = tmp_path / "requirements.txt"
    requirements.touch()
    info = runtime_info.inspect_runtime("gpu", "rocm", requirements, tmp_path / "venv", "0.11.0")
    assert info["status"] == "unavailable"
    assert info["reason"] == "No /dev/kfd access"


def test_rocm_archive_is_verified_and_only_extracts_wheels(tmp_path, monkeypatch):
    archive_path = tmp_path / "wheelhouse_legacy_rocm7.2.0.zip"
    with zipfile.ZipFile(archive_path, "w") as archive:
        archive.writestr("../unrelated.txt", "not extracted")
        archive.writestr("directory/jax_rocm7_pjrt-0.11.0-py3-none.whl", "wheel")
        for python_tag in ("cp312", "cp313", "cp314"):
            archive.writestr(f"directory/jax_rocm7_plugin-0.11.0-{python_tag}.whl", "wheel")
    monkeypatch.setattr(rocm, "ROCM_WHEELHOUSE_SHA256", hashlib.sha256(archive_path.read_bytes()).hexdigest())
    wheelhouse = rocm.prepare_wheelhouse(tmp_path)
    assert len(list(wheelhouse.glob("*.whl"))) == 4
    assert not (tmp_path.parent / "unrelated.txt").exists()
    assert rocm.prepare_wheelhouse(tmp_path) == wheelhouse


def test_rocm_archive_with_bad_checksum_is_not_extracted(tmp_path):
    (tmp_path / "wheelhouse_legacy_rocm7.2.0.zip").write_bytes(b"wrong archive")
    with pytest.raises(run_jax.LauncherError, match="checksum mismatch"):
        rocm.prepare_wheelhouse(tmp_path)
    assert not (tmp_path / "wheelhouse-legacy-rocm7.2.0").exists()


def test_rocm_environment_reuses_a_local_runtime_without_modifying_the_host(tmp_path):
    local_root = tmp_path / "rocm72-root" / "opt" / "rocm-7.2.0"
    local_root.mkdir(parents=True)
    env = {"LD_LIBRARY_PATH": "/existing/lib"}
    rocm.configure_environment(env, tmp_path)
    assert env["ROCM_PATH"] == str(local_root)
    assert env["LD_LIBRARY_PATH"] == f"{local_root}/lib:{local_root}/lib64:/existing/lib"
    assert env["XLA_PYTHON_CLIENT_PREALLOCATE"] == "false"


def test_rocm_environment_preserves_an_explicit_runtime_and_memory_setting(tmp_path):
    env = {"ROCM_PATH": "/custom/rocm", "XLA_PYTHON_CLIENT_PREALLOCATE": "true"}
    rocm.configure_environment(env, tmp_path)
    assert env["ROCM_PATH"] == "/custom/rocm"
    assert env["XLA_PYTHON_CLIENT_PREALLOCATE"] == "true"
    assert env["LD_LIBRARY_PATH"] == "/custom/rocm/lib:/custom/rocm/lib64"


def test_rocm_environment_prefers_the_newer_local_maintenance_runtime(tmp_path):
    older_root = tmp_path / "rocm72-root" / "opt" / "rocm-7.2.0"
    newer_root = tmp_path / "rocm724-root" / "opt" / "rocm-7.2.4"
    older_root.mkdir(parents=True)
    newer_root.mkdir(parents=True)
    env = {}
    rocm.configure_environment(env, tmp_path)
    assert env["ROCM_PATH"] == str(newer_root)
    assert env["LD_LIBRARY_PATH"] == f"{newer_root}/lib:{newer_root}/lib64"


def test_rocm_environment_can_replace_an_incomplete_system_runtime(tmp_path, monkeypatch):
    local_root = tmp_path / "rocm72-root" / "opt" / "rocm-7.2.0"
    local_root.mkdir(parents=True)
    monkeypatch.setattr(rocm.Path, "exists", lambda _: False)
    env = {"ROCM_PATH": "/opt/rocm"}
    rocm.configure_environment(env, tmp_path)
    assert env["ROCM_PATH"] == str(local_root)


def test_gfx1151_disables_unsupported_autotuning_kernels_but_preserves_other_flags(tmp_path, monkeypatch):
    monkeypatch.setattr(rocm, "has_gpu", lambda target=None: target == 110501)
    env = {"XLA_FLAGS": "--xla_gpu_enable_fast_min_max=false"}
    rocm.configure_environment(env, tmp_path)
    assert env["XLA_FLAGS"] == (
        "--xla_gpu_enable_fast_min_max=false --xla_gpu_autotune_level=0 "
        "--xla_gpu_enable_command_buffer="
    )
    assert env["CLUBB_JAX_PORTABLE_CHOLESKY"] == "1"


def test_rocm_preserves_an_explicit_autotuning_override(tmp_path, monkeypatch):
    monkeypatch.setattr(rocm, "has_gpu", lambda target=None: True)
    env = {"XLA_FLAGS": "--xla_gpu_autotune_level=4"}
    rocm.configure_environment(env, tmp_path)
    assert env["XLA_FLAGS"] == "--xla_gpu_autotune_level=4 --xla_gpu_enable_command_buffer="


def test_rocm_preserves_an_explicit_command_buffer_override(tmp_path, monkeypatch):
    monkeypatch.setattr(rocm, "has_gpu", lambda target=None: True)
    env = {"XLA_FLAGS": "--xla_gpu_enable_command_buffer=FUSION"}
    rocm.configure_environment(env, tmp_path)
    assert env["XLA_FLAGS"] == (
        "--xla_gpu_enable_command_buffer=FUSION --xla_gpu_autotune_level=0"
    )


def test_other_amd_targets_do_not_receive_the_gfx1151_workaround(tmp_path, monkeypatch):
    monkeypatch.setattr(rocm, "has_gpu", lambda target=None: False)
    env = {}
    rocm.configure_environment(env, tmp_path)
    assert "XLA_FLAGS" not in env
    assert "CLUBB_JAX_PORTABLE_CHOLESKY" not in env


def test_rocm_preserves_an_explicit_cholesky_override(tmp_path, monkeypatch):
    monkeypatch.setattr(rocm, "has_gpu", lambda target=None: True)
    env = {"CLUBB_JAX_PORTABLE_CHOLESKY": "0"}
    rocm.configure_environment(env, tmp_path)
    assert env["CLUBB_JAX_PORTABLE_CHOLESKY"] == "0"


@pytest.mark.parametrize("missing_component", [None, "_linalg", "_solver", "_sparse"])
def test_rocm_backend_verification_requires_loadable_math_libraries(missing_component, monkeypatch):
    device = SimpleNamespace(platform="gpu", id=0, device_kind="Radeon")
    fake_jax = SimpleNamespace(default_backend=lambda: "gpu", devices=lambda *args: [device],
                               __version__="0.11.0")
    monkeypatch.setitem(sys.modules, "jax", fake_jax)
    monkeypatch.setitem(sys.modules, "jaxlib", SimpleNamespace(__version__="0.11.0"))
    monkeypatch.setattr(sys, "argv", ["verify", "rocm"])
    original_import = importlib.import_module
    imports = []

    def import_component(name, *args):
        if name.startswith("jax_rocm7_plugin."):
            imports.append(name)
            if name == "jax_rocm7_plugin." + str(missing_component):
                raise ImportError("missing shared library")
            return SimpleNamespace()
        return original_import(name, *args)

    monkeypatch.setattr(importlib, "import_module", import_component)
    monkeypatch.setattr(run_jax.subprocess, "run", lambda command, **kwargs: exec(command[2], {}))
    if missing_component:
        with pytest.raises(RuntimeError, match="complete ROCm HIP SDK"):
            run_jax._verify_backend("python", "rocm", {})
    else:
        run_jax._verify_backend("python", "rocm", {})
        assert imports == ["jax_rocm7_plugin." + component
                           for component in ("_linalg", "_solver", "_sparse")]
