"""Saved tuner selections reach the JAX launcher without changing the parent environment."""

import os
import json
import runpy
import sys
from types import SimpleNamespace

import pytest

from clubb_jax.run_jax import runtime_arguments
from tuner.job_runtime import TunerJob, tuner_runtime_settings, tuner_worker_env


@pytest.mark.parametrize("prealloc", [False, True])
def test_tuner_child_honors_saved_gpu_and_preallocation(monkeypatch, prealloc):
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "parent-device")
    monkeypatch.setenv("XLA_PYTHON_CLIENT_PREALLOCATE", "parent-setting")
    request = {
        "backend": "jax", "jax_options": "gpu",
        "jax_gpu": "GPU-aaaaaaaa-bbbb-cccc-dddd-eeeeeeeeeeee",
        "jax_xla_prealloc": prealloc,
    }
    env = tuner_worker_env()
    settings = tuner_runtime_settings(request)
    args = runtime_arguments(settings["jax_profile"], device=settings["jax_gpu"], prealloc_gpu_mem=settings["jax_xla_prealloc"])
    assert env["CUDA_VISIBLE_DEVICES"] == "parent-device"
    assert args == [f"-options=gpu,device={request['jax_gpu']},prealloc_gpu_mem={str(prealloc).lower()}"]
    assert env["XLA_PYTHON_CLIENT_PREALLOCATE"] == "parent-setting"
    assert os.environ["CUDA_VISIBLE_DEVICES"] == "parent-device"
    assert os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] == "parent-setting"


def test_cpu_tuner_rejects_explicit_gpu_selection():
    with pytest.raises(ValueError, match="GPU profile"):
        tuner_runtime_settings({
            "backend": "jax", "jax_options": "cpu",
            "jax_gpu": "GPU-aaaaaaaa-bbbb-cccc-dddd-eeeeeeeeeeee",
        })


def test_tuner_rejects_conflicting_preallocation_settings():
    with pytest.raises(ValueError, match="conflicts"):
        tuner_runtime_settings({
            "backend": "jax", "jax_options": "gpu,xla_prealloc",
            "jax_xla_prealloc": False,
        })


@pytest.mark.parametrize("entry", ["job", "direct"])
@pytest.mark.parametrize("prealloc_gpu_mem", [False, True])
def test_saved_device_selection_is_forwarded_at_each_tuner_entry(tmp_path, monkeypatch, entry, prealloc_gpu_mem):
    gpu = "GPU-aaaaaaaa-bbbb-cccc-dddd-eeeeeeeeeeee"
    request = {
        "backend": "jax", "jax_options": "gpu",
        "jax_gpu": gpu, "jax_xla_prealloc": prealloc_gpu_mem,
    }
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "parent-device")
    monkeypatch.setenv("XLA_PYTHON_CLIENT_PREALLOCATE", "parent-setting")
    captured = {}
    if entry == "job":
        job = TunerJob.create(request, job_dir=tmp_path / "job")

        def launch(command, **kwargs):
            captured.update(command=command, env=kwargs["env"])
            return SimpleNamespace(pid=123)

        monkeypatch.setattr("tuner.job_runtime.subprocess.Popen", launch)
        job.start()
    else:
        (tmp_path / "request.json").write_text(json.dumps(request))
        monkeypatch.setattr(sys, "argv", ["tuner.tune_clubb", "-job_dir", str(tmp_path)])
        monkeypatch.delenv("_CLUBB_JAX_ENVIRONMENT_PYTHON", raising=False)

        class LaunchObserved(BaseException):
            pass

        def launch(python, command, env):
            captured.update(command=command, env=env)
            raise LaunchObserved

        monkeypatch.setattr(os, "execve", launch)
        with pytest.raises(LaunchObserved):
            runpy.run_module("tuner.tune_clubb", run_name="__main__")
    assert f"-options=gpu,device={gpu},prealloc_gpu_mem={str(prealloc_gpu_mem).lower()}" in captured["command"]
    assert "-module=tuner.tune_clubb" in captured["command"]
    assert captured["env"]["CUDA_VISIBLE_DEVICES"] == "parent-device"
    assert captured["env"]["XLA_PYTHON_CLIENT_PREALLOCATE"] == "parent-setting"


def test_tuner_cli_preserves_the_opaque_selection_value():
    from run_scripts.run_tuner_job import build_request, parse_args

    gpu = "GPU-aaaaaaaa-bbbb-cccc-dddd-eeeeeeeeeeee"
    value = f"gpu,device={gpu},prealloc_gpu_mem=false"
    args = parse_args([
        "-jax=" + value, "-cases", "bomex", "-fields", "cloud_frac",
        "-param_ranges", "C8:0.2:0.8", "-run_top", "never",
    ])
    request = build_request(args)
    assert request["jax_options"] == value
    settings = tuner_runtime_settings(request)
    assert settings["jax_profile"] == "gpu"
    assert settings["jax_gpu"] == gpu
    assert settings["jax_xla_prealloc"] is False
