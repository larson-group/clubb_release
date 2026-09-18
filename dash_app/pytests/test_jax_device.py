"""GPU choice stays physical, explicit and isolated across dashboard jobs."""

import json
import subprocess
from types import SimpleNamespace

import pytest
from dash import Dash

from dash_app.compile_tab import callbacks, build_selector
from dash_app.shared.jax_device import jax_device_env, normalize_jax_gpu
from dash_app.services.models import ScmRunBatchRequest
from dash_app.shared import actions


GPU_A = "GPU-aaaaaaaa-1111-2222-3333-000000000001"
GPU_B = "GPU-bbbbbbbb-1111-2222-3333-000000000002"


def gpu_settings(uuid=GPU_B):
    return {"implementation": "jax", "jax_profile": "gpu", "jax_gpu": uuid}


def inventory():
    return {"gpu": {"status": "ready", "selectable": True, "hardware": {"gpus": [
        {"index": 0, "uuid": GPU_A, "name": "Test GPU A", "memory_mib": 8192},
        {"index": 1, "uuid": GPU_B, "name": "Test GPU B", "memory_mib": 16384},
    ]}}}


def test_environments_are_per_job_and_default_preserves_server_selection(monkeypatch):
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "1")
    base = {"CUDA_VISIBLE_DEVICES": "1", "KEEP": "yes"}
    assert jax_device_env(gpu_settings(), base) == {**base, "CUDA_VISIBLE_DEVICES": GPU_B}
    assert jax_device_env(gpu_settings(GPU_A), base)["CUDA_VISIBLE_DEVICES"] == GPU_A
    assert jax_device_env(gpu_settings(""), base) == base
    assert jax_device_env({"implementation": "jax", "jax_profile": "cpu"}, base) == base
    assert base["CUDA_VISIBLE_DEVICES"] == "1"
    assert jax_device_env({})["CUDA_VISIBLE_DEVICES"] == "1"


@pytest.mark.parametrize("value", ["1", "GPU-short", GPU_A + "," + GPU_B, "MIG-example", "bad\nvalue"])
def test_only_full_single_gpu_uuids_are_accepted(value):
    with pytest.raises(ValueError):
        normalize_jax_gpu(value)


def test_explicit_gpu_requires_gpu_backend():
    with pytest.raises(ValueError, match="requires"):
        jax_device_env({**gpu_settings(), "jax_profile": "cpu"})


def test_cpu_launch_ignores_remembered_gpu(monkeypatch):
    monkeypatch.setattr(build_selector, "build_implementation_capability", lambda *args, **kwargs: (True, ""))
    cpu = build_selector.selected_launch_target("jax", jax_profile="cpu", jax_gpu=GPU_B)
    gpu = build_selector.selected_launch_target("jax", jax_profile="gpu", jax_gpu=GPU_B)
    assert "jax_gpu" not in cpu
    assert gpu["jax_gpu"] == GPU_B


def compute_callback():
    app = Dash(__name__, suppress_callback_exceptions=True)
    callbacks.register_compile_callbacks(app)
    entry = next(entry for key, entry in app.callback_map.items()
                 if "compile-run-jax-profile.data" in key)
    # Metadata changes never silently switch a selected backend/device.
    assert len(entry["inputs"]) == 1
    return entry["callback"].__wrapped__


def test_gpu_tile_reselect_enables_gpu_after_cpu_without_changing_uuid(monkeypatch):
    select = compute_callback()
    monkeypatch.setattr(callbacks, "clicked_trigger_id", lambda: {"index": GPU_B})
    assert select([1], inventory()) == ("gpu", GPU_B)
    monkeypatch.setattr(callbacks, "clicked_trigger_id", lambda: {"index": "cpu"})
    assert select([1], inventory()) == ("cpu", callbacks.no_update)
    monkeypatch.setattr(callbacks, "clicked_trigger_id", lambda: {"index": GPU_B})
    assert select([2], inventory()) == ("gpu", GPU_B)


def test_device_choices_display_physical_indices_and_store_uuids():
    menu = callbacks.render_compact_build_selector({}, implementation="jax", jax_profile="gpu",
        jax_runtime_info=inventory(), jax_gpu=GPU_B)
    assert len(menu) == 2  # Implementation and one flat compute tile grid.
    cpu, first, second = menu[1].children[1].children
    assert cpu.children.children[0].children[0].children == "CPU"
    assert first.children.children[0].children[0].children == "GPU 0"
    assert first.children.children[1].children == "Test GPU A · 8 GiB"
    assert second.children.children[0].children[0].children == "GPU 1"
    assert second.children.children[1].children == "Test GPU B · 16 GiB"
    assert all("compile-profile-card" in tile.className for tile in [cpu, first, second])
    assert second.id["index"] == GPU_B
    assert second.to_plotly_json()["props"]["aria-pressed"] == "true"
    missing = callbacks.render_compact_build_selector({}, implementation="jax", jax_profile="gpu",
        jax_runtime_info={}, jax_gpu=GPU_B)
    missing_tile = missing[1].children[1].children[-1]
    assert missing_tile.children.children[0].children[0].children == "Saved GPU"
    assert missing_tile.disabled is True


def test_click_callback_ignores_mount_and_unknown_device(monkeypatch):
    select = compute_callback()
    monkeypatch.setattr(callbacks, "clicked_trigger_id", lambda: None)
    assert select([0], inventory()) == (callbacks.no_update, callbacks.no_update)
    for value, expected in [(GPU_B, ("gpu", GPU_B)), ("gpu", ("gpu", "")),
                            ("unknown", (callbacks.no_update, callbacks.no_update))]:
        monkeypatch.setattr(callbacks, "clicked_trigger_id", lambda: {"index": value})
        assert select([1], inventory()) == expected


def test_cpu_selected_has_no_remembered_gpu_highlight_and_no_device_submenu():
    menu = callbacks.render_compact_build_selector({}, implementation="jax", jax_profile="cpu",
        jax_runtime_info=inventory(), jax_gpu=GPU_B)
    tiles = menu[1].children[1].children
    assert len(menu) == 2
    assert len(tiles) == 3
    assert [tile.to_plotly_json()["props"]["aria-pressed"] for tile in tiles] == ["true", "false", "false"]


def test_incompatible_gpu_does_not_disable_other_gpu_tile():
    info = inventory()
    info["gpu"].update(status="unavailable", selectable=False, reason="Test incompatibility")
    info["gpu"]["hardware"]["selected_gpu"] = info["gpu"]["hardware"]["gpus"][0]
    menu = callbacks.render_compact_build_selector({}, implementation="jax", jax_profile="gpu",
        jax_runtime_info=info, jax_gpu=GPU_A)
    cpu, first, second = menu[1].children[1].children
    assert first.disabled is True
    assert first.title == "Test incompatibility"
    assert second.disabled is False
    assert second.children.children[0].children[1].children == "Check on selection"


def test_preflight_uses_same_uuid_as_launch_without_changing_cpu(tmp_path, monkeypatch):
    wrapper = tmp_path / "clubb_jax" / "run_jax.py"
    wrapper.parent.mkdir()
    wrapper.touch(mode=0o755)
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "0")
    seen = {}

    def inspect(command, **kwargs):
        profile = command[1].split("=")[1]
        seen[profile] = kwargs["env"]["CUDA_VISIBLE_DEVICES"]
        return subprocess.CompletedProcess(command, 0, json.dumps({
            "schema_version": 1, "profile": profile, "status": "ready", "selectable": True,
        }), "")

    monkeypatch.setattr(build_selector.subprocess, "run", inspect)
    build_selector.inspect_jax_runtime_profiles(tmp_path, jax_gpu=GPU_B)
    assert seen == {"cpu": "0", "gpu": GPU_B}


def test_batch_children_and_provenance_retain_device():
    request = ScmRunBatchRequest(request_id="test-gpu-batch", cases=["arm", "bomex"], **gpu_settings())
    child = actions._batch_child_request(request, {"job_id": "batch"}, "arm")
    assert child.jax_gpu == GPU_B
    assert actions._normalize_dashboard_cli_options(gpu_settings())["jax_gpu"] == GPU_B
    assert actions._scm_build_identity(gpu_settings())["environment"]["CUDA_VISIBLE_DEVICES"] == GPU_B


@pytest.mark.parametrize("prealloc", [False, True])
def test_run_process_receives_uuid(monkeypatch, tmp_path, prealloc):
    from dash_app.run_tab import runtime

    captured = {}
    def launch(command, **kwargs):
        captured.update(kwargs)
        assert ("-jax=gpu,xla_prealloc" if prealloc else "-jax=gpu") in command
        return SimpleNamespace(pid=123)

    monkeypatch.setattr(runtime, "write_temp_namelist", lambda *args: None)
    monkeypatch.setattr(runtime, "run_child_env", lambda: {"CUDA_VISIBLE_DEVICES": "0", "XLA_PYTHON_CLIENT_PREALLOCATE": "true"})
    monkeypatch.setattr(runtime, "mark_case_started", lambda proc: None)
    monkeypatch.setattr(runtime.subprocess, "Popen", launch)
    monkeypatch.setattr(runtime.tempfile, "NamedTemporaryFile", lambda **kwargs: (tmp_path / "run.log").open("w"))
    runtime.start_case_process("arm", "none", {}, {**gpu_settings(), "jax_xla_prealloc": prealloc})
    assert captured["env"]["CUDA_VISIBLE_DEVICES"] == GPU_B
    assert captured["env"]["XLA_PYTHON_CLIENT_PREALLOCATE"] == str(prealloc).lower()


@pytest.mark.parametrize("prealloc", [False, True])
def test_profile_process_and_preview_receive_uuid(monkeypatch, tmp_path, prealloc):
    from dash_app.profile_tab import runtime
    from dash_app.pytests.test_profile_tab import settings

    selected = {**settings(tmp_path), **gpu_settings(), "jax_xla_prealloc": prealloc}
    captured = {}
    def launch(command, **kwargs):
        captured.update(kwargs)
        assert ("-jax=gpu,xla_prealloc" if prealloc else "-jax=gpu") in command
        return SimpleNamespace(pid=123)

    monkeypatch.setattr(runtime.subprocess, "Popen", launch)
    monkeypatch.setattr(runtime, "_process_started_at", lambda pid: 0)
    monkeypatch.setattr(runtime, "PROFILE_PROCESSES", {})
    monkeypatch.setattr(runtime.tempfile, "NamedTemporaryFile", lambda **kwargs: (tmp_path / "profile.log").open("wb"))
    job = runtime.start_profile_process(selected)
    assert captured["env"]["CUDA_VISIBLE_DEVICES"] == GPU_B
    assert captured["env"]["XLA_PYTHON_CLIENT_PREALLOCATE"] == str(prealloc).lower()
    assert job["settings"]["jax_gpu"] == GPU_B
    assert f"CUDA_VISIBLE_DEVICES={GPU_B} " in job["command_display"]
    assert f"XLA_PYTHON_CLIENT_PREALLOCATE={str(prealloc).lower()} " in job["command_display"]
    assert job["command_display"] == runtime.profile_command_display(selected)


@pytest.mark.parametrize("profile,status,selectable,selected_gpu,enabled", [
    ("gpu", "ready", True, GPU_B, True),
    ("gpu", "setup_required", True, GPU_B, True),
    ("cpu", "ready", True, GPU_B, False),
    ("gpu", "unavailable", False, GPU_B, False),
    ("gpu", "checking", True, GPU_B, False),
    ("gpu", "ready", True, GPU_A, False),
    ("gpu", "ready", True, None, False),
])
def test_preallocation_checkbox_requires_compatible_selected_gpu(profile, status, selectable, selected_gpu, enabled):
    info = inventory()
    info["gpu"].update(status=status, selectable=selectable)
    info["gpu"]["hardware"]["selected_gpu"] = {"uuid": selected_gpu}
    menu = callbacks.render_compact_build_selector({}, implementation="jax", jax_profile=profile,
        jax_runtime_info=info, jax_gpu=GPU_B, jax_xla_prealloc=True)
    checkbox = menu[1].children[2].children
    assert checkbox.value == ["enabled"]  # Remember preference even while disabled.
    assert checkbox.options[0]["disabled"] is not enabled


def test_preallocation_callback_rejects_disabled_changes_and_preserves_false():
    app = Dash(__name__, suppress_callback_exceptions=True)
    callbacks.register_compile_callbacks(app)
    change = app.callback_map["compile-run-jax-xla-prealloc.data"]["callback"].__wrapped__
    info = inventory()
    info["gpu"]["hardware"]["selected_gpu"] = {"uuid": GPU_B}
    assert change([["enabled"]], info, "cpu", GPU_B, False) is callbacks.no_update
    assert change([["enabled"]], info, "gpu", GPU_B, False) is True
    assert change([[]], info, "gpu", GPU_B, True) is False
    assert change([[]], info, "gpu", GPU_B, False) is callbacks.no_update
    assert change([], info, "gpu", GPU_B, True) is callbacks.no_update


def test_preallocation_is_validated_and_frozen_in_batch_children():
    selected = {**gpu_settings(), "jax_xla_prealloc": False}
    request = ScmRunBatchRequest(request_id="test-prealloc", cases=["arm"], **selected)
    assert actions._batch_child_request(request, {"job_id": "batch"}, "arm").jax_xla_prealloc is False
    assert actions._normalize_dashboard_cli_options(selected)["jax_xla_prealloc"] is False
    assert actions._scm_build_identity(selected)["environment"]["XLA_PYTHON_CLIENT_PREALLOCATE"] == "false"
    with pytest.raises(ValueError):
        actions._normalize_dashboard_cli_options({**selected, "jax_xla_prealloc": "false"})
    cpu = build_selector.selected_launch_target("jax", jax_profile="cpu", jax_xla_prealloc=True)
    assert "jax_xla_prealloc" not in cpu
