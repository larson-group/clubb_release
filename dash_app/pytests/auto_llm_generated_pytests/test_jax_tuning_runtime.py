"""Dash runtime selection must survive Tune launch, reload and result replay.

Exercise the ordinary Start callback, shared broker request and existing
runtime adapters without model runs. Browser/model verification is performed
through the running app; these checks isolate the lost-backend regression.
"""

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from dash import no_update
from dash_app.shared import actions, activity, broker_client
from dash_app.services import TuneRequest
from dash_app.tune_tab import callbacks_runs, runtime


class CallbackRegistry:
    def __init__(self):
        self.functions = {}

    def callback(self, *args, **kwargs):
        def register(function):
            self.functions[function.__name__] = function
            return function
        return register


def input_values():
    return dict(
        case_names=["bomex"], time_start_values=[10800], time_end_values=[11400],
        average_time_values=[600], altitude_min_values=[20], altitude_max_values=[2940],
        selected_fields=["cloud_frac"], min_values=[0.2], max_values=[0.8],
        batch_size=2, max_workers=1, strategy_mode="random", loss_mode="shape_first",
        aggregation_mode="quantile_weighted", aggregation_scope="overall",
        aggregation_weight_1=0.1, aggregation_weight_2=0.4,
        aggregation_weight_3=0.4, aggregation_weight_4=0.1,
        random_max_samples=4, resolve_spacing=0.1,
        simann_max_iters=2, simann_initial_temp=1.0, simann_final_temp=1e-12,
        adam_max_updates=2, adam_learning_rate_percent=1.0,
        adam_perturbation_percent=5.0, adam_spsa_pairs=2,
        selected_config="default", scm_override="",
    )


@pytest.mark.parametrize("implementation", ["fortran", "jax"])
def test_start_callback_sends_the_selected_backend(monkeypatch, implementation):
    registry = CallbackRegistry()
    callbacks_runs.register_run_callbacks(registry)
    monkeypatch.setattr(callbacks_runs, "build_validation_message", lambda *args, **kwargs: "")
    captured = {}

    def launch(action, payload, **kwargs):
        assert action == "launch_tuning_request"
        captured.update(payload["request"])
        return {"job": {"job_dir": "/tmp/callback-job", "pid": 1234}}

    monkeypatch.setattr(broker_client, "perform_action", launch)
    values = input_values()
    values.update(
        _n_clicks=1,
        member_ids=[{"type": "tune-range-member", "row": 0, "member": 0}],
        member_values=["C8"], range_row_order=[0], case_data={}, tunable_configs=["default"],
        tunable_names=["C8"], active_job={}, workspace_selection={"mode": "new"},
        displayed_status={}, run_implementation=implementation,
        jax_profile="cpu", jax_gpu="", jax_xla_prealloc=False,
    )
    response = registry.functions["start_tuning"](**values)
    assert response[0] is not no_update
    assert captured["backend"] == implementation
    if implementation == "jax":
        assert captured["jax_options"] == "cpu"
        assert captured["jax_xla_prealloc"] is False


@pytest.mark.parametrize("backend,options", [("fortran", "cpu"), ("jax", "cpu"), ("jax", "gpu,xla_prealloc")])
def test_broker_launch_retains_saved_runtime(tmp_path, monkeypatch, backend, options):
    monkeypatch.setattr(activity, "ACTIVITY_PATH", tmp_path / "activity.json")
    monkeypatch.setattr(activity, "LOCK_PATH", tmp_path / "activity.lock")
    activity.reset_activity()
    monkeypatch.setattr(actions, "recover_active_tuning_from_disk", lambda: None)
    monkeypatch.setattr(actions, "broker_jobs", lambda: {})
    monkeypatch.setattr(actions, "_background", lambda *args: None)
    captured = {}

    def start(request, **kwargs):
        captured.update(request)
        return {"job_dir": str(tmp_path), "status_path": str(tmp_path / "status.json")}

    monkeypatch.setattr(actions, "start_tuning_job", start)
    payload = {
        "backend": backend, "jax_options": options,
        "case_configs": [{"case_name": "bomex", "time_average_range": [10800, 11400], "num_time_windows": 1}],
        "parameter_ranges": [{"name": "C8", "min": 0.2, "max": 0.8}],
        "selected_fields": ["cloud_frac"], "strategy": {"name": "random", "options": {"max_samples": 4}},
        "batch_size": 2, "max_workers": 1,
    }
    actions.launch_tuning_request(payload)
    assert captured["backend"] == backend
    assert captured["jax_options"] == options
    controls = callbacks_runs.agent_request_to_tune_controls(captured, {"bomex": {"clubb_fields": ["cloud_frac"]}})
    assert controls["runtime"]["implementation"] == backend
    assert controls["runtime"]["jax_profile"] == options.split(",")[0]
    assert controls["runtime"]["jax_xla_prealloc"] is (True if "xla_prealloc" in options else None)


def test_typed_tune_request_retains_jax_runtime():
    request = TuneRequest(
        request_id="jax-runtime-probe", backend="jax", jax_options="gpu",
        jax_gpu="GPU-aaaaaaaa-bbbb-cccc-dddd-eeeeeeeeeeee", jax_xla_prealloc=False,
        cases=["bomex"], parameter_ranges=[{"name": "C8", "min": 0.2, "max": 0.8}],
    )
    payload = request.model_dump()
    assert payload["backend"] == "jax"
    assert payload["jax_gpu"] == request.jax_gpu
    assert payload["jax_xla_prealloc"] is False


@pytest.mark.parametrize("mode", ["window", "complete"])
def test_saved_jax_replay_uses_launcher_and_winning_columns(tmp_path, monkeypatch, mode):
    captured = []
    monkeypatch.setattr(runtime, "TUNE_RESULT_OUTPUT_ROOT", tmp_path / "outputs")

    def launch(command, **kwargs):
        captured.append((command, kwargs))
        return SimpleNamespace(pid=123, poll=lambda: 0)

    monkeypatch.setattr(runtime.subprocess, "Popen", launch)
    result = runtime.start_loss_run(
        ["bomex"], ["cloud_frac"], [{"C8": 0.7}, {"C8": 0.5}],
        run_mode=mode, runtime_request={"backend": "jax", "jax_options": "cpu"},
        override="C8=0.45,C11=0.4", workspace_name="runtime-check",
    )
    command, options = captured[0]
    assert "-jax=cpu" in command
    assert "-python" not in command
    saved = json.loads(Path(result["request_path"]).read_text())
    assert saved["runtime"]["implementation"] == "jax"
    assert saved["runtime"]["jax_profile"] == "cpu"
    assert json.loads(saved["override"])["bomex"] == {"C11": "0.4"}
    assert "C8 = 0.7, 0.5" in Path(result["params_path"]).read_text()


@pytest.mark.parametrize("mode", ["window", "complete"])
@pytest.mark.parametrize("preallocation", [False, True])
def test_gpu_replay_forwards_saved_selection_without_preparing_device_environment(tmp_path, monkeypatch, mode, preallocation):
    gpu = "GPU-aaaaaaaa-bbbb-cccc-dddd-eeeeeeeeeeee"
    monkeypatch.setattr(runtime, "TUNE_RESULT_OUTPUT_ROOT", tmp_path / "outputs")
    monkeypatch.setattr(runtime, "worker_env", lambda: {"CUDA_VISIBLE_DEVICES": "parent-device"})
    captured = []

    def launch(command, **kwargs):
        captured.append((command, kwargs))
        return SimpleNamespace(pid=123, poll=lambda: 0)

    monkeypatch.setattr(runtime.subprocess, "Popen", launch)
    runtime.start_loss_run(
        ["bomex"], ["cloud_frac"], [{"C8": 0.7}], run_mode=mode,
        workspace_name="device-forwarding", runtime_request={
            "backend": "jax", "jax_options": "gpu", "jax_gpu": gpu,
            "jax_xla_prealloc": preallocation,
        },
    )
    command, kwargs = captured[0]
    assert f"-jax=gpu,device={gpu},prealloc_gpu_mem={str(preallocation).lower()}" in command
    assert kwargs["env"]["CUDA_VISIBLE_DEVICES"] == "parent-device"


def test_selected_revision_replay_ignores_other_broker_job(tmp_path, monkeypatch):
    monkeypatch.setattr(activity, "ACTIVITY_PATH", tmp_path / "activity.json")
    monkeypatch.setattr(activity, "LOCK_PATH", tmp_path / "activity.lock")
    activity.reset_activity()
    request = {"backend": "jax", "jax_options": "cpu", "cases": ["bomex"], "selected_fields": ["cloud_frac"], "case_configs": [{"case_name": "bomex"}]}
    monkeypatch.setattr(actions, "load_tune_workspace_execution", lambda *args: {"job": {"status_path": "/tmp/status"}, "request": request})
    monkeypatch.setattr(actions, "read_tuning_status", lambda path: {"top_results": [{"params": {"C8": 0.7}}]})
    monkeypatch.setattr(actions, "_background", lambda *args: None)
    captured = {}

    def replay(*args, **kwargs):
        captured.update(kwargs)
        return {"run_id": "selected-replay", "pid": 123, "state": "running"}

    monkeypatch.setattr(actions, "start_loss_run", replay)
    actions.run_tuning_loss("window", 1, workspace_id="chosen", revision_id="rev1")
    assert captured["runtime_request"] == request
    assert captured["workspace_id"] == "chosen"
    assert captured["revision_id"] == "rev1"


def test_result_buttons_use_the_selected_revision_in_the_broker(monkeypatch):
    registry = CallbackRegistry()
    callbacks_runs.register_run_callbacks(registry)
    monkeypatch.setattr(callbacks_runs, "ctx", SimpleNamespace(
        triggered_id={"type": "tune-loss-run-button", "action": "window"},
        triggered=[{"value": [123]}],
    ))
    captured = {}

    def launch(action, payload, **kwargs):
        assert action == "run_tuning_loss"
        captured.update(payload)
        return {"run": {"run_id": "broker-replay", "pid": 123, "log_path": "/tmp/replay.log"}}

    monkeypatch.setattr(broker_client, "perform_action", launch)
    result = registry.functions["run_result_loss"](
        [123], [{"params": {"C8": 0.7}}], [], {},
        {"workspace_id": "selected", "revision_id": "rev2"},
    )
    assert captured["workspace_id"] == "selected"
    assert captured["revision_id"] == "rev2"
    assert result[0]["window"]["broker_managed"] is True


def test_result_poll_uses_broker_terminal_state(monkeypatch):
    registry = CallbackRegistry()
    callbacks_runs.register_run_callbacks(registry)
    monkeypatch.setattr(broker_client, "perform_action", lambda *args, **kwargs: {
        "loss_runs": {"saved-replay": {"state": "success", "returncode": 0}},
    })
    runs, disabled = registry.functions["poll_result_loss_runs"](
        1, {"window": {"broker_managed": True, "run_id": "saved-replay", "state": "running"}},
    )
    assert runs["window"]["state"] == "success"
    assert disabled is True
