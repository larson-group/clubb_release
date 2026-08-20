"""Focused checks for the structured agent-to-Tune handoff."""

from pathlib import Path

from dash_app.shared import actions, activity
from dash_app.tune_tab.callbacks_runs import agent_request_to_tune_controls


def test_agent_tuning_launches_a_canonical_visible_request(tmp_path, monkeypatch):
    """The action validates a small request before any background worker starts."""
    monkeypatch.setattr(activity, "ACTIVITY_PATH", Path(tmp_path) / "activity.json")
    monkeypatch.setattr(activity, "LOCK_PATH", Path(tmp_path) / "activity.lock")
    activity.reset_activity()
    monkeypatch.setattr(actions, "available_tunable_configs", lambda: [{"value": "default"}])
    monkeypatch.setattr(actions, "load_tunable_names", lambda _config: ["C4", "C8"])
    monkeypatch.setattr(
        actions,
        "load_case_defaults",
        lambda: {"arm": {"default_clubb_fields": ["wp3"], "clubb_fields": ["wp3", "wprcp"]}},
    )
    monkeypatch.setattr(actions, "available_fields_for_cases", lambda _cases, _data: ["wp3", "wprcp"])
    monkeypatch.setattr(
        actions,
        "read_case_tuner_defaults",
        lambda _case, overrides=None: {
            "time_average_range": [41400, 45000],
            "altitude_comparison_range": [0.0, 3000.0],
            "average_time_seconds": 3600,
            "num_time_windows": 1,
        },
    )
    monkeypatch.setattr(
        actions,
        "start_tuning_job",
        lambda request, **_kwargs: {
            "pid": 1234,
            "job_dir": "/tmp/tune-job",
            "request_path": "/tmp/tune-job/request.json",
            "status_path": "/tmp/tune-job/status.json",
            "results_path": "/tmp/tune-job/results.json",
            "log_path": "/tmp/tune-job/worker.log",
        },
    )
    monkeypatch.setattr(actions, "_background", lambda *_args: None)

    result = actions.launch_tuning(
        [{"case_name": "arm", "time_average_range": [41400, 45000], "average_time_seconds": 3600}],
        [{"name": "C4", "min": 1.0, "max": 1.2}],
        ["wp3"],
        max_samples=2,
        overrides="-override iiPDF_type=1",
    )

    assert result["status"] == "started"
    assert result["request"]["override"] == "iiPDF_type=1"
    assert result["request"]["strategy"] == {"name": "random", "options": {"max_samples": 2}}
    assert result["request"]["case_configs"][0]["time_average_range"] == [41400, 45000]
    snapshot = activity.read_activity()
    assert snapshot["ui_request"]["type"] == "tune"
    assert snapshot["ui_request"]["request"] == result["request"]


def test_agent_can_launch_a_windowed_loss_run_from_the_broker_leaderboard(tmp_path, monkeypatch):
    monkeypatch.setattr(activity, "ACTIVITY_PATH", Path(tmp_path) / "activity.json")
    monkeypatch.setattr(activity, "LOCK_PATH", Path(tmp_path) / "activity.lock")
    activity.reset_activity()
    activity.set_broker_job(
        "tune",
        {
            "state": "finished",
            "job": {"status_path": "/tmp/status.json", "request_path": "/tmp/request.json"},
            "request": {
                "config": "default",
                "case_configs": [{"case_name": "arm", "time_average_range": [41400, 45000], "average_time_seconds": 3600}],
                "selected_fields": ["wp3"],
            },
        },
    )
    monkeypatch.setattr(
        actions,
        "read_tuning_status",
        lambda _path: {"top_results": [{"params": {"C4": 1.1}}, {"params": {"C4": 1.2}}]},
    )
    captured = {}

    def fake_start(cases, fields, params, **kwargs):
        captured.update({"cases": cases, "fields": fields, "params": params, **kwargs})
        return {"run_id": "loss-1", "pid": 4321, "state": "running", "log_path": "/tmp/loss.log"}

    monkeypatch.setattr(actions, "start_loss_run", fake_start)
    monkeypatch.setattr(actions, "_background", lambda *_args: None)

    result = actions.run_tuning_loss("window", max_results=1)

    assert result["status"] == "started"
    assert result["result_count"] == 1
    assert captured["cases"] == ["arm"]
    assert captured["fields"] == ["wp3"]
    assert captured["params"] == [{"C4": 1.1}]
    assert captured["run_mode"] == "window"


def test_agent_windowed_loss_run_keeps_every_selected_case(tmp_path, monkeypatch):
    """A multi-case leaderboard replay is not silently reduced to its first case."""
    monkeypatch.setattr(activity, "ACTIVITY_PATH", Path(tmp_path) / "activity.json")
    monkeypatch.setattr(activity, "LOCK_PATH", Path(tmp_path) / "activity.lock")
    activity.reset_activity()
    case_configs = [
        {"case_name": "arm", "time_average_range": [54000, 55800], "average_time_seconds": 1800},
        {"case_name": "bomex", "time_average_range": [10800, 12600], "average_time_seconds": 1800},
        {"case_name": "dycoms2_rf01", "time_average_range": [25260, 27060], "average_time_seconds": 1800},
    ]
    activity.set_broker_job(
        "tune",
        {
            "state": "finished",
            "job": {"status_path": "/tmp/status.json", "request_path": "/tmp/request.json"},
            "request": {"config": "default", "case_configs": case_configs, "selected_fields": ["wp3"]},
        },
    )
    monkeypatch.setattr(actions, "read_tuning_status", lambda _path: {"top_results": [{"params": {"C4": 1.1}}]})
    captured = {}

    def fake_start(cases, fields, params, **kwargs):
        captured.update({"cases": cases, "fields": fields, "params": params, **kwargs})
        return {"run_id": "loss-3case", "pid": 4321, "state": "running", "log_path": "/tmp/loss.log"}

    monkeypatch.setattr(actions, "start_loss_run", fake_start)
    monkeypatch.setattr(actions, "_background", lambda *_args: None)

    actions.run_tuning_loss("window", max_results=1)

    assert captured["cases"] == ["arm", "bomex", "dycoms2_rf01"]
    assert captured["case_configs"] == case_configs
    assert captured["run_mode"] == "window"


def test_agent_request_is_rendered_as_native_tune_controls(monkeypatch):
    """Exact agent ranges and windows populate the normal editable Tune widgets."""
    # The helper resolves discovery through its own module namespace.
    from dash_app.tune_tab import callbacks_runs

    monkeypatch.setattr(callbacks_runs, "load_tunable_names", lambda _config: ["C4"])
    monkeypatch.setattr(callbacks_runs, "load_tunable_default_ranges", lambda _config: {"C4": {"min": "0.5", "max": "2"}})
    monkeypatch.setattr(callbacks_runs, "available_fields_for_cases", lambda _cases, _data: ["wp3"])
    controls = agent_request_to_tune_controls(
        {
            "config": "default",
            "case_configs": [
                {
                    "case_name": "arm",
                    "time_average_range": [41400, 45000],
                    "altitude_comparison_range": [0.0, 3000.0],
                    "average_time_seconds": 3600,
                }
            ],
            "selected_fields": ["wp3"],
            "parameter_ranges": [{"name": "C4", "min": 1.0, "max": 1.2}],
            "strategy": {"name": "random", "options": {"max_samples": 2}},
            "batch_size": 1,
            "max_workers": 1,
            "override": "iiPDF_type=1",
        },
        {"arm": {"clubb_fields": ["wp3"]}},
    )

    assert controls["case_rows"][0]["case_name"] == "arm"
    assert controls["parameter_rows"] == [
        {"id": 0, "param": "C4", "targets": ["C4"], "min": "1.0", "max": "1.2"}
    ]
    assert controls["fields"] == ["wp3"]


def test_stale_agent_parameter_name_is_migrated_using_current_config(monkeypatch):
    """A durable request follows the selected config, not an embedded name list."""
    from dash_app.tune_tab import callbacks_runs

    monkeypatch.setattr(callbacks_runs, "load_tunable_names", lambda _config: ["C_uu_shr"])
    monkeypatch.setattr(callbacks_runs, "load_tunable_default_ranges", lambda _config: {})
    monkeypatch.setattr(callbacks_runs, "available_fields_for_cases", lambda _cases, _data: [])
    controls = agent_request_to_tune_controls(
        {
            "config": "default",
            "parameter_ranges": [{"name": "C_uu_shr", "min": 0.1, "max": 0.9}],
        },
        {},
    )

    assert controls["parameter_rows"][0]["param"] == "C_uu_shr"


def test_native_tune_request_is_relaunched_by_the_broker_without_agent_cap(monkeypatch):
    """The ordinary Tune button can keep its larger visible random budget."""
    captured = {}

    def fake_launch(cases, ranges, fields, **kwargs):
        captured.update({"cases": cases, "ranges": ranges, "fields": fields, **kwargs})
        return {"status": "started", "job": {"pid": 42}}

    monkeypatch.setattr(actions, "launch_tuning", fake_launch)
    result = actions.launch_tuning_request(
        {
            "config": "default",
            "case_configs": [{"case_name": "arm", "time_average_range": [41400, 45000], "altitude_comparison_range": [0, 3000], "average_time_seconds": 3600}],
            "parameter_ranges": [{"name": "C4", "min": 1.0, "max": 1.2}],
            "selected_fields": ["wp3"],
            "strategy": {"name": "random", "options": {"max_samples": 10000}},
            "batch_size": 8,
            "max_workers": 4,
            "loss_mode": "shape_first",
            "aggregation_mode": "mean_max",
            "override": "iiPDF_type=1",
        }
    )

    assert result["job"]["pid"] == 42
    assert captured["cases"][0]["case_name"] == "arm"
    assert captured["max_samples"] == 10000
    assert captured["_max_samples_limit"] is None


def test_native_adam_request_keeps_only_its_own_strategy_options(monkeypatch):
    captured = {}

    def fake_launch(cases, ranges, fields, **kwargs):
        captured.update(kwargs)
        return {"status": "started", "job": {"pid": 42}}

    monkeypatch.setattr(actions, "launch_tuning", fake_launch)
    actions.launch_tuning_request(
        {
            "case_configs": [{"case_name": "arm"}],
            "parameter_ranges": [{"name": "C4", "min": 1.0, "max": 1.2}],
            "selected_fields": ["wp3"],
            "strategy": {
                "name": "adam",
                "options": {
                    "max_updates": 17,
                    "learning_rate": 0.025,
                    "perturbation": 0.08,
                    "spsa_pairs": 3,
                },
            },
            "batch_size": 12,
            "max_workers": 4,
        }
    )

    assert captured["strategy"] == "adam"
    assert captured["adam_max_updates"] == 17
    assert captured["adam_learning_rate"] == 0.025
    assert captured["adam_perturbation"] == 0.08
    assert captured["adam_spsa_pairs"] == 3


def test_shared_tune_strategy_normalizer_builds_adam_options():
    strategy = actions._normalize_tune_strategy(
        "adam",
        [{"name": "C4", "min": 1.0, "max": 1.2}],
        max_samples=12,
        resolve_spacing=0.1,
        simann_max_iters=200,
        simann_initial_temp=1.0,
        simann_final_temp=1.0e-12,
        adam_max_updates=17,
        adam_learning_rate=0.025,
        adam_perturbation=0.08,
        adam_spsa_pairs=3,
        batch_size=12,
        max_workers=4,
    )

    assert strategy == {
        "name": "adam",
        "options": {
            "max_updates": 17,
            "learning_rate": 0.025,
            "perturbation": 0.08,
            "spsa_pairs": 3,
        },
    }


def test_broker_accepts_a_preset_sized_linked_range_family(monkeypatch):
    """The native Tune button can launch all 13 coordinates in ``wpxp``.

    This guards against reintroducing the old 12-range typed-agent cap, which
    made valid dashboard presets look like a broken Start button.
    """
    monkeypatch.setattr(actions, "load_tunable_names", lambda _config: [f"P{index}" for index in range(14)])
    ranges = [
        {"name": f"P{index}", "targets": [f"P{index}"], "min": 0.0, "max": 1.0}
        for index in range(13)
    ]
    ranges[0]["targets"] = ["P0", "P13"]

    normalized = actions._normalize_tune_ranges(ranges, "default")

    assert len(normalized) == 13
    assert normalized[0]["targets"] == ["P0", "P13"]


def test_broker_recovers_latest_running_tuning_job_from_disk(tmp_path, monkeypatch):
    monkeypatch.setattr(activity, "ACTIVITY_PATH", Path(tmp_path) / "activity.json")
    monkeypatch.setattr(activity, "LOCK_PATH", Path(tmp_path) / "activity.lock")
    activity.reset_activity()
    output_root = Path(tmp_path) / "output_tuner"
    job_dir = output_root / "arm_live"
    job_dir.mkdir(parents=True)
    (job_dir / "request.json").write_text('{"cases": ["arm"], "strategy": {"name": "random"}}')
    (job_dir / "status.json").write_text('{"state": "running", "job_dir": "ignored"}')
    (job_dir / "results.json").write_text('{}')
    (job_dir / "control.json").write_text('{}')
    (job_dir / "worker.log").write_text("worker still running\n")
    monkeypatch.setattr(actions, "REPO_ROOT", Path(tmp_path))
    monkeypatch.setattr(actions, "read_tuning_status", lambda _path: {"state": "running", "samples_evaluated": 3})
    resumed = []
    monkeypatch.setattr(actions, "_background", lambda target, *args: resumed.append((target, args)))

    record = actions.recover_active_tuning_from_disk()

    assert record["recovered"] is True
    assert record["job"]["pid"] is None
    assert record["job"]["job_dir"] == str(job_dir.resolve())
    assert activity.broker_jobs()["tune"]["state"] == "running"
    assert resumed and resumed[0][0] is actions._watch_tuning


def test_broker_discards_pidless_recovery_record_with_idle_status(tmp_path, monkeypatch):
    """An abandoned recovery record must not permanently block a new Tune job."""
    monkeypatch.setattr(activity, "ACTIVITY_PATH", Path(tmp_path) / "activity.json")
    monkeypatch.setattr(activity, "LOCK_PATH", Path(tmp_path) / "activity.lock")
    activity.reset_activity()
    activity.set_broker_job(
        "tune",
        {
            "state": "running",
            "recovered": True,
            "job": {"pid": None, "job_dir": "/tmp/old-tune", "status_path": "/tmp/old-tune/status.json"},
            "request": {"cases": ["arm"]},
        },
    )
    monkeypatch.setattr(actions, "REPO_ROOT", Path(tmp_path))
    monkeypatch.setattr(actions, "read_tuning_status", lambda _path: {"state": "idle", "job_dir": None})

    assert actions.recover_active_tuning_from_disk() is None
    record = activity.broker_jobs()["tune"]
    assert record["state"] == "error"
    assert "no worker PID" in record["recovery_error"]


def test_broker_discards_recorded_tune_job_when_its_pid_has_exited(tmp_path, monkeypatch):
    """A cancelled worker cannot leave its broker record blocking the next start."""
    monkeypatch.setattr(activity, "ACTIVITY_PATH", Path(tmp_path) / "activity.json")
    monkeypatch.setattr(activity, "LOCK_PATH", Path(tmp_path) / "activity.lock")
    activity.reset_activity()
    activity.set_broker_job(
        "tune",
        {
            "state": "stopping",
            "job": {"pid": 1234, "job_dir": "/tmp/cancelled-tune", "status_path": "/tmp/cancelled-tune/status.json"},
            "request": {"cases": ["arm"]},
        },
    )
    monkeypatch.setattr(actions, "REPO_ROOT", Path(tmp_path))
    monkeypatch.setattr(actions, "read_tuning_status", lambda _path: {"state": "draft", "job_dir": "/tmp/cancelled-tune"})
    monkeypatch.setattr(actions, "active_job_exited", lambda _job: True)

    assert actions.recover_active_tuning_from_disk() is None
    record = activity.broker_jobs()["tune"]
    assert record["state"] == "error"
    assert "PID 1234 has exited" in record["recovery_error"]
