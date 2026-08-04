"""Focused regression tests for typed MCP service primitives."""

import json
import io
import multiprocessing
import os
import select
import stat
import subprocess
import sys
from pathlib import Path

import pytest
from pydantic import ValidationError

from dash_app.services.jobs import ArtifactStore, JobStore, SubmissionConflict
from dash_app.services.models import CompileRequest, LeaderboardRerunRequest, ParameterRange, ScmRunBatchRequest, ScmRunRequest, TuneRequest
from dash_app.services import profiles as profile_service
from dash_app.shared.runtime import atomic_write_json, private_runtime_dir, restrict_existing
from dash_app.shared.components import styled_dropdown


def _concurrent_job_submit(path, start, results):
    """Process target used to prove the persistent request lock, not threads."""
    store = JobStore(Path(path))
    start.wait(5)
    record, created = store.submit("scm", "concurrent-request-123", {"case": "arm"})
    results.put((record["job_id"], created))


def test_job_request_id_is_idempotent_and_conflicts(tmp_path):
    store = JobStore(tmp_path / "jobs.json")
    first, created = store.submit("scm", "request-123", {"case": "arm"})
    again, repeated = store.submit("scm", "request-123", {"case": "arm"})

    assert created is True
    assert repeated is False
    assert again["job_id"] == first["job_id"]
    with pytest.raises(SubmissionConflict):
        store.submit("scm", "request-123", {"case": "bomex"})


def test_job_request_id_lock_serializes_independent_processes(tmp_path):
    # CPython 3.14 defaults to forkserver here, but the test sandbox forbids
    # binding its private Unix socket.  The lock behavior itself is POSIX
    # process-level, so an explicit fork context tests the relevant boundary.
    context = multiprocessing.get_context("fork")
    start = context.Event()
    results = context.Queue()
    workers = [
        context.Process(target=_concurrent_job_submit, args=(str(tmp_path / "jobs.json"), start, results))
        for _ in range(2)
    ]
    for worker in workers:
        worker.start()
    start.set()
    received = [results.get(timeout=5) for _ in workers]
    for worker in workers:
        worker.join(timeout=5)
        assert worker.exitcode == 0
    assert {job_id for job_id, _created in received}
    assert len({job_id for job_id, _created in received}) == 1
    assert sum(created for _job_id, created in received) == 1


def test_strict_public_scm_request_rejects_paths_and_unknown_fields():
    with pytest.raises(ValidationError):
        ScmRunRequest(request_id="request-123", case="arm", stats_file="../../secret")
    with pytest.raises(ValidationError):
        ScmRunRequest(request_id="request-123", case="arm", unexpected=True)
    with pytest.raises(ValidationError):
        ScmRunRequest(
            request_id="request-123",
            case="arm",
            run_options={"max_iters": 0},
        )


def test_strict_public_scm_batch_request_rejects_duplicate_cases_and_paths():
    with pytest.raises(ValidationError):
        ScmRunBatchRequest(request_id="batch-request-123", cases=["arm", "arm"])
    with pytest.raises(ValidationError):
        ScmRunBatchRequest(request_id="batch-request-123", cases=["arm"], stats_file="../secret")


def test_mcp_output_dir_stays_below_repository_output_root(tmp_path, monkeypatch):
    from dash_app.shared import actions

    monkeypatch.setattr(actions, "REPO_ROOT", tmp_path)

    assert actions._resolve_mcp_output_dir("dash_default") == tmp_path / "output" / "dash_default"
    assert actions._resolve_mcp_output_dir("output/nested/run") == tmp_path / "output" / "nested" / "run"
    with pytest.raises(ValueError, match="inside the repository output"):
        actions._resolve_mcp_output_dir("../outside")
    with pytest.raises(ValueError, match="inside the repository output"):
        actions._resolve_mcp_output_dir(str(tmp_path / "outside"))


def test_scm_batch_submission_uses_one_flat_output_and_is_idempotent(tmp_path, monkeypatch):
    from dash_app.shared import actions

    store = JobStore(tmp_path / "jobs.json")
    artifacts = ArtifactStore(tmp_path / "agent_artifacts")
    calls = []

    def fake_run(case, overrides, config, stats, *, cli_options, output_dir, job_id, run_id):
        assert (tmp_path / "agent_artifacts" / job_id / "manifest.json").is_file()
        Path(output_dir, f"{case}_stats.nc").write_bytes(b"CDF\x01batch-run")
        calls.append((case, output_dir, job_id, run_id))
        return {
            "status": "started",
            "case": case,
            "log": str(tmp_path / f"{case}.log"),
            "output_directory": output_dir,
            "proc_data": {"pid": 100 + len(calls), "log": str(tmp_path / f"{case}.log")},
        }

    monkeypatch.setattr(actions, "_public_scm_output_root", lambda: tmp_path / "output" / "mcp_runs")
    monkeypatch.setattr(actions, "_JOB_STORE", store)
    monkeypatch.setattr(actions, "_ARTIFACT_STORE", artifacts)
    monkeypatch.setattr(actions, "run_scm", fake_run)
    monkeypatch.setattr(actions, "is_case_active", lambda _case: False)
    monkeypatch.setattr(actions, "broker_jobs", lambda: {"runs": {}})
    monkeypatch.setattr(actions, "snapshot_active_cases", lambda: {})

    request = ScmRunBatchRequest(request_id="batch-request-123", cases=["arm", "bomex"], max_workers=2)
    first = actions.submit_scm_batch(request)
    retry = actions.submit_scm_batch(request)

    assert first["status"] == "started"
    assert retry["status"] == "existing"
    assert len(calls) == 2
    assert len({path for _case, path, _job, _run in calls}) == 1
    output_directory = Path(calls[0][1])
    assert output_directory == tmp_path / "output" / "mcp_runs" / first["batch_id"]
    assert not (output_directory / "arm").is_dir()
    assert {path.name for path in output_directory.glob("*_stats.nc")} == {"arm_stats.nc", "bomex_stats.nc"}
    manifest = json.loads((output_directory / "batch_manifest.json").read_text())
    assert manifest["batch_id"] == first["batch_id"]
    assert {child["case"] for child in manifest["children"]} == {"arm", "bomex"}
    with pytest.raises(ValidationError):
        TuneRequest(
            request_id="request-123",
            cases=[{"case_name": "arm", "unrecognized": True}],
            parameter_ranges=[{"name": "c_K10", "min": 0.0, "max": 1.0}],
        )


def test_native_scm_batch_honors_selected_output_directory(tmp_path, monkeypatch):
    from dash_app.shared import actions

    store = JobStore(tmp_path / "jobs.json")
    artifacts = ArtifactStore(tmp_path / "agent_artifacts")
    selected_output = tmp_path / "output" / "dash_default"
    calls = []

    def fake_launch(case, stats, config, overrides, cli_options, *, job_id, run_id):
        calls.append((case, stats, config, overrides, cli_options, job_id, run_id))
        Path(cli_options["out_dir"]).mkdir(parents=True, exist_ok=True)
        Path(cli_options["out_dir"], f"{case}_stats.nc").write_bytes(b"CDF\x01native-run")
        return {
            "status": "started",
            "case": case,
            "output_directory": cli_options["out_dir"],
            "proc_data": {"pid": 100 + len(calls), "log": str(tmp_path / f"{case}.log")},
        }

    monkeypatch.setattr(actions, "_JOB_STORE", store)
    monkeypatch.setattr(actions, "_ARTIFACT_STORE", artifacts)
    monkeypatch.setattr(actions, "resolve_output_dir", lambda value: tmp_path / "output" / str(value))
    monkeypatch.setattr(actions, "launch_dashboard_scm_request", fake_launch)
    monkeypatch.setattr(actions, "is_case_active", lambda _case: False)
    monkeypatch.setattr(actions, "broker_jobs", lambda: {"runs": {}})
    monkeypatch.setattr(actions, "snapshot_active_cases", lambda: {})

    request = ScmRunBatchRequest(request_id="native-output-request", cases=["arm"], max_workers=1)
    result = actions.submit_scm_batch(
        request,
        native_overrides={},
        native_cli_options={"out_dir": "dash_default"},
    )

    assert result["status"] == "started"
    assert calls[0][4]["out_dir"] == str(selected_output)
    assert Path(result["output_directory"]) == selected_output
    assert (selected_output / "arm_stats.nc").is_file()
    assert (selected_output / "batch_manifest.json").is_file()


def test_mcp_scm_batch_honors_requested_output_directory(tmp_path, monkeypatch):
    from dash_app.shared import actions

    store = JobStore(tmp_path / "jobs.json")
    artifacts = ArtifactStore(tmp_path / "agent_artifacts")
    selected_output = tmp_path / "output" / "mcp_named"

    monkeypatch.setattr(actions, "_JOB_STORE", store)
    monkeypatch.setattr(actions, "_ARTIFACT_STORE", artifacts)
    monkeypatch.setattr(actions, "_resolve_mcp_output_dir", lambda value: selected_output)
    monkeypatch.setattr(actions, "_start_queued_scm_batch_children", lambda _job_id: None)
    monkeypatch.setattr(actions, "_ensure_scm_batch_watcher", lambda _job_id: None)

    request = ScmRunBatchRequest(
        request_id="mcp-output-request",
        cases=["arm", "bomex"],
        out_dir="mcp_named",
    )
    result = actions.submit_scm_batch(request)

    assert result["status"] == "queued"
    assert Path(result["output_directory"]) == selected_output
    assert (selected_output / "batch_manifest.json").is_file()


def test_typed_compile_request_exposes_validated_gptl_option():
    request = CompileRequest(request_id="compile-request-123", gptl=True)

    assert request.gptl is True
    assert request.model_dump()["gptl"] is True


def test_compile_wrapper_passes_gptl_through_shared_launcher(monkeypatch):
    from dash_app.shared import actions

    captured = {}

    def fake_launch(options, **kwargs):
        captured.update(options=options, kwargs=kwargs)
        return {"status": "planned"}

    monkeypatch.setattr(actions, "launch_compile_request", fake_launch)

    result = actions.compile_clubb(gptl=True, job_id="compile-job-123")

    assert result == {"status": "planned"}
    assert captured["options"]["gptl"] is True
    assert captured["kwargs"] == {"job_id": "compile-job-123"}


def test_scm_rejects_a_config_not_exposed_by_the_dashboard():
    from dash_app.shared import actions

    with pytest.raises(ValueError, match="unknown tunable configuration"):
        actions.run_scm("arm", config="not_a_checked_in_config")


def test_typed_scm_run_options_reuse_the_native_process_launcher(monkeypatch):
    from dash_app.shared import actions

    captured = {}

    def fake_launch(case, stats_file, config, overrides, cli_options, **kwargs):
        captured.update(
            case=case,
            stats_file=stats_file,
            config=config,
            overrides=overrides,
            cli_options=cli_options,
            kwargs=kwargs,
        )
        return {"status": "started"}

    monkeypatch.setattr(actions, "launch_dashboard_scm_request", fake_launch)

    result = actions.run_scm(
        "bomex",
        cli_options={"max_iters": 2, "dt_main": 30.0, "tout": 60.0},
    )

    assert result["status"] == "started"
    assert captured["cli_options"] == {
        "max_iters": "2",
        "dt_main": "30.0",
        "tout": "60.0",
    }


def test_mcp_scm_run_accepts_any_directory_below_output(monkeypatch, tmp_path):
    from dash_app.shared import actions

    captured = {}
    fake_repo = tmp_path / "repo"
    selected_output = fake_repo / "output" / "mcp_named"

    def fake_launch(case, stats_file, config, overrides, cli_options, **kwargs):
        captured.update(cli_options=cli_options)
        return {"status": "started"}

    monkeypatch.setattr(actions, "REPO_ROOT", fake_repo)
    monkeypatch.setattr(actions, "_validate_case", lambda value: str(value))
    monkeypatch.setattr(actions, "_validated_stats_file", lambda value: str(value))
    monkeypatch.setattr(actions, "launch_dashboard_scm_request", fake_launch)

    result = actions.run_scm("arm", output_dir=selected_output)

    assert result["status"] == "started"
    assert captured["cli_options"]["out_dir"] == str(selected_output)


def test_typed_run_and_tune_requests_share_settings_consequences():
    """MCP submissions must use the same implications as the browser tabs."""
    from dash_app.shared import actions

    scm, scm_resolution = actions._canonical_scm_request(
        ScmRunRequest(
            request_id="typed-scm-settings-123",
            case="arm",
            overrides={"iiPDF_type": 1, "l_predict_upwp_vpwp": True},
        )
    )
    assert scm.overrides["l_predict_upwp_vpwp"] is True
    assert scm_resolution["parameter_states"]["C_uu_shr"]["state"] == "active"
    assert actions._typed_override_text(scm.overrides).count(",") >= 1

    tune = actions._canonical_tune_request(
        TuneRequest(
            request_id="typed-tune-settings-123",
            cases=["arm"],
            fields=["cloud_frac"],
            parameter_ranges=[ParameterRange(name="C6rt", min=0.0, max=1.0)],
            overrides={"iiPDF_type": 1, "l_predict_upwp_vpwp": False},
        )
    )
    assert tune.parameter_ranges[0].targets == ["C6rt", "C6thl"]

    with pytest.raises(ValueError, match="inactive parameter"):
        actions._canonical_tune_request(
            TuneRequest(
                request_id="typed-tune-inactive-123",
                cases=["arm"],
                fields=["cloud_frac"],
                parameter_ranges=[ParameterRange(name="slope_coef_spread_DG_means_w", min=0.1, max=0.2)],
            )
        )


def test_tune_parameter_manifest_and_range_normalization_follow_selected_config(monkeypatch):
    """Agents discover names from the config rather than carrying a static list."""
    from dash_app.shared import actions

    monkeypatch.setattr(actions, "available_tunable_configs", lambda: [{"value": "default"}])
    monkeypatch.setattr(actions, "load_tunable_names", lambda _config: ["C_uu_shr"])
    monkeypatch.setattr(
        actions,
        "load_tunable_default_ranges",
        lambda _config: {"C_uu_shr": {"default": "0.5", "min": "0.0", "max": "1.0"}},
    )

    manifest = actions.list_tunable_parameters("default")
    canonical = actions._canonical_tune_request(
        TuneRequest(
            request_id="canonical-tune-123",
            cases=["arm"],
            parameter_ranges=[{"name": "C_uu_shr", "min": 0.1, "max": 0.9}],
        )
    )
    ranges = actions._normalize_tune_ranges(
        [{"name": "C_uu_shr", "min": 0.1, "max": 0.9}],
        "default",
    )

    assert manifest == {
        "config": "default",
        "parameters": [{"name": "C_uu_shr", "default": "0.5", "suggested_min": "0.0", "suggested_max": "1.0"}],
    }
    assert canonical.config == "default"
    assert canonical.parameter_ranges[0].name == "C_uu_shr"
    assert ranges == [{"name": "C_uu_shr", "targets": ["C_uu_shr"], "min": 0.1, "max": 0.9}]


def test_private_runtime_json_and_artifact_manifest_are_bounded(tmp_path, monkeypatch):
    runtime = private_runtime_dir(tmp_path)
    assert stat.S_IMODE(runtime.stat().st_mode) == 0o700
    payload_path = runtime / "state.json"
    atomic_write_json(payload_path, {"ok": True})
    assert stat.S_IMODE(payload_path.stat().st_mode) == 0o600
    assert json.loads(payload_path.read_text()) == {"ok": True}
    legacy = tmp_path / "legacy.json"
    legacy.write_text("{}")
    legacy.chmod(0o644)
    restrict_existing(legacy)
    assert stat.S_IMODE(legacy.stat().st_mode) == 0o600

    artifacts = ArtifactStore(tmp_path / "agent_artifacts")
    manifest = artifacts.create_manifest("artifact_123", {"result": {"state": "finished"}})
    assert stat.S_IMODE((tmp_path / "agent_artifacts").stat().st_mode) == 0o700
    assert stat.S_IMODE(manifest.parent.stat().st_mode) == 0o700
    assert stat.S_IMODE(manifest.stat().st_mode) == 0o600
    assert artifacts.get_manifest("artifact_123")["result"]["state"] == "finished"
    with pytest.raises(FileExistsError):
        artifacts.create_manifest("artifact_123", {})
    assert manifest.is_file()


def test_shared_profile_window_matches_plot_time_grid_and_dropdown_contract():
    from dash_app.plot_tab.plot_types.shared import snap_start_time_to_step

    case_data = {
        "time_slider_duration_min_minutes": 1,
        "time_slider_duration_max_minutes": 60,
        "time_slider_duration_step_minutes": 1,
        "time_slider_start_min_seconds": 0,
        "time_slider_start_max_seconds": 6600,
        "time_slider_final_end_seconds": 7200,
        "default_time_start_seconds": 0,
    }
    selected = profile_service.resolve_time_window(case_data, start_seconds=173, average_minutes=10)
    assert selected["time_start_seconds"] == snap_start_time_to_step(case_data, 173, 10)
    dropdown = styled_dropdown(id="shared-dropdown", options=[{"label": "A", "value": "a"}], value="a", clearable=False)
    assert dropdown.id == "shared-dropdown"
    assert dropdown.value == "a"
    assert dropdown.clearable is False
    assert "clubb-dropdown" in dropdown.className


def test_artifact_retention_prunes_only_completed_manifest_bundles(tmp_path, monkeypatch):
    monkeypatch.setenv("CLUBB_AGENT_ARTIFACT_RETENTION_COUNT", "1")
    artifacts = ArtifactStore(tmp_path / "agent_artifacts")
    artifacts.create_manifest("first", {"state": "finished"})
    artifacts.create_manifest("second", {"state": "finished"})
    assert not (tmp_path / "agent_artifacts" / "first").exists()
    assert (tmp_path / "agent_artifacts" / "second" / "manifest.json").is_file()


def test_active_artifact_bundle_survives_retention_until_released(tmp_path, monkeypatch):
    monkeypatch.setenv("CLUBB_AGENT_ARTIFACT_RETENTION_COUNT", "1")
    artifacts = ArtifactStore(tmp_path / "agent_artifacts")
    artifacts.create_manifest("active", {"state": "running"})
    artifacts.activate("active")
    artifacts.create_manifest("completed", {"state": "finished"})
    artifacts.create_manifest("latest", {"state": "finished"})

    assert (tmp_path / "agent_artifacts" / "active" / ".active").is_file()
    assert (tmp_path / "agent_artifacts" / "active" / "manifest.json").is_file()
    assert not (tmp_path / "agent_artifacts" / "completed").exists()
    assert (tmp_path / "agent_artifacts" / "latest").is_dir()

    artifacts.release("active")
    artifacts.prune()
    assert not (tmp_path / "agent_artifacts" / "active").exists()


def test_artifact_root_is_recreated_after_disposable_cleanup(tmp_path):
    artifacts = ArtifactStore(tmp_path / "agent_artifacts")
    artifacts.create_manifest("first", {})
    import shutil

    shutil.rmtree(tmp_path / "agent_artifacts")
    manifest = artifacts.create_manifest("second", {})
    assert manifest.is_file()
    assert artifacts.get_manifest("second")["lifecycle"]["storage"] == "ephemeral_staging"


def test_artifact_resource_rejects_traversal(monkeypatch, tmp_path):
    from dash_app.shared import actions

    artifacts = ArtifactStore(tmp_path / "agent_artifacts")
    artifacts.create_manifest("artifact_123", {})
    monkeypatch.setattr(actions, "_ARTIFACT_STORE", artifacts)
    with pytest.raises(ValueError):
        actions.read_artifact_file("artifact_123", "../manifest.json")
    with pytest.raises(ValueError):
        actions.get_artifact("../outside")


def test_typed_scm_submission_isolated_and_retries_without_relaunch(tmp_path, monkeypatch):
    from dash_app.shared import actions

    store = JobStore(tmp_path / "jobs.json")
    artifacts = ArtifactStore(tmp_path / "agent_artifacts")
    calls = []

    def fake_run(case, overrides, config, stats, *, cli_options, output_dir, job_id, run_id):
        assert (tmp_path / "agent_artifacts" / job_id / "manifest.json").is_file()
        assert cli_options == {}
        Path(output_dir, f"{case}_stats.nc").write_bytes(b"CDF\x01public-run")
        calls.append((case, output_dir, job_id))
        return {"status": "started", "case": case, "log": str(tmp_path / "run.log"), "output_directory": output_dir}

    monkeypatch.setattr(actions, "_public_scm_output_root", lambda: tmp_path / "output" / "mcp_runs")
    monkeypatch.setattr(actions, "_JOB_STORE", store)
    monkeypatch.setattr(actions, "_ARTIFACT_STORE", artifacts)
    monkeypatch.setattr(actions, "run_scm", fake_run)
    monkeypatch.setattr(actions, "is_case_active", lambda _case: False)
    monkeypatch.setattr(actions, "broker_jobs", lambda: {"runs": {}})
    monkeypatch.setattr(actions, "snapshot_active_cases", lambda: {})
    monkeypatch.setattr(
        actions,
        "exclusive_file_lock",
        lambda _path: (_ for _ in ()).throw(
            AssertionError("typed submission must let the shared launcher own the admission lock")
        ),
    )
    request = ScmRunRequest(request_id="request-typed-scm", case="arm")

    first = actions.submit_scm_run(request)
    retry = actions.submit_scm_run(request)

    assert first["status"] == "started"
    assert retry["status"] == "existing"
    assert len(calls) == 1
    batch_directory = Path(calls[0][1])
    assert batch_directory.name == first["batch_id"]
    assert batch_directory.parent == tmp_path / "output" / "mcp_runs"
    manifest = json.loads((tmp_path / "agent_artifacts" / first["job_id"] / "manifest.json").read_text())
    assert manifest["run_id"] == first["run_id"]
    assert manifest["inputs"]["case_setup"]["sha256"]
    assert manifest["execution"]["build_identity"]["precision"] == "double"
    assert manifest["execution"]["state"] == "planned"
    assert manifest["service"]["version"] == 3
    assert manifest["output_checksums"]["state"] == "pending"
    assert (tmp_path / "agent_artifacts" / first["job_id"] / "submission_result.json").is_file()
    public_manifest = json.loads((batch_directory / "batch_manifest.json").read_text())
    assert public_manifest["batch_id"] == first["batch_id"]
    assert public_manifest["children"][0]["run_id"] == first["run_id"]
    discovered = profile_service.discover_output_directories(tmp_path / "output")
    assert [item["path"] for item in discovered] == [calls[0][1]]


def test_same_case_submissions_get_distinct_durable_batch_roots(tmp_path, monkeypatch):
    from dash_app.shared import actions

    calls = []
    monkeypatch.setattr(actions, "_public_scm_output_root", lambda: tmp_path / "output" / "mcp_runs")
    monkeypatch.setattr(actions, "_JOB_STORE", JobStore(tmp_path / "jobs.json"))
    monkeypatch.setattr(actions, "_ARTIFACT_STORE", ArtifactStore(tmp_path / "agent_artifacts"))
    monkeypatch.setattr(actions, "is_case_active", lambda _case: False)
    monkeypatch.setattr(actions, "broker_jobs", lambda: {"runs": {}})
    monkeypatch.setattr(actions, "snapshot_active_cases", lambda: {})
    monkeypatch.setattr(
        actions,
        "run_scm",
        lambda _case, _overrides, _config, _stats, *, cli_options, output_dir, job_id, run_id: calls.append((job_id, output_dir)) or {"status": "started", "output_directory": output_dir},
    )
    first = actions.submit_scm_run(ScmRunRequest(request_id="case-isolation-one", case="arm"))
    second = actions.submit_scm_run(ScmRunRequest(request_id="case-isolation-two", case="arm"))
    assert first["run_id"] != second["run_id"]
    assert calls[0][1] != calls[1][1]
    assert all(Path(path).name.startswith("batch_") for _job_id, path in calls)
    assert all(Path(path).parent == tmp_path / "output" / "mcp_runs" for _job_id, path in calls)


def test_typed_scm_submission_respects_global_run_concurrency(monkeypatch, tmp_path):
    from dash_app.shared import actions

    monkeypatch.setattr(actions, "_JOB_STORE", JobStore(tmp_path / "jobs.json"))
    monkeypatch.setattr(actions, "is_case_active", lambda _case: False)
    monkeypatch.setattr(actions, "snapshot_active_cases", lambda: {})
    monkeypatch.setattr(actions, "broker_jobs", lambda: {"runs": {"bomex": {"state": "running"}}})
    monkeypatch.setattr(actions, "MAX_RUN_PROCS", 1)
    result = actions.submit_scm_run(ScmRunRequest(request_id="request-cap-123", case="arm"))
    assert result["status"] == "queued"


def test_scm_concurrency_counts_distinct_native_and_broker_processes(monkeypatch, tmp_path):
    from dash_app.shared import actions

    monkeypatch.setattr(actions, "_JOB_STORE", JobStore(tmp_path / "jobs.json"))
    monkeypatch.setattr(actions, "is_case_active", lambda _case: False)
    monkeypatch.setattr(
        actions,
        "broker_jobs",
        lambda: {"runs": {"arm": {"state": "running", "proc_data": {"pid": 1001}}}},
    )
    monkeypatch.setattr(actions, "snapshot_active_cases", lambda: {"bomex": {"pid": 1002}})
    monkeypatch.setattr(actions, "MAX_RUN_PROCS", 2)
    result = actions.submit_scm_run(ScmRunRequest(request_id="request-cap-union-123", case="rico"))
    assert result["status"] == "queued"


def test_native_run_and_rebuild_use_the_durable_broker_lifecycle(monkeypatch):
    """Native UI launches may retain rich controls, but not separate ownership."""
    from dash_app.shared import actions

    events = []
    broker_compile = []
    broker_runs = []
    watched = []
    monkeypatch.setattr(actions, "broker_jobs", lambda: {"compile": None, "runs": {}})
    monkeypatch.setattr(actions, "snapshot_active_cases", lambda: {})
    monkeypatch.setattr(actions, "is_case_active", lambda _case: False)
    monkeypatch.setattr(actions, "publish_event", lambda *args, **kwargs: events.append((args, kwargs)) or {"id": 1})
    monkeypatch.setattr(actions, "publish_run_request", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(actions, "set_broker_job", lambda kind, payload: broker_compile.append((kind, payload)))
    monkeypatch.setattr(actions, "set_broker_run", lambda case, payload: broker_runs.append((case, payload)))
    monkeypatch.setattr(actions, "_background", lambda target, *args: watched.append((target.__name__, args)))
    monkeypatch.setattr(
        actions,
        "start_rebuild_job",
        lambda builds, discovery, label: {"pid": 91, "log": "/tmp/rebuild.log", "command": f"rebuild {label}", "status": "running"},
    )
    rebuild = actions.launch_rebuild_request([{"path": "/tmp/build", "name": "debug"}], {}, "debug")
    assert rebuild["status"] == "started"
    assert broker_compile[-1][0] == "compile"
    assert broker_compile[-1][1]["operation"] == "rebuild"

    monkeypatch.setattr(
        actions,
        "start_case_process",
        lambda *_args: {"pid": 92, "log": "/tmp/run.log", "start_time": 0.0, "temp_files": []},
    )
    launched = actions.launch_dashboard_scm_request(
        "arm",
        "standard_stats.in",
        "default",
        {"flags": {"l_test": ".true."}, "tunable": {}, "silhs": {}},
        {"max_iters": "2", "extra_args": ["-override", "iiPDF_type=1"]},
    )
    assert launched["proc_data"]["pid"] == 92
    assert broker_runs[-1][0] == "arm"
    assert broker_runs[-1][1]["broker_managed"] is True
    assert {name for name, _args in watched} == {"_watch_compile", "_watch_run"}


def test_mcp_exposes_closed_world_typed_tools_and_annotations():
    from dash_app.agent_integration.mcp_server import create_server

    import asyncio

    tools = asyncio.run(create_server().list_tools())
    by_name = {tool.name: tool for tool in tools}
    assert {"list_tunable_parameters", "submit_compile", "submit_scm_run", "submit_scm_batch", "submit_tune", "submit_leaderboard_rerun", "create_profile_artifact", "cancel_job", "list_plots", "remove_plot"} <= set(by_name)
    assert by_name["get_job"].annotations.readOnlyHint is True
    assert by_name["cancel_job"].annotations.destructiveHint is True
    assert by_name["submit_scm_run"].annotations.idempotentHint is True
    assert by_name["list_tunable_parameters"].annotations.readOnlyHint is True
    assert by_name["connect_to_dashboard"].annotations.readOnlyHint is True
    assert by_name["list_plots"].annotations.readOnlyHint is True
    schema = by_name["submit_scm_run"].inputSchema
    assert schema["$defs"]["ScmRunRequest"]["additionalProperties"] is False
    batch_schema = by_name["submit_scm_batch"].inputSchema
    assert batch_schema["$defs"]["ScmRunBatchRequest"]["additionalProperties"] is False
    compile_schema = by_name["submit_compile"].inputSchema["$defs"]["CompileRequest"]
    assert compile_schema["properties"]["gptl"] == {"default": False, "title": "Gptl", "type": "boolean"}


def test_plot_mcp_schema_exposes_output_directory_selector():
    from dash_app.agent_integration.mcp_server import create_server

    import asyncio

    tools = asyncio.run(create_server().list_tools())
    plot_schema = next(tool for tool in tools if tool.name == "add_budget_plot").inputSchema
    assert "output_dir" in plot_schema["properties"]


def test_mcp_mutations_are_forwarded_to_the_detached_broker(monkeypatch):
    from dash_app.agent_integration import mcp_server

    captured = {}

    def fake_action(action, payload, **kwargs):
        captured.update(action=action, payload=payload, kwargs=kwargs)
        return {"status": "started", "job_id": "scm_detached"}

    monkeypatch.setattr(mcp_server, "perform_action", fake_action)

    result = mcp_server._broker_domain_result(
        "submit_scm_run",
        {"request": {"request_id": "detached-owner-123", "case": "arm"}},
    )

    assert result["job_id"] == "scm_detached"
    assert captured["action"] == "domain_submit_scm_run"
    assert captured["payload"]["request"]["case"] == "arm"
    assert captured["kwargs"] == {"internal": True, "timeout_seconds": 120.0}


def test_broker_domain_dispatch_revalidates_typed_requests(monkeypatch):
    from dash_app.shared import actions

    captured = {}

    def fake_submit(request):
        captured["request"] = request
        return {"status": "started"}

    monkeypatch.setattr(actions, "submit_scm_run", fake_submit)

    result = actions.dispatch(
        "domain_submit_scm_run",
        {"request": {"request_id": "broker-validation-123", "case": "bomex"}},
    )

    assert result["status"] == "started"
    assert isinstance(captured["request"], ScmRunRequest)
    assert captured["request"].case == "bomex"
    with pytest.raises(ValidationError):
        actions.dispatch(
            "domain_submit_scm_run",
            {
                "request": {
                    "request_id": "broker-validation-bad-123",
                    "case": "bomex",
                    "unexpected": True,
                }
            },
        )


def test_broker_domain_dispatch_revalidates_batch_requests(monkeypatch):
    from dash_app.shared import actions

    captured = {}

    def fake_submit(request, **kwargs):
        captured["request"] = request
        captured["kwargs"] = kwargs
        return {"status": "queued"}

    monkeypatch.setattr(actions, "submit_scm_batch", fake_submit)
    result = actions.dispatch(
        "domain_submit_scm_batch",
        {
            "request": {
                "request_id": "batch-dispatch-123",
                "cases": ["arm", "bomex"],
            }
        },
    )
    assert result["status"] == "queued"
    assert isinstance(captured["request"], ScmRunBatchRequest)
    assert captured["request"].cases == ["arm", "bomex"]
    assert captured["kwargs"] == {"native_overrides": None, "native_cli_options": None}
    with pytest.raises(ValidationError):
        actions.dispatch(
            "domain_submit_scm_batch",
            {"request": {"request_id": "batch-dispatch-123", "cases": ["arm", "arm"]}},
        )


def test_leaderboard_rerun_is_idempotent_and_cancellable(monkeypatch, tmp_path):
    from dash_app.shared import actions

    store = JobStore(tmp_path / "jobs.json")
    artifacts = ArtifactStore(tmp_path / "agent_artifacts")
    launched = []
    stopped = []
    monkeypatch.setattr(actions, "_JOB_STORE", store)
    monkeypatch.setattr(actions, "_ARTIFACT_STORE", artifacts)
    monkeypatch.setattr(
        actions,
        "run_tuning_loss",
        lambda mode, max_results, *, job_id: launched.append((mode, max_results, job_id)) or {"status": "started", "run": {"run_id": "loss_123", "log_path": str(tmp_path / "loss.log")}},
    )
    request = LeaderboardRerunRequest(request_id="leaderboard-request-123", mode="complete", max_results=2)
    first = actions.submit_leaderboard_rerun(request)
    retry = actions.submit_leaderboard_rerun(request)
    assert first["status"] == "started"
    assert retry["status"] == "existing"
    assert len(launched) == 1

    monkeypatch.setattr(actions, "broker_jobs", lambda: {"loss_runs": {"loss_123": {"run_id": "loss_123", "state": "running"}}})
    monkeypatch.setattr(actions, "stop_loss_run", lambda run: stopped.append(run) or {**run, "state": "stopping"})
    monkeypatch.setattr(actions, "update_broker_loss_run", lambda *_args, **_kwargs: None)
    cancelled = actions.cancel_job(first["job_id"])
    assert cancelled["status"] == "stop_requested"
    assert stopped and stopped[0]["run_id"] == "loss_123"


def test_bounded_log_reader_uses_nested_leaderboard_runtime(tmp_path, monkeypatch):
    from dash_app.shared import actions

    log = tmp_path / "leaderboard.log"
    log.write_text("leaderboard output", encoding="utf-8")
    store = JobStore(tmp_path / "jobs.json")
    record, _created = store.submit("leaderboard", "leaderboard-log-123", {"request_id": "leaderboard-log-123"})
    store.update(record["job_id"], runtime={"run": {"log_path": str(log)}})
    monkeypatch.setattr(actions, "_JOB_STORE", store)
    read = actions.read_job_log(record["job_id"], max_bytes=8)
    assert read["text"] == "leaderbo"
    assert read["next_cursor"] == 8


def test_mcp_stdio_discovery_and_structured_error_round_trip():
    root = Path(__file__).resolve().parents[2]
    server_path = root / "dash_app" / "agent_integration" / "mcp_server.py"
    process = subprocess.Popen(
        [sys.executable, str(server_path)],
        cwd=root,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        env={**os.environ, "PYTHONPATH": str(root) + os.pathsep + os.environ.get("PYTHONPATH", "")},
    )

    def call(request_id, method, params):
        assert process.stdin is not None and process.stdout is not None
        process.stdin.write(json.dumps({"jsonrpc": "2.0", "id": request_id, "method": method, "params": params}) + "\n")
        process.stdin.flush()
        ready, _, _ = select.select([process.stdout], [], [], 5.0)
        assert ready, f"MCP server did not answer {method}"
        return json.loads(process.stdout.readline())

    try:
        pre_init = call(1, "tools/list", {})
        assert pre_init["error"]["data"]["code"] == "NOT_INITIALIZED"
        initialized = call(2, "initialize", {"protocolVersion": "2025-06-18", "capabilities": {}, "clientInfo": {"name": "test", "version": "1"}})
        assert initialized["result"]["protocolVersion"] == "2025-06-18"
        assert initialized["result"]["capabilities"]["tools"] == {}
        tools = call(3, "tools/list", {})["result"]["tools"]
        assert any(tool["name"] == "submit_scm_run" for tool in tools)
        assert any(tool["name"] == "list_tunable_parameters" for tool in tools)
        assert any(tool["name"] == "resolve_clubb_settings" for tool in tools)
        # A temporarily unavailable browser/broker is an ordinary structured
        # tool result, not a ToolError that poisons the persistent MCP
        # transport.  Whether a developer has Dash open or not, ping must
        # still work immediately after this call.
        connected = call(
            30,
            "tools/call",
            {"name": "connect_to_dashboard", "arguments": {}},
        )
        assert "error" not in connected
        assert connected["result"]["content"]
        assert call(301, "ping", {})["result"] == {}
        resolved = call(
            31,
            "tools/call",
            {
                "name": "resolve_clubb_settings",
                "arguments": {
                    "flags": {"iiPDF_type": 1, "l_predict_upwp_vpwp": True},
                    "auto_correct": True,
                },
            },
        )
        resolution = json.loads(resolved["result"]["content"][0]["text"])["resolution"]
        assert resolution["forced_flags"] == {}
        assert resolution["parameter_states"]["C_uu_shr"]["state"] == "active"
        bad_job = call(4, "tools/call", {"name": "get_job", "arguments": {"job_id": "missing"}})
        assert "INVALID_REQUEST" in bad_job["result"]["content"][0]["text"]
        bad_schema = call(
            5,
            "tools/call",
            {"name": "submit_scm_run", "arguments": {"request": {"request_id": "schema-audit-123", "case": "arm", "unexpected": True}}},
        )
        assert bad_schema["error"]["code"] == -32602
        assert bad_schema["error"]["data"]["code"] == "INVALID_REQUEST"
        templates = call(6, "resources/templates/list", {})["result"]["resourceTemplates"]
        assert any(item["uriTemplate"].startswith("clubb-artifact://") for item in templates)
    finally:
        process.terminate()
        process.wait(timeout=5)


def test_mcp_stdio_keeps_tool_prints_out_of_protocol_stdout(monkeypatch):
    """Legacy launch diagnostics must not corrupt JSON-RPC framing."""
    from mcp.server.fastmcp import FastMCP
    from dash_app.agent_integration import mcp_server

    server = FastMCP("noisy-test")

    @server.tool()
    def noisy_tool() -> dict[str, bool]:
        print("ordinary launcher diagnostic")
        return {"ok": True}

    requests = [
        {
            "jsonrpc": "2.0",
            "id": 1,
            "method": "initialize",
            "params": {
                "protocolVersion": "2025-06-18",
                "capabilities": {},
                "clientInfo": {"name": "test", "version": "1"},
            },
        },
        {
            "jsonrpc": "2.0",
            "id": 2,
            "method": "tools/call",
            "params": {"name": "noisy_tool", "arguments": {}},
        },
    ]
    protocol_output = io.StringIO()
    diagnostics = io.StringIO()
    monkeypatch.setattr(
        mcp_server.sys,
        "stdin",
        io.StringIO("".join(json.dumps(item) + "\n" for item in requests)),
    )
    monkeypatch.setattr(mcp_server.sys, "stdout", protocol_output)
    monkeypatch.setattr(mcp_server.sys, "stderr", diagnostics)

    mcp_server._serve_stdio(server)

    responses = [
        json.loads(line)
        for line in protocol_output.getvalue().splitlines()
    ]
    assert [item["id"] for item in responses] == [1, 2]
    assert responses[1]["result"]["structuredContent"] == {"ok": True}
    assert "ordinary launcher diagnostic" in diagnostics.getvalue()


def test_scm_recovery_restarts_monitoring_from_durable_broker_state(monkeypatch):
    from dash_app.shared import actions

    launched = []
    monkeypatch.setattr(
        actions,
        "broker_jobs",
        lambda: {"runs": {"arm": {"state": "running", "case": "arm", "job_id": "scm_abc", "proc_data": {"pid": 1234, "log": "/tmp/arm.log"}}}},
    )
    monkeypatch.setattr(actions, "_pid_is_alive", lambda pid: pid == 1234)
    monkeypatch.setattr(actions, "_background", lambda target, *args: launched.append((target, args)))
    monkeypatch.setattr(actions, "publish_event", lambda *args, **kwargs: None)

    recovered = actions.recover_active_runs_from_state()

    assert recovered == [{"case": "arm", "pid": 1234, "job_id": "scm_abc"}]
    assert launched[0][0] is actions._watch_run


def test_scm_recovery_seals_a_dead_record_instead_of_leaving_it_running(monkeypatch, tmp_path):
    from dash_app.shared import actions

    store = JobStore(tmp_path / "jobs.json")
    record, _created = store.submit(
        "scm",
        "dead-recovery-123",
        {"case": "arm"},
    )
    artifacts = ArtifactStore(tmp_path / "agent_artifacts")
    artifacts.create_manifest(record["job_id"], {"job_id": record["job_id"]})
    artifacts.activate(record["job_id"])
    updates = []
    monkeypatch.setattr(actions, "_JOB_STORE", store)
    monkeypatch.setattr(actions, "_ARTIFACT_STORE", artifacts)
    monkeypatch.setattr(
        actions,
        "broker_jobs",
        lambda: {
            "runs": {
                "arm": {
                    "state": "running",
                    "case": "arm",
                    "job_id": record["job_id"],
                    "proc_data": {"pid": 1234},
                }
            }
        },
    )
    monkeypatch.setattr(actions, "_pid_is_alive", lambda _pid: False)
    monkeypatch.setattr(
        actions,
        "update_broker_run",
        lambda case, **values: updates.append((case, values)),
    )

    assert actions.recover_active_runs_from_state() == []
    assert updates[0][1]["state"] == "error"
    assert store.get(record["job_id"])["error"]["code"] == "SCM_RECOVERY_EXIT_UNKNOWN"
    assert not (artifacts.bundle_dir(record["job_id"]) / ".active").exists()


def test_scm_watcher_keeps_log_progress_out_of_global_activity(monkeypatch, tmp_path):
    """Raw SCM output remains available to Run, not global activity history."""
    from dash_app.shared import actions

    events = []

    class FinishedProcess:
        def poll(self):
            return 0

    monkeypatch.setattr(actions, "read_run_log_increment", lambda _path, _offset: ("iteration: 2 / 3\n", 16))
    monkeypatch.setattr(actions, "get_proc", lambda _pid: FinishedProcess())
    monkeypatch.setattr(actions, "record_case_finish", lambda *_args: None)
    monkeypatch.setattr(actions, "cleanup_temp_files", lambda *_args: None)
    monkeypatch.setattr(actions, "update_broker_run", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(actions, "publish_event", lambda *_args, **_kwargs: events.append(_args))

    actions._watch_run("arm", {"pid": 123, "log": str(tmp_path / "arm.log"), "temp_files": []})

    assert [event[1] for event in events] == ["arm finished"]


def test_scm_watcher_records_an_explicit_stop_as_cancelled(monkeypatch, tmp_path):
    from dash_app.shared import actions

    store = JobStore(tmp_path / "jobs.json")
    record, _created = store.submit("scm", "cancelled-scm-123", {"case": "arm"})

    class StoppedProcess:
        def poll(self):
            return -15

    broker_updates = []
    monkeypatch.setattr(actions, "_JOB_STORE", store)
    monkeypatch.setattr(actions, "_ARTIFACT_STORE", ArtifactStore(tmp_path / "artifacts"))
    monkeypatch.setattr(actions, "broker_jobs", lambda: {"runs": {"arm": {"state": "stopping"}}})
    monkeypatch.setattr(actions, "read_run_log_increment", lambda _path, _offset: ("", 0))
    monkeypatch.setattr(actions, "get_proc", lambda _pid: StoppedProcess())
    monkeypatch.setattr(actions, "record_case_finish", lambda *_args: None)
    monkeypatch.setattr(actions, "cleanup_temp_files", lambda *_args: None)
    monkeypatch.setattr(actions, "_write_run_result_summary", lambda *_args: None)
    monkeypatch.setattr(actions, "_read_log_tail", lambda *_args: "")
    monkeypatch.setattr(actions, "update_broker_run", lambda case, **values: broker_updates.append((case, values)))
    monkeypatch.setattr(actions, "publish_event", lambda *_args, **_kwargs: None)

    actions._watch_run(
        "arm",
        {"pid": 123, "log": str(tmp_path / "arm.log"), "temp_files": []},
        record["job_id"],
    )

    assert broker_updates[-1][1]["state"] == "stopped"
    assert store.get(record["job_id"])["state"] == "cancelled"


def test_compile_watcher_records_an_explicit_stop_as_cancelled(monkeypatch, tmp_path):
    from dash_app.shared import actions

    store = JobStore(tmp_path / "jobs.json")
    record, _created = store.submit("compile", "cancelled-compile-123", {"debug": True})
    broker_updates = []
    monkeypatch.setattr(actions, "_JOB_STORE", store)
    monkeypatch.setattr(actions, "_ARTIFACT_STORE", ArtifactStore(tmp_path / "artifacts"))
    monkeypatch.setattr(actions, "broker_jobs", lambda: {"compile": {"state": "stopping"}})
    monkeypatch.setattr(actions, "read_compile_log_increment", lambda _path, _offset: ("", 0))
    monkeypatch.setattr(actions, "poll_compile_job", lambda _job: -15)
    monkeypatch.setattr(actions, "finish_compile_job", lambda job, _returncode: job)
    monkeypatch.setattr(actions, "_read_log_tail", lambda *_args: "")
    monkeypatch.setattr(actions, "update_broker_job", lambda kind, **values: broker_updates.append((kind, values)))
    monkeypatch.setattr(actions, "publish_event", lambda *_args, **_kwargs: None)

    actions._watch_compile({"log": str(tmp_path / "compile.log"), "command": "compile"}, record["job_id"])

    assert broker_updates[-1][1]["state"] == "stopped"
    assert store.get(record["job_id"])["state"] == "cancelled"


def test_recovered_tune_blocks_duplicate_typed_submission(monkeypatch, tmp_path):
    from dash_app.shared import actions

    store = JobStore(tmp_path / "jobs.json")
    monkeypatch.setattr(actions, "_JOB_STORE", store)
    monkeypatch.setattr(actions, "recover_active_tuning_from_disk", lambda: {"state": "running"})
    monkeypatch.setattr(actions, "broker_jobs", lambda: {"tune": {"state": "running"}})
    request = TuneRequest(
        request_id="request-tune-123",
        cases=["arm"],
        parameter_ranges=[{"name": "c_K10", "min": 0.0, "max": 1.0}],
    )
    with pytest.raises(ValueError, match="already active"):
        actions.submit_tune(request)
