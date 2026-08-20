import base64
import threading
from pathlib import Path

from dash import Dash, html
import pytest

from dash_app import app as dashboard_app
from dash_app.services.jobs import ArtifactStore, JobStore
from dash_app.services.models import ScmRunBatchRequest
from dash_app.shared import actions, activity


def _isolated_activity(tmp_path, monkeypatch):
    monkeypatch.setattr(activity, "ACTIVITY_PATH", tmp_path / "activity.json")
    monkeypatch.setattr(activity, "LOCK_PATH", tmp_path / "activity.lock")
    activity.reset_activity()


def _isolated_batch_services(tmp_path, monkeypatch):
    store = JobStore(tmp_path / "jobs.json")
    monkeypatch.setattr(actions, "_JOB_STORE", store)
    monkeypatch.setattr(actions, "_ARTIFACT_STORE", ArtifactStore(tmp_path / "artifacts"))
    monkeypatch.setattr(actions, "_public_scm_output_root", lambda: tmp_path / "output")
    monkeypatch.setattr(actions, "_ensure_scm_batch_watcher", lambda _job_id: None)
    monkeypatch.setattr(
        dashboard_app,
        "run_telemetry",
        lambda revision, cursors: actions.run_telemetry(revision, cursors),
    )
    return store


def _set_runtime(store, child, *, pid, log):
    store.update(
        child["job_id"],
        state="running",
        runtime={"proc_data": {"pid": pid, "log": str(log), "start_time": 10}},
    )


def _native_output_resolver(root):
    def resolve(value):
        requested = Path(str(value or "output"))
        if requested.is_absolute():
            return requested
        parts = requested.parts[1:] if requested.parts[:1] == ("output",) else requested.parts
        return root / "output" / Path(*parts)

    return resolve


def test_batch_admission_is_scoped_to_canonical_case_and_output(
    tmp_path, monkeypatch
):
    store = _isolated_batch_services(tmp_path, monkeypatch)
    monkeypatch.setattr(
        actions, "resolve_output_dir", _native_output_resolver(tmp_path)
    )
    first = actions.submit_scm_batch(
        ScmRunBatchRequest(request_id="output-scope-one", cases=["arm"]),
        native_overrides={},
        native_cli_options={"out_dir": "default"},
    )
    alias = actions.submit_scm_batch(
        ScmRunBatchRequest(request_id="output-scope-two", cases=["arm"]),
        native_overrides={},
        native_cli_options={"out_dir": "output/default"},
    )
    other = actions.submit_scm_batch(
        ScmRunBatchRequest(request_id="output-scope-three", cases=["arm"]),
        native_overrides={},
        native_cli_options={"out_dir": "new"},
    )
    blank = actions.submit_scm_batch(
        ScmRunBatchRequest(request_id="output-scope-blank", cases=["rico"]),
        native_overrides={},
        native_cli_options={},
    )

    assert first["accepted_cases"] == ["arm"]
    assert alias["status"] == "no_op"
    assert alias["accepted_cases"] == []
    assert alias["skipped_cases"] == [
        {
            "case": "arm",
            "code": "CASE_OUTPUT_ACTIVE",
            "output_directory": str((tmp_path / "output" / "default").resolve()),
            "conflicting_job_id": first["children"][0]["job_id"],
            "conflicting_run_id": first["children"][0]["run_id"],
        }
    ]
    assert alias["state"] == "finished"
    assert alias["children"] == []
    assert other["accepted_cases"] == ["arm"]
    assert blank["output_directory"] == str((tmp_path / "output").resolve())
    assert len(store.list_kind("scm")) == 3


def test_persisted_queued_child_reserves_case_output_after_store_restart(
    tmp_path, monkeypatch
):
    store = _isolated_batch_services(tmp_path, monkeypatch)
    monkeypatch.setattr(
        actions, "resolve_output_dir", _native_output_resolver(tmp_path)
    )
    first = actions.submit_scm_batch(
        ScmRunBatchRequest(request_id="restart-reservation-one", cases=["arm"]),
        native_overrides={},
        native_cli_options={"out_dir": "default"},
    )

    restarted_store = JobStore(store.path)
    monkeypatch.setattr(actions, "_JOB_STORE", restarted_store)
    duplicate = actions.submit_scm_batch(
        ScmRunBatchRequest(request_id="restart-reservation-two", cases=["arm"]),
        native_overrides={},
        native_cli_options={"out_dir": "output/default"},
    )

    assert duplicate["status"] == "no_op"
    assert duplicate["skipped_cases"][0]["conflicting_job_id"] == (
        first["children"][0]["job_id"]
    )


def test_batch_admission_launches_only_cases_not_reserved_in_target(
    tmp_path, monkeypatch
):
    store = _isolated_batch_services(tmp_path, monkeypatch)
    monkeypatch.setattr(
        actions, "resolve_output_dir", _native_output_resolver(tmp_path)
    )
    first = actions.submit_scm_batch(
        ScmRunBatchRequest(
            request_id="partial-scope-one", cases=["arm", "bomex"]
        ),
        native_overrides={},
        native_cli_options={"out_dir": "default"},
    )
    arm, bomex = first["children"]
    store.update(
        bomex["job_id"],
        state="finished",
        returncode=0,
        finished_at_unix_seconds=20,
    )

    second_request = ScmRunBatchRequest(
        request_id="partial-scope-two", cases=["arm", "bomex"]
    )
    second = actions.submit_scm_batch(
        second_request,
        native_overrides={},
        native_cli_options={"out_dir": "default"},
    )
    retry = actions.submit_scm_batch(
        second_request,
        native_overrides={},
        native_cli_options={"out_dir": "default"},
    )

    assert second["status"] == "queued"
    assert second["accepted_cases"] == ["bomex"]
    assert [item["case"] for item in second["skipped_cases"]] == ["arm"]
    assert [child["case"] for child in second["children"]] == ["bomex"]
    assert retry["status"] == "existing"
    assert retry["accepted_cases"] == ["bomex"]
    assert retry["skipped_cases"] == second["skipped_cases"]
    assert store.get(arm["job_id"])["state"] == "queued"


def test_child_setup_failure_does_not_strand_later_queued_cases(
    tmp_path, monkeypatch
):
    store = _isolated_batch_services(tmp_path, monkeypatch)
    batch = actions.submit_scm_batch(
        ScmRunBatchRequest(
            request_id="continue-after-setup-failure",
            cases=["arm", "atex", "bomex"],
            max_workers=3,
        )
    )
    attempted = []

    def start_child(request, record, *_args, **_kwargs):
        attempted.append(request.case)
        if request.case == "arm":
            raise FileNotFoundError("retention cleanup raced another launcher")
        store.update(
            record["job_id"],
            state="running",
            runtime={"proc_data": {"pid": 9000 + len(attempted)}},
        )

    monkeypatch.setattr(actions, "_start_scm_submission", start_child)

    result = actions._start_queued_scm_batch_children(batch["job_id"])
    children = {
        (record.get("request") or {}).get("case"): record
        for record in store.list_kind("scm")
        if record.get("batch_job_id") == batch["job_id"]
    }

    assert attempted == ["arm", "atex", "bomex"]
    assert children["arm"]["state"] == "error"
    assert children["arm"]["error"]["code"] == "SCM_BATCH_CHILD_FAILED"
    assert children["atex"]["state"] == "running"
    assert children["bomex"]["state"] == "running"
    assert result["state"] == "running"


def test_run_completion_wakes_the_registered_watcher_without_draining_inline(
    tmp_path, monkeypatch
):
    store = _isolated_batch_services(tmp_path, monkeypatch)
    parent, _ = store.submit(
        "scm_batch",
        "completion-refill-parent",
        {"request_id": "completion-refill-parent", "cases": ["arm", "bomex"]},
    )
    arm, _ = store.submit("scm", "completion-refill-arm", {"case": "arm"})
    bomex, _ = store.submit(
        "scm", "completion-refill-bomex", {"case": "bomex"}
    )
    store.update(
        parent["job_id"],
        state="running",
        child_job_ids=[arm["job_id"], bomex["job_id"]],
    )
    store.update(
        arm["job_id"], state="running", batch_job_id=parent["job_id"]
    )
    store.update(
        bomex["job_id"], state="queued", batch_job_id=parent["job_id"]
    )
    awakened = []

    class FinishedProcess:
        def poll(self):
            return 0

    monkeypatch.setattr(actions, "get_proc", lambda _pid: FinishedProcess())
    monkeypatch.setattr(actions, "record_case_finish", lambda *_args: None)
    monkeypatch.setattr(actions, "cleanup_temp_files", lambda *_args: None)
    monkeypatch.setattr(actions, "publish_event", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(
        actions, "_ensure_scm_batch_watcher", lambda job_id: awakened.append(job_id)
    )
    monkeypatch.setattr(
        actions,
        "_start_queued_scm_batch_children",
        lambda _job_id: pytest.fail("completion thread drained the queue inline"),
    )

    actions._watch_run(
        "arm",
        {"pid": 123, "log": str(tmp_path / "arm.log"), "temp_files": []},
        arm["job_id"],
    )

    assert awakened == [parent["job_id"]]
    assert store.get(arm["job_id"])["state"] == "finished"
    assert store.get(bomex["job_id"])["state"] == "queued"


def test_default_new_scenario_keeps_four_reservations_and_restarts_finished_case(
    tmp_path, monkeypatch
):
    store = _isolated_batch_services(tmp_path, monkeypatch)
    monkeypatch.setattr(
        actions, "resolve_output_dir", _native_output_resolver(tmp_path)
    )
    default = actions.submit_scm_batch(
        ScmRunBatchRequest(
            request_id="four-run-default", cases=["arm", "bomex"]
        ),
        native_overrides={},
        native_cli_options={"out_dir": "default"},
    )
    blocked = actions.submit_scm_batch(
        ScmRunBatchRequest(
            request_id="four-run-default-blocked", cases=["arm", "bomex"]
        ),
        native_overrides={},
        native_cli_options={"out_dir": "output/default"},
    )
    new = actions.submit_scm_batch(
        ScmRunBatchRequest(request_id="four-run-new", cases=["arm", "bomex"]),
        native_overrides={},
        native_cli_options={"out_dir": "new"},
    )

    assert blocked["status"] == "no_op"
    assert default["accepted_cases"] == ["arm", "bomex"]
    assert new["accepted_cases"] == ["arm", "bomex"]
    assert len(
        [
            record
            for record in store.list_kind("scm")
            if record["state"] == "queued"
        ]
    ) == 4

    default_bomex = next(
        child for child in default["children"] if child["case"] == "bomex"
    )
    store.update(default_bomex["job_id"], state="finished", returncode=0)
    resumed = actions.submit_scm_batch(
        ScmRunBatchRequest(
            request_id="four-run-default-resume", cases=["arm", "bomex"]
        ),
        native_overrides={},
        native_cli_options={"out_dir": "default"},
    )

    assert resumed["accepted_cases"] == ["bomex"]
    assert [item["case"] for item in resumed["skipped_cases"]] == ["arm"]


def test_batch_cancel_targets_one_of_two_same_case_processes_by_job(
    tmp_path, monkeypatch
):
    store = _isolated_batch_services(tmp_path, monkeypatch)
    monkeypatch.setattr(
        actions, "resolve_output_dir", _native_output_resolver(tmp_path)
    )
    first = actions.submit_scm_batch(
        ScmRunBatchRequest(request_id="exact-cancel-one", cases=["arm"]),
        native_overrides={},
        native_cli_options={"out_dir": "default"},
    )
    second = actions.submit_scm_batch(
        ScmRunBatchRequest(request_id="exact-cancel-two", cases=["arm"]),
        native_overrides={},
        native_cli_options={"out_dir": "new"},
    )
    _set_runtime(store, first["children"][0], pid=4101, log="/tmp/arm-default.log")
    _set_runtime(store, second["children"][0], pid=4102, log="/tmp/arm-new.log")
    signalled = []
    monkeypatch.setattr(
        actions,
        "_signal_scm_record",
        lambda record: signalled.append(
            ((record.get("runtime") or {}).get("proc_data") or {}).get("pid")
        )
        or {"status": "stop_requested"},
    )

    with pytest.raises(ValueError, match="multiple arm runs are active"):
        actions.stop_run("arm")

    actions.cancel_scm_batch(first["job_id"])

    assert signalled == [4101]
    assert store.get(first["children"][0]["job_id"])["state"] == "stopping"
    assert store.get(second["children"][0]["job_id"])["state"] == "running"

    actions.cancel_job(second["children"][0]["job_id"])
    assert signalled == [4101, 4102]
    assert store.get(second["children"][0]["job_id"])["state"] == "stopping"


def test_batch_create_refresh_and_cancel_each_write_registry_once(tmp_path, monkeypatch):
    store = _isolated_batch_services(tmp_path, monkeypatch)
    writes = []
    original_write = store._write

    def counted_write(payload):
        writes.append(1)
        original_write(payload)

    monkeypatch.setattr(store, "_write", counted_write)
    batch = actions.submit_scm_batch(
        ScmRunBatchRequest(
            request_id="single-write-batch",
            cases=["arm", "bomex", "rico"],
            max_workers=2,
        )
    )
    assert len(writes) == 1
    stored_batch = store.get(batch["job_id"])
    assert "children" not in stored_batch
    assert stored_batch["child_job_ids"] == [
        child["job_id"] for child in batch["children"]
    ]

    first_child = batch["children"][0]["job_id"]
    store.update(first_child, state="running")
    writes.clear()
    refreshed = actions._refresh_scm_batch(batch["job_id"], first_child)
    assert refreshed["state"] == "running"
    assert len(writes) == 1

    monkeypatch.setattr(actions, "stop_run", lambda _case: {"status": "stop_requested"})
    store.update(first_child, runtime={"proc_data": {"pid": 123}})
    writes.clear()
    cancelled = actions.cancel_scm_batch(batch["job_id"])
    assert cancelled["state"] == "stopping"
    assert len(writes) == 1


def test_scm_records_never_embed_console_text(tmp_path, monkeypatch):
    _isolated_activity(tmp_path, monkeypatch)
    store = _isolated_batch_services(tmp_path, monkeypatch)
    batch = actions.submit_scm_batch(
        ScmRunBatchRequest(
            request_id="compact-batch-records",
            cases=["arm", "bomex"],
        )
    )
    assert "log_tail" not in str(store.snapshot())
    assert "log_offset" not in str(store.snapshot())


def test_telemetry_reports_all_user_runs_once_per_lifecycle_revision(
    tmp_path, monkeypatch
):
    store = _isolated_batch_services(tmp_path, monkeypatch)
    dash_run, _ = store.submit("scm", "dash-run", {"case": "arm"})
    mcp_run, _ = store.submit("scm", "mcp-run", {"case": "bomex"})
    internal_run, _ = store.submit("scm", "internal-run", {"case": "rico"})
    store.update(
        dash_run["job_id"],
        state="running",
        origin="dash",
        visibility="user",
        output_directory=str(tmp_path / "output" / "default"),
    )
    store.update(mcp_run["job_id"], state="queued", origin="mcp", visibility="user")
    store.update(
        internal_run["job_id"], state="running", origin="profile", visibility="internal"
    )

    first = actions.run_telemetry()
    assert {(run["case"], run["origin"]) for run in first["runs"]} == {
        ("arm", "dash"),
        ("bomex", "mcp"),
    }
    assert first["active"] is True
    arm = next(run for run in first["runs"] if run["case"] == "arm")
    assert arm["output_directory"] == str(
        (tmp_path / "output" / "default").resolve()
    )

    unchanged = actions.run_telemetry(first["revision"])
    assert "runs" not in unchanged


def test_direct_log_endpoint_opens_at_bounded_tail_then_streams_new_text(tmp_path, monkeypatch):
    store = _isolated_batch_services(tmp_path, monkeypatch)
    lines = [f"line {index:04d} " + "x" * 80 + "\n" for index in range(6000)]
    log_path = tmp_path / "wrapper.log"
    log_path.write_text("".join(lines), encoding="utf-8")
    record, _ = store.submit("scm", "complete-log", {"case": "arm"})
    store.update(
        record["job_id"],
        state="finished",
        runtime={"proc_data": {"log": str(log_path), "start_time": 10}},
        visibility="user",
    )
    run_id = record["run_id"]

    app = Dash(__name__)
    app.layout = html.Div()
    dashboard_app._register_run_telemetry_route(app)
    client = app.server.test_client()
    response = client.post(
        "/_clubb-run-telemetry",
        json={"log_cursors": {run_id: None}},
    )
    payload = response.get_json()["logs"][run_id]
    assert base64.b64decode(payload["data"]).decode().splitlines() == [
        line.rstrip("\n") for line in lines[-5000:]
    ]
    cursor = payload["next_cursor"]
    assert cursor == log_path.stat().st_size

    with log_path.open("a", encoding="utf-8") as handle:
        handle.write("new output\n")
    response = client.post(
        "/_clubb-run-telemetry",
        json={"known_revision": response.get_json()["revision"], "log_cursors": {run_id: cursor}},
    )
    payload = response.get_json()["logs"][run_id]
    assert base64.b64decode(payload["data"]) == b"new output\n"
    assert client.post(
        "/_clubb-run-telemetry",
        json={"log_cursors": {str(tmp_path / "other.log"): 0}},
    ).get_json()["logs"] == {}


def test_direct_log_endpoint_returns_bounded_progress_without_opening_console(
    tmp_path, monkeypatch
):
    store = _isolated_batch_services(tmp_path, monkeypatch)
    log_path = tmp_path / "wrapper.log"
    log_path.write_text(
        "iteration: 7 / 20 -- time = 420s / 1200s\n", encoding="utf-8"
    )
    record, _ = store.submit("scm", "progress-log", {"case": "arm"})
    _set_runtime(store, record, pid=123, log=log_path)

    app = Dash(__name__)
    app.layout = html.Div()
    dashboard_app._register_run_telemetry_route(app)
    response = app.server.test_client().post(
        "/_clubb-run-telemetry",
        json={},
    )

    assert response.get_json()["progress"] == {
        record["run_id"]: {"iteration": 7, "total": 20}
    }


def test_clear_deletes_only_terminal_wrapper_logs(tmp_path, monkeypatch):
    _isolated_activity(tmp_path, monkeypatch)
    store = JobStore(tmp_path / "jobs.json")
    monkeypatch.setattr(actions, "_JOB_STORE", store)
    terminal_log = tmp_path / "terminal.log"
    active_log = tmp_path / "active.log"
    terminal_log.write_text("done")
    active_log.write_text("working")
    terminal, _ = store.submit("scm", "terminal-job", {"case": "arm"})
    active, _ = store.submit("scm", "active-job", {"case": "bomex"})
    store.update(terminal["job_id"], state="finished", runtime={"proc_data": {"log": str(terminal_log)}})
    store.update(active["job_id"], state="running", runtime={"proc_data": {"log": str(active_log), "pid": 42}})

    actions.clear_terminal_scm_session()

    assert not terminal_log.exists()
    assert active_log.exists()
    assert store.get(terminal["job_id"]) is None
    assert store.get(active["job_id"])["state"] == "running"


def test_stop_run_terminates_the_process_group_then_escalates(tmp_path, monkeypatch):
    signals = []
    store = _isolated_batch_services(tmp_path, monkeypatch)
    record, _ = store.submit("scm", "stop-arm", {"case": "arm"})
    _set_runtime(store, record, pid=4242, log="/tmp/arm.log")
    monkeypatch.setattr(actions.os, "killpg", lambda pid, sig: signals.append((pid, sig)))
    monkeypatch.setattr(actions, "_pid_is_alive", lambda _pid: True)
    monkeypatch.setattr(actions.time, "sleep", lambda _seconds: None)
    monkeypatch.setattr(actions, "_background", lambda target, *args: target(*args))
    monkeypatch.setattr(actions, "publish_event", lambda *_args, **_kwargs: None)

    actions.stop_run("arm")

    assert signals == [(4242, actions.signal.SIGTERM), (4242, actions.signal.SIGKILL)]


def test_cancel_all_scm_signals_every_process_and_cancels_the_queue(
    tmp_path, monkeypatch
):
    _isolated_activity(tmp_path, monkeypatch)
    store = _isolated_batch_services(tmp_path, monkeypatch)
    batch = actions.submit_scm_batch(
        ScmRunBatchRequest(
            request_id="cancel-all-runs",
            cases=["arm", "bomex", "rico"],
            max_workers=2,
        )
    )
    arm, bomex, rico = batch["children"]
    _set_runtime(store, arm, pid=4101, log="/tmp/arm.log")
    _set_runtime(store, bomex, pid=4102, log="/tmp/bomex.log")

    signals = []
    writes = []
    original_write = store._write

    def counted_write(payload):
        writes.append(1)
        original_write(payload)

    monkeypatch.setattr(store, "_write", counted_write)
    monkeypatch.setattr(
        actions.os, "killpg", lambda pid, sig: signals.append((pid, sig))
    )
    monkeypatch.setattr(actions, "_background", lambda *_args: None)
    monkeypatch.setattr(
        actions,
        "stop_run",
        lambda _case: (_ for _ in ()).throw(AssertionError("per-case stop used")),
    )

    result = actions.cancel_all_scm_runs()

    assert writes == [1]
    assert signals == [
        (4101, actions.signal.SIGTERM),
        (4102, actions.signal.SIGTERM),
    ]
    assert store.get(arm["job_id"])["state"] == "stopping"
    assert store.get(bomex["job_id"])["state"] == "stopping"
    assert store.get(rico["job_id"])["state"] == "cancelled"
    assert store.get(batch["job_id"])["state"] == "stopping"
    assert result["stopped_cases"] == ["arm", "bomex"]
    assert result["summary"]["stopping"] == 2


def test_cancel_all_scm_cancels_a_fully_queued_batch(tmp_path, monkeypatch):
    _isolated_activity(tmp_path, monkeypatch)
    store = _isolated_batch_services(tmp_path, monkeypatch)
    batch = actions.submit_scm_batch(
        ScmRunBatchRequest(
            request_id="cancel-all-queued", cases=["arm", "bomex"]
        )
    )
    monkeypatch.setattr(actions, "_background", lambda *_args: None)

    result = actions.cancel_all_scm_runs()

    assert result["status"] == "cancelled"
    assert result["stopped_cases"] == []
    assert store.get(batch["job_id"])["state"] == "cancelled"
    assert all(
        store.get(child["job_id"])["state"] == "cancelled"
        for child in batch["children"]
    )


def test_cancel_all_does_not_wait_for_child_setup_or_allow_a_late_launch(
    tmp_path, monkeypatch
):
    _isolated_activity(tmp_path, monkeypatch)
    store = _isolated_batch_services(tmp_path, monkeypatch)
    batch = actions.submit_scm_batch(
        ScmRunBatchRequest(request_id="cancel-during-setup", cases=["arm"])
    )
    setup_started = threading.Event()
    release_setup = threading.Event()
    cancel_finished = threading.Event()
    launched = threading.Event()
    original_create_manifest = actions._ARTIFACT_STORE.create_manifest

    def slow_create_manifest(*args, **kwargs):
        setup_started.set()
        assert release_setup.wait(timeout=2)
        return original_create_manifest(*args, **kwargs)

    def unexpected_launch(*_args, **_kwargs):
        launched.set()
        raise AssertionError("cancelled child was launched")

    monkeypatch.setattr(
        actions._ARTIFACT_STORE, "create_manifest", slow_create_manifest
    )
    monkeypatch.setattr(actions, "_launch_scm_process", unexpected_launch)
    monkeypatch.setattr(actions, "_background", lambda *_args: None)

    starter = threading.Thread(
        target=actions._start_one_queued_scm_batch_child,
        args=(batch["job_id"],),
    )
    starter.start()
    assert setup_started.wait(timeout=1)

    def cancel():
        actions.cancel_all_scm_runs()
        cancel_finished.set()

    canceller = threading.Thread(target=cancel)
    canceller.start()
    try:
        assert cancel_finished.wait(timeout=1)
    finally:
        release_setup.set()
        starter.join(timeout=2)
        canceller.join(timeout=2)

    child = batch["children"][0]
    assert not starter.is_alive()
    assert not canceller.is_alive()
    assert launched.is_set() is False
    assert store.get(child["job_id"])["state"] == "cancelled"
    assert store.get(batch["job_id"])["state"] == "cancelled"


def test_dispatch_exposes_broker_owned_cancel_all(monkeypatch):
    monkeypatch.setattr(
        actions, "cancel_all_scm_runs", lambda: {"status": "stop_requested"}
    )

    assert actions.dispatch("domain_cancel_all_scm", {}) == {
        "status": "stop_requested"
    }
