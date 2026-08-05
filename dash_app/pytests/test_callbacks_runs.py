from dash_app.run_tab.callbacks_runs import (
    build_multicol_spec,
    discard_terminal_broker_runs,
    expand_linked_parameter_values,
    launch_broker_batch,
    fresh_batch_request_id,
)


def test_multicol_spec_preserves_current_parameter_names():
    spec = build_multicol_spec(
        ["C_uu_shr", "C_uu_buoy"],
        ["0", "0.1"],
        ["3", "1.2"],
        ["6", "4"],
        ["C_uu_shr", "C_uu_buoy"],
    )

    assert spec == "C_uu_shr/0:3/6,C_uu_buoy/0.1:1.2/4"


def test_fresh_batch_request_ids_do_not_reuse_browser_click_state():
    first = fresh_batch_request_id("same selected cases and settings")
    second = fresh_batch_request_id("same selected cases and settings")

    assert first != second
    assert first.startswith("dash-run-batch-")
    assert first.rsplit("-", 1)[1] == second.rsplit("-", 1)[1]


def test_multicol_spec_rejects_unknown_complete_parameter_rows():
    try:
        build_multicol_spec(["not_a_parameter"], ["0"], ["1"], ["4"], ["C1"])
    except ValueError as exc:
        assert "not available" in str(exc)
    else:
        raise AssertionError("unknown multicol parameter was accepted")


def test_linked_control_expands_to_every_physical_parameter_member():
    ids = [
        {"type": "run-param", "file": "tunable", "name": "C6thl"},
        {"type": "run-param", "file": "tunable", "name": "C1"},
        {"type": "run-param", "file": "tunable", "name": "C6rt"},
    ]

    expanded = expand_linked_parameter_values(
        ids,
        ["old-follower", "unchanged", "old-master"],
        [{"type": "run-linked-param", "group": "C6rt=C6thl"}],
        ["2.5"],
    )

    assert expanded == ["2.5", "unchanged", "2.5"]


def test_multicol_linked_member_becomes_one_equal_value_coordinate():
    spec = build_multicol_spec(
        ["C6rtb"],
        ["0"],
        ["4"],
        ["5"],
        ["C6rtb", "C6thlb"],
        [("C6rtb", "C6thlb")],
    )

    assert spec == "C6rtb=C6thlb/0:4/5"


def test_multicol_rejects_two_rows_for_one_linked_coordinate():
    try:
        build_multicol_spec(
            ["C6rtb", "C6thlb"],
            ["0", "0"],
            ["4", "4"],
            ["5", "5"],
            ["C6rtb", "C6thlb"],
            [("C6rtb", "C6thlb")],
        )
    except ValueError as exc:
        assert "already selected" in str(exc)
    else:
        raise AssertionError("duplicate linked hypergrid coordinate was accepted")


def test_clear_discards_only_terminal_broker_run_records():
    removed = []

    discarded = discard_terminal_broker_runs(
        {
            "arm": {"case": "arm", "state": "finished"},
            "bomex": {"case": "bomex", "state": "failed"},
            "dycoms2_rf01": {"case": "dycoms2_rf01", "state": "running"},
            "atex": {"case": "atex", "state": "stopping"},
        },
        lambda case, payload: removed.append((case, payload)),
    )

    assert discarded == ["arm", "bomex"]
    assert removed == [("arm", None), ("bomex", None)]


def test_native_multi_case_selection_submits_one_shared_batch(monkeypatch):
    from dash_app.shared import broker_client

    captured = {}

    def fake_action(action, payload, **kwargs):
        captured.update(action=action, payload=payload, kwargs=kwargs)
        return {
            "status": "started",
            "batch_id": "batch_native_test",
            "children": [
                {
                    "case": "arm",
                    "state": "running",
                    "runtime": {"proc_data": {"pid": 101, "log": "/tmp/arm.log"}},
                },
                {"case": "bomex", "state": "queued"},
            ],
        }

    monkeypatch.setattr(broker_client, "perform_action", fake_action)
    running = {}
    logs = {}
    queued, started, failures = launch_broker_batch(
        running,
        [
            {"case": "arm", "stats": "standard_stats.in", "config": "default", "overrides": {}, "cli_options": {}},
            {"case": "bomex", "stats": "standard_stats.in", "config": "default", "overrides": {}, "cli_options": {}},
        ],
        logs,
        2,
        "dash-batch-test-123",
    )

    assert queued == []
    assert started is True
    assert failures == []
    assert running["arm"]["broker_managed"] is True
    assert captured["action"] == "domain_submit_scm_batch"
    assert captured["payload"]["request"]["cases"] == ["arm", "bomex"]
    assert captured["payload"]["native_cli_options"] == {}
