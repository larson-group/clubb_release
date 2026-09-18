from dash import Dash

from dash_app.run_tab.callbacks_runs import (
    build_multicol_spec,
    complete_run_overrides,
    expand_linked_parameter_values,
    fresh_batch_request_id,
    prepared_run_with_output,
    run_output_rename_available,
    submit_prepared_run,
)
from dash_app.run_tab.tab import build_tab
from dash_app.run_tab.runtime import build_case_command, output_directory_details


def test_complete_run_overrides_freezes_effective_values_not_only_deltas():
    overrides = complete_run_overrides(
        {
            "normalized_flags": {
                "l_predict_upwp_vpwp": True,
                "l_diag_Lscale_from_tau": False,
                "iiPDF_type": "1",
            },
            "normalized_parameters": {
                "tunable": {"C8": "0.5"},
                "silhs": {"l_lh_straight_mc": ".false."},
            },
            "overrides": {},
        }
    )

    assert overrides == {
        "flags": {
            "l_predict_upwp_vpwp": ".true.",
            "l_diag_Lscale_from_tau": ".false.",
            "iiPDF_type": "1",
        },
        "tunable": {"C8": "0.5"},
        "silhs": {"l_lh_straight_mc": ".false."},
    }


def test_cancel_command_does_not_share_the_polling_callback():
    app = Dash(__name__, suppress_callback_exceptions=True)
    app.layout = build_tab(app)

    action_callback = next(
        entry
        for key, entry in app.callback_map.items()
        if "run-action-result.data" in key
    )
    action_inputs = {item["id"] for item in action_callback["inputs"]}
    assert {
        "run-button",
        "run-cancel",
        "run-clear",
        "run-overwrite-button",
        "run-rename-button",
        "run-overwrite-cancel-button",
    } <= action_inputs
    assert "run-sync-interval" not in action_inputs
    assert not any("run-snapshot.data" in key for key in app.callback_map)

    assert not any("run-console-container.children" in key for key in app.callback_map)
    assert not any("run-case-button" in key and ".style" in key for key in app.callback_map)
    assert any("run-ui-render-signal.data" in key for key in app.callback_map)
    render_callback = next(
        entry
        for key, entry in app.callback_map.items()
        if "run-ui-render-signal.data" in key
    )
    assert {item["id"] for item in render_callback["inputs"]} == {
        "run-selected-cases",
        "run-action-result",
        "run-resolved-output-dir",
    }


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


def test_run_command_uses_snapshotted_python_or_jax_build():
    python_command = build_case_command(
        "arm",
        "standard_stats.in",
        {"implementation": "python", "install_dir": "/tmp/build one"},
    )
    jax_command = build_case_command(
        "bomex",
        "standard_stats.in",
        {
            "implementation": "jax",
            "jax_profile": "gpu",
            "install_dir": "/tmp/build-two",
        },
    )

    assert "-python -install_dir '/tmp/build one' arm" in python_command
    assert "-jax=gpu bomex" in jax_command
    assert "-install_dir" not in jax_command


def test_output_directory_details_count_only_stats_cases(tmp_path):
    output = tmp_path / "results"
    output.mkdir()
    (output / "arm_stats.nc").write_bytes(b"CDF")
    (output / "bomex_stats.nc").write_bytes(b"CDF")
    (output / "run.log").write_text("done", encoding="utf-8")
    (output / "nested").mkdir()

    details = output_directory_details(output)

    assert details["path"] == str(output.resolve())
    assert details["nonempty"] is True
    assert details["case_count"] == 2
    assert details["created"] != "Not created yet"
    assert details["last_edited"] != "Not created yet"


def test_output_rename_changes_only_the_frozen_output_target(tmp_path):
    current = tmp_path / "current"
    renamed = tmp_path / "renamed"
    pending = {
        "cases": ["arm"],
        "output_dir": str(current),
        "cli_options": {"out_dir": str(current), "debug": "0"},
    }

    updated = prepared_run_with_output(pending, renamed)

    assert run_output_rename_available(renamed, pending) is True
    assert run_output_rename_available(current, pending) is False
    assert updated["output_dir"] == str(renamed)
    assert updated["cli_options"] == {"out_dir": str(renamed), "debug": "0"}
    assert pending["cli_options"]["out_dir"] == str(current)


def test_prepared_run_submission_preserves_frozen_settings():
    calls = []
    gpu = "GPU-aaaaaaaa-1111-2222-3333-000000000001"

    def perform_action(action, payload, *, internal):
        calls.append((action, payload, internal))
        return {"job_id": "batch-job"}

    result = submit_prepared_run(
        {
            "cases": ["arm", "bomex"],
            "stats": "standard_stats.in",
            "config": "default",
            "overrides": {"flags": {"l_uv_nudge": ".true."}},
            "typed_overrides": {"l_uv_nudge": ".true."},
            "cli_options": {
                "implementation": "jax",
                "jax_profile": "gpu",
                "jax_gpu": gpu,
                "jax_xla_prealloc": False,
                "out_dir": "renamed",
            },
            "typed_options": {"max_iters": 10},
            "max_workers": 2,
            "output_dir": "renamed",
            "implementation": "jax",
            "jax_profile": "gpu",
            "jax_gpu": gpu,
            "jax_xla_prealloc": False,
        },
        perform_action,
    )

    assert result["job_id"] == "batch-job"
    action, payload, internal = calls[0]
    assert action == "domain_submit_scm_batch"
    assert internal is True
    assert payload["request"]["cases"] == ["arm", "bomex"]
    assert payload["request"]["implementation"] == "jax"
    assert payload["request"]["jax_profile"] == "gpu"
    assert payload["request"]["jax_gpu"] == gpu
    assert payload["request"]["jax_xla_prealloc"] is False
    assert payload["request"]["max_workers"] == 2
    assert payload["native_cli_options"]["implementation"] == "jax"
    assert payload["native_cli_options"]["jax_profile"] == "gpu"
    assert payload["native_cli_options"]["jax_gpu"] == gpu
    assert payload["native_cli_options"]["jax_xla_prealloc"] is False
    assert payload["native_cli_options"]["out_dir"] == "renamed"
