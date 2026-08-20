from dash import Dash

from dash_app.run_tab.callbacks_runs import (
    build_multicol_spec,
    complete_run_overrides,
    expand_linked_parameter_values,
    fresh_batch_request_id,
)
from dash_app.run_tab.tab import build_tab
from dash_app.run_tab.runtime import build_case_command


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

    action_inputs = {
        item["id"] for item in app.callback_map["run-action-result.data"]["inputs"]
    }
    assert {"run-button", "run-cancel", "run-clear"} <= action_inputs
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
        {"implementation": "jax", "install_dir": "/tmp/build-two"},
    )

    assert "-python -install_dir '/tmp/build one' arm" in python_command
    assert "-jax -install_dir /tmp/build-two bomex" in jax_command
