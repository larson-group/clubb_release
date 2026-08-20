"""Regression checks for Tune's pattern-matching callback return values."""

from dash import no_update

from dash_app.tune_tab.callbacks_settings import (
    parameter_group_specs,
    parameter_options_by_row,
    parameter_rows_from_controls,
    parameter_targets_by_row,
    removal_click_is_real,
    tune_request_for_save,
    tuned_result_options,
    wildcard_no_update,
)
from dash_app.tune_tab.callbacks_display import (
    _runtime_best_loss_points,
    adam_layout_summary,
    build_results_table,
    mode_options_ready,
    rendered_config_names,
)
from dash_app.tune_tab.callbacks_runs import agent_request_to_tune_controls, build_request_payload
from dash_app.tune_tab.layout import build_param_range_row, build_tuned_config_save_dialog


def component_ids(component):
    ids = []
    if component is None:
        return ids
    component_id = getattr(component, "id", None)
    if component_id is not None:
        ids.append(component_id)
    children = getattr(component, "children", None)
    if not isinstance(children, (list, tuple)):
        children = [] if children is None else [children]
    for child in children:
        if hasattr(child, "children") or hasattr(child, "id"):
            ids.extend(component_ids(child))
    return ids


def test_wildcard_no_update_matches_every_dynamic_component():
    """Dash requires an element per ALL-pattern output, even for no_update."""
    values = wildcard_no_update(["arm", "bomex"])

    assert values == [no_update, no_update]


def test_newly_hydrated_remove_buttons_do_not_delete_their_rows():
    """Dynamic ALL inputs arrive at zero clicks during workspace hydration."""
    row_order = [0, 1]

    assert removal_click_is_real([0, 0], row_order, 0) is False
    assert removal_click_is_real([1, 0], row_order, 1) is False
    assert removal_click_is_real([1, 1], row_order, 1) is True


def test_tune_config_wildcard_lengths_follow_rendered_buttons_not_catalog():
    rendered_ids = [
        {"type": "tune-config-button", "name": "default"},
        {"type": "tune-config-button", "name": "adg2"},
    ]

    assert rendered_config_names(rendered_ids) == ["default", "adg2"]


def test_tuned_config_dialog_selects_only_the_retained_top_16(tmp_path):
    results = [
        {"rank": rank, "total_loss": rank / 10, "params": {"C4": rank}}
        for rank in range(1, 20)
    ]
    options = tuned_result_options(results)
    dialog_ids = component_ids(build_tuned_config_save_dialog())

    assert len(options) == 16
    assert options[0]["value"] == "0"
    assert "Rank 1" in options[0]["label"]
    assert {
        "tune-config-save-result",
        "tune-config-save-name",
        "tune-config-save-submit",
        "tune-config-save-cancel",
    } <= set(dialog_ids)

    job_dir = tmp_path / "job"
    job_dir.mkdir()
    (job_dir / "request.json").write_text(
        '{"config": "default", "override": "C4=2.5"}',
        encoding="utf-8",
    )
    request = tune_request_for_save(
        {"job_dir": str(job_dir)},
        output_roots=[tmp_path],
    )
    assert request == {"config": "default", "override": "C4=2.5"}


def test_runtime_best_loss_ignores_failed_penalty_points():
    points, _signature = _runtime_best_loss_points(
        {
            "best_loss_history": [
                {"sample_count": 1, "loss": 1.0e30},
                {"sample_count": 2, "loss": 4.0},
                {"sample_count": 3, "loss": float("nan")},
                {"sample_count": 4, "loss": 2.0},
            ]
        }
    )

    assert points == [(2, 4.0), (4, 2.0)]


def test_inactive_tune_parameters_remain_visible_with_reason_tooltip():
    states = {
        "C8": {"state": "active", "reason": "Available."},
        "omicron": {
            "state": "inactive-mode",
            "reason": "Only used when l_use_precip_frac is enabled.",
        },
    }

    options = parameter_options_by_row(
        [None],
        ["C8", "omicron"],
        states,
    )[0]

    assert options[0] == {"label": "C8", "value": "C8"}
    assert options[1] == {
        "label": "omicron",
        "value": "omicron",
        "disabled": True,
        "title": "Only used when l_use_precip_frac is enabled.",
    }


def test_selected_inactive_parameter_stays_visible_and_disabled_in_every_row():
    states = {
        "omicron": {
            "state": "inactive-mode",
            "reason": "Only used when l_use_precip_frac is enabled.",
        },
    }

    options = parameter_options_by_row(
        ["omicron", None],
        ["C8", "omicron"],
        states,
    )

    assert [option["value"] for option in options[0]] == ["C8", "omicron"]
    assert [option["value"] for option in options[1]] == ["C8", "omicron"]
    assert options[1][1]["disabled"] is True
    assert "l_use_precip_frac" in options[1][1]["title"]


def test_selected_parameter_remains_visible_but_disabled_elsewhere():
    options = parameter_options_by_row(["C8", None], ["C4", "C8"], {})

    assert options[0][1] == {"label": "C8", "value": "C8"}
    assert options[1][1]["disabled"] is True
    assert "another parameter row" in options[1][1]["title"]


def test_locked_group_helpers_preserve_member_and_row_order():
    member_ids = [
        {"row": 5, "member": 1},
        {"row": 2, "member": 0},
        {"row": 5, "member": 0},
    ]
    member_values = ["C2rtthl", "C8", "C2rt"]

    targets = parameter_targets_by_row(member_ids, member_values, [2, 5])
    rows = parameter_rows_from_controls(member_ids, member_values, [2, 5], ["1", "2"], ["3", "4"])

    assert targets == [["C8"], ["C2rt", "C2rtthl"]]
    assert rows[0]["locked"] is False
    assert rows[1] == {
        "id": 5,
        "targets": ["C2rt", "C2rtthl"],
        "min": "2",
        "max": "4",
        "locked": True,
    }
    assert parameter_group_specs(member_ids, member_values, [2, 5]) == [
        {"label": "C8", "targets": ["C8"]},
        {"label": "C2rt = C2rtthl", "targets": ["C2rt", "C2rtthl"]},
    ]


def test_locked_group_layout_has_members_bracket_action_and_one_shared_range():
    group = build_param_range_row(
        {"id": 7, "targets": ["C2rt", "C2thl"], "min": "0.1", "max": "2"},
        ["C2rt", "C2thl"],
    )
    ids = component_ids(group)

    assert {"type": "tune-range-member", "row": 7, "member": 0} in ids
    assert {"type": "tune-range-member", "row": 7, "member": 1} in ids
    assert {"type": "tune-range-member-add", "index": 7} in ids
    assert ids.count({"type": "tune-range-min", "index": 7}) == 1
    assert ids.count({"type": "tune-range-max", "index": 7}) == 1

    ordinary_ids = component_ids(
        build_param_range_row({"id": 7, "targets": ["C2thl"], "min": "0.1", "max": "2"}, ["C2thl"])
    )
    assert {"type": "tune-range-member-add", "index": 7} not in ordinary_ids


def test_adam_dash_controls_convert_percentages_and_show_derived_cost():
    request = build_request_payload(
        ["arm"], [0], [3600], [3600], [0], [1000],
        ["cloud_frac"], ["C8"], ["0.1"], ["0.2"],
        8, 5, "adam", "shape_first", "quantile_weighted",
        100, 0.1, 200, 1.0, 1.0e-12, "default", "",
        adam_max_updates=10,
        adam_learning_rate_percent=1.5,
        adam_perturbation_percent=4.0,
        adam_spsa_pairs=2,
    )

    assert request["strategy"] == {
        "name": "adam",
        "options": {
            "max_updates": 10,
            "learning_rate": 0.015,
            "perturbation": 0.04,
            "spsa_pairs": 2,
        },
    }
    assert mode_options_ready("adam", None, None, None, None, None, 10, 1.5, 4.0, 2, 8)
    assert "6 total chain(s)" in adam_layout_summary(8, 5, ["arm", "bomex"], 10, 2)
    assert "252 candidate evaluations" in adam_layout_summary(8, 5, ["arm", "bomex"], 10, 2)
    assert "4 baseline evaluation(s) separately" in adam_layout_summary(
        8, 5, ["arm", "bomex"], 10, 2, "l_diag_Lscale_from_tau=.true."
    )


def test_workspace_hydration_restores_locked_groups_and_adam_percentages():
    controls = agent_request_to_tune_controls(
        {
            "config": "default",
            "cases": ["arm"],
            "case_configs": [
                {
                    "case_name": "arm",
                    "time_average_range": [0, 3600],
                    "average_time_seconds": 3600,
                    "altitude_comparison_range": [0, 1000],
                }
            ],
            "selected_fields": ["cloud_frac"],
            "parameter_ranges": [
                {"name": "C2rt", "targets": ["C2rt", "C2thl"], "min": 0.1, "max": 2.0}
            ],
            "strategy": {
                "name": "adam",
                "options": {
                    "max_updates": 12,
                    "learning_rate": 0.025,
                    "perturbation": 0.08,
                    "spsa_pairs": 3,
                },
            },
        },
        {
            "arm": {
                "clubb_fields": ["cloud_frac"],
                "time_average_range": [0, 3600],
                "average_time_seconds": 3600,
                "altitude_comparison_range": [0, 1000],
            }
        },
    )

    assert controls["parameter_rows"][0]["targets"] == ["C2rt", "C2thl"]
    assert controls["adam_max_updates"] == 12
    assert controls["adam_learning_rate_percent"] == 2.5
    assert controls["adam_perturbation_percent"] == 8.0
    assert controls["adam_spsa_pairs"] == 3


def test_locked_group_uses_one_leaderboard_column():
    table = build_results_table(
        [{"total_loss": 1.0, "scaled_rmse_sum": 2.0, "params": {"C2rt": 0.7, "C2thl": 0.7}}],
        [{"label": "C2rt = C2thl", "targets": ["C2rt", "C2thl"]}],
    )

    assert "C2rt = C2thl" in str(table)
    assert str(table).count("0.7") == 1
