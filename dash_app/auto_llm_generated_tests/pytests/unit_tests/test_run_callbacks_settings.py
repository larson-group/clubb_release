from dash import Dash, no_update

from dash_app.run_tab import callbacks_settings
from dash_app.run_tab.callbacks_settings import (
    default_multicol_row,
    linked_flag_updates,
    multicol_rows_from_values,
    output_directory_warning,
    register_settings_callbacks,
    remove_multicol_row,
    resolved_output_directory,
    synchronized_flag_values,
    unchanged_multicol_bounds,
)
from dash_app.run_tab.layout import (
    build_console_shell,
    build_flag_controls,
    build_linked_tunable_input,
    build_multicol_section,
)


def _callback(app, name):
    return next(
        entry["callback"].__wrapped__
        for entry in app.callback_map.values()
        if entry.get("callback") is not None
        and entry["callback"].__name__ == name
    )


def test_multicol_section_starts_without_a_parameter_row():
    rows = next(
        component
        for component in build_multicol_section(["C1"])
        if getattr(component, "id", None) == "run-multicol-rows"
    )

    assert rows.children == []


def test_config_selection_resets_multicol_to_no_rows():
    app = Dash(__name__, suppress_callback_exceptions=True)
    register_settings_callbacks(app)

    result = _callback(app, "select_tunable_config")(
        {"name": "e3sm_maint32", "nonce": 1},
        callbacks_settings.available_tunable_configs(),
        "default",
        1,
        [0],
        [{"id": 0, "param": "C1", "min": "0", "max": "1", "npoints": "4"}],
        ["C1"],
        ["0"],
        ["1"],
        ["4"],
    )

    assert result[9:12] == (0, [], [])


def test_add_parameter_still_initializes_the_first_available_parameter():
    assert default_multicol_row(
        0,
        ["C1", "C4"],
        {"C1": {"min": "0.1", "max": "2"}},
        [],
    ) == {"id": 0, "param": "C1", "min": "0.1", "max": "2", "npoints": "4"}


def test_tunable_disabled_states_are_computed_clientside():
    app = Dash(__name__, suppress_callback_exceptions=True)
    register_settings_callbacks(app)
    callback = next(
        entry
        for key, entry in app.callback_map.items()
        if key.endswith('"type":"run-param"}.disabled')
    )

    assert "callback" not in callback
    assert [(item["id"], item["property"]) for item in callback["state"]] == [
        ('{"file":"tunable","name":["ALL"],"type":"run-param"}', "id")
    ]


def test_unchanged_multicol_bounds_matches_wildcard_output_cardinality():
    minimums, maximums = unchanged_multicol_bounds(["C_uu_shr", ""])

    assert minimums == [no_update, no_update]
    assert maximums == [no_update, no_update]


def test_final_multicol_row_can_be_removed_without_reappearing():
    row = {"id": 0, "param": "C8", "min": "0.1", "max": "1", "npoints": "4"}

    assert remove_multicol_row(0, [0], [row]) == (0, [], [])
    assert multicol_rows_from_values([], [], [], [], [], [], fallback_rows=[row]) == []


def test_output_directory_warning_is_hidden_for_missing_or_empty_paths(tmp_path):
    missing_text, missing_class = output_directory_warning(tmp_path / "missing")
    empty = tmp_path / "empty"
    empty.mkdir()
    empty_text, empty_class = output_directory_warning(empty)

    assert missing_text == empty_text == ""
    assert missing_class.endswith("--hidden")
    assert empty_class.endswith("--hidden")


def test_output_directory_warning_reports_existing_contents(tmp_path):
    destination = tmp_path / "existing"
    destination.mkdir()
    (destination / "old_stats.nc").write_text("old")

    text, class_name = output_directory_warning(destination)

    assert "1 existing item" in text
    assert "overwritten" in text
    assert class_name.endswith("--warning")


def test_resolved_output_directory_canonicalizes_aliases_and_rejects_files(
    tmp_path, monkeypatch
):
    monkeypatch.setattr(
        callbacks_settings,
        "resolve_output_dir",
        lambda value: tmp_path / "output" / str(value).removeprefix("output/"),
    )
    expected = str((tmp_path / "output" / "default").resolve())

    assert resolved_output_directory("default") == expected
    assert resolved_output_directory("output/default") == expected

    blocked = tmp_path / "output" / "blocked"
    blocked.parent.mkdir(parents=True)
    blocked.write_text("not a directory")
    assert resolved_output_directory("blocked") is None


def test_linked_tunable_has_one_visible_value_and_physical_member_ids():
    control = build_linked_tunable_input(
        [
            {"name": "C6rt", "value": "1.0"},
            {"name": "C6thl", "value": "1.0"},
        ],
        ("C6rt", "C6thl"),
        65,
        str,
    )
    children = list(control.children)
    visible = children[2]
    physical_rows = children[3:]

    assert visible.id == {"type": "run-linked-param", "group": "C6rt=C6thl"}
    assert visible.value == "1.0"
    assert [row.children.id["name"] for row in physical_rows] == ["C6rt", "C6thl"]
    assert all(row.style == {"display": "none"} for row in physical_rows)
    assert control.style == {"gridRow": "span 2"}


def test_related_flags_remain_individual_checkboxes_inside_one_group():
    controls = build_flag_controls(
        [
            {"name": "l_min_xp2_from_corr_wx", "value": ".true."},
            {"name": "l_enable_relaxed_clipping", "value": ".false."},
        ],
        lambda value: value == ".true.",
        [
            {
                "members": ["l_min_xp2_from_corr_wx", "l_enable_relaxed_clipping"],
                "label": "exclusive",
                "description": "Exactly one mode is enabled.",
            }
        ],
    )

    assert len(controls) == 1
    heading, members = controls[0].children
    assert heading.children[0].children == "exclusive"
    assert [row.children[0].id["name"] for row in members.children] == [
        "l_min_xp2_from_corr_wx",
        "l_enable_relaxed_clipping",
    ]
    assert controls[0].style == {"gridRow": "span 3"}


def test_flag_relationships_only_change_declared_checkbox_companions():
    assert linked_flag_updates("l_min_xp2_from_corr_wx", True) == {"l_enable_relaxed_clipping": False}
    assert linked_flag_updates("l_enable_relaxed_clipping", True) == {"l_min_xp2_from_corr_wx": False}
    assert linked_flag_updates("l_damp_wp2_using_em", True) == {"l_stability_correct_tau_zm": False}
    assert linked_flag_updates("l_stability_correct_tau_zm", True) == {"l_damp_wp2_using_em": False}
    assert linked_flag_updates("l_predict_upwp_vpwp", True) == {}


def test_exclusive_pair_evaluation_matches_the_companion_checkbox_update():
    flags = synchronized_flag_values(
        [
            {"type": "run-flag", "name": "l_min_xp2_from_corr_wx"},
            {"type": "run-flag", "name": "l_enable_relaxed_clipping"},
        ],
        [["on"], ["on"]],
        {"type": "run-flag", "name": "l_enable_relaxed_clipping"},
    )

    assert flags == {
        "l_min_xp2_from_corr_wx": False,
        "l_enable_relaxed_clipping": True,
    }


def test_console_shell_is_owned_by_broker_discovery_in_the_browser():
    shell = build_console_shell(["arm", "bomex"])

    assert shell.id == "run-console-container"
    assert len(shell.children) == 1
    assert shell.children[0].children == "No runs yet."
