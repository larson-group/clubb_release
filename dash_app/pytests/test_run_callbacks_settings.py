from dash import no_update

from dash_app.run_tab.callbacks_settings import (
    linked_flag_updates,
    output_directory_warning,
    synchronized_flag_values,
    unchanged_multicol_bounds,
)
from dash_app.run_tab.layout import build_flag_controls, build_linked_tunable_input


def test_unchanged_multicol_bounds_matches_wildcard_output_cardinality():
    minimums, maximums = unchanged_multicol_bounds(["C_uu_shr", ""])

    assert minimums == [no_update, no_update]
    assert maximums == [no_update, no_update]


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
