"""Tests for the compact Plot output-directory picker labels."""

from dash_app.plot_tab.layout import output_directory_options, selected_output_directory_chips
from dash_app.plot_tab.state import DEFAULT_OUTPUT_DIR


def test_output_directory_options_hide_selected_default_from_add_menu():
    options = output_directory_options([], [DEFAULT_OUTPUT_DIR])

    assert options == []


def test_output_directory_options_summarize_discovered_cases():
    options = output_directory_options(
        [
            {
                "path": "/tmp/run-one",
                "relative_path": "output/run-one",
                "case_names": ["arm", "bomex", "dycoms2_rf01", "rico"],
            }
        ]
    )

    assert options == [
        {
            "label": "run-one · arm, bomex, dycoms2_rf01 +1",
            "value": "/tmp/run-one",
        }
    ]


def test_selected_directory_chips_use_short_name_without_case_summary():
    chips = selected_output_directory_chips(
        [{"path": "/tmp/addg1", "relative_path": "output/addg1", "case_names": ["arm", "bomex"]}],
        ["/tmp/addg1"],
    )

    assert chips[0].children[0].children == "addg1"
