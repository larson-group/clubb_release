from pathlib import Path

from dash import html

from dash_app.plot_tab.callbacks_params import (
    _append_read_only_sections,
    _reconcile_compare_params,
)
from dash_app.plot_tab.plot_types.shared import (
    load_compare_flag_values,
    load_flag_values,
)


def _write_run_input(stats_path: Path, body: str, suffix=".in"):
    prefix = str(stats_path)[: -len("_stats.nc")]
    Path(prefix + suffix).write_text(body, encoding="utf-8")


def test_saved_configurable_flags_are_loaded_and_normalized(tmp_path):
    stats_path = tmp_path / "arm_stats.nc"
    _write_run_input(
        stats_path,
        """
&configurable_clubb_flags_nl
iiPDF_type = 2,
l_predict_upwp_vpwp = .false.,
saturation_formula = 3
/
""",
    )

    flags, available = load_flag_values([str(stats_path)])

    assert available
    assert flags == {
        "iiPDF_type": "2",
        "l_predict_upwp_vpwp": "False",
        "saturation_formula": "3",
    }


def test_setup_text_is_used_when_merged_input_is_missing(tmp_path):
    stats_path = tmp_path / "bomex_stats.nc"
    _write_run_input(
        stats_path,
        """
 &configurable_clubb_flags_nl
 --------------------------------------------------
 l_tke_aniso =  T
 iiPDF_type = 1
 --------------------------------------------------
 &next_group
""",
        suffix="_setup.txt",
    )

    flags, available = load_flag_values([str(stats_path)])

    assert available
    assert flags == {"l_tke_aniso": "True", "iiPDF_type": "1"}


def test_compared_flags_separate_matches_and_mismatches(tmp_path):
    first = tmp_path / "first_stats.nc"
    second = tmp_path / "second_stats.nc"
    _write_run_input(
        first,
        "&configurable_clubb_flags_nl\nl_same=.true.\niiPDF_type=1\n/\n",
    )
    _write_run_input(
        second,
        "&configurable_clubb_flags_nl\nl_same=T\niiPDF_type=2\n/\n",
    )

    matched, mismatched, available = load_compare_flag_values(
        [str(first), str(second)]
    )

    assert available
    assert matched == {"l_same": "True"}
    assert mismatched == [("iiPDF_type", "Mismatch across compared outputs")]


def test_read_only_sections_are_collapsed_and_follow_varying_controls():
    varying = html.Div("Displayed varying parameters", id="varying")
    children = _append_read_only_sections(
        [varying],
        [("C1", 1.0), ("beta", 2.4)],
        {
            "flags": {"iiPDF_type": "2", "l_predict_upwp_vpwp": "False"},
            "flag_mismatches": [],
            "has_flags": True,
        },
    )

    assert children[0].id == "varying"
    disclosures = children[1:]
    assert len(disclosures) == 2
    assert all(item.__class__.__name__ == "Details" for item in disclosures)
    assert all(item.open is False for item in disclosures)
    assert disclosures[0].children[0].children == "Constant parameters (2)"
    assert disclosures[1].children[0].children == "Configurable flags (2)"


def test_single_column_difference_does_not_hide_ensemble_axis():
    params, differences, conflicts = _reconcile_compare_params(
        [{"C1": [1.0]}, {"C1": [1.0, 2.0, 3.0]}],
        ["default", "ensemble"],
        3,
    )

    assert params == {"C1": [1.0, 2.0, 3.0]}
    assert conflicts == []
    assert differences == [
        {
            "name": "C1",
            "values": [
                {"source": "default", "value": "1"},
                {"source": "ensemble", "value": "1-3 [1, 2, 3]"},
            ],
            "conflict": False,
        }
    ]


def test_constant_multicol_difference_is_informational():
    params, differences, conflicts = _reconcile_compare_params(
        [{"C1": [2.0, 2.0, 2.0, 2.0]}, {"C1": [3.0]}],
        ["first", "second"],
        4,
    )

    assert params == {}
    assert conflicts == []
    assert differences[0]["conflict"] is False


def test_disagreeing_varying_multicol_axes_disable_parameter_selection():
    params, differences, conflicts = _reconcile_compare_params(
        [
            {"C1": [2.0, 3.0, 4.0, 5.0]},
            {"C1": [2.0, 4.0, 6.0, 8.0]},
        ],
        ["first", "second"],
        4,
    )

    assert params == {}
    assert conflicts == ["C1"]
    assert differences[0]["conflict"] is True
