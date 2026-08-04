from pathlib import Path

from dash import html

from dash_app.plot_tab.callbacks_params import _append_read_only_sections
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
