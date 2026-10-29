"""Provisional stats selection contracts owned by the shared Python unit stage."""

import json
import re
from pathlib import Path

import pytest

from input.scripts.stats_json_to_namelist import select_stats_entries


def test_presets_and_categories_are_independent_deduplicated_selections(tmp_path):
    """A preset alone adds no baseline; unions add exactly the requested categories."""
    path = tmp_path / "stats.jsonc"
    entries = [{"name": name, "grid": "zt", "units": "m", "long_name": name}
               for name in ("x", "y", "z")]
    path.write_text(json.dumps({"all": {"standard": {"core": entries[:2]}, "extra": entries[2:]},
                                "presets": {"small": ["z", "x", "z"]}}))
    assert select_stats_entries(None, path) == entries[:2]
    assert select_stats_entries("small", path) == [entries[2], entries[0]]
    assert select_stats_entries("small+core+small", path) == entries
    assert select_stats_entries("core,small", path) == entries
    custom = tmp_path / "custom.json"
    custom.write_text('["z","x","z"]')
    assert select_stats_entries(custom, path) == [entries[2], entries[0]]
    with pytest.raises(ValueError, match="Unknown stats category"):
        select_stats_entries("removed", path)


@pytest.mark.parametrize("category, names", [
    ("silhs", {
        "lh_rtp2_mc", "lh_thlp2_mc", "lh_wprtp_mc", "lh_wpthlp_mc", "lh_rtpthlp_mc",
        "lh_samp_frac_category", "silhs_variance_category", "lh_rrm_mc_nonadj",
        "lh_mixt_frac_unweighted", "lh_cloud_frac_unweighted", "lh_precip_frac_unweighted",
    }),
    ("cloud_diagnostics", {
        "w_up_in_cloud", "w_down_in_cloud", "cld_updr_frac", "cld_downdr_frac",
    }),
    ("microphysics/common_microphysics", {"rfrzm"}),
    ("scalars", {
        "wp2sclrp", "wpsclrp2", "wpsclrprtp", "wpsclrpthlp", "sclrm", "sclrm_f",
        "sclrprtp", "sclrp2", "sclrpthvp", "sclrpthlp", "sclrprcp", "wpsclrp",
        "edsclrm", "edsclrm_f", "wpedsclrp",
    }),
])
def test_optional_diagnostics_require_explicit_selection(category, names):
    """Diagnostics absent from the legacy standard list require explicit selection."""
    standard = {entry["name"] for entry in select_stats_entries()}
    optional = {entry["name"] for entry in select_stats_entries(f"non-standard/{category}")}
    combined = {entry["name"] for entry in select_stats_entries(f"standard+non-standard/{category}")}
    assert not standard & names
    assert names <= optional
    assert combined == standard | optional


def test_shared_config_cannot_change_the_stats_selection(tmp_path):
    """The normal input builder uses the same stats for baseline and e3sm settings."""
    from input.scripts.create_case_namelist import available_configs, create_case_settings

    expected = [" | ".join(entry[key] for key in ("name", "grid", "units", "long_name"))
                for entry in select_stats_entries("standard")]
    for config in (None, "e3sm"):
        settings = create_case_settings("bomex", tmp_path, config=config)
        assert settings["clubb_stats_nl"]["entry"] == expected
    assert "all_tuning_stats" not in available_configs()


def test_standard_contains_native_restart_fields():
    """Default output includes the catalog fields requested by native restart reads."""
    root = Path(__file__).resolve().parents[4]
    driver = (root / "src/clubb_driver.F90").read_text()
    restart = driver.split("  subroutine restart_clubb &", 1)[1].split("end subroutine restart_clubb", 1)[0]
    enabled = {name.lower() for name in re.findall(r"\b(l_input_\w+)\s*=\s*\.true\.", restart, re.I)}
    source = (root / "src/Input_fields/input_fields.F90").read_text()
    source = source.split("case ( clubb_input_type )", 1)[1].split("case ( coamps_input_type )", 1)[0]
    reads = re.findall(r'get_clubb_variable_interpolated\s*&?\s*\(\s*(l_input_\w+)\s*,\s*'
                       r'stat_files\([^)]*\)\s*,\s*"([^"]+)"', source, re.I)
    assert enabled and reads, "Native restart flags and reads must be discoverable"
    required = {name for flag, name in reads if flag.lower() in enabled}
    selected = {entry["name"] for entry in select_stats_entries("standard")}
    assert not required - selected, f"Native restart fields missing from standard: {sorted(required - selected)}"
