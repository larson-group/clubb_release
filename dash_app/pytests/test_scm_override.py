"""Focused checks for Tune-tab generic SCM overrides."""

import pytest

from dash_app.tune_tab.callbacks_runs import build_request_payload
from dash_app.tune_tab.runtime import write_loss_params_file
from tuner.taylor_metrics import DEFAULT_AGGREGATION_MODE, DEFAULT_LOSS_MODE
from utilities.create_case_namelist import normalize_override_string, override_value


def test_tune_request_keeps_a_generic_scm_override():
    """The free-form Config field reaches the persisted worker request."""
    request = build_request_payload(
        ["arm"],
        [41400],
        [41460],
        [60],
        [20.0],
        [100.0],
        ["cloud_frac"],
        ["C8"],
        ["0.1"],
        ["0.9"],
        1,
        1,
        "random",
        DEFAULT_LOSS_MODE,
        DEFAULT_AGGREGATION_MODE,
        1,
        0.1,
        1,
        1.0,
        1.0e-6,
        "default",
        "-override iiPDF_type=1",
        [],
    )

    assert request["override"] == "iiPDF_type=1"


def test_tune_request_migrates_a_documented_stale_parameter_name(monkeypatch):
    """The UI payload gets current names from its selected configuration."""
    from dash_app.tune_tab import callbacks_runs

    monkeypatch.setattr(callbacks_runs, "load_tunable_names", lambda _config: ["C_uu_shr"])
    request = build_request_payload(
        ["arm"], [41400], [41460], [60], [20.0], [100.0], ["cloud_frac"],
        ["C_uu_shr"], ["0.1"], ["0.9"], 1, 1, "random",
        DEFAULT_LOSS_MODE, DEFAULT_AGGREGATION_MODE, 1, 0.1, 1, 1.0, 1.0e-6,
        "default", "", [],
    )

    assert request["parameter_ranges"] == [
        {"name": "C_uu_shr", "targets": ["C_uu_shr"], "min": 0.1, "max": 0.9}
    ]


def test_tune_request_keeps_linked_physical_targets(monkeypatch):
    from dash_app.tune_tab import callbacks_runs

    monkeypatch.setattr(callbacks_runs, "load_tunable_names", lambda _config: ["C6rt", "C6thl"])
    request = build_request_payload(
        ["arm"], [41400], [41460], [60], [20.0], [100.0], ["cloud_frac"],
        ["C6rt"], ["0"], ["4"], 1, 1, "random", DEFAULT_LOSS_MODE,
        DEFAULT_AGGREGATION_MODE, 1, 0.1, 1, 1.0, 1.0e-6, "default", "", [],
        range_targets=[["C6rt", "C6thl"]],
    )
    assert request["parameter_ranges"] == [
        {"name": "C6rt", "targets": ["C6rt", "C6thl"], "min": 0.0, "max": 4.0}
    ]


def test_top_result_namelist_preserves_all_linked_physical_assignments(tmp_path):
    path = write_loss_params_file(tmp_path, [{"C6rt": 2.5, "C6thl": 2.5}])
    rendered = path.read_text(encoding="utf-8")
    assert "C6rt = 2.5" in rendered
    assert "C6thl = 2.5" in rendered


def test_override_value_accepts_the_optional_cli_prefix():
    """Pasting the complete run_scm fragment remains convenient and safe."""
    assert normalize_override_string(" -override iiPDF_type=1 ") == "iiPDF_type=1"
    assert normalize_override_string("iiPDF_type=1") == "iiPDF_type=1"


def test_override_value_rejects_a_setting_outside_all_namelists():
    with pytest.raises(ValueError, match="not present"):
        override_value("missing_flag=.true.", "&configurable_clubb_flags_nl\n/\n")
