"""Focused contract tests for shared Tune presets and linked ranges."""

from collections import deque
import json

import pytest

from run_scripts.run_tuner_job import parse_param_spec
from tuner.enhanced_simann_strategy import EnhancedSimulatedAnnealingStrategy
from tuner.presets import apply_preset, get_preset, list_presets
from tuner.request import load_request
from tuner.tuning_strategy import RandomUniformStrategy, ResolveGridStrategy


def _linked_spec():
    return {"name": "C6rt", "targets": ["C6rt", "C6thl"], "min": 0.0, "max": 4.0}


def _assert_linked(sample):
    assert sample["selected_params"]["C6rt"] == sample["selected_params"]["C6thl"]
    assert sample["param_row"][0] == sample["param_row"][1]


def test_preset_catalog_has_the_reviewed_experiments():
    catalog = {item["name"] for item in list_presets()}
    assert catalog == {"wpxp", "xp2_xpyp", "adg1_geometry", "wp23"}
    wpxp = get_preset("wpxp")
    assert wpxp["parameter_ranges"][0]["targets"] == ["C6rt", "C6thl"]
    wp23 = get_preset("wp23")
    assert wp23["selected_fields"] == ["wp2", "wp3", "Skw_zt"]
    assert len(wp23["parameter_ranges"]) == 22
    assert "min" not in wpxp["parameter_ranges"][0]


def test_presets_define_hourly_subwindows_for_their_standard_cases():
    """A preset must not silently reduce its time range to one average."""
    for preset in (get_preset(item["name"]) for item in list_presets()):
        assert [item["average_time_seconds"] for item in preset["case_configs"]] == [3600, 3600, 3600]


def test_preset_request_retains_explicit_hourly_case_windows():
    request = apply_preset({"preset": "wpxp"})
    assert [item["average_time_seconds"] for item in request["case_configs"]] == [3600, 3600, 3600]


def test_preset_baseline_keeps_explicit_request_pieces():
    request = apply_preset(
        {
            "preset": "wpxp",
            "selected_fields": ["cloud_frac"],
            "parameter_ranges": [{"name": "C8", "min": 0.1, "max": 0.2}],
            "batch_size": 6,
            "max_workers": 10,
        }
    )
    assert request["selected_fields"] == ["cloud_frac"]
    assert request["parameter_ranges"] == [{"name": "C8", "min": 0.1, "max": 0.2}]
    assert request["cases"] == ["arm", "bomex", "dycoms2_rf01"]
    assert request["batch_size"] == 6
    assert request["max_workers"] == 10


def test_linked_cli_parameter_syntax():
    assert parse_param_spec("C6rt=C6thl:0:4") == {
        "name": "C6rt", "targets": ["C6rt", "C6thl"], "min": 0.0, "max": 4.0
    }
    assert parse_param_spec("C7:0:1")["targets"] == ["C7"]


@pytest.mark.parametrize("strategy_class,kwargs", [
    (RandomUniformStrategy, {"max_samples": 1, "seed": 1}),
    (ResolveGridStrategy, {"spacing": 4.0}),
    (EnhancedSimulatedAnnealingStrategy, {"max_iters": 1, "chain_count": 1, "seed": 1}),
])
def test_linked_range_is_one_coordinate_but_two_physical_assignments(strategy_class, kwargs):
    strategy = strategy_class(
        param_names=["C6rt", "C6thl", "C8"],
        default_params_row=[1.0, 1.0, 2.0],
        parameter_ranges=[_linked_spec()],
        **kwargs,
    )
    pending = deque()
    strategy.fill(pending, capacity=4)
    assert pending
    for sample in pending:
        _assert_linked(sample)
        assert set(sample["selected_params"]) == {"C6rt", "C6thl"}


def test_request_expands_preset_then_normalizes_linked_range(tmp_path, monkeypatch):
    path = tmp_path / "request.json"
    path.write_text(
        json.dumps(
            {
                "preset": "wpxp",
                "batch_size": 1,
                "max_workers": 1,
                "strategy": {"name": "random", "options": {"max_samples": 1}},
            }
        ),
        encoding="utf-8",
    )
    monkeypatch.setattr("tuner.request.supported_fields", lambda: ["wprtp", "wpthlp", "rtp2", "thlp2", "rtpthlp"])
    monkeypatch.setattr(
        "tuner.presets.default_parameter_ranges",
        lambda _config: {
            "C6rt": {"default": 2.0, "min": 0.0, "max": 8.0},
            "C6thl": {"default": 2.0, "min": 0.0, "max": 8.0},
            "C6rtb": {"default": 2.0, "min": 0.0, "max": 8.0},
            "C6thlb": {"default": 2.0, "min": 0.0, "max": 8.0},
            "C6rtc": {"default": 1.0, "min": 0.0, "max": 4.0},
            "C6thlc": {"default": 1.0, "min": 0.0, "max": 4.0},
            "C7": {"default": 0.5, "min": 0.0, "max": 1.0},
            "C7b": {"default": 0.5, "min": 0.0, "max": 1.0},
            "C7c": {"default": 0.5, "min": 0.0, "max": 2.0},
            "C6rt_Lscale0": {"default": 5.0, "min": 0.0, "max": 20.0},
            "C6thl_Lscale0": {"default": 5.0, "min": 0.0, "max": 20.0},
            "C7_Lscale0": {"default": 0.2, "min": 0.0, "max": 0.8},
            "wpxp_L_thresh": {"default": 50.0, "min": 0.0, "max": 200.0},
            "c_K6": {"default": 0.2, "min": 0.0, "max": 0.8},
            "nu6": {"default": 4.0, "min": 0.0, "max": 16.0},
            "C_uu_shr": {"default": 0.4, "min": 0.0, "max": 1.0},
            "altitude_threshold": {"default": 75.0, "min": 0.0, "max": 300.0},
        },
    )
    monkeypatch.setattr(
        "tuner.request.read_case_tuner_defaults",
        lambda _case, overrides=None: {
            "les_stats_file": "benchmark.nc", "altitude_comparison_range": [0.0, 1000.0],
            "time_average_range": [0, 3600], "num_time_windows": 1,
        },
    )
    request = load_request(path)
    first = request["parameter_ranges"][0]
    assert request["preset"] == "wpxp"
    assert first["name"] == "C6rt"
    assert first["targets"] == ["C6rt", "C6thl"]


def test_request_rejects_overlapping_linked_targets(tmp_path, monkeypatch):
    path = tmp_path / "request.json"
    path.write_text(
        json.dumps({
            "case_name": "arm", "selected_fields": ["cloud_frac"], "batch_size": 1, "max_workers": 1,
            "strategy": {"name": "random", "options": {"max_samples": 1}},
            "parameter_ranges": [
                {"name": "C6rt", "targets": ["C6rt", "C6thl"], "min": 0, "max": 1},
                {"name": "C6thl", "targets": ["C6thl"], "min": 0, "max": 1},
            ],
        }), encoding="utf-8")
    monkeypatch.setattr("tuner.request.supported_fields", lambda: ["cloud_frac"])
    with pytest.raises(RuntimeError, match="assigned by more than one range"):
        load_request(path)


def test_request_rejects_unknown_linked_target(tmp_path, monkeypatch):
    path = tmp_path / "request.json"
    path.write_text(
        json.dumps({
            "case_name": "arm", "selected_fields": ["cloud_frac"], "batch_size": 1, "max_workers": 1,
            "strategy": {"name": "random", "options": {"max_samples": 1}},
            "parameter_ranges": [{"name": "C6rt", "targets": ["C6rt", "not_a_param"], "min": 0, "max": 1}],
        }), encoding="utf-8")
    monkeypatch.setattr("tuner.request.supported_fields", lambda: ["cloud_frac"])
    with pytest.raises(RuntimeError, match="Unknown tunable parameter target"):
        load_request(path)
