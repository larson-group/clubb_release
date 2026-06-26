from collections import deque
import json
import math
from pathlib import Path

from tuner.tuning_scheduler import TuningScheduler
from tuner.taylor_metrics import (
    DEFAULT_AGGREGATION_MODE,
    DEFAULT_LOSS_MODE,
    INVALID_CORRELATION_PENALTY,
    INVALID_LOG_STD_RATIO_PENALTY,
    LOSS_METRIC_NAMES,
    LOSS_POLICY_VERSION,
    aggregate_losses,
    blended_mean_max_loss,
    compute_field_loss_diagnostics,
    field_smart_loss,
)
from tuner.request import load_request
from tuner.job_runtime import TunerJob
from tuner.status import renew_keepalive, should_stop, stop_reason, write_control
from utilities.create_case_namelist import build_tuner_namelist
from tuner.tuning_strategy import (
    RandomUniformStrategy,
    ResolveGridStrategy,
    normalize_strategy_config,
)


def test_field_smart_loss_combines_centered_rmse_and_bias():
    metrics = {"centered_rmse_norm": 0.25, "bias_norm": -0.75}

    assert field_smart_loss(metrics) == 1.0


def test_compute_field_loss_diagnostics_all_named_modes():
    metrics = {
        "scaled_rmse": 1.5,
        "centered_rmse_norm": 0.25,
        "bias_norm": -0.5,
        "correlation": 0.8,
        "std_ratio": 2.0,
    }

    diagnostics = compute_field_loss_diagnostics(metrics)
    losses = diagnostics["per_field_losses"]

    assert math.isclose(losses["scaled_rmse"], 1.5)
    assert math.isclose(losses["centered_rmse_bias"], 0.75)
    assert math.isclose(
        losses["taylor_components"],
        0.25 + 0.5 * 0.5 + 0.25 * 0.2 + 0.25 * abs(math.log(2.0)),
    )
    assert math.isclose(
        losses["taylor_components_squared"],
        0.25**2 + 0.5 * 0.5**2 + 0.25 * 0.2**2 + 0.25 * abs(math.log(2.0)) ** 2,
    )
    assert math.isclose(
        losses["shape_first"],
        0.5 * 0.25 + 0.5 * 0.5 + 1.0 * 0.2 + 0.5 * abs(math.log(2.0)),
    )
    assert math.isclose(
        losses["bias_light_taylor"],
        0.25 + 0.25 * 0.5 + 0.5 * 0.2 + 0.5 * abs(math.log(2.0)),
    )
    assert math.isclose(
        losses["decomposed_taylor"],
        0.5 + 0.2 + abs(math.log(2.0)),
    )
    assert math.isclose(
        diagnostics["loss_components_by_mode"]["taylor_components"]["correlation_penalty"],
        0.25 * 0.2,
    )
    assert "centered_rmse_norm" not in diagnostics["loss_components_by_mode"]["decomposed_taylor"]


def test_invalid_correlation_uses_named_penalty_and_flag():
    diagnostics = compute_field_loss_diagnostics(
        {
            "centered_rmse_norm": 0.0,
            "bias_norm": 0.0,
            "correlation": float("nan"),
            "std_ratio": 1.0,
        }
    )

    assert diagnostics["correlation_invalid"] is True
    assert diagnostics["possible_degenerate_taylor_metrics"] is True
    assert diagnostics["correlation_penalty"] == INVALID_CORRELATION_PENALTY


def test_small_correlation_roundoff_is_clamped_without_invalid_flag():
    diagnostics = compute_field_loss_diagnostics(
        {
            "centered_rmse_norm": 0.0,
            "bias_norm": 0.0,
            "correlation": 1.0 + 1.0e-8,
            "std_ratio": 1.0,
        }
    )

    assert diagnostics["correlation_clamped"] is True
    assert diagnostics["correlation_invalid"] is False
    assert diagnostics["sanitized_correlation"] == 1.0
    assert diagnostics["correlation_penalty"] == 0.0


def test_invalid_std_ratio_uses_named_log_penalty_and_flag():
    diagnostics = compute_field_loss_diagnostics(
        {
            "centered_rmse_norm": 0.0,
            "bias_norm": 0.0,
            "correlation": 1.0,
            "std_ratio": 0.0,
        }
    )

    assert diagnostics["std_ratio_invalid"] is True
    assert diagnostics["std_ratio_nonpositive"] is True
    assert diagnostics["possible_degenerate_taylor_metrics"] is True
    assert diagnostics["abs_log_std_ratio"] == INVALID_LOG_STD_RATIO_PENALTY


def test_blended_mean_max_loss_uses_normalized_mean_and_worst_loss():
    loss = blended_mean_max_loss([1.0, 3.0], [1.0, 1.0])

    assert loss == 0.6 * 2.0 + 0.4 * 3.0


def test_mean_max_aggregation_emits_diagnostics():
    diagnostics = aggregate_losses([1.0, 3.0, 5.0], [1.0, 1.0, 0.0], "mean_max")

    assert diagnostics["loss"] == 0.6 * 2.0 + 0.4 * 3.0
    assert diagnostics["weighted_mean"] == 2.0
    assert diagnostics["max_loss"] == 3.0
    assert diagnostics["active_item_count"] == 2
    assert diagnostics["excluded_item_count"] == 1


def test_mean_worst_quantile_selects_worst_raw_losses_then_weights_subset():
    diagnostics = aggregate_losses([1.0, 10.0, 8.0, 2.0], [1.0, 9.0, 1.0, 1.0], "mean_worst_quantile")
    weighted_mean = (1.0 * 1.0 + 10.0 * 9.0 + 8.0 * 1.0 + 2.0 * 1.0) / 12.0

    assert diagnostics["worst_quantile_item_count"] == 1
    assert diagnostics["worst_quantile_mean"] == 10.0
    assert diagnostics["loss"] == 0.7 * weighted_mean + 0.3 * 10.0


def test_normalize_strategy_accepts_legacy_max_samples():
    strategy = normalize_strategy_config({"max_samples": "12"})

    assert strategy == {"name": "random", "options": {"max_samples": 12}}


def test_tuner_job_create_initializes_keepalive_control(tmp_path):
    job = TunerJob.create({"cases": ["bomex"]}, job_dir=tmp_path / "job")

    assert should_stop(job.control_path) is False
    assert stop_reason(job.control_path) is None
    control = json.loads(job.control_path.read_text(encoding="utf-8"))
    assert control["keepalive_required"] is True
    assert control["keepalive_timeout_seconds"] == 300.0
    assert control["keepalive_action"] == "stop"
    assert control["keepalive_updated_at"]


def test_keepalive_expiration_requests_stop(tmp_path):
    control_path = tmp_path / "control.json"
    write_control(
        control_path,
        stop_requested=False,
        keepalive_required=True,
        keepalive_timeout_seconds=1,
        keepalive_updated_at="2000-01-01T00:00:00Z",
        keepalive_action="stop",
    )

    assert should_stop(control_path) is True
    assert stop_reason(control_path) == "keepalive_expired"


def test_keepalive_renewal_preserves_explicit_stop(tmp_path):
    control_path = tmp_path / "control.json"
    write_control(control_path, stop_requested=True)
    renew_keepalive(control_path)

    assert should_stop(control_path) is True
    assert stop_reason(control_path) == "stop_requested"


def test_load_request_defaults_loss_policy_metadata(tmp_path, monkeypatch):
    request_path = tmp_path / "request.json"
    request_path.write_text(
        json.dumps(
            {
                "case_name": "atex",
                "selected_fields": ["cloud_frac"],
                "parameter_ranges": [{"name": "C8", "min": 0.1, "max": 0.2}],
                "batch_size": 1,
                "max_workers": 1,
                "strategy": {"name": "random", "options": {"max_samples": 1}},
            }
        ),
        encoding="utf-8",
    )
    monkeypatch.setattr(
        "tuner.request.read_case_tuner_defaults",
        lambda _case, overrides=None: {
            "les_stats_file": "benchmark.nc",
            "altitude_comparison_range": [0.0, 1000.0],
            "time_average_range": [0, 3600],
            "num_time_windows": 1,
        },
    )
    monkeypatch.setattr("tuner.request.supported_fields", lambda: ["cloud_frac"])

    request = load_request(request_path)

    assert request["loss_mode"] == DEFAULT_LOSS_MODE
    assert request["aggregation_mode"] == DEFAULT_AGGREGATION_MODE
    assert request["loss_policy_version"] == LOSS_POLICY_VERSION
    assert request["case_configs"] == [
        {
            "case_name": "atex",
            "altitude_comparison_range": [0.0, 1000.0],
            "time_average_range": [0, 3600],
            "num_time_windows": 1,
        }
    ]
    assert "time_window_mode" not in request
    assert "invalid_correlation_penalty" in request["loss_policy_constants"]


def test_load_request_accepts_per_case_window_counts(tmp_path, monkeypatch):
    request_path = tmp_path / "request.json"

    def fake_defaults(_case, overrides=None):
        defaults = {
            "les_stats_file": "benchmark.nc",
            "altitude_comparison_range": [0.0, 1000.0],
            "time_average_range": [0, 3600],
            "num_time_windows": 1,
        }
        defaults.update(overrides or {})
        return defaults

    request_path.write_text(
        json.dumps(
            {
                "case_configs": [
                    {
                        "case_name": "atex",
                        "altitude_comparison_range": [10.0, 900.0],
                        "time_average_range": [600, 3600],
                        "num_time_windows": 1,
                    },
                    {
                        "case_name": "bomex",
                        "altitude_comparison_range": [20.0, 1200.0],
                        "time_average_range": [1200, 4800],
                        "num_time_windows": 4,
                    },
                ],
                "selected_fields": ["cloud_frac"],
                "parameter_ranges": [{"name": "C8", "min": 0.1, "max": 0.2}],
                "batch_size": 1,
                "max_workers": 1,
                "strategy": {"name": "random", "options": {"max_samples": 1}},
            }
        ),
        encoding="utf-8",
    )
    monkeypatch.setattr("tuner.request.read_case_tuner_defaults", fake_defaults)
    monkeypatch.setattr("tuner.request.supported_fields", lambda: ["cloud_frac"])

    request = load_request(request_path)

    assert request["cases"] == ["atex", "bomex"]
    assert request["case_configs"][0]["num_time_windows"] == 1
    assert request["case_configs"][1]["num_time_windows"] == 4
    assert request["case_defaults"]["bomex"]["time_average_range"] == [1200, 4800]


def test_load_request_rejects_les_stats_file_case_config_override(tmp_path, monkeypatch):
    request_path = tmp_path / "request.json"
    request_path.write_text(
        json.dumps(
            {
                "case_configs": [{"case_name": "atex", "les_stats_file": "other.nc"}],
                "selected_fields": ["cloud_frac"],
                "parameter_ranges": [{"name": "C8", "min": 0.1, "max": 0.2}],
                "batch_size": 1,
                "max_workers": 1,
                "strategy": {"name": "random", "options": {"max_samples": 1}},
            }
        ),
        encoding="utf-8",
    )
    monkeypatch.setattr("tuner.request.supported_fields", lambda: ["cloud_frac"])

    try:
        load_request(request_path)
    except RuntimeError as exc:
        assert "les_stats_file" in str(exc)
    else:
        raise AssertionError("expected les_stats_file override to be rejected")


def test_build_tuner_namelist_uses_num_time_windows_only():
    text = build_tuner_namelist(
        {
            "les_stats_file": "benchmark.nc",
            "altitude_comparison_range": [0.0, 1000.0],
            "time_average_range": [0, 3600],
        },
        ["cloud_frac"],
        num_time_windows=1,
    )

    assert "time_window_mode" not in text
    assert "num_time_windows = 1" in text


def test_random_uniform_strategy_fills_pending_queue_in_bounds():
    strategy = RandomUniformStrategy(
        param_names=["C1", "C8"],
        default_params_row=[1.0, 2.0],
        parameter_ranges=[{"name": "C8", "min": 10.0, "max": 11.0}],
        max_samples=2,
        seed=1,
    )
    pending = deque()

    strategy.fill(pending, capacity=4)

    assert len(pending) == 2
    assert strategy.is_exhausted()
    assert all(10.0 <= sample["param_row"][1] <= 11.0 for sample in pending)
    assert all(sample["param_row"][0] == 1.0 for sample in pending)


def test_resolve_grid_strategy_includes_endpoints_and_uses_smaller_even_spacing():
    strategy = ResolveGridStrategy(
        param_names=["C4"],
        default_params_row=[0.0],
        parameter_ranges=[{"name": "C4", "min": 0.0, "max": 1.0}],
        spacing=0.3,
    )
    pending = deque()

    strategy.fill(pending, capacity=10)

    assert [sample["selected_params"]["C4"] for sample in pending] == [0.0, 0.25, 0.5, 0.75, 1.0]
    assert strategy.is_exhausted()
    assert strategy.estimated_sample_count() == 5


def test_scheduler_packs_partial_pending_samples_with_default_rows(tmp_path):
    scheduler = TuningScheduler(
        request={
            "cases": ["atex"],
            "selected_fields": ["cloud_frac"],
            "batch_size": 4,
            "max_workers": 1,
            "strategy": {"name": "random", "options": {"max_samples": 2}},
            "parameter_ranges": [{"name": "C8", "min": 10.0, "max": 11.0}],
        },
        job_dir=tmp_path,
        control_path=Path(tmp_path) / "control.json",
        status_path=Path(tmp_path) / "status.json",
        results_path=Path(tmp_path) / "results.json",
    )
    scheduler.default_params_row = [1.0, 2.0]
    scheduler.strategy = RandomUniformStrategy(
        param_names=["C1", "C8"],
        default_params_row=[1.0, 2.0],
        parameter_ranges=scheduler.request["parameter_ranges"],
        max_samples=2,
        seed=2,
    )

    scheduler._fill_pending_samples()
    scheduler._pack_pending_batches()

    assert len(scheduler.pending_samples) == 0
    assert len(scheduler.batches) == 1
    batch = next(iter(scheduler.batches.values()))
    assert batch["active_sample_count"] == 2
    assert len(batch["sample_ids"]) == 2
    assert len(batch["params_batch"]) == 4
    assert batch["params_batch"][2:] == [[1.0, 2.0], [1.0, 2.0]]


def test_scheduler_stores_selected_loss_and_aggregation_diagnostics(tmp_path):
    scheduler = TuningScheduler(
        request={
            "cases": ["atex"],
            "selected_fields": ["cloud_frac"],
            "batch_size": 1,
            "max_workers": 1,
            "strategy": {"name": "random", "options": {"max_samples": 1}},
            "parameter_ranges": [{"name": "C8", "min": 10.0, "max": 11.0}],
            "loss_mode": "taylor_components",
            "aggregation_mode": "mean_worst_quantile",
        },
        job_dir=tmp_path,
        control_path=Path(tmp_path) / "control.json",
        status_path=Path(tmp_path) / "status.json",
        results_path=Path(tmp_path) / "results.json",
    )
    batch = {
        "active_sample_count": 1,
        "sample_ids": [7],
        "selected_params_by_sample": [{"C8": 10.0}],
        "all_params_by_sample": [{"C8": 10.0}],
        "case_loss_metrics": {
            "atex": {
                "cloud_frac": {
                    "scaled_rmse": [[9.0]],
                    "correlation": [[0.8]],
                    "std_ratio": [[1.0]],
                    "centered_rmse_norm": [[0.25]],
                    "bias_norm": [[-0.5]],
                }
            }
        },
    }
    assert set(batch["case_loss_metrics"]["atex"]["cloud_frac"]) == set(LOSS_METRIC_NAMES)

    scheduler.batches[3] = batch
    scheduler._complete_batch(3, batch)

    result = scheduler.completed_samples[0]
    metrics = result["field_metrics"]["atex"]["cloud_frac"]
    assert result["loss_mode"] == "taylor_components"
    assert result["aggregation_mode"] == "mean_worst_quantile"
    assert result["case_loss_diagnostics"]["atex"]["aggregation_mode"] == "mean_worst_quantile"
    assert result["total_loss_diagnostics"]["aggregation_mode"] == "mean_worst_quantile"
    assert metrics["scaled_rmse"] == 9.0
    assert "centered_rmse_bias" in metrics["per_field_losses"]
    assert metrics["smart_loss"] == metrics["per_field_losses"]["taylor_components"]
    assert metrics["loss_components_by_mode"]["taylor_components"]["abs_bias_norm"] == 0.25


def test_scheduler_handles_different_window_counts_by_case(tmp_path):
    scheduler = TuningScheduler(
        request={
            "cases": ["atex", "bomex"],
            "case_configs": [
                {
                    "case_name": "atex",
                    "altitude_comparison_range": [0.0, 1000.0],
                    "time_average_range": [0, 3600],
                    "num_time_windows": 1,
                },
                {
                    "case_name": "bomex",
                    "altitude_comparison_range": [0.0, 1000.0],
                    "time_average_range": [0, 3600],
                    "num_time_windows": 2,
                },
            ],
            "case_defaults": {
                "atex": {
                    "les_stats_file": "atex.nc",
                    "altitude_comparison_range": [0.0, 1000.0],
                    "time_average_range": [0, 3600],
                    "num_time_windows": 1,
                },
                "bomex": {
                    "les_stats_file": "bomex.nc",
                    "altitude_comparison_range": [0.0, 1000.0],
                    "time_average_range": [0, 3600],
                    "num_time_windows": 2,
                },
            },
            "selected_fields": ["cloud_frac"],
            "batch_size": 1,
            "max_workers": 1,
            "strategy": {"name": "random", "options": {"max_samples": 1}},
            "parameter_ranges": [{"name": "C8", "min": 10.0, "max": 11.0}],
        },
        job_dir=tmp_path,
        control_path=Path(tmp_path) / "control.json",
        status_path=Path(tmp_path) / "status.json",
        results_path=Path(tmp_path) / "results.json",
    )
    batch = {
        "active_sample_count": 1,
        "sample_ids": [8],
        "selected_params_by_sample": [{"C8": 10.0}],
        "all_params_by_sample": [{"C8": 10.0}],
        "case_loss_metrics": {
            "atex": {
                "cloud_frac": {
                    "scaled_rmse": [[1.0]],
                    "correlation": [[1.0]],
                    "std_ratio": [[1.0]],
                    "centered_rmse_norm": [[0.1]],
                    "bias_norm": [[0.1]],
                }
            },
            "bomex": {
                "cloud_frac": {
                    "scaled_rmse": [[2.0], [3.0]],
                    "correlation": [[1.0], [0.9]],
                    "std_ratio": [[1.0], [1.0]],
                    "centered_rmse_norm": [[0.2], [0.3]],
                    "bias_norm": [[0.1], [0.2]],
                }
            },
        },
    }

    scheduler.batches[4] = batch
    scheduler._complete_batch(4, batch)

    result = scheduler.completed_samples[0]
    assert result["case_window_counts"] == {"atex": 1, "bomex": 2}
    assert result["field_metrics"]["atex"]["cloud_frac"]["num_time_windows"] == 1
    assert result["field_metrics"]["bomex"]["cloud_frac"]["num_time_windows"] == 2
    assert len(result["field_metrics"]["bomex"]["cloud_frac"]["subwindows"]) == 2
