from collections import deque
import json
import math
from pathlib import Path
import threading
from types import SimpleNamespace

import numpy as np
import pytest

from tuner.sample_history import SAMPLE_HISTORY_DIR, SampleHistoryWriter, sample_history_paths
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
    normalize_aggregation_weights,
)
from tuner.request import load_request
from tuner.job_runtime import TunerJob
from tuner import status as status_module
from tuner.status import atomic_write_json, renew_keepalive, should_stop, stop_reason, write_control
from utilities.create_case_namelist import build_tuner_namelist
from tuner.enhanced_simann_strategy import EnhancedSimulatedAnnealingStrategy
from tuner.tuning_strategy import (
    RandomUniformStrategy,
    ResolveGridStrategy,
    normalize_strategy_config,
)


def test_field_smart_loss_combines_centered_rmse_and_bias():
    metrics = {"centered_rmse_norm": 0.25, "bias_norm": -0.75}

    assert field_smart_loss(metrics, "centered_rmse_bias") == 1.0


def test_default_loss_mode_prefers_shape_first():
    assert DEFAULT_LOSS_MODE == "shape_first"


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


def test_quantile_weighted_aggregation_uses_best_to_worst_bin_weights():
    diagnostics = aggregate_losses(
        [1.0, 2.0, 3.0, 4.0],
        aggregation_mode="quantile_weighted",
        aggregation_weights=[0.1, 0.4, 0.4, 0.1],
    )

    assert diagnostics["loss"] == pytest.approx(2.5)
    assert [item["mean"] for item in diagnostics["quantile_bin_means"]] == [1.0, 2.0, 3.0, 4.0]


def test_quantile_weights_normalize_and_require_four_nonnegative_entries():
    assert normalize_aggregation_weights([1.0, 4.0, 4.0, 1.0]) == [0.1, 0.4, 0.4, 0.1]
    with pytest.raises(ValueError, match="exactly 4"):
        normalize_aggregation_weights([0.1, 0.9])
    with pytest.raises(ValueError, match="nonnegative"):
        normalize_aggregation_weights([0.1, -0.1, 0.4, 0.6])


def test_normalize_strategy_accepts_legacy_max_samples():
    strategy = normalize_strategy_config({"max_samples": "12"})

    assert strategy == {"name": "random", "options": {"max_samples": 12}}


def test_normalize_strategy_accepts_simann_defaults():
    strategy = normalize_strategy_config({"strategy": "simann"})

    assert strategy["name"] == "simann"
    assert strategy["options"]["max_iters"] == 200
    assert strategy["options"]["initial_temp"] == 1.0
    assert strategy["options"]["max_final_temp"] == 1.0e-12


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


def test_atomic_json_writers_do_not_share_a_temporary_path(tmp_path, monkeypatch):
    """Overlapping control writers must not delete each other's staging file."""
    control_path = tmp_path / "control.json"
    original_replace = status_module.os.replace
    first_replace_reached = threading.Event()
    release_first_replace = threading.Event()
    call_lock = threading.Lock()
    replace_call_count = 0
    failures = []

    def pause_the_first_replace(source, destination):
        nonlocal replace_call_count
        with call_lock:
            replace_call_count += 1
            is_first = replace_call_count == 1
        if is_first:
            first_replace_reached.set()
            assert release_first_replace.wait(timeout=2)
        return original_replace(source, destination)

    monkeypatch.setattr(status_module.os, "replace", pause_the_first_replace)

    def first_writer():
        try:
            atomic_write_json(control_path, {"writer": "first"})
        except Exception as exc:  # pragma: no cover - assertion below reports it
            failures.append(exc)

    worker = threading.Thread(target=first_writer)
    worker.start()
    assert first_replace_reached.wait(timeout=2)
    atomic_write_json(control_path, {"writer": "second"})
    release_first_replace.set()
    worker.join(timeout=2)

    assert not worker.is_alive()
    assert failures == []
    assert json.loads(control_path.read_text(encoding="utf-8"))["writer"] == "first"


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
    assert request["aggregation_weights"] == [0.1, 0.4, 0.4, 0.1]
    assert request["time_window_aggregation_scope"] == "overall"
    assert request["config"] == "default"
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


def test_load_request_accepts_named_tunable_config(tmp_path, monkeypatch):
    request_path = tmp_path / "request.json"
    request_path.write_text(
        json.dumps(
            {
                "config": "compatible_r8029",
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

    assert request["config"] == "compatible_r8029"


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


def test_load_request_normalizes_simann_strategy_and_chain_count(tmp_path, monkeypatch):
    request_path = tmp_path / "request.json"
    request_path.write_text(
        json.dumps(
            {
                "case_name": "atex",
                "selected_fields": ["cloud_frac"],
                "parameter_ranges": [{"name": "C8", "min": 0.1, "max": 0.2}],
                "batch_size": 8,
                "max_workers": 1,
                "strategy": {"name": "simann", "options": {"max_iters": 12}},
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

    assert request["strategy"]["name"] == "simann"
    assert request["strategy"]["options"]["max_iters"] == 12
    assert request["strategy"]["options"]["chain_count"] == 8
    assert request["batch_size"] == 8
    assert request["total_samples"] == 96


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


def test_enhanced_simann_strategy_evaluates_initial_point_then_reacts_to_loss():
    strategy = EnhancedSimulatedAnnealingStrategy(
        param_names=["C1", "C8"],
        default_params_row=[1.0, 2.0],
        parameter_ranges=[{"name": "C8", "min": 10.0, "max": 11.0}],
        max_iters=3,
        chain_count=1,
        seed=1,
    )
    pending = deque()

    strategy.fill(pending, capacity=4)
    strategy.fill(pending, capacity=4)

    assert len(pending) == 1
    initial = pending.popleft()
    assert 10.0 <= initial["selected_params"]["C8"] <= 11.0
    assert initial["param_row"] == [1.0, initial["selected_params"]["C8"]]
    assert initial["chain_id"] == 0

    strategy.tell([{"sample_id": initial["sample_id"], "total_loss": 4.0}])
    strategy.fill(pending, capacity=4)

    assert len(pending) == 1
    moved = pending.popleft()
    assert moved["sample_id"] == 1
    assert 10.0 <= moved["selected_params"]["C8"] <= 11.0
    strategy.tell([{"sample_id": moved["sample_id"], "total_loss": 3.0}])

    assert strategy.nrgy_opt == 3.0
    assert strategy.xopt == [moved["selected_params"]["C8"]]
    assert strategy.is_exhausted() is False


def test_enhanced_simann_strategy_starts_chains_with_latin_hypercube_points():
    strategy = EnhancedSimulatedAnnealingStrategy(
        param_names=["C8"],
        default_params_row=[0.5],
        parameter_ranges=[{"name": "C8", "min": 0.0, "max": 1.0}],
        max_iters=3,
        chain_count=4,
        seed=3,
    )
    pending = deque()

    strategy.fill(pending, capacity=4)
    strategy.fill(pending, capacity=4)

    assert len(pending) == 4
    assert {sample["chain_id"] for sample in pending} == {0, 1, 2, 3}
    values = sorted(sample["selected_params"]["C8"] for sample in pending)
    for idx, value in enumerate(values):
        assert idx / 4 <= value < (idx + 1) / 4


def test_enhanced_simann_strategy_stops_at_max_iters():
    strategy = EnhancedSimulatedAnnealingStrategy(
        param_names=["C8"],
        default_params_row=[0.5],
        parameter_ranges=[{"name": "C8", "min": 0.0, "max": 1.0}],
        max_iters=2,
        chain_count=1,
        seed=2,
    )
    pending = deque()

    strategy.fill(pending, capacity=1)
    initial = pending.popleft()
    strategy.tell([{"sample_id": initial["sample_id"], "total_loss": 1.0}])
    loss = 0.99
    while not strategy.is_exhausted():
        strategy.fill(pending, capacity=1)
        moved = pending.popleft()
        strategy.tell([{"sample_id": moved["sample_id"], "total_loss": loss}])
        loss -= 0.01

    assert strategy.is_exhausted() is True
    assert strategy.iter >= strategy.max_iters


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


def test_scheduler_stop_during_worker_initialization_is_graceful(tmp_path, monkeypatch):
    scheduler = TuningScheduler(
        request={
            "cases": ["bomex"],
            "selected_fields": ["cloud_frac"],
            "batch_size": 1,
            "max_workers": 1,
            "strategy": {"name": "random", "options": {"max_samples": 2}},
            "parameter_ranges": [{"name": "C8", "min": 0.1, "max": 0.2}],
        },
        job_dir=tmp_path,
        control_path=tmp_path / "control.json",
        status_path=tmp_path / "status.json",
        results_path=tmp_path / "results.json",
    )
    monkeypatch.setattr(scheduler, "_initialize_case_barrier", lambda _started_at: False)

    assert scheduler.run() == 0
    assert json.loads(scheduler.status_path.read_text())["state"] == "stopped"
    assert json.loads(scheduler.results_path.read_text())["state"] == "stopped"


def test_scheduler_rebalances_only_an_idle_extra_worker_to_the_slow_case(tmp_path, monkeypatch):
    """A short-case warm actor can be safely repurposed at the next idle point."""
    scheduler = TuningScheduler(
        request={
            "cases": ["short", "long"],
            "selected_fields": ["cloud_frac"],
            "batch_size": 1,
            "max_workers": 3,
            "strategy": {"name": "random", "options": {"max_samples": 4}},
            "parameter_ranges": [{"name": "C8", "min": 0.1, "max": 0.2}],
        },
        job_dir=tmp_path,
        control_path=tmp_path / "control.json",
        status_path=tmp_path / "status.json",
        results_path=tmp_path / "results.json",
    )
    sent = []
    scheduler.workers = [
        SimpleNamespace(worker_id=0, case_name="short", state="idle", conn=SimpleNamespace(send=lambda value: sent.append(value))),
        SimpleNamespace(worker_id=1, case_name="short", state="busy", conn=SimpleNamespace(send=lambda value: sent.append(value))),
        SimpleNamespace(worker_id=2, case_name="long", state="busy", conn=SimpleNamespace(send=lambda value: sent.append(value))),
    ]
    scheduler.case_evaluation_seconds.update({"short": 1.0, "long": 10.0})
    scheduler.batches = {4: {"batch_id": 4}}
    scheduler.case_jobs = deque([
        {"batch_id": 4, "case_name": "long"},
        {"batch_id": 4, "case_name": "long"},
        {"batch_id": 4, "case_name": "long"},
    ])
    scheduler.next_rebalance_monotonic = 0.0

    scheduler._rebalance_warm_workers_if_due()

    assert scheduler.workers[0].state == "stopping"
    assert scheduler.pending_replacements == deque([{"from_case": "short", "to_case": "long"}])
    assert scheduler.last_rebalance["from_case"] == "short"
    assert scheduler.last_rebalance["to_case"] == "long"
    assert scheduler.rebalance_interval_index == 0
    assert sent == [{"type": "stop"}]


def test_scheduler_reassigns_one_idle_worker_when_all_case_drains_are_under_five_seconds(tmp_path):
    """An idle worker can help the largest short queue without batch churn."""
    scheduler = TuningScheduler(
        request={
            "cases": ["short", "long"],
            "selected_fields": ["cloud_frac"],
            "batch_size": 1,
            "max_workers": 2,
            "strategy": {"name": "random", "options": {"max_samples": 4}},
            "parameter_ranges": [{"name": "C8", "min": 0.1, "max": 0.2}],
        },
        job_dir=tmp_path,
        control_path=tmp_path / "control.json",
        status_path=tmp_path / "status.json",
        results_path=tmp_path / "results.json",
    )
    sent = []
    scheduler.workers = [
        SimpleNamespace(
            worker_id=0,
            case_name="short",
            state="idle",
            conn=SimpleNamespace(send=lambda value: sent.append(value)),
        ),
        SimpleNamespace(worker_id=1, case_name="long", state="busy"),
    ]
    scheduler.case_evaluation_seconds.update({"short": 0.2, "long": 1.0})
    scheduler.batches = {4: {"batch_id": 4}}
    scheduler.case_jobs = deque([{"batch_id": 4, "case_name": "long"}] * 4)
    scheduler.next_rebalance_monotonic = 0.0

    scheduler._rebalance_warm_workers_if_due()

    assert scheduler.workers[0].state == "stopping"
    assert scheduler.pending_replacements == deque([{"from_case": "short", "to_case": "long"}])
    assert scheduler.last_rebalance["reason"] == "idle_short_drain"
    assert sent == [{"type": "stop"}]


def test_scheduler_never_rebalances_a_case_away_from_its_only_warm_worker(tmp_path):
    scheduler = TuningScheduler(
        request={
            "cases": ["short", "long"],
            "selected_fields": ["cloud_frac"],
            "batch_size": 1,
            "max_workers": 2,
            "strategy": {"name": "random", "options": {"max_samples": 2}},
            "parameter_ranges": [{"name": "C8", "min": 0.1, "max": 0.2}],
        },
        job_dir=tmp_path,
        control_path=tmp_path / "control.json",
        status_path=tmp_path / "status.json",
        results_path=tmp_path / "results.json",
    )
    scheduler.workers = [
        SimpleNamespace(worker_id=0, case_name="short", state="idle"),
        SimpleNamespace(worker_id=1, case_name="long", state="busy"),
    ]
    scheduler.case_evaluation_seconds.update({"short": 1.0, "long": 10.0})
    scheduler.batches = {4: {"batch_id": 4}}
    scheduler.case_jobs = deque([{"batch_id": 4, "case_name": "long"}] * 4)
    scheduler.next_rebalance_monotonic = 0.0

    scheduler._rebalance_warm_workers_if_due()

    assert scheduler.workers[0].state == "idle"
    assert not scheduler.pending_replacements


def test_scheduler_marks_busy_donor_to_retire_after_current_batch(tmp_path):
    """A rebalanced busy donor finishes safely instead of waiting for an idle poll."""
    scheduler = TuningScheduler(
        request={
            "cases": ["short", "long"],
            "selected_fields": ["cloud_frac"],
            "batch_size": 1,
            "max_workers": 3,
            "strategy": {"name": "random", "options": {"max_samples": 4}},
            "parameter_ranges": [{"name": "C8", "min": 0.1, "max": 0.2}],
        },
        job_dir=tmp_path,
        control_path=tmp_path / "control.json",
        status_path=tmp_path / "status.json",
        results_path=tmp_path / "results.json",
    )
    sent = []
    scheduler.workers = [
        SimpleNamespace(worker_id=0, case_name="short", state="busy", conn=SimpleNamespace(send=lambda value: sent.append(value))),
        SimpleNamespace(worker_id=1, case_name="short", state="busy", conn=SimpleNamespace(send=lambda value: sent.append(value))),
        SimpleNamespace(worker_id=2, case_name="long", state="busy", conn=SimpleNamespace(send=lambda value: sent.append(value))),
    ]
    scheduler.case_evaluation_seconds.update({"short": 1.0, "long": 10.0})
    scheduler.batches = {4: {"batch_id": 4}}
    scheduler.case_jobs = deque([{"batch_id": 4, "case_name": "long"}] * 3)
    scheduler.next_rebalance_monotonic = 0.0

    scheduler._rebalance_warm_workers_if_due()

    assert scheduler.workers[0].state == "busy"
    assert scheduler.workers[0].retire_after_batch is True
    assert scheduler.pending_replacements == deque([{"from_case": "short", "to_case": "long"}])
    assert sent == []


def test_scheduler_reports_stale_compiled_python_api_parameter_names(tmp_path):
    scheduler = TuningScheduler(
        request={
            "cases": ["atex"],
            "selected_fields": ["cloud_frac"],
            "batch_size": 1,
            "max_workers": 1,
            "strategy": {"name": "random", "options": {"max_samples": 1}},
            "parameter_ranges": [{"name": "C_uu_shr", "min": 0.0, "max": 1.0}],
        },
        job_dir=tmp_path,
        control_path=Path(tmp_path) / "control.json",
        status_path=Path(tmp_path) / "status.json",
        results_path=Path(tmp_path) / "results.json",
    )
    scheduler.param_names = ["C8"]

    with pytest.raises(RuntimeError, match="compiled CLUBB Python API.*compile.py -python"):
        scheduler._assert_compiled_parameter_compatibility()


def test_scheduler_rejects_ranges_outside_compiled_fortran_hard_bounds(tmp_path):
    scheduler = TuningScheduler(
        request={
            "cases": ["atex"],
            "selected_fields": ["cloud_frac"],
            "batch_size": 1,
            "max_workers": 1,
            "strategy": {"name": "random", "options": {"max_samples": 1}},
            "parameter_ranges": [{"name": "C8", "min": -0.1, "max": 1.0}],
        },
        job_dir=tmp_path,
        control_path=Path(tmp_path) / "control.json",
        status_path=Path(tmp_path) / "status.json",
        results_path=Path(tmp_path) / "results.json",
    )
    scheduler.param_names = ["C8"]
    scheduler.parameter_hard_bounds = {
        "C8": {"name": "C8", "min": 0.0, "max": 1.0}
    }

    with pytest.raises(RuntimeError, match="violate compiled CLUBB hard bounds.*C8"):
        scheduler._assert_compiled_parameter_compatibility()


def test_scheduler_honors_interior_compiled_fortran_hard_bounds(tmp_path):
    scheduler = TuningScheduler(
        request={
            "cases": ["atex"],
            "selected_fields": ["cloud_frac"],
            "batch_size": 1,
            "max_workers": 1,
            "strategy": {"name": "random", "options": {"max_samples": 1}},
            "parameter_ranges": [{"name": "C_uu_shr", "min": 0.0, "max": 0.5}],
        },
        job_dir=tmp_path,
        control_path=Path(tmp_path) / "control.json",
        status_path=Path(tmp_path) / "status.json",
        results_path=Path(tmp_path) / "results.json",
    )
    scheduler.param_names = ["C_uu_shr"]
    scheduler.parameter_hard_bounds = {
        "C_uu_shr": {
            "name": "C_uu_shr", "min": 1.0e-12, "max": 1.0 - 1.0e-12,
        }
    }

    with pytest.raises(RuntimeError, match="requested min 0 must be >= 1e-12"):
        scheduler._assert_compiled_parameter_compatibility()


def test_scheduler_defaults_simann_chain_count_from_parallel_capacity(tmp_path):
    scheduler = TuningScheduler(
        request={
            "cases": ["atex"],
            "selected_fields": ["cloud_frac"],
            "batch_size": 8,
            "max_workers": 4,
            "strategy": {"name": "simann", "options": {"max_iters": 4}},
            "parameter_ranges": [{"name": "C8", "min": 10.0, "max": 11.0}],
        },
        job_dir=tmp_path,
        control_path=Path(tmp_path) / "control.json",
        status_path=Path(tmp_path) / "status.json",
        results_path=Path(tmp_path) / "results.json",
    )

    assert scheduler.batch_size == 8
    assert scheduler.request["batch_size"] == 8
    assert scheduler.request["strategy"]["options"]["chain_count"] == 32
    assert scheduler.max_pending_batches == 4
    assert scheduler.max_pending_samples == 32


def test_scheduler_bounds_each_case_queue_to_twice_the_worker_count(tmp_path):
    scheduler = TuningScheduler(
        request={
            "cases": ["arm", "bomex", "dycoms2_rf01"],
            "selected_fields": ["cloud_frac"],
            "batch_size": 2,
            "max_workers": 4,
            "strategy": {"name": "random", "options": {"max_samples": 20}},
            "parameter_ranges": [{"name": "C8", "min": 10.0, "max": 11.0}],
        },
        job_dir=tmp_path,
        control_path=Path(tmp_path) / "control.json",
        status_path=Path(tmp_path) / "status.json",
        results_path=Path(tmp_path) / "results.json",
    )

    assert scheduler.max_pending_case_jobs_per_case == 8
    assert scheduler.max_pending_batches == 8
    assert scheduler.max_pending_batches == scheduler.max_pending_case_jobs_per_case
    assert scheduler.max_pending_samples == 16


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
            "loss_mode": "centered_rmse_bias",
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
    assert result["time_window_aggregation_scope"] == "overall"
    # The default policy ranks all three active time-window losses together:
    # 0.2, 0.3, and 0.5.  The empty worst bin is renormalized away.
    assert result["total_loss"] == pytest.approx((0.1 * 0.2 + 0.4 * 0.3 + 0.4 * 0.5) / 0.9)


def test_sample_history_writer_stores_observation_axis_without_window_padding(tmp_path):
    writer = SampleHistoryWriter(
        job_dir=tmp_path,
        param_names=["C1", "C8"],
        case_names=["atex", "bomex"],
        field_names=["cloud_frac", "rcm"],
        metric_names=list(LOSS_METRIC_NAMES),
        case_configs=[
            {
                "case_name": "atex",
                "time_average_range": [0, 3600],
                "num_time_windows": 1,
            },
            {
                "case_name": "bomex",
                "time_average_range": [0, 3600],
                "num_time_windows": 2,
            },
        ],
        case_window_counts={"atex": 1, "bomex": 2},
        chunk_size=2,
    )

    def field_entry(base):
        return {
            "subwindows": [
                {
                    "scaled_rmse": base,
                    "correlation": base + 0.1,
                    "std_ratio": base + 0.2,
                    "centered_rmse_norm": base + 0.3,
                    "bias_norm": base + 0.4,
                },
                {
                    "scaled_rmse": base + 1.0,
                    "correlation": base + 1.1,
                    "std_ratio": base + 1.2,
                    "centered_rmse_norm": base + 1.3,
                    "bias_norm": base + 1.4,
                },
            ]
        }

    entries = [
        {
            "sample_id": 0,
            "batch_id": 5,
            "all_params": {"C1": 1.0, "C8": 0.5},
            "field_metrics": {
                "atex": {
                    "cloud_frac": {"subwindows": field_entry(10.0)["subwindows"][:1]},
                    "rcm": {"subwindows": field_entry(20.0)["subwindows"][:1]},
                },
                "bomex": {
                    "cloud_frac": field_entry(30.0),
                    "rcm": field_entry(40.0),
                },
            },
        },
        {
            "sample_id": 1,
            "batch_id": 6,
            "all_params": {"C1": 2.0, "C8": 0.7},
            "field_metrics": {
                "atex": {
                    "cloud_frac": {"subwindows": field_entry(50.0)["subwindows"][:1]},
                    "rcm": {"subwindows": field_entry(60.0)["subwindows"][:1]},
                },
                "bomex": {
                    "cloud_frac": field_entry(70.0),
                    "rcm": field_entry(80.0),
                },
            },
        },
    ]

    writer.append(entries)
    writer.close()

    paths = sample_history_paths(tmp_path)
    assert len(paths) == 1
    assert paths[0].parent == tmp_path / SAMPLE_HISTORY_DIR
    history = np.load(paths[0])
    assert history["all_params"].shape == (2, 2)
    assert history["loss_metrics"].shape == (2, 6, len(LOSS_METRIC_NAMES))
    assert list(history["param_names"]) == ["C1", "C8"]
    assert list(history["case_names"]) == ["atex", "bomex"]
    assert list(history["field_names"]) == ["cloud_frac", "rcm"]
    assert list(history["metric_names"]) == list(LOSS_METRIC_NAMES)
    assert list(history["obs_case"]) == [0, 0, 1, 1, 1, 1]
    assert list(history["obs_field"]) == [0, 1, 0, 0, 1, 1]
    assert list(history["obs_window"]) == [0, 0, 0, 1, 0, 1]
    np.testing.assert_allclose(history["all_params"], [[1.0, 0.5], [2.0, 0.7]])
    assert history["loss_metrics"][0, 0, list(LOSS_METRIC_NAMES).index("scaled_rmse")] == 10.0
    assert history["loss_metrics"][0, 3, list(LOSS_METRIC_NAMES).index("bias_norm")] == 31.4
    assert history["loss_metrics"][1, 5, list(LOSS_METRIC_NAMES).index("correlation")] == 81.1
