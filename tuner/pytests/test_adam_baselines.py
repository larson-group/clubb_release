from collections import deque
import json
from pathlib import Path

import numpy as np
import pytest

from run_scripts.run_tuner_job import estimate_total_samples, parse_strategy_spec
from tuner.adam_spsa_strategy import AdamSPSAStrategy
from tuner.request import load_request
from tuner.tuning_scheduler import TuningScheduler
from tuner.tuning_strategy import normalize_strategy_config


def _case_defaults(_case, overrides=None):
    values = {
        "les_stats_file": "benchmark.nc",
        "altitude_comparison_range": [0.0, 1000.0],
        "time_average_range": [0, 3600],
        "num_time_windows": 1,
    }
    values.update(overrides or {})
    return values


def _adam_request(tmp_path, monkeypatch, *, batch_size=8, max_workers=5, pairs=2):
    request_path = tmp_path / "request.json"
    request_path.write_text(
        json.dumps(
            {
                "cases": ["atex", "bomex"],
                "selected_fields": ["cloud_frac"],
                "parameter_ranges": [{"name": "C8", "min": 0.1, "max": 0.2}],
                "batch_size": batch_size,
                "max_workers": max_workers,
                "strategy": {
                    "name": "adam",
                    "options": {
                        "max_updates": 3,
                        "learning_rate": 0.01,
                        "perturbation": 0.05,
                        "spsa_pairs": pairs,
                    },
                },
            }
        ),
        encoding="utf-8",
    )
    monkeypatch.setattr("tuner.request.read_case_tuner_defaults", _case_defaults)
    monkeypatch.setattr("tuner.request.supported_fields", lambda: ["cloud_frac"])
    return request_path


def _adam_strategy(*, max_updates=3, pairs=2, chains_per_batch=2, concurrent_batches=1, seed=7):
    return AdamSPSAStrategy(
        param_names=["A", "B"],
        default_params_row=[0.25, 0.75],
        parameter_ranges=[
            {"name": "A", "targets": ["A"], "min": 0.0, "max": 1.0},
            {"name": "B", "targets": ["B"], "min": 0.0, "max": 1.0},
        ],
        max_updates=max_updates,
        learning_rate=0.04,
        perturbation=0.05,
        spsa_pairs=pairs,
        chains_per_batch=chains_per_batch,
        concurrent_batches=concurrent_batches,
        chain_count=chains_per_batch * concurrent_batches,
        seed=seed,
    )


def test_adam_config_and_cli_use_normalized_fractions():
    expected = {
        "name": "adam",
        "options": {
            "max_updates": 20,
            "learning_rate": 0.01,
            "perturbation": 0.05,
            "spsa_pairs": 3,
        },
    }

    assert normalize_strategy_config({"strategy": expected}) == expected
    assert parse_strategy_spec("adam:20:0.01:0.05:3") == expected
    assert estimate_total_samples(
        {
            "cases": ["arm", "bomex"],
            "batch_size": 12,
            "max_workers": 5,
            "strategy": expected,
        }
    ) == str(6 * (2 + 2 * 3 * 20))


def test_adam_request_derives_parallel_chain_layout(tmp_path, monkeypatch):
    request = load_request(_adam_request(tmp_path, monkeypatch))
    options = request["strategy"]["options"]

    assert options["chains_per_batch"] == 2
    assert options["concurrent_batches"] == 3
    assert options["chain_count"] == 6
    assert request["total_samples"] == 6 * (2 + 2 * 2 * 3)


def test_adam_request_requires_complete_spsa_pairs_per_batch(tmp_path, monkeypatch):
    with pytest.raises(RuntimeError, match=r"divisible by 2 \* spsa_pairs"):
        load_request(_adam_request(tmp_path, monkeypatch, batch_size=6, pairs=2))


def test_adam_starts_at_defaults_then_emits_one_exact_probe_batch():
    strategy = _adam_strategy(max_updates=1)
    pending = deque()

    strategy.fill(pending, capacity=20)
    centers = list(pending)
    pending.clear()
    assert len(centers) == 2
    assert centers[0]["selected_params"] == {"A": 0.25, "B": 0.75}
    strategy.tell([
        {"sample_id": sample["sample_id"], "total_loss": 1.0}
        for sample in centers
    ])

    strategy.fill(pending, capacity=8)
    assert len(pending) == 8
    assert all(
        0.0 <= value <= 1.0
        for sample in pending
        for value in sample["selected_params"].values()
    )


def test_adam_nondefault_starts_form_a_seeded_latin_hypercube():
    strategy = _adam_strategy(chains_per_batch=5, seed=19)
    repeated = _adam_strategy(chains_per_batch=5, seed=19)
    starts = np.stack([chain.x for chain in strategy.chains[1:]])

    assert np.array_equal(starts, np.stack([chain.x for chain in repeated.chains[1:]]))
    for dimension in range(starts.shape[1]):
        occupied_strata = np.floor(starts[:, dimension] * len(starts)).astype(int)
        assert sorted(occupied_strata) == list(range(len(starts)))


def test_adam_locked_group_is_one_coordinate_with_equal_physical_targets():
    strategy = AdamSPSAStrategy(
        param_names=["C2rt", "C2thl", "fixed"],
        default_params_row=[0.4, 0.4, 9.0],
        parameter_ranges=[
            {"name": "C2rt", "targets": ["C2rt", "C2thl"], "min": 0.0, "max": 1.0},
        ],
        max_updates=1,
        learning_rate=0.01,
        perturbation=0.05,
        spsa_pairs=1,
        chains_per_batch=1,
        concurrent_batches=1,
        chain_count=1,
        seed=2,
    )
    pending = deque()

    strategy.fill(pending, capacity=1)
    sample = pending.popleft()

    assert sample["selected_params"]["C2rt"] == sample["selected_params"]["C2thl"]
    assert sample["param_row"] == [0.4, 0.4, 9.0]


def test_adam_is_seeded_and_converges_on_two_dimensional_quadratic():
    strategy = _adam_strategy(max_updates=80, pairs=4, chains_per_batch=1, seed=11)
    pending = deque()
    evaluated = []

    while not strategy.is_exhausted():
        strategy.fill(pending, capacity=64)
        assert pending
        completed = []
        while pending:
            sample = pending.popleft()
            x = sample["selected_params"]["A"]
            y = sample["selected_params"]["B"]
            completed.append(
                {"sample_id": sample["sample_id"], "total_loss": (x - 0.2) ** 2 + (y - 0.8) ** 2}
            )
            evaluated.append((x, y))
        strategy.tell(completed)

    assert strategy.estimated_sample_count() == 2 + 2 * 4 * 80
    assert np.linalg.norm(strategy.chains[0].x - np.asarray([0.2, 0.8])) < 0.08
    assert len(evaluated) == strategy.estimated_sample_count()


def test_adam_fixed_coordinates_evaluate_once():
    strategy = AdamSPSAStrategy(
        param_names=["A"],
        default_params_row=[2.0],
        parameter_ranges=[{"name": "A", "min": 2.0, "max": 2.0}],
        max_updates=100,
        learning_rate=0.01,
        perturbation=0.05,
        spsa_pairs=2,
        chains_per_batch=2,
        concurrent_batches=3,
        chain_count=6,
        seed=1,
    )
    pending = deque()
    strategy.fill(pending, capacity=20)
    sample = pending.popleft()
    strategy.tell([{"sample_id": sample["sample_id"], "total_loss": 1.0}])

    assert strategy.is_exhausted()
    assert strategy.estimated_sample_count() == 1


def test_adam_checkpoint_preserves_incomplete_generation(tmp_path):
    request = {
        "cases": ["atex", "bomex"],
        "selected_fields": ["cloud_frac"],
        "batch_size": 4,
        "max_workers": 2,
        "strategy": {"name": "random", "options": {"max_samples": 1}},
        "parameter_ranges": [{"name": "A", "min": 0.0, "max": 1.0}],
    }
    scheduler = TuningScheduler(
        request=request,
        job_dir=tmp_path,
        control_path=tmp_path / "control.json",
        status_path=tmp_path / "status.json",
        results_path=tmp_path / "results.json",
    )
    strategy = _adam_strategy(max_updates=1, pairs=1, chains_per_batch=1)
    pending = deque()
    strategy.fill(pending, capacity=1)
    center = pending.popleft()
    strategy.tell([{"sample_id": center["sample_id"], "total_loss": 1.0}])
    strategy.fill(pending, capacity=2)
    scheduler.strategy = strategy
    scheduler.pending_samples = pending
    scheduler.batches = {
        4: {
            "case_loss_metrics": {"atex": {"already": "complete"}},
        }
    }
    scheduler.next_batch_id = 5
    scheduler._write_resume_checkpoint()

    restored = TuningScheduler(
        request=request,
        job_dir=tmp_path,
        control_path=tmp_path / "control.json",
        status_path=tmp_path / "status.json",
        results_path=tmp_path / "results.json",
    )
    restored._restore_resume_checkpoint()

    assert len(restored.pending_samples) == 2
    assert restored.strategy.chains[0].phase == "probe_waiting"
    assert len(restored.strategy._pending_context) == 2
    assert list(restored.case_jobs) == [{"batch_id": 4, "case_name": "bomex"}]


def test_improvement_scores_are_signed_and_unavailable_for_zero_baseline(tmp_path):
    scheduler = TuningScheduler(
        request={
            "cases": ["atex"],
            "selected_fields": ["cloud_frac"],
            "batch_size": 1,
            "max_workers": 1,
            "strategy": {"name": "random", "options": {"max_samples": 1}},
            "parameter_ranges": [{"name": "C8", "min": 0.0, "max": 1.0}],
        },
        job_dir=tmp_path,
        control_path=Path(tmp_path) / "control.json",
        status_path=Path(tmp_path) / "status.json",
        results_path=Path(tmp_path) / "results.json",
    )
    scheduler.baselines = {
        "clubb_default": {"total_loss": 10.0},
        "override_defaults": {"total_loss": 0.0},
    }

    assert scheduler._improvement_scores(8.0) == {
        "improvement_vs_clubb_default_percent": 20.0,
        "improvement_vs_override_defaults_percent": None,
    }
    assert scheduler._improvement_scores(12.0)["improvement_vs_clubb_default_percent"] == -20.0
    scheduler.baselines["override_defaults"]["total_loss"] = float("nan")
    assert scheduler._improvement_scores(8.0)["improvement_vs_override_defaults_percent"] is None
    assert scheduler._improvement_scores(1.0e30) == {
        "improvement_vs_clubb_default_percent": None,
        "improvement_vs_override_defaults_percent": None,
    }


def test_baseline_and_candidates_share_the_exact_aggregation_path(tmp_path):
    scheduler = TuningScheduler(
        request={
            "cases": ["atex"],
            "selected_fields": ["cloud_frac"],
            "batch_size": 1,
            "max_workers": 1,
            "strategy": {"name": "random", "options": {"max_samples": 2}},
            "parameter_ranges": [{"name": "C8", "min": 0.0, "max": 1.0}],
        },
        job_dir=tmp_path,
        control_path=tmp_path / "control.json",
        status_path=tmp_path / "status.json",
        results_path=tmp_path / "results.json",
    )

    def metrics(centered, bias):
        return {
            "atex": {
                "cloud_frac": {
                    "scaled_rmse": [[centered + bias]],
                    "correlation": [[1.0]],
                    "std_ratio": [[1.0]],
                    "centered_rmse_norm": [[centered]],
                    "bias_norm": [[bias]],
                }
            }
        }

    scheduler.baseline_case_metrics["clubb_default"] = metrics(0.4, 0.2)
    scheduler._finalize_baselines()
    baseline_loss = scheduler.baselines["clubb_default"]["total_loss"]

    for batch_id, case_metrics in ((1, metrics(0.4, 0.2)), (2, metrics(0.2, 0.1))):
        batch = {
            "active_sample_count": 1,
            "sample_ids": [batch_id],
            "selected_params_by_sample": [{"C8": 0.5}],
            "all_params_by_sample": [{"C8": 0.5}],
            "case_loss_metrics": case_metrics,
        }
        scheduler.batches[batch_id] = batch
        scheduler._complete_batch(batch_id, batch)

    same, better = list(scheduler.completed_samples)
    assert same["total_loss"] == baseline_loss
    assert same["improvement_vs_clubb_default_percent"] == pytest.approx(0.0)
    assert better["improvement_vs_clubb_default_percent"] == pytest.approx(50.0)
    assert "override_defaults" not in scheduler.baselines


def test_rejected_baseline_is_unavailable_without_stopping_tune(tmp_path):
    scheduler = TuningScheduler(
        request={
            "cases": ["lba"],
            "selected_fields": ["cloud_frac"],
            "batch_size": 1,
            "max_workers": 1,
            "override": "l_diag_Lscale_from_tau=.true.",
            "strategy": {"name": "random", "options": {"max_samples": 1}},
            "parameter_ranges": [{"name": "C8", "min": 0.0, "max": 1.0}],
        },
        job_dir=tmp_path,
        control_path=tmp_path / "control.json",
        status_path=tmp_path / "status.json",
        results_path=tmp_path / "results.json",
    )
    scheduler.baseline_case_metrics["clubb_default"] = {
        "lba": {"__baseline_unavailable__": "CLUBB rejected this baseline parameter set"}
    }
    scheduler.baseline_case_metrics["override_defaults"] = {
        "lba": {"__baseline_unavailable__": "CLUBB rejected this baseline parameter set"}
    }

    scheduler._finalize_baselines()

    assert scheduler.baselines["clubb_default"] == {
        "available": False,
        "total_loss": None,
        "case_loss": {},
        "unavailable_cases": {"lba": "CLUBB rejected this baseline parameter set"},
    }
    assert scheduler._improvement_scores(1.0) == {
        "improvement_vs_clubb_default_percent": None,
        "improvement_vs_override_defaults_percent": None,
    }
