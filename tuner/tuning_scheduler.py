"""Evaluate proposed parameter sets across cases and return one score to the optimizer.

For example, a candidate C8/C11 row must run in both BOMEX and ATEX before the
optimizer can use its combined loss or gradient. TuningScheduler owns that
coordination: it starts reusable case workers, packs candidates into model
column batches, waits for all selected cases, combines their prepared losses
and optional derivatives, then sends completed feedback to the algorithm.
The algorithm selects new points; the case loss driver computes each case's
objective. This scheduler does neither numerical search nor profile-error math.

Usage
-----
run_scheduler(request, job_dir=..., control_path=..., status_path=...,
results_path=...) runs one validated experiment and returns an exit code.
TuningScheduler provides the same lifecycle through construction followed by
run(). WorkerHandle records a worker process, case ownership and outstanding
work. Worker setup, baseline scoring, result ranking, full sample history and
public progress files are managed during that lifecycle.

The request fixes cases, per-case windows, fields, physical parameter ranges,
algorithm options and backend. Every case receives the same candidate rows;
workers remain warm across batches. A graceful stop drains active evaluations
and saves pending proposals, algorithm moments/random state and scheduling
progress in resume_checkpoint.pkl. A checkpoint is a saved Python state file
used to continue the stopped search without discarding its progress. It is
separate from JAX timestep checkpoints, which manage differentiation memory.
Current-format state loads through pickle; obsolete formats are unsupported.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass
import math
import numpy as np
import multiprocessing
import os
from pathlib import Path
import pickle
import tempfile
import time
import traceback

from tuner.sample_history import SampleHistoryWriter
from tuner.tuning_strategy import build_tuning_strategy
from tuner.tuning_worker import worker_main
from tuner.taylor_metrics import (
    AGGREGATION_MODE_NAMES,
    DEFAULT_AGGREGATION_WEIGHTS,
    DEFAULT_AGGREGATION_MODE,
    LOSS_METRIC_NAMES,
    LOSS_POLICY_CONSTANTS,
    LOSS_POLICY_VERSION,
    WORST_QUANTILE_FRACTION,
    aggregate_losses,
    normalize_aggregation_weights,
)
from utilities.loss_metrics import normalize_loss_weights, normalize_named_weights
from tuner.status import (
    RESULTS_FILE_LIMIT,
    should_stop,
    stop_reason,
    utc_now_iso,
    write_results,
    write_status,
)


POLL_INTERVAL_SECONDS = 0.05
STATUS_HEARTBEAT_SECONDS = 1.0
# Rebalancing intentionally starts responsively, then becomes cheap for long
# jobs.  These are elapsed intervals between checks, not a busy polling loop.
REBALANCE_INTERVALS_SECONDS = (5.0, 10.0, 15.0)
REBALANCE_IMBALANCE_RATIO = 2.5
# Do not churn warm F2PY actors when every case queue will clear essentially
# immediately.  Relative imbalance is useful only once there is enough work
# to amortize an actor restart.
REBALANCE_MIN_DRAIN_SECONDS = 5.0
REBALANCE_MIN_WARM_WORKERS_PER_CASE = 1
EVALUATION_TIME_EWMA_ALPHA = 0.30
BEST_LOSS_HISTORY_LIMIT = 300


@dataclass
class WorkerHandle:
    """Track one worker process and its current scheduler state."""

    worker_id: int
    case_name: str
    process: multiprocessing.Process
    conn: object
    worker_dir: Path
    state: str = "starting"
    current_batch_id: int | None = None
    current_job_started_monotonic: float | None = None
    retire_after_batch: bool = False


class TuningScheduler:
    """Coordinate workers, parameter batches, aggregation, and public JSON writes."""

    def __init__(
        self,
        *,
        request: dict,
        job_dir: Path,
        control_path: Path,
        status_path: Path,
        results_path: Path,
    ):
        self.request = request
        self.job_dir = Path(job_dir)
        self.control_path = Path(control_path)
        self.status_path = Path(status_path)
        self.results_path = Path(results_path)
        self.cases = list(request["cases"])
        self.case_configs = {
            str(config["case_name"]): dict(config)
            for config in request.get("case_configs", [])
        }
        self.case_defaults = dict(request.get("case_defaults", {}))
        self.config = str(request.get("config") or "default").strip() or "default"
        self.override = str(request.get("override") or "").strip()
        self.selected_fields = list(request["selected_fields"])
        self.batch_size = int(request["batch_size"])
        self.max_workers = int(request["max_workers"])
        from tuner.algorithm_registry import prepare_strategy
        prepare_strategy(request)
        self.strategy_name = request["strategy"]["name"]
        from utilities.clubb_loss_runner import validate_gradient_backend
        self.grad = request.get("grad", False)
        validate_gradient_backend(request.get("backend", "fortran"), self.grad)
        self.gradient_parameters = list(dict.fromkeys(target for spec in request["parameter_ranges"]
                                                     for target in spec.get("targets", [spec["name"]])))
        self.case_weights = dict(request.get("case_weights", {}))
        self.field_weights = normalize_named_weights(request.get("field_weights"), self.selected_fields, "field")
        self.loss_weights = normalize_loss_weights(request.get("loss_weights"))
        self.aggregation_mode = request.get("aggregation_mode", DEFAULT_AGGREGATION_MODE)
        self.aggregation_weights = normalize_aggregation_weights(
            request.get("aggregation_weights", DEFAULT_AGGREGATION_WEIGHTS)
        )
        self.case_num_time_windows = {
            case_name: int(
                self.case_configs.get(case_name, {}).get(
                    "num_time_windows",
                    self.case_defaults.get(case_name, {}).get("num_time_windows", 1),
                )
            )
            for case_name in self.cases
        }
        if self.aggregation_mode not in AGGREGATION_MODE_NAMES:
            raise RuntimeError(f"Unknown aggregation mode: {self.aggregation_mode}")
        invalid_window_cases = [
            case_name for case_name, count in self.case_num_time_windows.items()
            if count < 1
        ]
        if invalid_window_cases:
            raise RuntimeError("num_time_windows must be >= 1 for: " + ", ".join(invalid_window_cases))
        self.loss_policy_version = request.get("loss_policy_version", LOSS_POLICY_VERSION)
        if self.loss_policy_version != LOSS_POLICY_VERSION:
            raise RuntimeError(
                f"Unsupported loss_policy_version {self.loss_policy_version}; expected {LOSS_POLICY_VERSION}"
            )
        self.loss_policy_constants = dict(LOSS_POLICY_CONSTANTS)
        self.aggregation_options = {
            "worst_quantile_fraction": WORST_QUANTILE_FRACTION,
            "quantile_weights": list(self.aggregation_weights),
        }
        self.request["loss_weights"] = dict(self.loss_weights)
        self.request["config"] = self.config
        self.request["override"] = self.override
        self.request["aggregation_mode"] = self.aggregation_mode
        self.request["case_configs"] = [
            {
                "case_name": case_name,
                "altitude_comparison_range": list(self.case_defaults.get(case_name, {}).get("altitude_comparison_range", [])),
                "time_average_range": list(self.case_defaults.get(case_name, {}).get("time_average_range", [])),
                "num_time_windows": int(self.case_num_time_windows[case_name]),
                **(
                    {"average_time_seconds": int(self.case_defaults[case_name]["average_time_seconds"])}
                    if self.case_defaults.get(case_name, {}).get("average_time_seconds") is not None
                    else {}
                ),
            }
            for case_name in self.cases
        ]
        self.request["case_window_counts"] = dict(self.case_num_time_windows)
        self.request["aggregation_weights"] = list(self.aggregation_weights)
        self.request["loss_policy_version"] = self.loss_policy_version
        self.request["loss_policy_constants"] = dict(self.loss_policy_constants)
        self.request["aggregation_options"] = dict(self.aggregation_options)
        # Each batch contributes one case job to every selected case. Keep up
        # to two workers' worth of jobs queued *per case*, so a slower case
        # retains enough visible backlog for the drain-time balancer to act.
        self.max_pending_case_jobs_per_case = 2 * self.max_workers
        case_job_limited_batches = self.max_pending_case_jobs_per_case
        self.max_pending_batches = case_job_limited_batches
        self.max_pending_samples = self.max_pending_batches * self.batch_size
        self.worker_cap = max(self.max_workers, len(self.cases))
        self.results_file_limit = RESULTS_FILE_LIMIT

        self.ctx = multiprocessing.get_context("spawn")
        self.workers: list[WorkerHandle] = []
        self.next_worker_id = 0
        self.next_batch_id = 0
        self.samples_evaluated = 0
        self.completed_batches = 0
        self.best_results: list[dict] = []
        self.best_results_by_case: dict[str, list[dict]] = {case_name: [] for case_name in self.cases}
        self.best_loss_history: list[dict[str, float | int]] = []
        self.pending_samples = deque()
        self.completed_samples = deque()
        self.batches: dict[int, dict] = {}
        self.case_jobs = deque()
        self.started_monotonic = time.monotonic()
        # A CLUBB/F2PY loss session owns global Fortran state and therefore
        # remains bound to one case.  We can still adapt the *pool* safely by
        # retiring only idle warm workers and spawning their replacement for a
        # backlogged case.  No in-flight evaluation is ever interrupted.
        self.case_evaluation_seconds: dict[str, float | None] = {
            case_name: None for case_name in self.cases
        }
        self.case_completed_evaluations: dict[str, int] = {
            case_name: 0 for case_name in self.cases
        }
        self.pending_replacements: deque[dict[str, str]] = deque()
        self.rebalance_interval_index = 0
        self.next_rebalance_monotonic = self.started_monotonic + REBALANCE_INTERVALS_SECONDS[0]
        self.last_rebalance: dict[str, object] | None = None
        self.strategy = None
        self.sample_history_writer: SampleHistoryWriter | None = None
        self.param_names: list[str] | None = None
        self.parameter_hard_bounds: dict[str, dict] | None = None
        self.default_params_row: list[float] | None = None
        self.baseline_case_results: dict[str, dict[str, dict]] = {
            "clubb_default": {},
            **({"override_defaults": {}} if self.override else {}),
        }
        self.baselines: dict[str, dict] = {}
        self.error_message = ""
        self.elapsed_before_resume_seconds = 0.0
        self.last_status_write = 0.0
        self.checkpoint_path = self.job_dir / "resume_checkpoint.pkl"

    def run(self) -> int:
        """Run the full tuning lifecycle and return a process exit code."""
        started_at = utc_now_iso()
        self._write_public_state("initializing", started_at, started_at, force=True)

        try:
            if not self._initialize_case_barrier(started_at):
                finished_at = utc_now_iso()
                write_status(
                    self.status_path,
                    state="stopped",
                    job_dir=self.job_dir,
                    samples_evaluated=self.samples_evaluated,
                    elapsed_seconds=self._elapsed_seconds(),
                    best_results=self.best_results,
                    metrics=self._metrics(),
                )
                write_results(
                    self.results_path,
                    state="stopped",
                    job_dir=self.job_dir,
                    request=self.request,
                    samples_evaluated=self.samples_evaluated,
                    best_results=self.best_results,
                    best_results_by_case=self.best_results_by_case,
                    baselines=self.baselines,
                    started_at=started_at,
                    updated_at=finished_at,
                    finished_at=finished_at,
                )
                return 0
            self._assert_compiled_parameter_compatibility()
            self._finalize_baselines()
            self._initialize_sample_history_writer()
            seed = self.request.get("seed")
            self.strategy = build_tuning_strategy(
                self.request["strategy"],
                param_names=self.param_names or [],
                default_params_row=self.default_params_row or [],
                parameter_ranges=self.request["parameter_ranges"],
                seed=None if seed is None else int(seed),
            )
            pending_limit = self.strategy.pending_batch_limit(self.batch_size)
            if pending_limit is not None:
                self.max_pending_batches = min(
                    self.max_pending_case_jobs_per_case,
                    max(1, int(pending_limit)),
                )
                self.max_pending_samples = self.max_pending_batches * self.batch_size
            self._restore_resume_checkpoint()
            self.request["total_samples"] = self.strategy.estimated_sample_count()
            self._write_public_state("running", started_at, utc_now_iso(), force=True)

            stopping = False
            while True:
                self._poll_workers()
                self._start_pending_replacements()
                self._drain_completed_samples()
                if self.sample_history_writer is not None:
                    self.sample_history_writer.flush_if_due()
                if self.error_message:
                    raise RuntimeError(self.error_message)

                if not stopping and self._should_stop():
                    stopping = True
                    self.case_jobs.clear()
                    self._stop_all_workers()

                if stopping:
                    if self._active_evaluations() == 0:
                        finished_at = utc_now_iso()
                        self._close_sample_history()
                        self._write_resume_checkpoint()
                        write_status(
                            self.status_path,
                            state="stopped",
                            job_dir=self.job_dir,
                            samples_evaluated=self.samples_evaluated,
                            elapsed_seconds=self._elapsed_seconds(),
                            best_results=self.best_results,
                            metrics=self._metrics(),
                        )
                        write_results(
                            self.results_path,
                            state="stopped",
                            job_dir=self.job_dir,
                            request=self.request,
                            samples_evaluated=self.samples_evaluated,
                            best_results=self.best_results,
                            best_results_by_case=self.best_results_by_case,
                            baselines=self.baselines,
                            started_at=started_at,
                            updated_at=finished_at,
                            finished_at=finished_at,
                        )
                        return 0
                    self._write_public_state("stopping", started_at, utc_now_iso())
                    time.sleep(POLL_INTERVAL_SECONDS)
                    continue

                self._fill_pending_samples()
                self._pack_pending_batches()
                self._rebalance_warm_workers_if_due()
                self._dispatch_jobs()

                if self._is_finished():
                    self._stop_all_workers()
                    finished_at = utc_now_iso()
                    self._close_sample_history()
                    self._remove_resume_checkpoint()
                    write_status(
                        self.status_path,
                        state="finished",
                        job_dir=self.job_dir,
                        samples_evaluated=self.samples_evaluated,
                        elapsed_seconds=self._elapsed_seconds(),
                        best_results=self.best_results,
                        metrics=self._metrics(),
                    )
                    write_results(
                        self.results_path,
                        state="finished",
                        job_dir=self.job_dir,
                        request=self.request,
                        samples_evaluated=self.samples_evaluated,
                        best_results=self.best_results,
                        best_results_by_case=self.best_results_by_case,
                        baselines=self.baselines,
                        started_at=started_at,
                        updated_at=finished_at,
                        finished_at=finished_at,
                    )
                    return 0

                self._write_public_state("running", started_at, utc_now_iso())
                time.sleep(POLL_INTERVAL_SECONDS)

        except Exception as exc:
            self.error_message = f"{exc}\n{traceback.format_exc(limit=10)}"
            self._stop_all_workers()
            self._close_sample_history(suppress_errors=True)
            self._remove_resume_checkpoint()
            finished_at = utc_now_iso()
            write_status(
                self.status_path,
                state="error",
                job_dir=self.job_dir,
                samples_evaluated=self.samples_evaluated,
                elapsed_seconds=self._elapsed_seconds(),
                best_results=self.best_results,
                error_message=self.error_message,
                metrics=self._metrics(),
            )
            write_results(
                self.results_path,
                state="error",
                job_dir=self.job_dir,
                request=self.request,
                samples_evaluated=self.samples_evaluated,
                best_results=self.best_results,
                best_results_by_case=self.best_results_by_case,
                baselines=self.baselines,
                started_at=started_at,
                updated_at=finished_at,
                finished_at=finished_at,
                error_message=self.error_message,
            )
            return 1

    def _assert_compiled_parameter_compatibility(self) -> None:
        """Fail clearly when the Python F2PY build predates the selected config.

        Worker metadata comes from ``clubb_api.get_param_names()``, i.e. the
        actual compiled API that will evaluate every candidate.  Comparing it
        with the request here keeps a namelist/API rename mismatch from being
        mistaken for a strategy error or allowing a partial Tune launch.
        """
        available = set(self.param_names or [])
        requested = {
            str(target).strip()
            for spec in self.request.get("parameter_ranges", [])
            for target in (
                spec.get("targets", [spec.get("name")])
                if not isinstance(spec.get("targets", [spec.get("name")]), str)
                else [spec.get("targets")]
            )
        }
        missing = sorted(name for name in requested if name and name not in available)
        if missing:
            raise RuntimeError(
                "The selected Tune configuration requests parameter(s) absent from "
                "the compiled CLUBB Python API: "
                + ", ".join(missing)
                + ". Rebuild the Python API from this checkout (./compile.py -python) "
                "before running Tune."
            )

        self._assert_requested_ranges_within_hard_bounds()

    def _assert_requested_ranges_within_hard_bounds(self) -> None:
        """Reject a Tune envelope outside the compiled Fortran contract.

        Request parsing validates syntax and checked-in configuration names.
        This check deliberately runs after worker initialization because the
        worker's F2PY module is the exact CLUBB build that will score samples.
        """
        if self.parameter_hard_bounds is None:
            raise RuntimeError(
                "The compiled CLUBB Python API did not return parameter hard bounds. "
                "Rebuild it from this checkout (./compile.py -python)."
            )

        violations: list[str] = []
        for spec in self.request.get("parameter_ranges", []):
            name = str(spec.get("name", "")).strip()
            requested_min = float(spec["min"])
            requested_max = float(spec["max"])
            raw_targets = spec.get("targets", [name])
            targets = [raw_targets] if isinstance(raw_targets, str) else list(raw_targets)
            for target in targets:
                target = str(target).strip()
                bound = self.parameter_hard_bounds.get(target)
                if bound is None:
                    violations.append(f"{target}: no compiled hard-bound metadata")
                    continue
                lower = bound.get("min")
                upper = bound.get("max")
                if lower is not None and requested_min < float(lower):
                    violations.append(
                        f"{target}: requested min {requested_min:g} must be >= {float(lower):g}"
                    )
                if upper is not None and requested_max > float(upper):
                    violations.append(
                        f"{target}: requested max {requested_max:g} must be <= {float(upper):g}"
                    )

        if violations:
            raise RuntimeError(
                "Requested Tune range(s) violate compiled CLUBB hard bounds: "
                + "; ".join(violations)
            )

    def _initialize_case_barrier(self, started_at: str) -> bool:
        """Start one worker per case and report whether initialization completed."""
        for case_name in self.cases:
            self._start_worker(case_name, evaluate_baselines=True)

        while True:
            self._poll_workers()
            if self.error_message:
                raise RuntimeError(self.error_message)
            if self._should_stop():
                self._stop_all_workers()
                return False
            initialized_cases = {
                worker.case_name
                for worker in self.workers
                if worker.state in {"idle", "busy"}
            }
            if all(case_name in initialized_cases for case_name in self.cases):
                return True
            self._write_public_state("initializing", started_at, utc_now_iso())
            time.sleep(POLL_INTERVAL_SECONDS)

    def _start_worker(self, case_name: str, *, evaluate_baselines: bool = False) -> WorkerHandle:
        """Spawn one worker process for a case."""
        worker_id = self.next_worker_id
        self.next_worker_id += 1
        worker_dir = self.job_dir / "workers" / f"{case_name}_{worker_id}"
        parent_conn, child_conn = self.ctx.Pipe()
        payload = {
            "backend": self.request.get("backend", "fortran"),
            "case_name": case_name,
            "selected_fields": self.selected_fields,
            "batch_size": self.batch_size,
            "worker_dir": str(worker_dir),
            "num_time_windows": self.case_num_time_windows.get(case_name, 1),
            "case_defaults": self.case_defaults.get(case_name, {}),
            "config": self.config,
            "override": self.override,
            "evaluate_baselines": evaluate_baselines,
            "grad": self.grad,
            "loss_weights": dict(self.loss_weights),
            "field_weights": dict(self.field_weights),
            "gradient_parameters": self.gradient_parameters,
            "gradient_mode": self.request.get('gradient_mode', 'auto'),
        }
        process = self.ctx.Process(target=worker_main, args=(child_conn, payload))
        process.start()
        child_conn.close()
        handle = WorkerHandle(
            worker_id=worker_id,
            case_name=case_name,
            process=process,
            conn=parent_conn,
            worker_dir=worker_dir,
        )
        self.workers.append(handle)
        return handle

    def _poll_workers(self) -> None:
        """Drain worker messages and surface failed processes."""
        for worker in list(self.workers):
            try:
                while worker.conn.poll():
                    message = worker.conn.recv()
                    self._handle_worker_message(worker, message)
            except EOFError:
                if worker.state not in {"stopping", "failed"}:
                    self.error_message = f"Worker for case {worker.case_name} exited unexpectedly"

            if (
                worker.state not in {"stopping", "failed"}
                and not worker.process.is_alive()
                and worker.process.exitcode not in (None, 0)
            ):
                self.error_message = (
                    f"Worker for case {worker.case_name} exited with code {worker.process.exitcode}"
                )
            if worker.state == "stopping" and not worker.process.is_alive():
                try:
                    worker.process.join(timeout=0)
                except Exception:
                    pass
                try:
                    worker.conn.close()
                except Exception:
                    pass
                if worker in self.workers:
                    self.workers.remove(worker)

    def _handle_worker_message(self, worker: WorkerHandle, message: dict) -> None:
        message_type = message.get("type")
        if message_type == "initialized":
            self._handle_initialized(worker, message)
            return
        if message_type == "result":
            self._record_worker_evaluation_time(worker)
            worker.current_batch_id = None
            worker.current_job_started_monotonic = None
            self._handle_result(message)
            if worker.retire_after_batch:
                worker.retire_after_batch = False
                self._stop_worker(worker)
            else:
                worker.state = "idle"
            return
        if message_type == "error":
            worker.state = "failed"
            self.error_message = message.get("error_message", "Worker failed")
            return
        self.error_message = f"Unknown worker response type: {message_type}"

    def _handle_initialized(self, worker: WorkerHandle, message: dict) -> None:
        field_names = list(message.get("field_names", []))
        if field_names != self.selected_fields:
            self.error_message = (
                f"Worker field order for {worker.case_name} did not match request: "
                f"{field_names} != {self.selected_fields}"
            )
            worker.state = "failed"
            return

        param_names = list(message.get("param_names", []))
        if self.param_names is None:
            self.param_names = param_names
        elif param_names != self.param_names:
            self.error_message = f"Worker parameter names for {worker.case_name} did not match canonical order"
            worker.state = "failed"
            return

        raw_bounds = list(message.get("hard_parameter_bounds", []))
        bounds_by_name = {
            str(item.get("name", "")).strip(): dict(item)
            for item in raw_bounds
            if isinstance(item, dict) and str(item.get("name", "")).strip()
        }
        if set(bounds_by_name) != set(param_names):
            self.error_message = (
                f"Worker for {worker.case_name} did not return hard-bound metadata "
                "for every compiled tunable parameter"
            )
            worker.state = "failed"
            return
        if self.parameter_hard_bounds is None:
            self.parameter_hard_bounds = bounds_by_name
        elif bounds_by_name != self.parameter_hard_bounds:
            self.error_message = (
                f"Worker hard-bound metadata for {worker.case_name} did not match canonical metadata"
            )
            worker.state = "failed"
            return

        default_params = message.get("default_params", [])
        if not default_params:
            self.error_message = f"Worker for {worker.case_name} did not return default parameters"
            worker.state = "failed"
            return
        if self.default_params_row is None:
            self.default_params_row = [float(value) for value in default_params[0]]

        for baseline_name, baseline_message in dict(message.get("baselines") or {}).items():
            if baseline_name not in self.baseline_case_results:
                self.error_message = f"Worker for {worker.case_name} returned unknown baseline {baseline_name}"
                worker.state = "failed"
                return
            if not baseline_message.get("available", True):
                self.baseline_case_results[baseline_name][worker.case_name] = {
                    "__baseline_unavailable__": str(
                        baseline_message.get("unavailable_reason")
                        or "CLUBB rejected this baseline parameter set"
                    )
                }
                continue
            case_metrics = self._case_result_from_message(baseline_message, worker.case_name)
            if case_metrics is None:
                worker.state = "failed"
                return
            self.baseline_case_results[baseline_name][worker.case_name] = case_metrics

        worker.state = "idle"

    def _initialize_sample_history_writer(self) -> None:
        """Create the raw sample-history writer once canonical axes are known."""
        if self.sample_history_writer is not None:
            return
        if self.param_names is None:
            raise RuntimeError("Cannot initialize sample history before parameter names are known")
        self.sample_history_writer = SampleHistoryWriter(
            job_dir=self.job_dir,
            param_names=self.param_names,
            case_names=self.cases,
            field_names=self.selected_fields,
            metric_names=list(LOSS_METRIC_NAMES),
            case_configs=list(self.request.get("case_configs", [])),
            case_window_counts=dict(self.case_num_time_windows),
            gradient_parameters=self.gradient_parameters if self.grad else [],
        )

    def _close_sample_history(self, *, suppress_errors: bool = False) -> None:
        """Flush buffered raw sample-history records."""
        if self.sample_history_writer is None:
            return
        try:
            self.sample_history_writer.close()
        except Exception:
            if not suppress_errors:
                raise

    def _handle_result(self, message: dict) -> None:
        batch_id = int(message["batch_id"])
        case_name = message["case_name"]
        batch = self.batches.get(batch_id)
        if batch is None:
            return

        case_metrics = self._case_result_from_message(message, case_name)
        if case_metrics is None:
            return

        if len(case_metrics['loss']) != len(batch['params_batch']):
            self.error_message = f'Result candidate count for {case_name} did not match request'
            return
        if self.grad:
            if message.get('gradient_parameters') != self.gradient_parameters:
                self.error_message = f'Gradient parameter order for {case_name} did not match request'
                return
            gradients = message.get('gradient', {})
            expected = (len(batch['params_batch']), len(self.gradient_parameters))
            if message.get('gradient_target') != 'loss' or set(gradients) != {'loss'} or any(
                np.asarray(values).shape != expected for values in gradients.values()
            ):
                self.error_message = f'Gradient shape for {case_name} did not match request'
                return
            batch.setdefault('case_gradients', {})[case_name] = gradients
        batch["case_results"][case_name] = case_metrics

        if all(case_name in batch["case_results"] for case_name in self.cases):
            self._complete_batch(batch_id, batch)

    def _case_result_from_message(self, message: dict, case_name: str) -> dict | None:
        """Validate one case result and return its canonical metric structure."""
        field_names = list(message.get("field_names", []))
        if field_names != self.selected_fields:
            self.error_message = f"Result field order for {case_name} did not match request"
            return None

        metric_names = list(message.get("loss_metric_names", []))
        if metric_names != list(LOSS_METRIC_NAMES):
            self.error_message = f"Result loss metric names for {case_name} did not match expected names"
            return None

        metric_arrays = message.get("loss_metrics_by_metric_window_field_and_column", {})
        if set(metric_arrays) != set(LOSS_METRIC_NAMES):
            self.error_message = f"Result loss metric keys for {case_name} did not match expected names"
            return None

        case_loss = np.asarray(message.get("case_loss", []), dtype=float)
        field_loss = np.asarray(message.get("field_loss", []), dtype=float)
        expected = (self.case_num_time_windows.get(case_name, 1), len(self.selected_fields), case_loss.size)
        if (case_loss.ndim != 1 or not case_loss.size
                or field_loss.shape != expected[1:]
                or any(np.asarray(values).shape != expected for values in metric_arrays.values())
                or not np.all(np.isfinite(case_loss)) or np.any(case_loss < 0)):
            self.error_message = f"Result loss shapes or values for {case_name} did not match request"
            return None
        return {"metrics": metric_arrays, "loss": case_loss.tolist(), "field_loss": field_loss.tolist()}

    def _complete_batch(self, batch_id: int, batch: dict) -> None:
        active_count = int(batch["active_sample_count"])
        entries = []
        for col_idx in range(active_count):
            field_loss = {}
            field_metrics = {}
            case_loss = {}
            case_loss_diagnostics = {}
            for case_name in self.cases:
                result = batch["case_results"][case_name]
                case_loss[case_name] = float(result["loss"][col_idx])
                case_loss_diagnostics[case_name] = {
                    "aggregation": "rms_fields_and_windows",
                    "loss_weights": dict(self.loss_weights),
                    "loss": case_loss[case_name],
                }
                field_loss[case_name] = {}
                field_metrics[case_name] = {}
                for field_idx, field_name in enumerate(self.selected_fields):
                    subwindows = [
                        {"window_index": window_idx + 1,
                         **{name: float(result["metrics"][name][window_idx][field_idx][col_idx])
                            for name in LOSS_METRIC_NAMES}}
                        for window_idx in range(self.case_num_time_windows[case_name])
                    ]
                    for window in subwindows:
                        window["loss"] = window["rmse"]
                    representative = max(subwindows, key=lambda item: item["rmse"])
                    value = float(result["field_loss"][field_idx][col_idx])
                    field_loss[case_name][field_name] = value
                    field_metrics[case_name][field_name] = {
                        **representative, "rmse": value, "loss": value,
                        "num_time_windows": len(subwindows),
                        "representative_window_index": representative["window_index"],
                        "subwindows": subwindows,
                    }
            total_loss_diagnostics = aggregate_losses(
                [case_loss[name] for name in self.cases],
                [self.case_weights.get(name, 1.0) for name in self.cases],
                self.aggregation_mode, self.aggregation_weights,
            )
            total_loss = float(total_loss_diagnostics["loss"])

            entry = {
                "sample_id": int(batch["sample_ids"][col_idx]),
                "batch_id": int(batch_id),
                "total_loss": float(total_loss),
                "loss_weights": dict(self.loss_weights),
                "aggregation_mode": self.aggregation_mode,
                "case_window_counts": dict(self.case_num_time_windows),
                "config": self.config,
                "case_configs": list(self.request.get("case_configs", [])),
                "aggregation_weights": list(self.aggregation_weights),
                "loss_policy_version": self.loss_policy_version,
                "loss_policy_constants": dict(self.loss_policy_constants),
                "aggregation_options": dict(self.aggregation_options),
                "selected_params": dict(batch["selected_params_by_sample"][col_idx]),
                "all_params": dict(batch["all_params_by_sample"][col_idx]),
                "case_loss": case_loss,
                "case_loss_diagnostics": case_loss_diagnostics,
                "total_loss_diagnostics": total_loss_diagnostics,
                "field_loss": field_loss,
                "field_metrics": field_metrics,
            }
            if not batch.get("baseline_name"):
                entry.update(self._improvement_scores(total_loss))
            if self.grad and not batch.get('baseline_name'):
                entry['gradient_parameters'] = self.gradient_parameters
                from tuner.gradients import combine_case_gradients
                entry['total_gradient'] = combine_case_gradients(
                    total_loss_diagnostics['loss_coefficients'],
                    [np.asarray(batch['case_gradients'][case]['loss'])[col_idx]
                     for case in self.cases],
                ).tolist()
                entry['gradient_target'] = 'loss'
                entry['gradient'] = {
                    case: {'loss': np.asarray(batch['case_gradients'][case]['loss'])[col_idx].tolist()}
                    for case in self.cases}
            entries.append(entry)

        if batch.get("baseline_name"):
            baseline = dict(entries[0])
            for key in ("sample_id", "batch_id", "selected_params", "all_params"):
                baseline.pop(key, None)
            self.baselines[str(batch["baseline_name"])] = baseline
        else:
            self.completed_samples.extend(entries)
            self.completed_batches += 1
            self.batches.pop(batch_id, None)

    def _finalize_baselines(self) -> None:
        """Aggregate the one-column initialization baselines like candidates."""
        for baseline_name, case_metrics in self.baseline_case_results.items():
            missing = [case_name for case_name in self.cases if case_name not in case_metrics]
            if missing:
                raise RuntimeError(
                    f"Baseline {baseline_name} was not evaluated for: {', '.join(missing)}"
                )
            unavailable_cases = {
                case_name: str(metrics["__baseline_unavailable__"])
                for case_name, metrics in case_metrics.items()
                if "__baseline_unavailable__" in metrics
            }
            if unavailable_cases:
                self.baselines[baseline_name] = {
                    "available": False,
                    "total_loss": None,
                    "case_loss": {},
                    "unavailable_cases": unavailable_cases,
                }
                continue
            self._complete_batch(
                -1,
                {
                    "baseline_name": baseline_name,
                    "active_sample_count": 1,
                    "sample_ids": [-1],
                    "selected_params_by_sample": [{}],
                    "all_params_by_sample": [{}],
                    "case_results": case_metrics,
                },
            )
            self.baselines[baseline_name]["available"] = True

    def _improvement_scores(self, candidate_loss: float) -> dict[str, float | None]:
        """Return signed percent improvement relative to available baselines."""
        scores = {}
        candidate_loss = float(candidate_loss)
        candidate_available = (
            math.isfinite(candidate_loss)
            and candidate_loss < float(LOSS_POLICY_CONSTANTS["invalid_rmse_penalty"])
        )
        for baseline_name, output_name in (
            ("clubb_default", "improvement_vs_clubb_default_percent"),
            ("override_defaults", "improvement_vs_override_defaults_percent"),
        ):
            baseline = self.baselines.get(baseline_name)
            if baseline is None:
                continue
            baseline_value = baseline.get("total_loss")
            baseline_loss = float("nan") if baseline_value is None else float(baseline_value)
            scores[output_name] = (
                100.0 * (baseline_loss - candidate_loss) / baseline_loss
                if candidate_available
                and baseline.get("available", True)
                and math.isfinite(baseline_loss)
                and baseline_loss > 0.0
                else None
            )
        return scores

    def _drain_completed_samples(self) -> list[dict]:
        completed = []
        while self.completed_samples:
            entry = self.completed_samples.popleft()
            previous_best = float(self.best_results[0]["total_loss"]) if self.best_results else float("inf")
            self.best_results = self._update_best_results(
                self.best_results,
                entry,
                key=lambda item: float(item["total_loss"]),
            )
            entry_loss = float(entry["total_loss"])
            if entry_loss < previous_best:
                self.best_loss_history.append(
                    {"sample_count": self.samples_evaluated + len(completed) + 1, "loss": entry_loss}
                )
                self.best_loss_history = self.best_loss_history[-BEST_LOSS_HISTORY_LIMIT:]
            for case_name in self.cases:
                self.best_results_by_case[case_name] = self._update_best_results(
                    self.best_results_by_case.get(case_name, []),
                    entry,
                    key=lambda item, name=case_name: float(item.get("case_loss", {}).get(name, item["total_loss"])),
                )
            completed.append(entry)

        if completed:
            if self.sample_history_writer is not None:
                self.sample_history_writer.append(completed)
            self.samples_evaluated += len(completed)
            if self.strategy is not None:
                self.strategy.tell(completed)
        return completed

    def _fill_pending_samples(self) -> None:
        if self.strategy is None:
            return
        packed_samples = sum(int(batch["active_sample_count"]) for batch in self.batches.values())
        capacity = max(0, self.max_pending_samples - packed_samples)
        self.strategy.fill(self.pending_samples, capacity)

    def _pack_pending_batches(self) -> None:
        if self.default_params_row is None:
            return
        while len(self.batches) < self.max_pending_batches:
            if not self.pending_samples:
                return

            samples = []
            while self.pending_samples and len(samples) < self.batch_size:
                samples.append(self.pending_samples.popleft())

            batch_id = self.next_batch_id
            self.next_batch_id += 1
            active_count = len(samples)
            params_batch = [list(sample["param_row"]) for sample in samples]
            for _ in range(self.batch_size - active_count):
                params_batch.append(list(self.default_params_row))

            self.batches[batch_id] = {
                "batch_id": batch_id,
                "active_sample_count": active_count,
                "sample_ids": [sample["sample_id"] for sample in samples],
                "params_batch": params_batch,
                "selected_params_by_sample": [sample["selected_params"] for sample in samples],
                "all_params_by_sample": [sample["all_params"] for sample in samples],
                "case_results": {},
            }
            for case_name in self.cases:
                self.case_jobs.append({"batch_id": batch_id, "case_name": case_name})

    def _dispatch_jobs(self) -> None:
        while self.case_jobs and self._active_evaluations() < self.max_workers:
            made_progress = False
            for _ in range(len(self.case_jobs)):
                if self._active_evaluations() >= self.max_workers:
                    break

                job = self.case_jobs.popleft()
                if job["batch_id"] not in self.batches:
                    made_progress = True
                    continue

                worker = self._idle_worker(job["case_name"])
                if worker is not None:
                    self._dispatch_to_worker(worker, job)
                    made_progress = True
                    continue

                if self._can_start_worker(job["case_name"]):
                    self._start_worker(job["case_name"])
                    made_progress = True

                self.case_jobs.append(job)

            if not made_progress:
                break

    def _dispatch_to_worker(self, worker: WorkerHandle, job: dict) -> None:
        batch = self.batches[job["batch_id"]]
        worker.conn.send(
            {
                "type": "evaluate_batch",
                "batch_id": int(job["batch_id"]),
                "params_batch": batch["params_batch"],
            }
        )
        worker.state = "busy"
        worker.current_batch_id = int(job["batch_id"])
        worker.current_job_started_monotonic = time.monotonic()

    def _idle_worker(self, case_name: str) -> WorkerHandle | None:
        for worker in self.workers:
            if worker.case_name == case_name and worker.state == "idle":
                return worker
        return None

    def _can_start_worker(self, case_name: str) -> bool:
        if len([worker for worker in self.workers if worker.state != "failed"]) >= self.worker_cap:
            return False
        return not any(worker.case_name == case_name and worker.state == "starting" for worker in self.workers)

    def _record_worker_evaluation_time(self, worker: WorkerHandle) -> None:
        """Update the case's stable batch-duration estimate after a result."""
        started = worker.current_job_started_monotonic
        if started is None:
            return
        duration = max(0.0, time.monotonic() - started)
        case_name = worker.case_name
        previous = self.case_evaluation_seconds.get(case_name)
        if previous is None:
            self.case_evaluation_seconds[case_name] = duration
        else:
            self.case_evaluation_seconds[case_name] = (
                EVALUATION_TIME_EWMA_ALPHA * duration
                + (1.0 - EVALUATION_TIME_EWMA_ALPHA) * previous
            )
        self.case_completed_evaluations[case_name] = (
            self.case_completed_evaluations.get(case_name, 0) + 1
        )

    def _queued_jobs_by_case(self) -> dict[str, int]:
        """Return outstanding valid case jobs, grouped by initialized case."""
        counts = {case_name: 0 for case_name in self.cases}
        for job in self.case_jobs:
            if job.get("batch_id") in self.batches and job.get("case_name") in counts:
                counts[job["case_name"]] += 1
        return counts

    def _warm_workers_by_case(self) -> dict[str, int]:
        """Count usable warm workers; stopping/failed actors are excluded."""
        counts = {case_name: 0 for case_name in self.cases}
        for worker in self.workers:
            if worker.case_name in counts and worker.state in {"idle", "busy", "starting"}:
                counts[worker.case_name] += 1
        return counts

    def _estimated_case_drain_seconds(
        self,
        case_name: str,
        *,
        queued_jobs: int,
        warm_workers: int,
    ) -> float | None:
        """Estimate queue drain time from observed batch timings.

        A case with no observed completion remains eligible for work, but is
        not used to make a timing-driven reassignment until at least one
        result has established its cost.  This avoids guessing during the
        expensive F2PY/Fortran initialization phase.
        """
        duration = self.case_evaluation_seconds.get(case_name)
        if queued_jobs <= 0:
            return 0.0
        if duration is None or warm_workers <= 0:
            return None
        return float(queued_jobs) * float(duration) / float(warm_workers)

    def _rebalance_warm_workers_if_due(self) -> None:
        """Move one warm actor toward the case with the longest queue.

        Workers cannot change cases in-process because ``init_clubb_loss``
        owns process-global CLUBB state.  An idle donor stops immediately; a
        busy donor is marked to retire after its current batch.  Both paths
        are followed by one replacement process for the target case.  We keep
        one warm worker for each active case and make at most one move at a
        progressively slower 5/10/15-second reassessment cadence.
        """
        now = time.monotonic()
        if now < self.next_rebalance_monotonic:
            return
        self.rebalance_interval_index = min(
            self.rebalance_interval_index + 1,
            len(REBALANCE_INTERVALS_SECONDS) - 1,
        )
        interval = REBALANCE_INTERVALS_SECONDS[self.rebalance_interval_index]
        self.next_rebalance_monotonic = now + interval

        if self.pending_replacements:
            return
        queued = self._queued_jobs_by_case()
        warm = self._warm_workers_by_case()
        drain = {
            case_name: self._estimated_case_drain_seconds(
                case_name,
                queued_jobs=queued[case_name],
                warm_workers=warm[case_name],
            )
            for case_name in self.cases
        }
        eligible_targets = [
            case_name for case_name in self.cases
            if queued[case_name] > 0 and drain[case_name] is not None
        ]
        if not eligible_targets:
            return
        target = max(eligible_targets, key=lambda case_name: float(drain[case_name] or 0.0))
        target_drain = float(drain[target] or 0.0)

        # Normally actor churn is not worthwhile when every measured queue
        # drains in under five seconds.  The one useful exception is an
        # already-idle actor: it has no in-flight work to protect and can be
        # repurposed immediately for the case with the largest short drain.
        # Make only one such transfer per check. Every case retains one warm
        # actor so later batches can never strand work in a workerless queue.
        all_drains_short = (
            all(value is not None for value in drain.values())
            and max((float(value or 0.0) for value in drain.values()), default=0.0)
            < REBALANCE_MIN_DRAIN_SECONDS
        )
        if all_drains_short:
            idle_donors = [
                worker for worker in self.workers
                if worker.state == "idle"
                and worker.case_name != target
                and not getattr(worker, "retire_after_batch", False)
            ]
            if not idle_donors:
                return
            donor = min(
                idle_donors,
                key=lambda worker: (
                    float(drain[worker.case_name] or 0.0),
                    worker.worker_id,
                ),
            )
            self._schedule_replacement(
                donor=donor,
                target=target,
                target_drain=target_drain,
                donor_drain=drain[donor.case_name],
                now=now,
                reason="idle_short_drain",
            )
            return

        donor_candidates = [
            worker for worker in self.workers
            if worker.state in {"idle", "busy"}
            and worker.case_name != target
            and warm.get(worker.case_name, 0) > REBALANCE_MIN_WARM_WORKERS_PER_CASE
            and not getattr(worker, "retire_after_batch", False)
        ]
        if not donor_candidates:
            return

        # Prefer an idle donor when one exists.  Otherwise select the worker
        # with the shortest projected drain and let it finish its one active
        # batch before retiring; this is event-driven, not another 5-second
        # reassessment wait.
        donor = min(
            donor_candidates,
            key=lambda worker: (
                0 if worker.state == "idle" else 1,
                float(drain[worker.case_name])
                if drain[worker.case_name] is not None else float("inf"),
                worker.worker_id,
            ),
        )
        donor_drain = drain[donor.case_name]
        if donor_drain is not None and target_drain < REBALANCE_IMBALANCE_RATIO * float(donor_drain):
            return

        self._schedule_replacement(
            donor=donor,
            target=target,
            target_drain=target_drain,
            donor_drain=donor_drain,
            now=now,
            reason="drain_imbalance",
        )

    def _schedule_replacement(
        self,
        *,
        donor: WorkerHandle,
        target: str,
        target_drain: float,
        donor_drain: float | None,
        now: float,
        reason: str,
    ) -> None:
        """Retire one donor and queue its replacement for ``target``."""
        warm_donors = sum(
            worker.case_name == donor.case_name
            and worker.state in {"idle", "busy", "starting"}
            for worker in self.workers
        )
        if warm_donors <= REBALANCE_MIN_WARM_WORKERS_PER_CASE:
            return

        self.pending_replacements.append({"from_case": donor.case_name, "to_case": target})
        self.last_rebalance = {
            "from_case": donor.case_name,
            "to_case": target,
            "target_drain_seconds": target_drain,
            "donor_drain_seconds": donor_drain,
            "reason": reason,
            "at_elapsed_seconds": self._elapsed_seconds(),
        }
        if donor.state == "idle":
            self._stop_worker(donor)
        else:
            donor.retire_after_batch = True
        # A successful move changes the queue geometry.  Reassess quickly so
        # a still-backlogged case can receive another actor without waiting
        # for the no-op backoff cadence.  Only successive checks that make no
        # move grow from 5 to 10 and finally 15 seconds.
        self.rebalance_interval_index = 0
        self.next_rebalance_monotonic = now + REBALANCE_INTERVALS_SECONDS[0]

    def _start_pending_replacements(self) -> None:
        """Start rebalanced actors only after their donor has fully exited."""
        while self.pending_replacements and len(self.workers) < self.worker_cap:
            replacement = self.pending_replacements.popleft()
            self._start_worker(str(replacement["to_case"]))

    def _stop_all_workers(self) -> None:
        for worker in list(self.workers):
            self._stop_worker(worker)
        for worker in list(self.workers):
            try:
                worker.process.join(timeout=5)
            except Exception:
                pass
            if worker.process.is_alive():
                try:
                    worker.process.terminate()
                    worker.process.join(timeout=2)
                except Exception:
                    pass

    def _stop_worker(self, worker: WorkerHandle) -> None:
        if worker.state == "stopping":
            return
        try:
            worker.conn.send({"type": "stop"})
        except Exception:
            pass
        worker.state = "stopping"

    def _should_stop(self) -> bool:
        return should_stop(self.control_path)

    def _is_finished(self) -> bool:
        if self.strategy is None or not self.strategy.is_exhausted():
            return False
        return (
            not self.pending_samples
            and not self.batches
            and not self.case_jobs
            and self._active_evaluations() == 0
        )

    def _write_resume_checkpoint(self) -> None:
        """Save the stopped search so Continue can resume its proposals and optimizer state.

        Write algorithm moments/RNG, pending candidates, partial case batches, ranking
        and elapsed progress to a temporary pickle file, then atomically replace
        resume_checkpoint.pkl. This preserves the search trajectory, not just results.
        """
        if self.strategy is None:
            return
        payload = {
            "version": 1,
            "strategy": self.strategy,
            "pending_samples": list(self.pending_samples),
            "batches": self.batches,
            "baselines": self.baselines,
            "samples_evaluated": self.samples_evaluated,
            "completed_batches": self.completed_batches,
            "next_batch_id": self.next_batch_id,
            "best_results": self.best_results,
            "best_results_by_case": self.best_results_by_case,
            "best_loss_history": self.best_loss_history,
            "elapsed_seconds": self._elapsed_seconds(),
        }
        descriptor, temp_name = tempfile.mkstemp(
            prefix=f".{self.checkpoint_path.name}.", suffix=".tmp", dir=self.job_dir
        )
        try:
            with os.fdopen(descriptor, "wb") as handle:
                pickle.dump(payload, handle, protocol=pickle.HIGHEST_PROTOCOL)
            os.replace(temp_name, self.checkpoint_path)
        finally:
            try:
                Path(temp_name).unlink()
            except FileNotFoundError:
                pass

    def _restore_resume_checkpoint(self) -> None:
        """Load saved current-format search state when continuing a stopped execution.

        Restore the algorithm and pending work, then rebuild case jobs only for case
        results still missing from partial batches. With no saved file there is
        nothing to restore; an unreadable/unsupported file raises a clear error.
        """
        if not self.checkpoint_path.is_file():
            return
        try:
            with self.checkpoint_path.open("rb") as handle:
                payload = pickle.load(handle)
            if int(payload.get("version", 0)) != 1 or payload.get("strategy") is None:
                raise ValueError("unsupported checkpoint")
            self.strategy = payload["strategy"]
            self.pending_samples = deque(payload.get("pending_samples", []))
            self.batches = dict(payload.get("batches", {}))
            self.baselines = dict(payload.get("baselines", self.baselines))
            self.case_jobs = deque(
                {"batch_id": int(batch_id), "case_name": case_name}
                for batch_id, batch in self.batches.items()
                for case_name in self.cases
                if case_name not in batch.get("case_results", {})
            )
            self.samples_evaluated = int(payload.get("samples_evaluated", 0))
            self.completed_batches = int(payload.get("completed_batches", 0))
            self.next_batch_id = int(payload.get("next_batch_id", 0))
            self.best_results = list(payload.get("best_results", []))
            self.best_results_by_case = {
                str(name): list(values)
                for name, values in dict(payload.get("best_results_by_case", {})).items()
            }
            for case_name in self.cases:
                self.best_results_by_case.setdefault(case_name, [])
            self.best_loss_history = [
                {"sample_count": int(item["sample_count"]), "loss": float(item["loss"])}
                for item in list(payload.get("best_loss_history", []))[-BEST_LOSS_HISTORY_LIMIT:]
                if isinstance(item, dict) and "sample_count" in item and "loss" in item
            ]
            self.elapsed_before_resume_seconds = float(payload.get("elapsed_seconds", 0.0))
        except Exception as exc:
            raise RuntimeError(f"Could not restore stopped tuning revision: {exc}") from exc

    def _remove_resume_checkpoint(self) -> None:
        try:
            self.checkpoint_path.unlink()
        except FileNotFoundError:
            pass

    def _active_evaluations(self) -> int:
        return sum(1 for worker in self.workers if worker.state == "busy")

    def _metrics(self) -> dict:
        queued = self._queued_jobs_by_case()
        warm = self._warm_workers_by_case()
        case_workers = {}
        for case_name in self.cases:
            active = sum(
                1 for worker in self.workers
                if worker.case_name == case_name and worker.state == "busy"
            )
            idle = sum(
                1 for worker in self.workers
                if worker.case_name == case_name and worker.state == "idle"
            )
            estimated = self._estimated_case_drain_seconds(
                case_name,
                queued_jobs=queued[case_name],
                warm_workers=warm[case_name],
            )
            case_workers[case_name] = {
                "warm_workers": warm[case_name],
                "active_workers": active,
                "idle_workers": idle,
                "queued_jobs": queued[case_name],
                "mean_evaluation_seconds": self.case_evaluation_seconds[case_name],
                "estimated_drain_seconds": estimated,
                "completed_evaluations": self.case_completed_evaluations[case_name],
            }
        return {
            "active_evaluations": self._active_evaluations(),
            "idle_workers": sum(1 for worker in self.workers if worker.state == "idle"),
            "initialized_workers": sum(1 for worker in self.workers if worker.state in {"idle", "busy"}),
            "queued_case_jobs": len(self.case_jobs),
            "case_worker_metrics": case_workers,
            "worker_rebalance": {
                "next_check_in_seconds": max(0.0, self.next_rebalance_monotonic - time.monotonic()),
                "next_interval_seconds": REBALANCE_INTERVALS_SECONDS[
                    min(self.rebalance_interval_index, len(REBALANCE_INTERVALS_SECONDS) - 1)
                ],
                "pending_replacements": len(self.pending_replacements),
                "last_move": dict(self.last_rebalance or {}),
            },
            "completed_batches": self.completed_batches,
            "best_loss_history": list(self.best_loss_history),
            "total_samples": None if self.strategy is None else self.strategy.estimated_sample_count(),
            "selected_fields": list(self.selected_fields),
            "loss_weights": dict(self.loss_weights),
            "aggregation_mode": self.aggregation_mode,
            "case_window_counts": dict(self.case_num_time_windows),
            "case_configs": list(self.request.get("case_configs", [])),
            "config": self.config,
            "aggregation_weights": list(self.aggregation_weights),
            "loss_policy_version": self.loss_policy_version,
            "loss_policy_constants": dict(self.loss_policy_constants),
            "aggregation_options": dict(self.aggregation_options),
            "control_stop_reason": stop_reason(self.control_path),
            "baselines": {
                name: {
                    "available": baseline.get("available", True),
                    "total_loss": baseline.get("total_loss"),
                    "case_loss": dict(baseline.get("case_loss", {})),
                    "unavailable_cases": dict(baseline.get("unavailable_cases", {})),
                }
                for name, baseline in self.baselines.items()
            },
        }

    def _elapsed_seconds(self) -> float:
        return self.elapsed_before_resume_seconds + time.monotonic() - self.started_monotonic

    def _write_public_state(self, state: str, started_at: str, updated_at: str, *, force: bool = False) -> None:
        now = time.monotonic()
        if not force and now - self.last_status_write < STATUS_HEARTBEAT_SECONDS:
            return
        self.last_status_write = now
        write_status(
            self.status_path,
            state=state,
            job_dir=self.job_dir,
            samples_evaluated=self.samples_evaluated,
            elapsed_seconds=self._elapsed_seconds(),
            best_results=self.best_results,
            metrics=self._metrics(),
        )
        write_results(
            self.results_path,
            state=state,
            job_dir=self.job_dir,
            request=self.request,
            samples_evaluated=self.samples_evaluated,
            best_results=self.best_results,
            best_results_by_case=self.best_results_by_case,
            baselines=self.baselines,
            started_at=started_at,
            updated_at=updated_at,
        )

    def _update_best_results(self, best_results: list[dict], entry: dict, *, key) -> list[dict]:
        updated = list(best_results)
        updated.append(entry)
        updated.sort(key=key)
        return updated[: self.results_file_limit]


def run_scheduler(
    request: dict,
    *,
    job_dir: Path,
    control_path: Path,
    status_path: Path,
    results_path: Path,
) -> int:
    """Run one validated experiment and return its scheduler exit code.

    The job directory and control/status/results paths belong to the caller.
    TuningScheduler owns workers, feedback, stop/continuation and file updates.
    """
    return TuningScheduler(
        request=request,
        job_dir=job_dir,
        control_path=control_path,
        status_path=status_path,
        results_path=results_path,
    ).run()
