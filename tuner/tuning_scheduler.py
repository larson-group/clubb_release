"""Master-side scheduler for multi-case CLUBB tuning."""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass
import multiprocessing
from pathlib import Path
import time
import traceback

from tuner.tuning_strategy import build_tuning_strategy
from tuner.tuning_worker import worker_main
from tuner.taylor_metrics import (
    AGGREGATION_MODE_NAMES,
    DEFAULT_AGGREGATION_MODE,
    DEFAULT_LOSS_MODE,
    LOSS_METRIC_NAMES,
    LOSS_MODE_NAMES,
    LOSS_POLICY_CONSTANTS,
    LOSS_POLICY_VERSION,
    WORST_QUANTILE_FRACTION,
    aggregate_losses,
    compute_field_loss_diagnostics,
)
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
        self.selected_fields = list(request["selected_fields"])
        self.batch_size = int(request["batch_size"])
        self.max_workers = int(request["max_workers"])
        self.case_weights = dict(request.get("case_weights", {}))
        self.field_weights = dict(request.get("field_weights", {}))
        self.loss_mode = request.get("loss_mode", DEFAULT_LOSS_MODE)
        self.aggregation_mode = request.get("aggregation_mode", DEFAULT_AGGREGATION_MODE)
        self.time_window_aggregation_mode = request.get("time_window_aggregation_mode", self.aggregation_mode)
        self.case_num_time_windows = {
            case_name: int(
                self.case_configs.get(case_name, {}).get(
                    "num_time_windows",
                    self.case_defaults.get(case_name, {}).get("num_time_windows", 1),
                )
            )
            for case_name in self.cases
        }
        if self.loss_mode not in LOSS_MODE_NAMES:
            raise RuntimeError(f"Unknown loss mode: {self.loss_mode}")
        if self.aggregation_mode not in AGGREGATION_MODE_NAMES:
            raise RuntimeError(f"Unknown aggregation mode: {self.aggregation_mode}")
        if self.time_window_aggregation_mode not in AGGREGATION_MODE_NAMES:
            raise RuntimeError(f"Unknown time window aggregation mode: {self.time_window_aggregation_mode}")
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
        }
        self.request["loss_mode"] = self.loss_mode
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
        self.request.pop("time_window_mode", None)
        self.request.pop("num_time_windows", None)
        self.request["time_window_aggregation_mode"] = self.time_window_aggregation_mode
        self.request["loss_policy_version"] = self.loss_policy_version
        self.request["loss_policy_constants"] = dict(self.loss_policy_constants)
        self.request["aggregation_options"] = dict(self.aggregation_options)
        self.max_pending_batches = max(1, 2 * self.max_workers)
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
        self.pending_samples = deque()
        self.completed_samples = deque()
        self.batches: dict[int, dict] = {}
        self.case_jobs = deque()
        self.strategy = None
        self.param_names: list[str] | None = None
        self.default_params_row: list[float] | None = None
        self.error_message = ""
        self.started_monotonic = time.monotonic()
        self.last_status_write = 0.0

    def run(self) -> int:
        """Run the full tuning lifecycle and return a process exit code."""
        started_at = utc_now_iso()
        self._write_public_state("initializing", started_at, started_at, force=True)

        try:
            self._initialize_case_barrier(started_at)
            seed = self.request.get("seed")
            self.strategy = build_tuning_strategy(
                self.request["strategy"],
                param_names=self.param_names or [],
                default_params_row=self.default_params_row or [],
                parameter_ranges=self.request["parameter_ranges"],
                seed=None if seed is None else int(seed),
            )
            self.request["total_samples"] = self.strategy.estimated_sample_count()
            self._write_public_state("running", started_at, utc_now_iso(), force=True)

            stopping = False
            while True:
                self._poll_workers()
                self._drain_completed_samples()
                if self.error_message:
                    raise RuntimeError(self.error_message)

                if not stopping and self._should_stop():
                    stopping = True
                    self.pending_samples.clear()
                    self.case_jobs.clear()
                    self._stop_all_workers()

                if stopping:
                    if self._active_evaluations() == 0:
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
                self._dispatch_jobs()

                if self._is_finished():
                    self._stop_all_workers()
                    finished_at = utc_now_iso()
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
                started_at=started_at,
                updated_at=finished_at,
                finished_at=finished_at,
                error_message=self.error_message,
            )
            return 1

    def _initialize_case_barrier(self, started_at: str) -> None:
        """Start one worker per case and wait until metadata are validated."""
        for case_name in self.cases:
            self._start_worker(case_name)

        while True:
            self._poll_workers()
            if self.error_message:
                raise RuntimeError(self.error_message)
            if self._should_stop():
                self._stop_all_workers()
                raise RuntimeError("Tuning stopped during worker initialization")
            initialized_cases = {
                worker.case_name
                for worker in self.workers
                if worker.state in {"idle", "busy"}
            }
            if all(case_name in initialized_cases for case_name in self.cases):
                return
            self._write_public_state("initializing", started_at, utc_now_iso())
            time.sleep(POLL_INTERVAL_SECONDS)

    def _start_worker(self, case_name: str) -> WorkerHandle:
        """Spawn one worker process for a case."""
        worker_id = self.next_worker_id
        self.next_worker_id += 1
        worker_dir = self.job_dir / "workers" / f"{case_name}_{worker_id}"
        parent_conn, child_conn = self.ctx.Pipe()
        payload = {
            "case_name": case_name,
            "selected_fields": self.selected_fields,
            "batch_size": self.batch_size,
            "worker_dir": str(worker_dir),
            "num_time_windows": self.case_num_time_windows.get(case_name, 1),
            "case_defaults": self.case_defaults.get(case_name, {}),
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
            worker.state = "idle"
            worker.current_batch_id = None
            self._handle_result(message)
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

        default_params = message.get("default_params", [])
        if not default_params:
            self.error_message = f"Worker for {worker.case_name} did not return default parameters"
            worker.state = "failed"
            return
        if self.default_params_row is None:
            self.default_params_row = [float(value) for value in default_params[0]]

        worker.state = "idle"

    def _handle_result(self, message: dict) -> None:
        batch_id = int(message["batch_id"])
        case_name = message["case_name"]
        batch = self.batches.get(batch_id)
        if batch is None:
            return

        field_names = list(message.get("field_names", []))
        if field_names != self.selected_fields:
            self.error_message = f"Result field order for {case_name} did not match request"
            return

        metric_names = list(message.get("loss_metric_names", []))
        if metric_names != list(LOSS_METRIC_NAMES):
            self.error_message = f"Result loss metric names for {case_name} did not match expected names"
            return

        metric_arrays = message.get("loss_metrics_by_metric_window_field_and_column", {})
        if set(metric_arrays) != set(LOSS_METRIC_NAMES):
            self.error_message = f"Result loss metric keys for {case_name} did not match expected names"
            return

        expected_window_count = self.case_num_time_windows.get(case_name, 1)
        for metric_name in LOSS_METRIC_NAMES:
            metric_matrix = metric_arrays.get(metric_name, [])
            if len(metric_matrix) != expected_window_count:
                self.error_message = f"Result {metric_name} shape for {case_name} did not match time-window count"
                return
            if any(len(window_matrix) != len(self.selected_fields) for window_matrix in metric_matrix):
                self.error_message = f"Result {metric_name} shape for {case_name} did not match selected fields"
                return

        batch["case_loss_metrics"][case_name] = {
            field_name: {
                metric_name: [
                    [float(value) for value in metric_arrays[metric_name][window_idx][field_idx]]
                    for window_idx in range(expected_window_count)
                ]
                for metric_name in LOSS_METRIC_NAMES
            }
            for field_idx, field_name in enumerate(self.selected_fields)
        }

        if all(case_name in batch["case_loss_metrics"] for case_name in self.cases):
            self._complete_batch(batch_id, batch)

    def _complete_batch(self, batch_id: int, batch: dict) -> None:
        active_count = int(batch["active_sample_count"])
        for col_idx in range(active_count):
            field_loss = {}
            scaled_rmse_by_field = {}
            field_metrics = {}
            case_loss = {}
            case_loss_diagnostics = {}
            scaled_rmse_case_sum = {}
            smart_losses = []
            smart_weights = []
            scaled_rmse_sum = 0.0
            for case_name in self.cases:
                case_window_count = self.case_num_time_windows.get(case_name, 1)
                case_smart_losses = []
                case_smart_weights = []
                field_metrics[case_name] = {}
                field_loss[case_name] = {}
                scaled_rmse_by_field[case_name] = {}
                for field_name in self.selected_fields:
                    field_metric_values = batch["case_loss_metrics"][case_name][field_name]
                    subwindows = []
                    subwindow_losses = []
                    subwindow_scaled_rmse_values = []
                    for window_idx in range(case_window_count):
                        subwindow_metrics = {
                            "window_index": window_idx + 1,
                        }
                        for metric_name in LOSS_METRIC_NAMES:
                            subwindow_metrics[metric_name] = float(field_metric_values[metric_name][window_idx][col_idx])
                        subwindow_diagnostics = compute_field_loss_diagnostics(subwindow_metrics)
                        subwindow_metrics.update(subwindow_diagnostics)
                        subwindow_smart_loss = float(subwindow_metrics["per_field_losses"][self.loss_mode])
                        subwindow_metrics["loss"] = subwindow_smart_loss
                        subwindow_metrics["smart_loss"] = subwindow_smart_loss
                        subwindow_metrics["loss_mode"] = self.loss_mode
                        subwindow_metrics["loss_policy_version"] = self.loss_policy_version
                        subwindows.append(subwindow_metrics)
                        subwindow_losses.append(subwindow_smart_loss)
                        subwindow_scaled_rmse_values.append(float(subwindow_metrics["scaled_rmse"]))

                    subwindow_aggregation = aggregate_losses(
                        subwindow_losses,
                        None,
                        self.time_window_aggregation_mode,
                    )
                    scaled_rmse_aggregation = aggregate_losses(
                        subwindow_scaled_rmse_values,
                        None,
                        self.time_window_aggregation_mode,
                    )
                    scaled_rmse = float(scaled_rmse_aggregation["loss"])
                    smart_loss = float(subwindow_aggregation["loss"])
                    worst_window = max(subwindows, key=lambda item: float(item["loss"]))
                    metrics = dict(worst_window)
                    metrics.update(
                        {
                            "scaled_rmse": scaled_rmse,
                            "scaled_rmse_aggregation": scaled_rmse_aggregation,
                            "loss": smart_loss,
                            "smart_loss": smart_loss,
                            "num_time_windows": case_window_count,
                            "time_window_aggregation_mode": self.time_window_aggregation_mode,
                            "time_window_aggregation": subwindow_aggregation,
                            "representative_window_index": int(worst_window["window_index"]),
                            "subwindows": subwindows,
                        }
                    )
                    metrics["loss_mode"] = self.loss_mode
                    metrics["loss_policy_version"] = self.loss_policy_version
                    field_loss[case_name][field_name] = smart_loss
                    field_metrics[case_name][field_name] = metrics
                    scaled_rmse_by_field[case_name][field_name] = scaled_rmse
                case_weight = float(self.case_weights.get(case_name, 1.0))
                for field_name, field_scaled_rmse in scaled_rmse_by_field[case_name].items():
                    field_weight = float(self.field_weights.get(field_name, 1.0))
                    combined_weight = case_weight * field_weight
                    scaled_rmse_sum += combined_weight * field_scaled_rmse
                    smart_loss = field_loss[case_name][field_name]
                    case_smart_losses.append(smart_loss)
                    case_smart_weights.append(field_weight)
                    smart_losses.append(smart_loss)
                    smart_weights.append(combined_weight)
                case_aggregation = aggregate_losses(
                    case_smart_losses,
                    case_smart_weights,
                    self.aggregation_mode,
                )
                case_loss[case_name] = float(case_aggregation["loss"])
                case_loss_diagnostics[case_name] = case_aggregation
                scaled_rmse_case_sum[case_name] = float(sum(scaled_rmse_by_field[case_name].values()))

            total_loss_diagnostics = aggregate_losses(
                smart_losses,
                smart_weights,
                self.aggregation_mode,
            )
            total_loss = float(total_loss_diagnostics["loss"])

            entry = {
                "sample_id": int(batch["sample_ids"][col_idx]),
                "batch_id": int(batch_id),
                "total_loss": float(total_loss),
                "loss_mode": self.loss_mode,
                "aggregation_mode": self.aggregation_mode,
                "case_window_counts": dict(self.case_num_time_windows),
                "case_configs": list(self.request.get("case_configs", [])),
                "time_window_aggregation_mode": self.time_window_aggregation_mode,
                "loss_policy_version": self.loss_policy_version,
                "loss_policy_constants": dict(self.loss_policy_constants),
                "aggregation_options": dict(self.aggregation_options),
                "scaled_rmse_sum": float(scaled_rmse_sum),
                "selected_params": dict(batch["selected_params_by_sample"][col_idx]),
                "all_params": dict(batch["all_params_by_sample"][col_idx]),
                "case_loss": case_loss,
                "case_loss_diagnostics": case_loss_diagnostics,
                "total_loss_diagnostics": total_loss_diagnostics,
                "field_loss": field_loss,
                "scaled_rmse_case_sum": scaled_rmse_case_sum,
                "scaled_rmse_by_field": scaled_rmse_by_field,
                "field_metrics": field_metrics,
            }
            self.completed_samples.append(entry)

        self.completed_batches += 1
        self.batches.pop(batch_id, None)

    def _drain_completed_samples(self) -> list[dict]:
        completed = []
        while self.completed_samples:
            entry = self.completed_samples.popleft()
            self.best_results = self._update_best_results(
                self.best_results,
                entry,
                key=lambda item: float(item["total_loss"]),
            )
            for case_name in self.cases:
                self.best_results_by_case[case_name] = self._update_best_results(
                    self.best_results_by_case.get(case_name, []),
                    entry,
                    key=lambda item, name=case_name: float(item.get("case_loss", {}).get(name, item["total_loss"])),
                )
            completed.append(entry)

        if completed:
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
                "case_loss_metrics": {},
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

    def _idle_worker(self, case_name: str) -> WorkerHandle | None:
        for worker in self.workers:
            if worker.case_name == case_name and worker.state == "idle":
                return worker
        return None

    def _can_start_worker(self, case_name: str) -> bool:
        if len([worker for worker in self.workers if worker.state != "failed"]) >= self.worker_cap:
            return False
        return not any(worker.case_name == case_name and worker.state == "starting" for worker in self.workers)

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

    def _active_evaluations(self) -> int:
        return sum(1 for worker in self.workers if worker.state == "busy")

    def _metrics(self) -> dict:
        return {
            "active_evaluations": self._active_evaluations(),
            "idle_workers": sum(1 for worker in self.workers if worker.state == "idle"),
            "initialized_workers": sum(1 for worker in self.workers if worker.state in {"idle", "busy"}),
            "queued_case_jobs": len(self.case_jobs),
            "completed_batches": self.completed_batches,
            "total_samples": None if self.strategy is None else self.strategy.estimated_sample_count(),
            "loss_mode": self.loss_mode,
            "aggregation_mode": self.aggregation_mode,
            "case_window_counts": dict(self.case_num_time_windows),
            "case_configs": list(self.request.get("case_configs", [])),
            "time_window_aggregation_mode": self.time_window_aggregation_mode,
            "loss_policy_version": self.loss_policy_version,
            "loss_policy_constants": dict(self.loss_policy_constants),
            "aggregation_options": dict(self.aggregation_options),
            "control_stop_reason": stop_reason(self.control_path),
        }

    def _elapsed_seconds(self) -> float:
        return time.monotonic() - self.started_monotonic

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
    """Run one master scheduler instance."""
    return TuningScheduler(
        request=request,
        job_dir=job_dir,
        control_path=control_path,
        status_path=status_path,
        results_path=results_path,
    ).run()
