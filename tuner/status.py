"""Status, control, and results-file helpers for tuning jobs."""

from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path


STATUS_RESULT_LIMIT = 16
RESULTS_FILE_LIMIT = 16
DEFAULT_KEEPALIVE_TIMEOUT_SECONDS = 300
DEFAULT_KEEPALIVE_ACTION = "stop"


def utc_now_iso() -> str:
    """Return the current UTC timestamp in ISO-8601 form."""
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def _utc_from_iso(value: str) -> datetime | None:
    """Parse the UTC timestamp format written by utc_now_iso."""
    try:
        return datetime.fromisoformat(str(value).replace("Z", "+00:00"))
    except (TypeError, ValueError):
        return None


def atomic_write_json(path: Path, payload: dict) -> None:
    """Write JSON atomically so readers never observe a partial file."""
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    tmp_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    tmp_path.replace(path)


def read_json_or_default(path: Path, default: dict) -> dict:
    """Read a JSON file when present, otherwise return the provided default."""
    if not path.is_file():
        return dict(default)
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return dict(default)


def read_control(control_path: Path) -> dict:
    """Read the control file, defaulting to stop_requested=false."""
    if not control_path.is_file():
        return {"stop_requested": False}
    try:
        payload = json.loads(control_path.read_text(encoding="utf-8"))
    except Exception:
        return {"stop_requested": False}
    if not isinstance(payload, dict):
        return {"stop_requested": False}
    payload["stop_requested"] = bool(payload.get("stop_requested", False))
    payload["keepalive_required"] = bool(payload.get("keepalive_required", False))
    if payload.get("keepalive_timeout_seconds") is not None:
        try:
            payload["keepalive_timeout_seconds"] = float(payload["keepalive_timeout_seconds"])
        except (TypeError, ValueError):
            payload["keepalive_timeout_seconds"] = DEFAULT_KEEPALIVE_TIMEOUT_SECONDS
    return payload


def write_control(
    control_path: Path,
    *,
    stop_requested: bool | None = None,
    keepalive_required: bool | None = None,
    keepalive_timeout_seconds: float | None = None,
    keepalive_updated_at: str | None = None,
    keepalive_action: str | None = None,
) -> None:
    """Update the canonical control file while preserving unrelated keys."""
    payload = read_control(control_path)
    if stop_requested is not None:
        payload["stop_requested"] = bool(stop_requested)
    if keepalive_required is not None:
        payload["keepalive_required"] = bool(keepalive_required)
    if keepalive_timeout_seconds is not None:
        payload["keepalive_timeout_seconds"] = float(keepalive_timeout_seconds)
    if keepalive_updated_at is not None:
        payload["keepalive_updated_at"] = str(keepalive_updated_at)
    if keepalive_action is not None:
        payload["keepalive_action"] = str(keepalive_action)
    payload.setdefault("stop_requested", False)
    atomic_write_json(control_path, payload)


def renew_keepalive(
    control_path: Path,
    *,
    timeout_seconds: float = DEFAULT_KEEPALIVE_TIMEOUT_SECONDS,
    action: str = DEFAULT_KEEPALIVE_ACTION,
) -> None:
    """Renew the controller lease for a running tuner job."""
    write_control(
        control_path,
        keepalive_required=True,
        keepalive_timeout_seconds=timeout_seconds,
        keepalive_updated_at=utc_now_iso(),
        keepalive_action=action,
    )


def keepalive_expired(control: dict, *, now: datetime | None = None) -> bool:
    """Return whether a required controller lease has expired."""
    if not control.get("keepalive_required", False):
        return False
    updated_at = _utc_from_iso(control.get("keepalive_updated_at"))
    if updated_at is None:
        return True
    timeout_seconds = float(control.get("keepalive_timeout_seconds", DEFAULT_KEEPALIVE_TIMEOUT_SECONDS))
    now = now or datetime.now(timezone.utc)
    return (now - updated_at).total_seconds() > timeout_seconds


def stop_reason(control_path: Path) -> str | None:
    """Return the current graceful-stop reason, if any."""
    control = read_control(control_path)
    if control.get("stop_requested", False):
        return "stop_requested"
    if control.get("keepalive_action", DEFAULT_KEEPALIVE_ACTION) == "stop" and keepalive_expired(control):
        return "keepalive_expired"
    return None


def should_stop(control_path: Path) -> bool:
    """Return whether the control file requests a graceful stop."""
    return stop_reason(control_path) is not None


def summarize_for_status(best_results: list[dict]) -> list[dict]:
    """Build the lightweight top-results summary for Dash polling."""
    summary = []
    for idx, result in enumerate(best_results[:STATUS_RESULT_LIMIT], start=1):
        entry = {
            "rank": idx,
            "total_loss": float(result["total_loss"]),
            "params": dict(result["selected_params"]),
            "loss_mode": result.get("loss_mode"),
            "aggregation_mode": result.get("aggregation_mode"),
            "case_window_counts": dict(result.get("case_window_counts", {})),
            "time_window_aggregation_mode": result.get("time_window_aggregation_mode"),
            "loss_policy_version": result.get("loss_policy_version"),
        }
        if "scaled_rmse_sum" in result:
            entry["scaled_rmse_sum"] = float(result["scaled_rmse_sum"])
        elif "simple_rms_sum" in result:
            entry["scaled_rmse_sum"] = float(result["simple_rms_sum"])
        if "case_loss" in result:
            entry["case_loss"] = dict(result["case_loss"])
        summary.append(entry)
    return summary


def write_status(
    status_path: Path,
    *,
    state: str,
    job_dir: Path,
    samples_evaluated: int = 0,
    elapsed_seconds: float = 0.0,
    best_results: list[dict] | None = None,
    error_message: str = "",
    metrics: dict | None = None,
) -> None:
    """Write the lightweight live status file used by Dash polling."""
    best_results = best_results or []
    metrics = metrics or {}
    payload = {
        "state": state,
        "job_dir": str(job_dir),
        "samples_evaluated": int(samples_evaluated),
        "total_samples": None if metrics.get("total_samples") is None else int(metrics.get("total_samples")),
        "elapsed_seconds": float(round(elapsed_seconds, 3)),
        "best_total_loss": None if not best_results else float(best_results[0]["total_loss"]),
        "top_results": summarize_for_status(best_results),
        "error_message": error_message,
        "active_evaluations": int(metrics.get("active_evaluations", 0)),
        "idle_workers": int(metrics.get("idle_workers", 0)),
        "initialized_workers": int(metrics.get("initialized_workers", 0)),
        "queued_case_jobs": int(metrics.get("queued_case_jobs", 0)),
        "completed_batches": int(metrics.get("completed_batches", 0)),
        "loss_mode": metrics.get("loss_mode"),
        "aggregation_mode": metrics.get("aggregation_mode"),
        "case_window_counts": dict(metrics.get("case_window_counts", {})),
        "case_configs": list(metrics.get("case_configs", [])),
        "time_window_aggregation_mode": metrics.get("time_window_aggregation_mode"),
        "loss_policy_version": metrics.get("loss_policy_version"),
        "loss_policy_constants": dict(metrics.get("loss_policy_constants", {})),
        "aggregation_options": dict(metrics.get("aggregation_options", {})),
        "control_stop_reason": metrics.get("control_stop_reason"),
    }
    atomic_write_json(status_path, payload)


def write_results(
    results_path: Path,
    *,
    state: str,
    job_dir: Path,
    request: dict | None,
    samples_evaluated: int,
    best_results: list[dict],
    best_results_by_case: dict[str, list[dict]] | None = None,
    started_at: str,
    updated_at: str,
    finished_at: str | None = None,
    error_message: str = "",
) -> None:
    """Write the retained results file with full parameter vectors."""
    payload = {
        "state": state,
        "job_dir": str(job_dir),
        "case_name": None if request is None else request.get("case_name"),
        "cases": [] if request is None else list(request.get("cases", [])),
        "selected_fields": [] if request is None else list(request.get("selected_fields", [])),
        "case_configs": [] if request is None else list(request.get("case_configs", [])),
        "case_window_counts": {} if request is None else dict(request.get("case_window_counts", {})),
        "loss_mode": None if request is None else request.get("loss_mode"),
        "aggregation_mode": None if request is None else request.get("aggregation_mode"),
        "time_window_aggregation_mode": None if request is None else request.get("time_window_aggregation_mode"),
        "loss_policy_version": None if request is None else request.get("loss_policy_version"),
        "loss_policy_constants": {} if request is None else dict(request.get("loss_policy_constants", {})),
        "aggregation_options": {} if request is None else dict(request.get("aggregation_options", {})),
        "samples_evaluated": int(samples_evaluated),
        "total_samples": None if request is None else request.get("total_samples"),
        "started_at": started_at,
        "updated_at": updated_at,
        "finished_at": finished_at,
        "error_message": error_message,
        "request": request,
        "best_results": [],
        "best_results_by_case": {},
    }
    for idx, result in enumerate(best_results[:RESULTS_FILE_LIMIT], start=1):
        payload["best_results"].append(
            {
                "rank": idx,
                "total_loss": float(result["total_loss"]),
                "loss_mode": result.get("loss_mode"),
                "aggregation_mode": result.get("aggregation_mode"),
                "case_window_counts": dict(result.get("case_window_counts", {})),
                "case_configs": list(result.get("case_configs", [])),
                "time_window_aggregation_mode": result.get("time_window_aggregation_mode"),
                "loss_policy_version": result.get("loss_policy_version"),
                "loss_policy_constants": dict(result.get("loss_policy_constants", {})),
                "aggregation_options": dict(result.get("aggregation_options", {})),
                "scaled_rmse_sum": float(result.get("scaled_rmse_sum", result.get("simple_rms_sum", result["total_loss"]))),
                "selected_params": dict(result["selected_params"]),
                "all_params": dict(result["all_params"]),
                "sample_id": result.get("sample_id"),
                "batch_id": result.get("batch_id"),
                "case_loss": dict(result.get("case_loss", {})),
                "case_loss_diagnostics": {
                    case_name: dict(diagnostics)
                    for case_name, diagnostics in result.get("case_loss_diagnostics", {}).items()
                },
                "total_loss_diagnostics": dict(result.get("total_loss_diagnostics", {})),
                "field_loss": {
                    case_name: dict(case_fields)
                    for case_name, case_fields in result.get("field_loss", {}).items()
                },
                "scaled_rmse_case_sum": dict(
                    result.get("scaled_rmse_case_sum", result.get("simple_rms_case_sum", {}))
                ),
                "scaled_rmse_by_field": {
                    case_name: dict(case_fields)
                    for case_name, case_fields in result.get(
                        "scaled_rmse_by_field",
                        result.get("simple_rms_by_field", {}),
                    ).items()
                },
                "field_metrics": {
                    case_name: {
                        field_name: dict(metrics)
                        for field_name, metrics in case_fields.items()
                    }
                    for case_name, case_fields in result.get("field_metrics", {}).items()
                },
            }
        )
    for case_name, case_results in (best_results_by_case or {}).items():
        payload["best_results_by_case"][case_name] = [
            _lightweight_case_result(case_name, result, idx)
            for idx, result in enumerate(case_results[:RESULTS_FILE_LIMIT], start=1)
        ]
    atomic_write_json(results_path, payload)


def _lightweight_case_result(case_name: str, result: dict, rank: int) -> dict:
    """Return a compact per-case top-result row for display-oriented uses."""
    raw_case_loss = result.get("case_loss", {})
    case_loss_value = raw_case_loss.get(case_name, result["total_loss"]) if isinstance(raw_case_loss, dict) else raw_case_loss
    raw_scaled_rmse_case_sum = result.get("scaled_rmse_case_sum", result.get("simple_rms_case_sum", {}))
    scaled_rmse_case_value = (
        raw_scaled_rmse_case_sum.get(case_name, 0.0)
        if isinstance(raw_scaled_rmse_case_sum, dict)
        else raw_scaled_rmse_case_sum
    )
    return {
        "rank": int(rank),
        "case_name": str(case_name),
        "sample_id": result.get("sample_id"),
        "batch_id": result.get("batch_id"),
        "total_loss": float(result["total_loss"]),
        "case_loss": float(case_loss_value),
        "scaled_rmse_sum": float(result.get("scaled_rmse_sum", result.get("simple_rms_sum", result["total_loss"]))),
        "scaled_rmse_case_sum": float(scaled_rmse_case_value),
        "selected_params": dict(result.get("selected_params", {})),
        "all_params": dict(result.get("all_params", {})),
        "loss_mode": result.get("loss_mode"),
        "aggregation_mode": result.get("aggregation_mode"),
        "case_window_counts": dict(result.get("case_window_counts", {})),
        "time_window_aggregation_mode": result.get("time_window_aggregation_mode"),
        "loss_policy_version": result.get("loss_policy_version"),
    }
