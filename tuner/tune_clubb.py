#!/usr/bin/env python3
"""Standalone tuner entrypoint that communicates through a job directory."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
import traceback

from tuner.request import load_request
from tuner.status import (
    read_json_or_default,
    utc_now_iso,
    write_control,
    write_results,
    write_status,
)
from tuner.tuning_scheduler import run_scheduler


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse the standalone tuner CLI."""
    parser = argparse.ArgumentParser(description="Run one CLUBB tuning job from a job directory.")
    parser.add_argument("--job-dir", required=True, help="Job directory containing request/control/status/results files.")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """CLI entrypoint for the standalone tuner."""
    args = parse_args(argv)
    job_dir = Path(args.job_dir).resolve()
    job_dir.mkdir(parents=True, exist_ok=True)

    request_path = job_dir / "request.json"
    control_path = job_dir / "control.json"
    status_path = job_dir / "status.json"
    results_path = job_dir / "results.json"

    if not control_path.exists():
        write_control(control_path, stop_requested=False)

    request = None
    try:
        request = load_request(request_path)
        return run_scheduler(
            request,
            job_dir=job_dir,
            control_path=control_path,
            status_path=status_path,
            results_path=results_path,
        )
    except Exception as exc:
        error_message = f"{exc}\n{traceback.format_exc(limit=10)}"
        finished_at = utc_now_iso()
        existing_status = read_json_or_default(status_path, {})
        existing_results = read_json_or_default(results_path, {})
        preserved_best_results = existing_results.get("best_results", [])
        preserved_best_results_by_case = existing_results.get("best_results_by_case", {})
        write_status(
            status_path,
            state="error",
            job_dir=job_dir,
            samples_evaluated=int(existing_status.get("samples_evaluated", 0)),
            elapsed_seconds=float(existing_status.get("elapsed_seconds", 0.0)),
            best_results=preserved_best_results,
            error_message=error_message,
        )
        write_results(
            results_path,
            state="error",
            job_dir=job_dir,
            request=request or existing_results.get("request"),
            samples_evaluated=int(existing_results.get("samples_evaluated", 0)),
            best_results=preserved_best_results,
            best_results_by_case=preserved_best_results_by_case,
            started_at=existing_results.get("started_at", finished_at),
            updated_at=finished_at,
            finished_at=finished_at,
            error_message=error_message,
        )
        return 1


if __name__ == "__main__":
    sys.exit(main())
