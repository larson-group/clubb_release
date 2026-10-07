#!/usr/bin/env python3
"""Standalone tuner entrypoint that communicates through a job directory."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys
import traceback

if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    # Choose the runtime before importing scheduler dependencies. Direct
    # job-directory calls use the same managed JAX launcher as TunerJob.
    bootstrap_parser = argparse.ArgumentParser(add_help=False)
    bootstrap_parser.add_argument("-job_dir")
    bootstrap_args, _ = bootstrap_parser.parse_known_args()
    try:
        raw_request = (
            json.loads((Path(bootstrap_args.job_dir) / "request.json").read_text())
            if bootstrap_args.job_dir else {}
        )
    except (OSError, ValueError):
        raw_request = {}  # Normal error reporting below retains the job results.
    if not isinstance(raw_request, dict):
        raw_request = {}  # Let load_request report the invalid saved JSON below.
    if str(raw_request.get("backend", "")).strip().lower() == "jax":
        if os.environ.get("_CLUBB_JAX_ENVIRONMENT_PYTHON") != sys.executable:
            from tuner.job_runtime import tuner_worker_env
            from clubb_jax.run_jax import runtime_arguments
            from tuner.status import write_job_error

            try:
                os.execve(sys.executable, [
                    sys.executable,
                    str(Path(__file__).resolve().parents[1] / "clubb_jax" / "run_jax.py"),
                    *runtime_arguments(
                        str(raw_request.get("jax_options", "cpu")),
                        device=raw_request.get("jax_gpu", ""),
                        prealloc_gpu_mem=raw_request.get("jax_xla_prealloc"),
                    ),
                    "-module=tuner.tune_clubb", *sys.argv[1:],
                ], tuner_worker_env())
            except Exception as exc:
                job_dir = Path(bootstrap_args.job_dir)
                write_job_error(job_dir / "status.json", job_dir / "results.json", str(exc))
                print(f"ERROR: {exc}", file=sys.stderr)
                raise SystemExit(1) from exc
    else:
        from utilities.setup_python_venv import ensure_python_venv

        ensure_python_venv()

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
    parser = argparse.ArgumentParser(description="Run one CLUBB tuning job from a job directory.", add_help=False, allow_abbrev=False)
    parser.add_argument("-h", "-help", action="help", help="Show this help and exit.")
    parser.add_argument('-job_dir', dest='job_dir', required=True, help="Job directory containing request/control/status/results files.")
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
