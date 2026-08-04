"""Shared constants and runtime state for the tuning tab."""

from __future__ import annotations

import os
import threading


REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
# New Tune lineages are stored under output/tuner.  ``output_tuner`` remains a
# read-only legacy location for old jobs; no new worker writes there.
OUTPUT_TUNER_DIR = os.path.join(REPO_ROOT, "output", "tuner")
TUNE_STATUS_TEMPLATE = {
    "state": "idle",
    "job_dir": None,
    "samples_evaluated": 0,
    "total_samples": None,
    "elapsed_seconds": 0.0,
    "best_total_loss": None,
    "top_results": [],
    "error_message": "",
    "config": "default",
    "active_evaluations": 0,
    "idle_workers": 0,
    "initialized_workers": 0,
    "queued_case_jobs": 0,
    "case_worker_metrics": {},
    "worker_rebalance": {},
    "completed_batches": 0,
    "best_loss_history": [],
    "selected_fields": [],
    "case_configs": [],
    "case_window_counts": {},
}
TUNE_ACTIVE = {}
TUNE_LOCK = threading.Lock()
LOSS_RUN_PROCS = {}
LOSS_RUN_LOCK = threading.Lock()
