"""Shared constants and runtime state for the tuning tab."""

from __future__ import annotations

import os
import threading


REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
OUTPUT_TUNER_DIR = os.path.join(REPO_ROOT, "output_tuner")
TUNE_STATUS_TEMPLATE = {
    "state": "idle",
    "job_dir": None,
    "samples_evaluated": 0,
    "total_samples": None,
    "elapsed_seconds": 0.0,
    "best_total_loss": None,
    "top_results": [],
    "error_message": "",
    "active_evaluations": 0,
    "idle_workers": 0,
    "initialized_workers": 0,
    "queued_case_jobs": 0,
    "completed_batches": 0,
    "case_configs": [],
    "case_window_counts": {},
}
TUNE_ACTIVE = {}
TUNE_LOCK = threading.Lock()
LOSS_RUN_PROCS = {}
LOSS_RUN_LOCK = threading.Lock()
