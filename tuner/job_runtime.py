"""Shared runtime wrapper for file-backed tuner jobs."""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
import os
from pathlib import Path
import subprocess
import sys

from tuner.paths import REPO_ROOT
from tuner.status import (
    DEFAULT_KEEPALIVE_ACTION,
    DEFAULT_KEEPALIVE_TIMEOUT_SECONDS,
    atomic_write_json,
    read_json_or_default,
    renew_keepalive,
    write_control,
)


TERMINAL_STATES = {"finished", "stopped", "error"}


def tuner_worker_env() -> dict:
    """Return the environment used by tuner subprocesses."""
    run_env = os.environ.copy()
    existing_pythonpath = run_env.get("PYTHONPATH", "")
    pythonpath_entries = [str(REPO_ROOT), str(REPO_ROOT / "clubb_python_api")]
    if existing_pythonpath:
        pythonpath_entries.append(existing_pythonpath)
    run_env["PYTHONPATH"] = os.pathsep.join(pythonpath_entries)
    return run_env


def _case_text_from_request(request: dict) -> str:
    cases = request.get("cases") or [request.get("case_name", "tune")]
    case_text = "_".join(str(case).strip() for case in cases if str(case).strip())
    return case_text or "tune"


def _make_default_status(job_dir: Path, state: str) -> dict:
    return {
        "state": state,
        "job_dir": str(job_dir),
        "samples_evaluated": 0,
        "total_samples": None,
        "elapsed_seconds": 0.0,
        "best_total_loss": None,
        "top_results": [],
        "error_message": "",
    }


def _make_default_results(job_dir: Path, request: dict, state: str) -> dict:
    return {
        "state": state,
        "job_dir": str(job_dir),
        "request": request,
        "samples_evaluated": 0,
        "best_results": [],
        "best_results_by_case": {},
    }


@dataclass
class TunerJob:
    """Thin object wrapper around the tuner job-directory contract."""

    job_dir: Path
    request_path: Path
    control_path: Path
    status_path: Path
    results_path: Path
    log_path: Path
    pid: int | None = None
    proc: subprocess.Popen | None = field(default=None, repr=False)

    @classmethod
    def from_dir(
        cls,
        job_dir: str | Path,
        *,
        pid: int | None = None,
        proc: subprocess.Popen | None = None,
    ) -> "TunerJob":
        job_dir = Path(job_dir).resolve()
        return cls(
            job_dir=job_dir,
            request_path=job_dir / "request.json",
            control_path=job_dir / "control.json",
            status_path=job_dir / "status.json",
            results_path=job_dir / "results.json",
            log_path=job_dir / "worker.log",
            pid=pid if pid is not None else (None if proc is None else proc.pid),
            proc=proc,
        )

    @classmethod
    def from_dict(
        cls,
        payload: dict,
        *,
        proc: subprocess.Popen | None = None,
    ) -> "TunerJob":
        job_dir = payload.get("job_dir") or payload.get("work_dir")
        if not job_dir:
            raise RuntimeError("Tuner job payload is missing job_dir/work_dir")
        job = cls.from_dir(job_dir, pid=payload.get("pid"), proc=proc)
        for attr, key in (
            ("request_path", "request_path"),
            ("control_path", "control_path"),
            ("status_path", "status_path"),
            ("results_path", "results_path"),
            ("log_path", "log_path"),
        ):
            if payload.get(key):
                setattr(job, attr, Path(payload[key]).resolve())
        return job

    @classmethod
    def create(
        cls,
        request: dict,
        *,
        output_root: str | Path = REPO_ROOT / "output_tuner",
        job_dir: str | Path | None = None,
        prefix: str = "",
        initial_status: dict | None = None,
        initial_state: str = "created",
        keepalive_required: bool = True,
        keepalive_timeout_seconds: float = DEFAULT_KEEPALIVE_TIMEOUT_SECONDS,
        keepalive_action: str = DEFAULT_KEEPALIVE_ACTION,
    ) -> "TunerJob":
        if job_dir is not None:
            job_path = Path(job_dir).resolve()
            if job_path.exists() and any(job_path.iterdir()):
                raise RuntimeError(f"Job directory already exists and is not empty: {job_path}")
            job_path.mkdir(parents=True, exist_ok=True)
        else:
            output_root = Path(output_root).resolve()
            output_root.mkdir(parents=True, exist_ok=True)
            case_text = _case_text_from_request(request)
            name_parts = [part for part in (prefix, case_text) if part]
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
            job_name = "_".join(name_parts + [timestamp, str(os.getpid())])
            job_path = output_root / job_name
            job_path.mkdir(parents=True, exist_ok=False)

        job = cls.from_dir(job_path)
        atomic_write_json(job.request_path, request)
        write_control(
            job.control_path,
            stop_requested=False,
            keepalive_required=keepalive_required,
            keepalive_timeout_seconds=keepalive_timeout_seconds,
            keepalive_action=keepalive_action,
        )
        if keepalive_required:
            job.heartbeat(
                timeout_seconds=keepalive_timeout_seconds,
                action=keepalive_action,
            )

        status_payload = _make_default_status(job.job_dir, initial_state)
        status_payload.update(initial_status or {})
        status_payload["job_dir"] = str(job.job_dir)
        atomic_write_json(job.status_path, status_payload)
        atomic_write_json(job.results_path, _make_default_results(job.job_dir, request, initial_state))
        return job

    def start(self) -> subprocess.Popen:
        """Launch ``python -m tuner.tune_clubb`` for this job."""
        if self.proc is not None:
            raise RuntimeError("This tuner job already owns a subprocess")
        log_handle = open(self.log_path, "w", encoding="utf-8")
        try:
            self.proc = subprocess.Popen(
                [sys.executable, "-m", "tuner.tune_clubb", "--job-dir", str(self.job_dir)],
                cwd=str(REPO_ROOT),
                env=tuner_worker_env(),
                stdout=log_handle,
                stderr=subprocess.STDOUT,
                text=True,
            )
        finally:
            log_handle.close()
        self.pid = self.proc.pid
        return self.proc

    def poll(self) -> int | None:
        if self.proc is None:
            return None
        return self.proc.poll()

    def request_stop(self) -> None:
        """Request graceful stop through the file-backed control API."""
        write_control(self.control_path, stop_requested=True)

    def heartbeat(
        self,
        *,
        timeout_seconds: float = DEFAULT_KEEPALIVE_TIMEOUT_SECONDS,
        action: str = DEFAULT_KEEPALIVE_ACTION,
    ) -> None:
        """Renew this job's controller lease."""
        renew_keepalive(
            self.control_path,
            timeout_seconds=timeout_seconds,
            action=action,
        )

    def terminate(self, *, timeout: float = 5.0) -> int | None:
        """Terminate a live process owned by this Python process."""
        if self.proc is None:
            return None
        self.proc.terminate()
        try:
            return self.proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            self.proc.kill()
            return self.proc.wait(timeout=timeout)

    def status(self, default: dict | None = None) -> dict:
        return read_json_or_default(self.status_path, default or {})

    def results(self, default: dict | None = None) -> dict:
        return read_json_or_default(self.results_path, default or {})

    def is_process_alive(self) -> bool:
        if self.proc is not None:
            return self.proc.poll() is None
        if not self.pid:
            return False
        try:
            os.kill(int(self.pid), 0)
        except ProcessLookupError:
            return False
        except PermissionError:
            return True
        except OSError:
            return False
        return True

    def has_exited(self) -> bool:
        return not self.is_process_alive()

    def to_dict(self) -> dict:
        pid = self.pid if self.pid is not None else (None if self.proc is None else self.proc.pid)
        return {
            "pid": pid,
            "job_dir": str(self.job_dir),
            "work_dir": str(self.job_dir),
            "request_path": str(self.request_path),
            "control_path": str(self.control_path),
            "status_path": str(self.status_path),
            "results_path": str(self.results_path),
            "log_path": str(self.log_path),
        }
