#!/usr/bin/env python3
"""Foreground supervisor for one checkout's broker and Dash processes."""

from __future__ import annotations

import fcntl
import json
import os
import signal
import subprocess
import sys
import threading
import time
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

import psutil

from dash_app.shared import dashboard_registry
from dash_app.shared.activity import ACTIVE_JOB_STATES, broker_jobs, count_active_jobs
from dash_app.shared.broker import connection_is_live, ensure_broker, reap_started_brokers, stop_broker
from dash_app.shared.gateway import read_connection
from dash_app.shared.manager_lease import (
    HEARTBEAT_INTERVAL_SECONDS,
    MANAGER_LOCK_PATH,
    MANAGER_REQUIRED_ENV,
    clear_manager_lease,
    manager_lease_is_live,
    process_started_at,
    read_manager_lease,
    write_manager_lease,
)
from dash_app.shared.provenance import runtime_source_fingerprint


REPO_ROOT = Path(__file__).resolve().parents[1]
DASH_APP = REPO_ROOT / "dash_app" / "app.py"
BROWSER_OPENED_ENV = "CLUBB_DASH_BROWSER_OPENED"
RESTART_INTERVAL_SECONDS = 10.0
RESTART_DEADLINE_SECONDS = 5.0 * 60.0
STARTUP_HEALTH_TIMEOUT_SECONDS = 45.0
RUNNING_HEALTH_TIMEOUT_SECONDS = 15.0
BROKER_RETRY_SECONDS = 10.0
JOB_STOP_TIMEOUT_SECONDS = 12.0
POLL_SECONDS = 0.5


def _active_job_labels(jobs: dict[str, Any] | None = None) -> list[str]:
    snapshot = broker_jobs() if jobs is None else jobs
    labels: list[str] = []
    compile_job = snapshot.get("compile") or {}
    if str(compile_job.get("state") or "") in ACTIVE_JOB_STATES:
        labels.append("compile")
    tune_job = snapshot.get("tune") or {}
    if str(tune_job.get("state") or "") in ACTIVE_JOB_STATES:
        labels.append("Tune")
    for case, record in (snapshot.get("runs") or {}).items():
        if str((record or {}).get("state") or "") in ACTIVE_JOB_STATES:
            labels.append(f"SCM:{case}")
    for run_id, record in (snapshot.get("loss_runs") or {}).items():
        if str((record or {}).get("state") or "") in ACTIVE_JOB_STATES:
            labels.append(f"Tune-result:{run_id}")
    return labels


def _report_active_jobs(prefix: str) -> None:
    labels = _active_job_labels()
    if labels:
        print(f"{prefix}: {', '.join(labels)}", file=sys.stderr, flush=True)
    else:
        print(f"{prefix}: none", file=sys.stderr, flush=True)


def _wait_for_pid_exit(pid: int, timeout: float) -> bool:
    deadline = time.monotonic() + max(0.0, float(timeout))
    while time.monotonic() < deadline:
        try:
            if not psutil.pid_exists(int(pid)):
                return True
        except (TypeError, ValueError):
            return True
        time.sleep(0.1)
    return not psutil.pid_exists(int(pid))


def _stop_orphaned_dashboard(status: dict[str, Any]) -> None:
    pid = int(status.get("pid") or 0)
    if pid < 1:
        return
    print(f"Adopting this checkout from a stopped manager; replacing orphaned Dash process {pid}.", flush=True)
    try:
        os.kill(pid, signal.SIGTERM)
    except ProcessLookupError:
        return
    except OSError as exc:
        raise RuntimeError(f"could not stop orphaned Dash process {pid}: {exc}") from exc
    if not _wait_for_pid_exit(pid, 5.0):
        raise RuntimeError(f"orphaned Dash process {pid} did not stop after SIGTERM")


def _post_broker_action(action: str, payload: dict[str, Any]) -> dict[str, Any]:
    """Call an existing broker without allowing this shutdown path to spawn one."""
    connection = read_connection()
    request = Request(
        str(connection["url"]).rstrip("/") + "/actions",
        data=json.dumps({"action": action, "payload": payload}).encode("utf-8"),
        headers={
            "Content-Type": "application/json",
            "X-CLUBB-Agent-Token": str(connection["token"]),
            "X-CLUBB-Dash-Internal": "1",
        },
        method="POST",
    )
    with urlopen(request, timeout=5.0) as response:  # nosec B310: private loopback record
        return json.loads(response.read().decode("utf-8"))


def _stop_dashboard(proc: subprocess.Popen[Any] | None) -> None:
    if proc is None or proc.poll() is not None:
        return
    try:
        os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
    except OSError:
        try:
            proc.terminate()
        except ProcessLookupError:
            return
    try:
        proc.wait(timeout=5.0)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
        except OSError:
            try:
                proc.kill()
            except ProcessLookupError:
                return
        proc.wait(timeout=3.0)


def _stop_managed_runtime(dash_proc: subprocess.Popen[Any] | None, reason: str) -> None:
    _stop_dashboard(dash_proc)
    _report_active_jobs("Active jobs at manager shutdown")
    try:
        _post_broker_action("stop_all_broker_work", {"reason": reason})
    except (HTTPError, URLError, OSError, RuntimeError, ValueError) as exc:
        print(f"Could not ask the broker to stop its workers: {exc}", file=sys.stderr, flush=True)
        # The durable records contain enough process identity for the same
        # idempotent stop implementation to work after a broker crash.
        try:
            from dash_app.shared.actions import stop_all_broker_work

            stop_all_broker_work(reason=reason)
        except Exception as fallback_exc:  # shutdown must continue to the broker
            print(f"Local worker cleanup also failed: {fallback_exc}", file=sys.stderr, flush=True)

    deadline = time.monotonic() + JOB_STOP_TIMEOUT_SECONDS
    while count_active_jobs(broker_jobs()) and time.monotonic() < deadline:
        time.sleep(0.25)
    try:
        stop_broker(timeout=5.0)
    except RuntimeError as exc:
        print(f"Could not stop the runtime broker cleanly: {exc}", file=sys.stderr, flush=True)
    reap_started_brokers()


class DashboardManager:
    def __init__(self, dash_args: list[str]) -> None:
        self.dash_args = list(dash_args)
        self.stop_event = threading.Event()
        self.stop_reason = "manager shutdown requested"
        self.exit_code = 0
        self.dash_proc: subprocess.Popen[Any] | None = None
        self.launched_once = False
        self.manager_pid = os.getpid()
        self.manager_started_at = process_started_at(self.manager_pid)

    def request_stop(self, reason: str, *, exit_code: int = 0) -> None:
        self.stop_reason = str(reason)
        self.exit_code = int(exit_code)
        self.stop_event.set()

    def _signal(self, signum: int, _frame: Any) -> None:
        name = signal.Signals(signum).name
        self.request_stop(f"received {name}")

    def _launch_dash(self) -> subprocess.Popen[Any]:
        env = os.environ.copy()
        env[MANAGER_REQUIRED_ENV] = "1"
        if self.launched_once:
            env[BROWSER_OPENED_ENV] = "1"
        command = [sys.executable, str(DASH_APP), *self.dash_args]
        # A separate process group lets the manager also stop a debug reloader
        # and any Dash background-callback children without signalling itself.
        proc = subprocess.Popen(  # noqa: S603
            command,
            cwd=str(REPO_ROOT),
            env=env,
            start_new_session=True,
        )
        self.launched_once = True
        print(f"Started CLUBB Dash (pid {proc.pid}).", flush=True)
        return proc

    def _broker_is_live(self) -> bool:
        try:
            return connection_is_live(read_connection())
        except RuntimeError:
            return False

    def run(self) -> int:
        os.environ[MANAGER_REQUIRED_ENV] = "1"
        write_manager_lease(pid=self.manager_pid, started_at=self.manager_started_at)
        signal.signal(signal.SIGINT, self._signal)
        signal.signal(signal.SIGTERM, self._signal)
        if hasattr(signal, "SIGHUP"):
            signal.signal(signal.SIGHUP, self._signal)

        ensure_broker(expected_fingerprint=runtime_source_fingerprint(REPO_ROOT))
        print("CLUBB dashboard manager is running; it owns Dash and the runtime broker.", flush=True)

        next_heartbeat = 0.0
        next_broker_attempt = 0.0
        next_dash_start = time.monotonic()
        dash_started = 0.0
        unhealthy_since: float | None = None
        run_was_healthy = False
        outage_started: float | None = None
        last_failure = "Dash did not start"

        while not self.stop_event.is_set():
            now = time.monotonic()
            if now >= next_heartbeat:
                write_manager_lease(pid=self.manager_pid, started_at=self.manager_started_at)
                next_heartbeat = now + HEARTBEAT_INTERVAL_SECONDS

            if now >= next_broker_attempt:
                reap_started_brokers()
                if not self._broker_is_live():
                    try:
                        ensure_broker(expected_fingerprint=runtime_source_fingerprint(REPO_ROOT))
                        print("Runtime broker restarted; durable job monitoring was recovered.", flush=True)
                    except RuntimeError as exc:
                        print(f"Runtime broker restart failed; retrying in 10 seconds: {exc}", file=sys.stderr, flush=True)
                next_broker_attempt = now + BROKER_RETRY_SECONDS

            if outage_started is not None and now - outage_started >= RESTART_DEADLINE_SECONDS:
                self.request_stop(
                    "Dash did not become healthy within 5 minutes; "
                    f"last failure: {last_failure}",
                    exit_code=1,
                )
                break

            if self.dash_proc is None:
                if now >= next_dash_start:
                    try:
                        self.dash_proc = self._launch_dash()
                        dash_started = now
                        unhealthy_since = now
                        run_was_healthy = False
                    except OSError as exc:
                        last_failure = f"could not launch Dash: {exc}"
                        outage_started = outage_started or now
                        next_dash_start = now + RESTART_INTERVAL_SECONDS
                        print(f"{last_failure}; retrying in 10 seconds.", file=sys.stderr, flush=True)
                self.stop_event.wait(POLL_SECONDS)
                continue

            returncode = self.dash_proc.poll()
            if returncode is not None:
                last_failure = f"Dash process {self.dash_proc.pid} exited with code {returncode}"
                outage_started = outage_started or now
                print(f"{last_failure}; retrying in 10 seconds.", file=sys.stderr, flush=True)
                _report_active_jobs("Broker-owned jobs remain visible while Dash is unavailable")
                self.dash_proc = None
                next_dash_start = now + RESTART_INTERVAL_SECONDS
                unhealthy_since = None
                continue

            if dashboard_registry.dashboard_is_live():
                if not run_was_healthy:
                    print(f"CLUBB Dash is healthy (pid {self.dash_proc.pid}).", flush=True)
                    if outage_started is not None:
                        print("Dash restart succeeded; restart deadline cleared.", flush=True)
                run_was_healthy = True
                outage_started = None
                unhealthy_since = None
            else:
                unhealthy_since = unhealthy_since or now
                health_timeout = (
                    RUNNING_HEALTH_TIMEOUT_SECONDS if run_was_healthy else STARTUP_HEALTH_TIMEOUT_SECONDS
                )
                if now - unhealthy_since >= health_timeout:
                    last_failure = (
                        f"Dash process {self.dash_proc.pid} stopped reporting a healthy dashboard "
                        f"for {health_timeout:.0f} seconds"
                    )
                    outage_started = outage_started or unhealthy_since
                    print(f"{last_failure}; retrying in 10 seconds.", file=sys.stderr, flush=True)
                    _report_active_jobs("Broker-owned jobs remain visible while Dash is unavailable")
                    _stop_dashboard(self.dash_proc)
                    self.dash_proc = None
                    next_dash_start = time.monotonic() + RESTART_INTERVAL_SECONDS
                    unhealthy_since = None

            self.stop_event.wait(POLL_SECONDS)

        print(f"CLUBB dashboard manager stopping: {self.stop_reason}", file=sys.stderr, flush=True)
        _stop_managed_runtime(self.dash_proc, self.stop_reason)
        return self.exit_code


def _run_locked(dash_args: list[str]) -> int:
    MANAGER_LOCK_PATH.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
    with MANAGER_LOCK_PATH.open("a+", encoding="utf-8") as lock_file:
        os.chmod(MANAGER_LOCK_PATH, 0o600)
        try:
            fcntl.flock(lock_file.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            status = dashboard_registry.dashboard_status()
            location = f" on port {status.get('port')}" if status.get("status") == "available" else ""
            print(f"A CLUBB dashboard manager is already running{location}.", file=sys.stderr)
            return 0 if location else 1

        prior_lease = read_manager_lease()
        dashboard = dashboard_registry.dashboard_status()
        if dashboard.get("status") == "available":
            if prior_lease and not manager_lease_is_live(prior_lease):
                _stop_orphaned_dashboard(dashboard)
            else:
                print(
                    f"An unmanaged CLUBB Dash is already running on port {dashboard.get('port')}; "
                    "stop it before switching this checkout to the manager.",
                    file=sys.stderr,
                )
                return 0

        manager = DashboardManager(dash_args)
        try:
            return manager.run()
        except Exception as exc:
            manager.request_stop(f"manager failure: {exc}", exit_code=1)
            print(f"CLUBB dashboard manager failed: {exc}", file=sys.stderr)
            _stop_managed_runtime(manager.dash_proc, manager.stop_reason)
            return 1
        finally:
            clear_manager_lease(pid=manager.manager_pid, started_at=manager.manager_started_at)
            fcntl.flock(lock_file.fileno(), fcntl.LOCK_UN)


def main() -> None:
    raise SystemExit(_run_locked(sys.argv[1:]))


if __name__ == "__main__":
    main()
