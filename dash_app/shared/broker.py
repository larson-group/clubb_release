#!/usr/bin/env python3
"""Durable local owner for dashboard runtime state and semantic jobs.

Dash is intentionally not a process supervisor.  This broker owns the local
gateway, activity stream, action watchers, and tuner keepalive loop.  The Dash
process can therefore reload or crash during development; when it returns, its
tabs read the same durable snapshot and re-adopt the jobs.
Agent adapters are clients of this service, not its owner.
"""

from __future__ import annotations

import argparse
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
from urllib.error import URLError
from urllib.request import Request, urlopen

from flask import Flask, jsonify
import psutil
from werkzeug.serving import make_server

from .activity import REPO_ROOT, broker_jobs, count_active_jobs, publish_event, set_broker_metadata, update_broker_job
from . import dashboard_registry
from dash_app.agent_integration.broker_endpoint import (
    process_started_at as broker_process_started_at,
    start_broker_endpoint,
    stop_broker_endpoint,
)
from .gateway import API_PREFIX, BROKER_LOCK_PATH, BROKER_LOG_PATH, CONNECTION_PATH, install_gateway_routes, read_connection, update_connection_endpoint, update_connection_logs, write_connection
from .broker_protocol import BROKER_PROTOCOL_VERSION
from .provenance import runtime_source_fingerprint
from .manager_lease import LEASE_TIMEOUT_SECONDS, MANAGER_REQUIRED_ENV, manager_lease_is_live


BROKER_START_TIMEOUT_SECONDS = 8.0
BROKER_WORK_STOP_TIMEOUT_SECONDS = 12.0
_STARTED_BROKER_PROCESSES: list[subprocess.Popen[Any]] = []


def _pid_is_alive(pid: Any) -> bool:
    try:
        process = psutil.Process(int(pid))
        return process.status() != psutil.STATUS_ZOMBIE
    except (psutil.Error, TypeError, ValueError, OSError):
        return False


def connection_is_live(connection: dict[str, Any], *, timeout: float = 0.5) -> bool:
    """Return whether an authenticated broker accepts a lightweight status request."""
    try:
        request = Request(
            str(connection["url"]).rstrip("/") + "/status",
            headers={"X-CLUBB-Agent-Token": str(connection["token"])},
            method="GET",
        )
        with urlopen(request, timeout=timeout) as response:  # nosec B310: mode-0600 loopback record
            return response.status == 200
    except (KeyError, URLError, OSError, ValueError):
        return False


def _existing_live_connection() -> dict[str, Any] | None:
    try:
        connection = read_connection()
    except RuntimeError:
        return None
    return connection if connection_is_live(connection) else None


def _connection_is_compatible(connection: dict[str, Any]) -> bool:
    """Whether a live broker can serve this Dash source revision's actions."""
    try:
        return int(connection.get("version") or 0) >= BROKER_PROTOCOL_VERSION
    except (TypeError, ValueError):
        return False


def _connection_is_stale(connection: dict[str, Any], expected_fingerprint: str | None) -> bool:
    """Return whether a live broker was started from different runtime code."""
    if os.environ.get(MANAGER_REQUIRED_ENV) == "1" and not bool(connection.get("manager_required")):
        return True
    if not expected_fingerprint:
        return False
    return str(connection.get("runtime_fingerprint") or "") != str(expected_fingerprint)


def _warn_stale_runtime(connection: dict[str, Any], expected_fingerprint: str) -> None:
    actual = str(connection.get("runtime_fingerprint") or "missing")
    print(
        "CLUBB Dash runtime broker code is stale; preserving it because it owns active work "
        f"(actual={actual[:12]}, expected={expected_fingerprint[:12]}). "
        "Restart Dash after the work finishes to replace the broker and MCP endpoint.",
        file=sys.stderr,
    )


def _stable_endpoint_is_published(connection: dict[str, Any] | None) -> bool:
    endpoint = dict((connection or {}).get("mcp_endpoint") or {})
    return bool(endpoint.get("endpoint_url") and endpoint.get("endpoint_token") and endpoint.get("instance_id"))


def _broker_has_active_work(connection: dict[str, Any]) -> bool:
    """Avoid replacing a compatible-unknown sidecar while it owns real work."""
    try:
        request = Request(
            str(connection["url"]).rstrip("/") + "/status",
            headers={"X-CLUBB-Agent-Token": str(connection["token"])},
            method="GET",
        )
        with urlopen(request, timeout=0.75) as response:  # nosec B310: local private connection record
            jobs = dict((json.loads(response.read().decode("utf-8")).get("jobs") or {}))
    except (KeyError, OSError, URLError, ValueError, TypeError):
        # A live connection that cannot be inspected is safer to preserve than
        # to terminate during development.
        return True
    active_states = {"queued", "submitting", "running", "stopping"}

    def contains_active_state(value: Any) -> bool:
        if isinstance(value, dict):
            if str(value.get("state") or "") in active_states:
                return True
            return any(contains_active_state(item) for item in value.values())
        if isinstance(value, (list, tuple)):
            return any(contains_active_state(item) for item in value)
        return False

    return contains_active_state(jobs)


def _stop_connection_process(connection: dict[str, Any], *, timeout: float = 3.0) -> bool:
    """Stop only the recorded idle sidecar and wait for its listener to exit."""
    try:
        pid = int(connection.get("pid"))
        if pid < 1:
            return False
        os.kill(pid, signal.SIGTERM)
    except (OSError, TypeError, ValueError):
        return False
    deadline = time.monotonic() + max(0.1, float(timeout))
    while _pid_is_alive(pid) and time.monotonic() < deadline:
        time.sleep(0.05)
    return not _pid_is_alive(pid)


def _start_broker_process() -> None:
    command = [sys.executable, "-m", "dash_app.shared.broker", "serve"]
    BROKER_LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    with BROKER_LOG_PATH.open("a", encoding="utf-8") as log_file:
        os.chmod(BROKER_LOG_PATH, 0o600)
        process = subprocess.Popen(  # noqa: S603 - fixed module command, never agent-provided text
            command,
            cwd=str(REPO_ROOT),
            stdin=subprocess.DEVNULL,
            stdout=log_file,
            stderr=subprocess.STDOUT,
            # Signal isolation keeps a terminal Ctrl-C ordered through the
            # manager. This is still its direct child, not an orphan: the
            # manager lease makes abrupt parent death terminal as well.
            start_new_session=True,
            close_fds=True,
        )
    _STARTED_BROKER_PROCESSES.append(process)


def reap_started_brokers() -> list[int]:
    """Reap broker children launched by this process and return their exit codes."""
    returncodes: list[int] = []
    survivors: list[subprocess.Popen[Any]] = []
    for process in _STARTED_BROKER_PROCESSES:
        returncode = process.poll()
        if returncode is None:
            survivors.append(process)
        else:
            returncodes.append(int(returncode))
    _STARTED_BROKER_PROCESSES[:] = survivors
    return returncodes


def ensure_broker(
    *,
    timeout: float = BROKER_START_TIMEOUT_SECONDS,
    expected_fingerprint: str | None = None,
    require_fresh: bool = False,
) -> Path:
    """Return the stable connection record, launching one detached broker if needed."""
    existing = _existing_live_connection()
    stale = existing is not None and _connection_is_stale(existing, expected_fingerprint)
    if existing is not None and _connection_is_compatible(existing) and not stale:
        if _stable_endpoint_is_published(existing):
            return CONNECTION_PATH
        deadline = time.monotonic() + max(0.5, float(timeout))
        while time.monotonic() < deadline:
            ready = _existing_live_connection()
            if ready is not None and _connection_is_compatible(ready) and not _connection_is_stale(ready, expected_fingerprint) and _stable_endpoint_is_published(ready):
                return CONNECTION_PATH
            time.sleep(0.1)
        return CONNECTION_PATH

    BROKER_LOCK_PATH.parent.mkdir(parents=True, exist_ok=True)
    with BROKER_LOCK_PATH.open("a+", encoding="utf-8") as lock_file:
        os.chmod(BROKER_LOCK_PATH, 0o600)
        fcntl.flock(lock_file.fileno(), fcntl.LOCK_EX)
        try:
            existing = _existing_live_connection()
            stale = existing is not None and _connection_is_stale(existing, expected_fingerprint)
            if existing is not None and _connection_is_compatible(existing) and not stale:
                if _stable_endpoint_is_published(existing):
                    return CONNECTION_PATH
                deadline = time.monotonic() + max(0.5, float(timeout))
                while time.monotonic() < deadline:
                    ready = _existing_live_connection()
                    if ready is not None and _connection_is_compatible(ready) and not _connection_is_stale(ready, expected_fingerprint) and _stable_endpoint_is_published(ready):
                        return CONNECTION_PATH
                    time.sleep(0.1)
                return CONNECTION_PATH

            stale = existing is not None and (
                not _connection_is_compatible(existing)
                or _connection_is_stale(existing, expected_fingerprint)
            )
            if existing is not None and stale:
                dashboard_live = dashboard_registry.dashboard_is_live()
                active_work = _broker_has_active_work(existing)
                if dashboard_live:
                    if require_fresh:
                        raise RuntimeError("cannot replace the stale runtime broker while a dashboard is live; close the dashboard and retry --restart-runtime")
                    if expected_fingerprint and _connection_is_stale(existing, expected_fingerprint):
                        _warn_stale_runtime(existing, expected_fingerprint)
                    return CONNECTION_PATH
                if active_work:
                    if require_fresh:
                        raise RuntimeError("cannot replace the stale runtime broker while broker-owned work is active; wait for completion and retry --restart-runtime")
                    if expected_fingerprint and _connection_is_stale(existing, expected_fingerprint):
                        _warn_stale_runtime(existing, expected_fingerprint)
                    return CONNECTION_PATH
                if not _stop_connection_process(existing):
                    raise RuntimeError("The dashboard runtime broker is out of date and could not be restarted automatically.")

            try:
                stale_connection = read_connection()
            except RuntimeError:
                stale_connection = None
            if stale_connection and _pid_is_alive(stale_connection.get("pid")):
                raise RuntimeError(
                    "The dashboard runtime broker process exists but is not responding; "
                    "stop that broker before starting another one."
                )

            _start_broker_process()
            deadline = time.monotonic() + max(0.5, float(timeout))
            while time.monotonic() < deadline:
                ready = _existing_live_connection()
                if ready is not None and _stable_endpoint_is_published(ready):
                    return CONNECTION_PATH
                time.sleep(0.1)
        finally:
            fcntl.flock(lock_file.fileno(), fcntl.LOCK_UN)
    raise RuntimeError(f"Timed out while starting the local dashboard runtime broker; see {BROKER_LOG_PATH}")


def _recover_tune_keepalive() -> None:
    """Resume the durable tuner lease if a broker itself was restarted mid-job."""
    tune = (broker_jobs().get("tune") or {})
    if str(tune.get("state") or "") not in {"running", "stopping"}:
        try:
            from . import actions

            actions.recover_active_tuning_from_disk()
        except (OSError, RuntimeError, ValueError) as exc:
            publish_event("broker", "Could not scan for a running Tune job", str(exc), status="error")
        return
    job = dict(tune.get("job") or {})
    if not job:
        return
    try:
        from . import actions

        actions._background(actions._watch_tuning, job, dict(tune.get("request") or {}))
        publish_event("broker", "Resumed tuner keepalive", str(job.get("job_dir") or ""), status="info")
    except (OSError, RuntimeError, ValueError) as exc:
        update_broker_job("tune", state="error", recovery_error=str(exc))


def _recover_scm_monitoring() -> None:
    """Resume log/status monitoring for SCM children that outlived a broker."""
    try:
        from . import actions

        actions.recover_active_runs_from_state()
        actions.recover_queued_scm_batches()
    except (OSError, RuntimeError, ValueError) as exc:
        publish_event("broker", "Could not recover SCM job monitoring", str(exc), status="error")


def _recover_compile_monitoring() -> None:
    """Resume monitoring of a compile that outlived the broker process."""
    try:
        from . import actions

        actions.recover_active_compile_from_state()
    except (OSError, RuntimeError, ValueError) as exc:
        publish_event("broker", "Could not recover compile monitoring", str(exc), status="error")


def stop_broker(*, timeout: float = 3.0) -> None:
    """Deliberately stop the detached broker, e.g. after editing broker code."""
    try:
        connection = read_connection()
    except RuntimeError as exc:
        raise RuntimeError("No local dashboard runtime broker is running") from exc
    pid = connection.get("pid")
    if not _pid_is_alive(pid):
        return
    try:
        os.kill(int(pid), signal.SIGTERM)
    except (OSError, TypeError, ValueError) as exc:
        raise RuntimeError("Could not stop the local dashboard runtime broker") from exc
    deadline = time.monotonic() + max(0.1, float(timeout))
    while _pid_is_alive(pid) and time.monotonic() < deadline:
        time.sleep(0.05)
    if _pid_is_alive(pid):
        raise RuntimeError("Broker did not stop within the requested timeout")


def create_broker_app(connection: dict[str, Any]) -> Flask:
    """Build the lightweight broker Flask app without importing Dash itself."""
    app = Flask("clubb_dashboard_runtime_broker")
    install_gateway_routes(app, connection)

    @app.get(f"{API_PREFIX}/health")
    def health():
        return jsonify({"status": "ok", "pid": os.getpid()})

    return app


def _stop_after_manager_lease_expires(server: Any, stop_event: threading.Event) -> None:
    """Stop owned work if an abruptly lost manager is not replaced in time."""
    while not stop_event.wait(1.0):
        if manager_lease_is_live():
            continue
        reason = (
            f"dashboard manager heartbeat expired after {LEASE_TIMEOUT_SECONDS:.0f} seconds; "
            "stopping broker-owned work and the broker"
        )
        print(f"CLUBB runtime broker stopping: {reason}", file=sys.stderr, flush=True)
        dashboard = dashboard_registry.dashboard_status()
        if dashboard.get("status") == "available":
            try:
                dashboard_pid = int(dashboard["pid"])
                os.killpg(os.getpgid(dashboard_pid), signal.SIGTERM)
                deadline = time.monotonic() + 5.0
                while (
                    dashboard_registry.process_identity_is_live(
                        dashboard_pid, dashboard.get("started_at")
                    )
                    and time.monotonic() < deadline
                ):
                    time.sleep(0.1)
                if dashboard_registry.process_identity_is_live(
                    dashboard_pid, dashboard.get("started_at")
                ):
                    os.killpg(os.getpgid(dashboard_pid), signal.SIGKILL)
            except (KeyError, OSError, TypeError, ValueError) as exc:
                print(f"Could not stop the orphaned Dash process: {exc}", file=sys.stderr, flush=True)
        try:
            from . import actions

            result = actions.stop_all_broker_work(reason=reason)
            if result.get("errors"):
                print(
                    "Some broker-owned workers could not be stopped cleanly: "
                    + "; ".join(str(item) for item in result["errors"]),
                    file=sys.stderr,
                    flush=True,
                )
        except Exception as exc:
            publish_event("broker", "Manager lease cleanup failed", str(exc), status="error")
            print(f"Broker worker cleanup failed: {exc}", file=sys.stderr, flush=True)

        deadline = time.monotonic() + BROKER_WORK_STOP_TIMEOUT_SECONDS
        while count_active_jobs(broker_jobs()) and time.monotonic() < deadline:
            if stop_event.wait(0.25):
                return
        stop_event.set()
        server.shutdown()
        return


def serve() -> None:
    """Run one loopback broker until its manager or an explicit stop ends it."""
    app = Flask("clubb_dashboard_runtime_broker")
    server = make_server("127.0.0.1", 0, app, threaded=True)
    managed = os.environ.get(MANAGER_REQUIRED_ENV) == "1"
    stop_event = threading.Event()

    def request_shutdown(_signum: int, _frame: Any) -> None:
        stop_event.set()
        threading.Thread(target=server.shutdown, name="clubb-broker-shutdown", daemon=True).start()

    previous_sigterm = signal.signal(signal.SIGTERM, request_shutdown)
    previous_sigint = None
    if managed:
        # Ctrl-C belongs to the foreground manager, which performs the ordered
        # Dash -> workers -> broker shutdown. The broker waits for that request.
        previous_sigint = signal.signal(signal.SIGINT, signal.SIG_IGN)
    # The stable MCP endpoint validates this timestamp against the broker's
    # actual OS process creation time.  ``time.time()`` here would measure how
    # long imports/socket setup took and can fail that identity check on slower
    # macOS starts.
    broker_started_at = broker_process_started_at(os.getpid())
    connection = write_connection(
        server.server_port,
        pid=os.getpid(),
        runtime_fingerprint=runtime_source_fingerprint(REPO_ROOT),
    )
    install_gateway_routes(app, connection)
    set_broker_metadata(pid=os.getpid(), url=connection["url"], started_at=broker_started_at, state="running")
    _recover_compile_monitoring()
    _recover_tune_keepalive()
    _recover_scm_monitoring()
    server_thread = threading.Thread(target=server.serve_forever, name="clubb-broker-http", daemon=True)
    server_thread.start()
    lease_thread = None
    if managed:
        lease_thread = threading.Thread(
            target=_stop_after_manager_lease_expires,
            args=(server, stop_event),
            name="clubb-manager-lease",
            daemon=True,
        )
        lease_thread.start()
    try:
        endpoint = start_broker_endpoint(
            connection,
            owner_pid=os.getpid(),
            owner_started_at=broker_started_at,
        )
        update_connection_endpoint(endpoint)
        log_paths = {"broker": str(BROKER_LOG_PATH)}
        if endpoint.get("endpoint_log_path"):
            log_paths["endpoint"] = str(endpoint["endpoint_log_path"])
        update_connection_logs(log_paths)
        server_thread.join()
    finally:
        stop_event.set()
        stop_broker_endpoint()
        server.shutdown()
        signal.signal(signal.SIGTERM, previous_sigterm)
        if previous_sigint is not None:
            signal.signal(signal.SIGINT, previous_sigint)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run or inspect the durable local CLUBB dashboard runtime broker.")
    parser.add_argument("command", choices=("serve", "status", "stop"), nargs="?", default="serve")
    args = parser.parse_args()
    if args.command == "status":
        connection = _existing_live_connection()
        if connection is None:
            print("broker: not running")
            raise SystemExit(1)
        print(f"broker: running ({connection.get('url')})")
        return
    if args.command == "stop":
        stop_broker()
        print("broker: stopped")
        return
    serve()


if __name__ == "__main__":
    main()
