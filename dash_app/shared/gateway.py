"""Authenticated gateway helpers for the durable local dashboard broker."""

from __future__ import annotations

import hashlib
import hmac
import json
import os
import secrets
from pathlib import Path
from typing import Any

from flask import Flask, jsonify, request
from werkzeug.exceptions import HTTPException

from . import actions
from .broker_protocol import BROKER_PROTOCOL_VERSION
from .provenance import RUNTIME_FINGERPRINT_SCOPE
from .activity import REPO_ROOT, broker_jobs, read_activity
from . import dashboard_registry
from dash_app.shared.runtime import atomic_write_json, private_path, readable_private_paths, restrict_existing
from dash_app.shared.manager_lease import MANAGER_REQUIRED_ENV


_TOKEN = hashlib.sha256(str(REPO_ROOT).encode("utf-8")).hexdigest()[:16]
_LEGACY_CONNECTION_PATH = Path("/tmp") / f"clubb_dash_agent_gateway_{_TOKEN}.json"
CONNECTION_PATH = private_path(REPO_ROOT, "connection.json")
BROKER_LOCK_PATH = private_path(REPO_ROOT, "broker.lock")
BROKER_LOG_PATH = private_path(REPO_ROOT, "broker.log")
restrict_existing(_LEGACY_CONNECTION_PATH)
API_PREFIX = "/api/agent-integration"


def write_connection(
    port: int,
    *,
    pid: int | None = None,
    runtime_fingerprint: str | None = None,
) -> dict[str, str | int]:
    """Publish a stable, broker-owned loopback connection record."""
    payload: dict[str, str | int] = {
        "version": BROKER_PROTOCOL_VERSION,
        "url": f"http://127.0.0.1:{int(port)}{API_PREFIX}",
        "token": secrets.token_urlsafe(32),
        "repository": str(REPO_ROOT),
        "manager_required": os.environ.get(MANAGER_REQUIRED_ENV) == "1",
    }
    if pid is not None:
        payload["pid"] = int(pid)
    if runtime_fingerprint:
        payload["runtime_fingerprint"] = str(runtime_fingerprint)
        payload["runtime_fingerprint_scope"] = RUNTIME_FINGERPRINT_SCOPE
    atomic_write_json(CONNECTION_PATH, payload)
    return payload


def update_connection_endpoint(endpoint: dict[str, Any]) -> dict[str, Any]:
    """Publish stable MCP details without changing the broker API identity."""
    connection = read_connection()
    connection["mcp_endpoint"] = {
        "instance_id": endpoint.get("instance_id"),
        "endpoint_url": endpoint.get("endpoint_url"),
        "endpoint_token": endpoint.get("endpoint_token"),
        "endpoint_port": endpoint.get("endpoint_port"),
    }
    atomic_write_json(CONNECTION_PATH, connection)
    return connection


def update_connection_logs(log_paths: dict[str, str]) -> dict[str, Any]:
    """Publish component log paths without publishing credentials or payloads."""

    connection = read_connection()
    connection["log_paths"] = {
        str(name): str(Path(path).resolve())
        for name, path in log_paths.items()
        if name and path
    }
    atomic_write_json(CONNECTION_PATH, connection)
    return connection


def read_connection() -> dict[str, Any]:
    """Read the broker-published local connection information.

    Check all private runtime roots before falling back to the pre-hardening
    path.  This keeps an editor-launched MCP adapter attached when Dash was
    launched from a shell with a different ``XDG_RUNTIME_DIR``.
    """

    candidates = [CONNECTION_PATH, *readable_private_paths(REPO_ROOT, "connection.json"), _LEGACY_CONNECTION_PATH]
    seen: set[Path] = set()
    last_error: Exception | None = None
    for path in candidates:
        if path in seen:
            continue
        seen.add(path)
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError, TypeError) as exc:
            last_error = exc
            continue
        if isinstance(payload, dict) and payload.get("url") and payload.get("token"):
            return payload
        last_error = RuntimeError("Dash published an invalid agent connection record")
    if isinstance(last_error, RuntimeError):
        raise last_error
    raise RuntimeError("The dashboard runtime broker is not running; start ./launch_dashboard.sh first") from last_error


def _local_request() -> bool:
    address = str(request.remote_addr or "")
    return address in {"127.0.0.1", "::1"}


def _authenticated(connection: dict[str, Any]) -> bool:
    supplied = str(request.headers.get("X-CLUBB-Agent-Token") or "")
    expected = str(connection.get("token") or "")
    return _local_request() and bool(expected) and hmac.compare_digest(supplied, expected)


def install_gateway_routes(server: Flask, connection: dict[str, Any]) -> None:
    """Install broker routes on one Flask server using an already-issued token."""
    if getattr(server, "_clubb_agent_gateway_configured", False):
        return
    server._clubb_agent_gateway_configured = True
    # Flask installs ``None`` by default, so ``setdefault`` would not enforce
    # the broker boundary.  This is intentionally a hard local API limit.
    server.config["MAX_CONTENT_LENGTH"] = 256 * 1024

    @server.errorhandler(HTTPException)
    def _http_error(error):
        return jsonify({"error": {"code": error.name.upper().replace(" ", "_"), "message": error.description}}), error.code

    @server.errorhandler(Exception)
    def _unexpected_error(_error):
        return jsonify({"error": {"code": "INTERNAL_ERROR", "message": "dashboard broker action failed"}}), 500

    def api_error(code: str, message: str, status: int):
        """Return the one stable error envelope used by every broker route."""
        return jsonify({"error": {"code": str(code), "message": str(message)[:500]}}), int(status)

    def require_auth():
        if not _authenticated(connection):
            return api_error("AUTH_REQUIRED", "local authenticated agent connection required", 403)
        return None

    def require_internal():
        denied = require_auth()
        if denied:
            return denied
        if request.headers.get("X-CLUBB-Dash-Internal") != "1":
            return api_error("INTERNAL_ONLY", "dashboard lifecycle requests are broker-internal", 403)
        return None

    @server.post(f"{API_PREFIX}/dashboard/register")
    def register_dashboard():
        denied = require_internal()
        if denied:
            return denied
        payload = request.get_json(silent=True) or {}
        try:
            result = dashboard_registry.register_dashboard(
                pid=int(payload.get("pid")),
                started_at=float(payload.get("started_at")),
                port=int(payload.get("port")),
            )
        except dashboard_registry.DashboardAlreadyRegistered as exc:
            return api_error("DASHBOARD_ALREADY_REGISTERED", str(exc), 409)
        except (TypeError, ValueError, OSError) as exc:
            return api_error("INVALID_DASHBOARD_REGISTRATION", str(exc), 400)
        return jsonify(result), 200

    @server.post(f"{API_PREFIX}/dashboard/heartbeat")
    def heartbeat_dashboard():
        denied = require_internal()
        if denied:
            return denied
        payload = request.get_json(silent=True) or {}
        try:
            result = dashboard_registry.heartbeat_dashboard(
                pid=int(payload.get("pid")), started_at=float(payload.get("started_at"))
            )
        except (TypeError, ValueError, OSError) as exc:
            return api_error("INVALID_DASHBOARD_HEARTBEAT", str(exc), 400)
        return jsonify(result), 200

    @server.post(f"{API_PREFIX}/dashboard/unregister")
    def unregister_dashboard():
        denied = require_internal()
        if denied:
            return denied
        payload = request.get_json(silent=True) or {}
        try:
            removed = dashboard_registry.unregister_dashboard(
                pid=int(payload.get("pid")), started_at=float(payload.get("started_at"))
            )
        except (TypeError, ValueError, OSError) as exc:
            return api_error("INVALID_DASHBOARD_UNREGISTER", str(exc), 400)
        return jsonify({"status": "unregistered" if removed else "not_current"}), 200

    @server.post(f"{API_PREFIX}/actions")
    def perform_action():
        denied = require_auth()
        if denied:
            return denied
        payload = request.get_json(silent=True) or {}
        # Native Dash callbacks use the same broker service but are marked as
        # internal so they may dispatch typed domain operations.
        internal = request.headers.get("X-CLUBB-Dash-Internal") == "1"
        action_name = str(payload.get("action") or "")
        if not internal and action_name not in actions.allowed_action_names():
            return api_error(
                "DEPRECATED_DIRECT_ACTION",
                "external actions are limited to inspect_dashboard/invoke_dashboard browser handoff; use typed MCP domain tools for scientific work",
                400,
            )
        try:
            result = actions.dispatch(action_name, dict(payload.get("payload") or {}))
        except (TypeError, ValueError) as exc:
            return api_error("INVALID_REQUEST", str(exc), 400)
        return jsonify(result), 202

    @server.get(f"{API_PREFIX}/status")
    def status():
        denied = require_auth()
        if denied:
            return denied
        snapshot = read_activity()
        return jsonify(
            {
                "jobs": broker_jobs(),
                "broker": snapshot.get("broker") or {},
                "dashboard": dashboard_registry.dashboard_status(),
                "latest_event_id": snapshot.get("next_id", 1) - 1,
            }
        )


def configure_gateway(server: Flask, port: int) -> Path:
    """Compatibility helper for focused tests; production uses ``broker.py``."""
    connection = write_connection(port, pid=os.getpid())
    install_gateway_routes(server, connection)

    return CONNECTION_PATH
