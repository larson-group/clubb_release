"""Tiny generic adapter client for a locally running CLUBB Dash gateway.

An agent can use this module directly, or an integration can copy the four
HTTP operations it uses.  It intentionally needs no model-specific package.
"""

from __future__ import annotations

import argparse
import json
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

from dash_app.shared.gateway import read_connection


class DashboardAlreadyRegisteredError(RuntimeError):
    """Raised when another live dashboard owns this checkout."""


def request_gateway(
    method: str,
    path: str,
    payload: dict[str, Any] | None = None,
    *,
    internal: bool = False,
    timeout_seconds: float = 10.0,
    connection: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Send one authenticated, loopback-only gateway request."""
    if internal:
        # Dash callbacks must not keep talking to a durable broker imported
        # from an older checkout after a debug reload.  In particular, Tune's
        # typed request schema evolves independently of the dashboard UI.
        # ``ensure_broker`` replaces only an *idle* incompatible broker and
        # preserves one with active compile/run/Tune work.
        from dash_app.shared.broker import ensure_broker

        ensure_broker()
    data = json.dumps(payload).encode("utf-8") if payload is not None else None
    def _send(connection: dict[str, Any]) -> dict[str, Any]:
        headers = {"Content-Type": "application/json", "X-CLUBB-Agent-Token": str(connection["token"])}
        if internal:
            headers["X-CLUBB-Dash-Internal"] = "1"
        request = Request(
            str(connection["url"]).rstrip("/") + "/" + path.lstrip("/"),
            data=data,
            method=method.upper(),
            headers=headers,
        )
        with urlopen(  # nosec B310: fixed loopback URL from mode-0600 record
            request,
            timeout=max(0.1, float(timeout_seconds)),
        ) as response:
            return json.loads(response.read().decode("utf-8"))

    # A Dash reload can race a broker replacement: the first request may use a
    # record whose port/token has just changed.  Re-read and retry once.  This
    # does not create a broker for an external agent, so a failed discovery is
    # never mistaken for permission to launch a dashboard.
    first_connection = dict(connection or read_connection())
    try:
        return _send(first_connection)
    except (HTTPError, URLError, OSError, ValueError) as first_error:
        try:
            retry_connection = dict(connection or read_connection())
        except RuntimeError:
            retry_connection = first_connection
        if retry_connection != first_connection:
            try:
                return _send(retry_connection)
            except (HTTPError, URLError, OSError, ValueError) as retry_error:
                first_error = retry_error
        if not isinstance(first_error, HTTPError):
            raise RuntimeError(
                "Dashboard broker connection failed while it may be reloading; "
                "retry the connection once Dash is ready."
            ) from None
        exc = first_error
        # Native Dash callbacks used to surface only ``HTTP Error 400``.  The
        # broker already returns a bounded structured envelope; preserve that
        # useful reason in the console while never leaking an implementation
        # traceback from the local service.
        try:
            payload = json.loads(exc.read().decode("utf-8", errors="replace"))
            raw_error = payload.get("error")
            if isinstance(raw_error, dict):
                code = str(raw_error.get("code") or f"HTTP_{exc.code}")
                message = str(raw_error.get("message") or exc.reason)
            else:
                # Version-2 brokers returned a bare string.  Name the real
                # remediation rather than hiding it behind a generic 400.
                message = str(raw_error or exc.reason)
                code = "BROKER_RESTART_REQUIRED" if "unknown action" in message.lower() else f"HTTP_{exc.code}"
                if code == "BROKER_RESTART_REQUIRED":
                    message = "the persistent dashboard broker predates this Dash reload; refresh Dash to replace its idle broker"
        except (OSError, ValueError, TypeError, json.JSONDecodeError):
            code = f"HTTP_{exc.code}"
            message = str(exc.reason)
        if code == "DASHBOARD_ALREADY_REGISTERED":
            raise DashboardAlreadyRegisteredError(message) from None
        raise RuntimeError(f"Dashboard broker rejected {path}: {code}: {message}") from None


def connect(*, connection: dict[str, Any] | None = None) -> dict[str, Any]:
    """Verify one short-lived, authenticated broker connection."""
    result = request_gateway("GET", "status", connection=connection)
    result["status"] = "connected"
    result["allowed_actions"] = ["inspect_dashboard", "invoke_dashboard"]
    return result


def perform_action(
    action: str,
    payload: dict[str, Any],
    *,
    internal: bool = False,
    timeout_seconds: float = 10.0,
    connection: dict[str, Any] | None = None,
) -> dict[str, Any]:
    return request_gateway(
        "POST",
        "actions",
        {"action": action, "payload": payload},
        internal=internal,
        timeout_seconds=timeout_seconds,
        connection=connection,
    )


def register_dashboard(*, pid: int, started_at: float, port: int) -> dict[str, Any]:
    return request_gateway(
        "POST",
        "dashboard/register",
        {"pid": int(pid), "started_at": float(started_at), "port": int(port)},
        internal=True,
        connection=None,
    )


def heartbeat_dashboard(*, pid: int, started_at: float) -> dict[str, Any]:
    return request_gateway(
        "POST",
        "dashboard/heartbeat",
        {"pid": int(pid), "started_at": float(started_at)},
        internal=True,
        connection=None,
    )


def unregister_dashboard(*, pid: int, started_at: float) -> dict[str, Any]:
    return request_gateway(
        "POST",
        "dashboard/unregister",
        {"pid": int(pid), "started_at": float(started_at)},
        internal=True,
        connection=None,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Use one transient authenticated request with CLUBB Dash.")
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("connect", help="verify the broker and list its current status")
    action = subparsers.add_parser("action", help="send one documented semantic action to Dash")
    action.add_argument("name")
    action.add_argument("--payload-json", default="{}", help="JSON object passed as the action payload")
    args = parser.parse_args()
    if args.command == "connect":
        result = connect()
    elif args.command == "action":
        try:
            payload = json.loads(args.payload_json)
        except json.JSONDecodeError as exc:
            parser.error(f"--payload-json must be a JSON object: {exc}")
        if not isinstance(payload, dict):
            parser.error("--payload-json must be a JSON object")
        result = perform_action(args.name, payload)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
