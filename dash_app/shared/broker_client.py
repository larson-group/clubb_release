"""Internal client for Dash callbacks talking to the durable runtime broker.

This client deliberately has no agent identity, session, or adapter concept.
Agent-facing clients live under :mod:`dash_app.agent_integration` and use the
same broker gateway through their own authentication path.
"""

from __future__ import annotations

import json
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

from .broker import ensure_broker
from .gateway import read_connection


def perform_action(
    action: str,
    payload: dict[str, Any],
    *,
    internal: bool = True,
    timeout_seconds: float = 10.0,
) -> dict[str, Any]:
    """Dispatch one Dash-owned action through the durable broker."""
    if not internal:
        raise ValueError("shared broker client only supports internal Dash actions")
    ensure_broker()
    connection = read_connection()
    request = Request(
        str(connection["url"]).rstrip("/") + "/actions",
        data=json.dumps({"action": action, "payload": payload}).encode("utf-8"),
        method="POST",
        headers={
            "Content-Type": "application/json",
            "X-CLUBB-Agent-Token": str(connection["token"]),
            "X-CLUBB-Dash-Internal": "1",
        },
    )
    try:
        with urlopen(request, timeout=max(0.1, float(timeout_seconds))) as response:  # nosec B310
            return json.loads(response.read().decode("utf-8"))
    except HTTPError as exc:
        try:
            error_payload = json.loads(exc.read().decode("utf-8", errors="replace"))
            error = error_payload.get("error") or {}
            message = str(error.get("message") or exc.reason)
        except (OSError, ValueError, TypeError, json.JSONDecodeError):
            message = str(exc.reason)
        raise RuntimeError(f"Dashboard broker rejected {action}: {message}") from None
    except (URLError, OSError) as exc:
        raise RuntimeError(f"Dashboard broker connection failed for {action}: {exc}") from None
