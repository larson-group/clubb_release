"""Dash-side registration for the broker's current browser consumer."""

from __future__ import annotations

import threading
import time
from dataclasses import dataclass, field

from dash_app.agent_integration import client
from dash_app.agent_integration.broker_endpoint import process_started_at


HEARTBEAT_INTERVAL_SECONDS = 1.5


@dataclass
class DashboardRegistration:
    pid: int
    started_at: float
    port: int
    _stop_event: threading.Event = field(default_factory=threading.Event, repr=False)
    _thread: threading.Thread | None = field(default=None, repr=False)

    @classmethod
    def register(cls, *, port: int) -> "DashboardRegistration":
        pid = __import__("os").getpid()
        started_at = process_started_at(pid)
        client.register_dashboard(pid=pid, started_at=started_at, port=port)
        handle = cls(pid=pid, started_at=started_at, port=int(port))
        handle._thread = threading.Thread(target=handle._heartbeat_loop, name="clubb-dashboard-heartbeat", daemon=True)
        handle._thread.start()
        return handle

    def _heartbeat_loop(self) -> None:
        while not self._stop_event.wait(HEARTBEAT_INTERVAL_SECONDS):
            try:
                result = client.heartbeat_dashboard(pid=self.pid, started_at=self.started_at)
                if result.get("status") == "unavailable":
                    # The broker may have restarted; retrying registration is
                    # safe because the registry rejects another live browser.
                    client.register_dashboard(pid=self.pid, started_at=self.started_at, port=self.port)
            except (OSError, RuntimeError, ValueError):
                continue

    def stop(self) -> None:
        if self._stop_event.is_set():
            return
        self._stop_event.set()
        try:
            client.unregister_dashboard(pid=self.pid, started_at=self.started_at)
        except (OSError, RuntimeError, ValueError):
            pass
