"""Regression tests for broker-owned endpoint and dashboard liveness state."""

import json
from pathlib import Path

import pytest

from dash_app.agent_integration import broker_endpoint
from dash_app.shared import dashboard_registry


def test_dashboard_registry_rejects_a_second_live_browser(tmp_path, monkeypatch):
    monkeypatch.setattr(dashboard_registry, "REGISTRY_PATH", tmp_path / "dashboard.json")
    monkeypatch.setattr(dashboard_registry, "REGISTRY_LOCK_PATH", tmp_path / "dashboard.lock")
    monkeypatch.setattr(
        dashboard_registry,
        "process_identity_is_live",
        lambda pid, started: (int(pid), float(started)) in {(101, 1.0), (202, 2.0)},
    )

    first = dashboard_registry.register_dashboard(pid=101, started_at=1.0, port=8060)
    assert first["status"] == "available"
    with pytest.raises(dashboard_registry.DashboardAlreadyRegistered):
        dashboard_registry.register_dashboard(pid=202, started_at=2.0, port=8061)

    assert dashboard_registry.heartbeat_dashboard(pid=101, started_at=1.0)["status"] == "available"
    assert dashboard_registry.unregister_dashboard(pid=202, started_at=2.0) is False
    assert dashboard_registry.unregister_dashboard(pid=101, started_at=1.0) is True


def test_dashboard_registry_reconcile_removes_stale_record(tmp_path, monkeypatch):
    monkeypatch.setattr(dashboard_registry, "REGISTRY_PATH", tmp_path / "dashboard.json")
    monkeypatch.setattr(dashboard_registry, "REGISTRY_LOCK_PATH", tmp_path / "dashboard.lock")
    dashboard_registry.REGISTRY_PATH.write_text(
        json.dumps({"pid": 7, "started_at": 1.0, "last_heartbeat": 1.0}), encoding="utf-8"
    )
    monkeypatch.setattr(dashboard_registry, "process_identity_is_live", lambda *_args: False)

    assert dashboard_registry.reconcile_dashboard()["status"] == "unavailable"
    assert not dashboard_registry.REGISTRY_PATH.exists()


def test_stable_endpoint_identity_survives_runtime_reconciliation(tmp_path, monkeypatch):
    root = tmp_path / "mcp-endpoint"
    root.mkdir()
    monkeypatch.setattr(broker_endpoint, "ENDPOINT_ROOT", root)
    monkeypatch.setattr(broker_endpoint, "IDENTITY_PATH", root / "identity.json")
    monkeypatch.setattr(broker_endpoint, "_local_port", lambda: 43210)

    first = broker_endpoint._identity()
    second = broker_endpoint._identity()

    assert first == second
    assert first["endpoint_port"] == 43210
    assert Path(broker_endpoint.IDENTITY_PATH).stat().st_mode & 0o777 == 0o600
