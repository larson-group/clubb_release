"""PyPlotGen export lifecycle and fixed-path serving tests."""

from pathlib import Path

from dash import Dash, html

from dash_app.plot_tab import pyplotgen_runtime, static
from dash_app.plot_tab.callbacks_pyplotgen import _current_job
from dash_app.plot_tab.layout import _directory_case_selector, _plots_stores
from dash_app.shared import actions


class _Process:
    pid = 4321


def _find(component, component_id):
    if getattr(component, "id", None) == component_id:
        return component
    children = getattr(component, "children", None)
    if not isinstance(children, (list, tuple)):
        children = [children] if children is not None else []
    for child in children:
        found = _find(child, component_id)
        if found is not None:
            return found
    return None


def test_runtime_uses_fixed_unique_output(tmp_path, monkeypatch):
    input_a = tmp_path / "one"
    input_b = tmp_path / "two"
    input_a.mkdir()
    input_b.mkdir()
    output_root = tmp_path / "pyplotgen"
    captured = {}

    def fake_popen(command, **kwargs):
        captured["command"] = command
        captured["kwargs"] = kwargs
        return _Process()

    monkeypatch.setattr(pyplotgen_runtime, "PYPLOTGEN_OUTPUT_ROOT", output_root)
    monkeypatch.setattr(pyplotgen_runtime.subprocess, "Popen", fake_popen)

    job, _process = pyplotgen_runtime.start_pyplotgen([str(input_a), str(input_b)])

    output_dir = Path(job["output_dir"])
    assert output_dir.parent == output_root
    assert job["html_path"] == str(output_dir / "index.html")
    assert job["html_url"].endswith(f"/{output_dir.name}/index.html")
    assert captured["command"][:4] == [
        pyplotgen_runtime.sys.executable,
        "-u",
        str(pyplotgen_runtime.PYPLOTGEN_SCRIPT),
        "-l",
    ]
    assert captured["command"][4:7] == ["-c", str(input_a), str(input_b)]
    assert captured["kwargs"]["start_new_session"] is True
    assert captured["kwargs"]["env"]["PYTHONUNBUFFERED"] == "1"
    pyplotgen_runtime.release_pyplotgen(_process)


def test_progress_reader_returns_latest_complete_counter(tmp_path):
    log = tmp_path / "pyplotgen.log"
    log.write_bytes(
        b"Progress:    7 of   30 total .png panels complete\r"
        b"Progress:   21 of   44 total .png panels complete\r"
    )

    assert pyplotgen_runtime.read_pyplotgen_progress(str(log)) == (21, 44)


def test_stop_request_marks_stopping_before_signalling_process_group(monkeypatch):
    record = {"state": "running", "pid": 4321, "run_id": "run", "output_dir": "/tmp/gallery"}
    calls = []
    monkeypatch.setattr(actions, "broker_jobs", lambda: {"pyplotgen": record})
    monkeypatch.setattr(
        actions,
        "update_broker_job",
        lambda kind, **updates: calls.append(("update", kind, updates)) or {**record, **updates},
    )
    monkeypatch.setattr(actions, "stop_pyplotgen", lambda pid: calls.append(("stop", pid)))
    monkeypatch.setattr(actions, "_background", lambda *_args: None)
    monkeypatch.setattr(actions, "publish_event", lambda *_args, **_kwargs: None)

    result = actions.stop_pyplotgen_request()

    assert result["job"]["state"] == "stopping"
    assert calls == [
        ("update", "pyplotgen", {"state": "stopping"}),
        ("stop", 4321),
    ]


def test_terminal_broker_state_replaces_local_stopping_snapshot():
    broker = {"pyplotgen": {"run_id": "same", "state": "stopped", "returncode": -15}}
    action = {"kind": "stopping", "job": {"run_id": "same", "state": "stopping"}}

    assert _current_job(broker, action)["state"] == "stopped"


def test_static_route_serves_only_completed_fixed_exports(tmp_path, monkeypatch):
    folder = "20260814_120000_abcdef12"
    export = tmp_path / folder
    export.mkdir()
    (export / "index.html").write_text("gallery", encoding="utf-8")
    (export / "plot.png").write_bytes(b"image")
    monkeypatch.setattr(static, "PYPLOTGEN_OUTPUT_ROOT", tmp_path)
    app = Dash(__name__)
    app.layout = html.Div()
    static.register_pyplotgen_routes(app)
    client = app.server.test_client()

    assert client.get(f"/_clubb-pyplotgen/{folder}/index.html").data == b"gallery"
    assert client.get(f"/_clubb-pyplotgen/{folder}/plot.png").data == b"image"
    assert client.get("/_clubb-pyplotgen/../../etc/passwd").status_code == 404


def test_plot_header_contains_export_control_and_browser_state():
    initial = {"case_data": {}, "enabled_benchmarks": []}
    selector = _directory_case_selector(initial)
    stores = _plots_stores(
        {
            **initial,
            "plot_order": [],
            "plot_state": {},
            "next_id": 1,
            "selected_column": 0,
        }
    )

    assert _find(selector, "plots-pyplotgen-run") is not None
    assert {store.id for store in stores} >= {
        "plots-pyplotgen-action",
        "plots-pyplotgen-popup",
        "plots-pyplotgen-opened-run",
        "plots-pyplotgen-progress-interval",
    }
