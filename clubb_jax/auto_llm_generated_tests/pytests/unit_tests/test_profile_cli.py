"""Check profile option forwarding, capture discovery and frontend presentation.

These focused checks use synthetic captures and model calls. Real CLI captures
and browser/server transfers are exercised separately as workflow validation.
The clubb_jax "Auto-LLM-generated JAX pytests" stage runs these checks.
"""

import json
import os
from pathlib import Path
import shlex
import socket
import threading
import urllib.error
import urllib.request

import pytest

from clubb_jax import run_jax, view_profile


@pytest.mark.parametrize("options,backend", [
    ("profile", None), ("cpu,profile", "cpu"), ("gpu,profile", "gpu"),
    ("gpu,profile,prealloc_gpu_mem=false", "gpu"),
])
def test_profile_survives_runtime_argument_round_trip(options, backend):
    """Keep profiling and backend selection through the shared encoding boundary."""
    selection = run_jax.runtime_selection(options)
    assert selection["trace"] is True
    assert selection["profile"] == backend
    args = run_jax.runtime_arguments(options)
    values, remaining = run_jax.parse_launcher_args([*args, "case.in"])
    assert values["trace"] is True
    assert values["profile"] == backend
    assert remaining == ["case.in"]
    assert run_jax.runtime_arguments(options, scm=True)[0].startswith("-jax=")


@pytest.mark.parametrize("options", [
    "profile,profile", "cpu,profile,profile", "gpu,profile=xprof",
    "profile=perfetto", "profile,", "profile,unknown",
])
def test_invalid_profile_modifiers_fail_without_preparing_environment(options, monkeypatch):
    """Reject duplicate/unsupported modifiers before launching or installing anything."""
    monkeypatch.setattr(run_jax, "_prepare_environment", lambda *a: pytest.fail("unexpected setup"))
    with pytest.raises((run_jax.LauncherError, ValueError)):
        run_jax.main([f"-options={options}", "case.in"])


def test_profile_uses_inherited_backend_and_enables_only_child_environment(monkeypatch):
    """A profiling-only request preserves bare -jax selection and overrides trace=0."""
    monkeypatch.setenv("CLUBB_JAX_ACCELERATOR", "cpu")
    monkeypatch.setenv("CLUBB_JAX_TRACE", "0")
    configuration = run_jax.runtime_configuration("profile")
    assert configuration["accelerator"] == "cpu"
    assert configuration["environment"]["CLUBB_JAX_TRACE"] == "1"
    assert os.environ["CLUBB_JAX_TRACE"] == "0"


def test_commands_quote_capture_paths_without_starting_viewers(tmp_path, monkeypatch, capsys):
    """Printed shell commands preserve spaces/metacharacters and launch no process."""
    capture = tmp_path / "a capture; $literal"
    monkeypatch.setattr(view_profile.subprocess, "Popen", lambda *a, **k: pytest.fail("viewer started"))
    monkeypatch.setattr(view_profile.webbrowser, "open", lambda *a: pytest.fail("browser started"))
    view_profile.print_profile_commands(None)
    assert capsys.readouterr().out == ""

    view_profile.print_profile_commands(capture)
    lines = capsys.readouterr().out.splitlines()
    commands = [shlex.split(line.strip()) for line in lines if line.startswith("  ")]
    assert len(commands) == 2
    assert [command[1:3] for command in commands] == [["-viewer", "xprof"], ["-viewer", "perfetto"]]
    assert all(Path(command[-1]).resolve() == capture for command in commands)


def test_latest_capture_ignores_unfinished_and_interrupted_runs(tmp_path, monkeypatch):
    """Choose the newest completed export, while allowing explicit interrupted captures."""
    monkeypatch.setattr(view_profile, "TRACE_ROOT", tmp_path)
    for name, interrupted, timestamp in [("first", False, 1), ("latest", False, 2), ("failed", True, 3)]:
        capture = tmp_path / name
        capture.mkdir()
        info = capture / "trace_info.json"
        info.write_text(json.dumps({"interrupted": interrupted}))
        os.utime(info, (timestamp, timestamp))
    (tmp_path / "still_exporting").mkdir()

    assert view_profile._find_capture(None) == tmp_path / "latest"
    assert view_profile._find_capture(tmp_path / "failed") == tmp_path / "failed"


def test_missing_capture_fails_with_actionable_message(tmp_path, monkeypatch, capsys):
    """An absent default capture fails without creating environments or servers."""
    monkeypatch.setattr(view_profile, "TRACE_ROOT", tmp_path)
    assert view_profile.main([]) == 1
    assert "No completed JAX profiles" in capsys.readouterr().err
    assert list(tmp_path.iterdir()) == []


def test_perfetto_serves_one_capture_and_releases_its_port(tmp_path, monkeypatch):
    """Fetch the handoff page and exact trace, reject other files, then finish serving."""
    trace = tmp_path / "perfetto_trace.json.gz"
    contents = b"synthetic compressed trace bytes"
    trace.write_bytes(contents)
    ready = threading.Event()
    servers, urls, outcomes = [], [], []
    http_server = view_profile.HTTPServer

    def create_server(*args):
        """Bound the test server's idle waits so failed assertions can clean it up."""
        server = http_server(*args)
        server.timeout = 0.1
        servers.append(server)
        return server

    def open_browser(url):
        """Record the real listener URL instead of opening a browser during pytest."""
        urls.append(url)
        ready.set()

    def serve():
        """Record the actual viewer helper's completion or exception."""
        try:
            outcomes.append(view_profile._view_perfetto(tmp_path, 0))
        except BaseException as exc:
            outcomes.append(exc)

    monkeypatch.setattr(view_profile, "HTTPServer", create_server)
    monkeypatch.setattr(view_profile.webbrowser, "open", open_browser)
    thread = threading.Thread(target=serve, daemon=True)
    thread.start()
    try:
        assert ready.wait(5), outcomes
        with urllib.request.urlopen(urls[0], timeout=5) as response:
            assert response.headers.get_content_type() == "text/html"
            assert b"Open in Perfetto" in response.read()
        with pytest.raises(urllib.error.HTTPError) as error:
            urllib.request.urlopen(urls[0] + "unrelated_file", timeout=5)
        assert error.value.code == 404
        with urllib.request.urlopen(urls[0] + trace.name, timeout=5) as response:
            assert response.read() == contents
        thread.join(5)
        assert not thread.is_alive()
        assert outcomes == [0]
        with socket.socket() as probe:
            probe.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
            probe.bind(servers[0].server_address)
    finally:
        if servers:
            servers[0].trace_delivered = True
        thread.join(5)


def test_frontend_prints_cleanup_capture_after_success(monkeypatch, tmp_path, capsys):
    """Standalone presents the capture returned by cleanup after reporting completion."""
    from clubb_jax.src import clubb_standalone

    capture = tmp_path / "trace_case"
    state = {"ngrdcol": 1, "total_param_sets": 1, "ifinal": 2}
    monkeypatch.setattr(clubb_standalone, "init_clubb_case", lambda *a: state)
    monkeypatch.setattr(clubb_standalone, "set_case_initial_conditions", lambda *a, **k: None)
    monkeypatch.setattr(clubb_standalone, "advance_clubb_to_end", lambda *a, **k: None)
    monkeypatch.setattr(clubb_standalone, "clean_up_clubb", lambda *a: capture)

    clubb_standalone.run_case("synthetic.in")
    text = capsys.readouterr().out
    assert text.index("Completed 2 timesteps") < text.index("JAX profile saved in:")
    assert "-viewer xprof" in text and "-viewer perfetto" in text


def test_cleanup_returns_capture_without_viewer_presentation(monkeypatch, tmp_path, capsys):
    """Driver cleanup passes the timer's result back without printing viewer commands."""
    from clubb_jax.src import clubb_case_initalization as initialization
    from clubb_jax.src.Microphys import microphys_init_cleanup

    capture = tmp_path / "trace_case"
    monkeypatch.setattr(initialization, "timer_start", lambda *a: None)
    monkeypatch.setattr(initialization, "timer_stop", lambda *a: None)
    monkeypatch.setattr(initialization, "timer_finalize", lambda *a, **k: capture)
    monkeypatch.setattr(microphys_init_cleanup, "cleanup_microphys", lambda: None)
    state = {"l_stats": False, "l_t_dependent": False, "output_file_prefix": str(tmp_path / "case")}

    assert initialization.clean_up_clubb(state) == capture
    assert "view_profile.py" not in capsys.readouterr().out


def test_one_shot_loss_preserves_metrics_and_presents_capture(monkeypatch, tmp_path, capsys):
    """Keep the public metric tuple while the loss driver presents its cleanup result."""
    from clubb_jax.src import clubb_loss_driver

    capture = tmp_path / "trace_loss"
    names, params, metrics = ["thlm"], object(), tuple(object() for _ in range(5))
    monkeypatch.setattr(clubb_loss_driver, "init_clubb_loss", lambda *a: (names, params))
    monkeypatch.setattr(clubb_loss_driver, "clubb_get_loss_for_params", lambda *a: metrics)
    monkeypatch.setattr(clubb_loss_driver, "finalize_clubb_loss", lambda: capture)

    result = clubb_loss_driver.clubb_get_loss("synthetic.in")
    assert result == (names, *metrics)
    assert "-viewer xprof" in capsys.readouterr().out


def test_reusable_loss_finalization_returns_capture_and_clears_state(monkeypatch, tmp_path, capsys):
    """A reusable loss API caller receives the path without presentation side effects."""
    from clubb_jax.src import clubb_loss_driver

    capture = tmp_path / "trace_loss"
    request = clubb_loss_driver.loss_request_type(state={"synthetic": True})
    monkeypatch.setattr(clubb_loss_driver, "active_request", request)
    monkeypatch.setattr(clubb_loss_driver, "clean_up_clubb", lambda *a: capture)

    assert clubb_loss_driver.finalize_clubb_loss() == capture
    assert not clubb_loss_driver.active_request.state
    assert clubb_loss_driver.finalize_clubb_loss() is None
    assert capsys.readouterr().out == ""
