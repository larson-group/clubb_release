from dash_app.run_tab import runtime
from dash_app.run_tab.runtime import refresh_running_runtimes


def test_refresh_running_runtimes_advances_quiet_broker_cases():
    values = refresh_running_runtimes(
        {"finished": 3.0},
        {
            "arm": {"start_time": 40.0, "broker_managed": True},
            "bomex": {"start_time": 55.0},
        },
        now=100.0,
    )

    assert values == {"finished": 3.0, "arm": 60.0, "bomex": 45.0}


def test_pid_is_alive_uses_cross_process_probe(monkeypatch):
    calls = []
    monkeypatch.setattr(runtime.os, "kill", lambda pid, signal: calls.append((pid, signal)))

    assert runtime.pid_is_alive("123")
    assert calls == [(123, 0)]
    assert not runtime.pid_is_alive(0)


def test_pid_is_alive_handles_a_finished_process(monkeypatch):
    def missing(_pid, _signal):
        raise ProcessLookupError

    monkeypatch.setattr(runtime.os, "kill", missing)
    assert not runtime.pid_is_alive(123)
