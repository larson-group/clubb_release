"""Regression checks for bounded multi-column Tune result replays."""

from pathlib import Path

from dash_app.tune_tab import runtime


def test_loss_replay_uses_sequential_safe_parameter_batches(tmp_path, monkeypatch):
    """A 16-row replay is two eight-column calls, never one unsafe 16-column call."""

    started = []

    class FinishedProcess:
        def __init__(self, command, **_kwargs):
            self.command = command
            self.pid = 1000 + len(started)
            self._terminated = False
            started.append(self)

        def poll(self):
            return 143 if self._terminated else 0

        def terminate(self):
            self._terminated = True

    monkeypatch.setattr(runtime, "TUNE_RESULT_OUTPUT_ROOT", Path(tmp_path) / "output" / "tuner")
    monkeypatch.setattr(runtime.subprocess, "Popen", FinishedProcess)
    runtime.LOSS_RUN_PROCS.clear()

    run = runtime.start_loss_run(
        ["dycoms2_rf01"],
        ["wp3"],
        [{"C1": 1.0}] * 16,
        run_mode="window",
        batch_size=8,
        workspace_id="arm-bomex-test",
        revision_id="rev2",
        workspace_name="Arm BOMEX test",
    )

    assert run["batch_count"] == 2
    assert [batch["param_count"] for batch in run["batches"]] == [8, 8]
    result_dir = Path(tmp_path) / "output" / "tuner" / "Arm_BOMEX_test_rev2_loss_window"
    assert Path(run["output_dir"]) == result_dir
    assert Path(run["batches"][0]["output_dir"]) == result_dir
    assert Path(run["batches"][1]["output_dir"]) == result_dir / "batch_002"
    assert len(started) == 1

    updated, any_running = runtime.poll_loss_runs({"window": run})
    assert any_running is True
    assert updated["window"]["active_batch"] == 2
    assert len(started) == 2

    updated, any_running = runtime.poll_loss_runs(updated)
    assert any_running is False
    assert updated["window"]["state"] == "success"
