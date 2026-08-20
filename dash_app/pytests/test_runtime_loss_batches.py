"""Regression checks for unified multi-column Tune result replays."""

from pathlib import Path

from dash_app.tune_tab import runtime


def test_loss_replay_uses_one_multicol_run_for_every_parameter_row(tmp_path, monkeypatch):

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
        workspace_id="arm-bomex-test",
        revision_id="rev2",
        workspace_name="Arm BOMEX test",
    )

    assert run["column_count"] == 16
    result_dir = Path(tmp_path) / "output" / "tuner" / "Arm_BOMEX_test_rev2_loss_window"
    assert Path(run["output_dir"]) == result_dir
    assert len(started) == 1
    params_text = Path(run["params_path"]).read_text(encoding="utf-8")
    assert "ngrdcol = 16" in params_text
    assert "batch_size = 16" in params_text
    assert "batch_002" not in params_text

    updated, any_running = runtime.poll_loss_runs({"window": run})
    assert any_running is False
    assert updated["window"]["state"] == "success"
