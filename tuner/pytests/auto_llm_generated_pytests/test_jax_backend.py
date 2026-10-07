"""Check backend selection through the existing managed tuner workflow.

Request, worker and launcher contracts follow the shared Python tuner path
used by the native checks in jenkins_tests/clubb_tuner/Jenkinsfile. The native
backend calls src/clubb_loss_driver.F90; the JAX backend calls its port through
the same scheduler. These focused tests isolate forwarding, validation and
launch-error reporting with mocked workers instead of advancing model cases.
"""

import importlib
import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import pytest

from run_scripts import run_tuner_job as cli
from tuner.job_runtime import TunerJob
from tuner.request import load_request
from tuner.tuning_scheduler import TuningScheduler


def request_args(extra=()):
    return cli.parse_args([
        "-cases", "bomex",
        "-fields", "wp2",
        "-param_ranges", "C8:0.3:0.5",
        "-strategy", "random:3",
        "-batch_size", "2",
        *extra,
    ])


@pytest.mark.parametrize(
    "extra,backend,options",
    [
        ([], "fortran", None),
        (["-jax"], "jax", "cpu"),
        (["-jax=gpu"], "jax", "gpu"),
    ],
)
def test_cli_and_request_backend(tmp_path, extra, backend, options):
    request = cli.build_request(request_args(extra))
    path = tmp_path / "request.json"
    path.write_text(json.dumps(request))
    normalized = load_request(path)
    assert normalized["backend"] == backend
    assert normalized.get("jax_options") == options


@pytest.mark.parametrize(
    "backend,options",
    [("unknown", "cpu"), ("jax", "other")],
)
def test_invalid_backend_configuration(tmp_path, backend, options):
    request = cli.build_request(request_args())
    request.update(backend=backend, jax_options=options)
    path = tmp_path / "request.json"
    path.write_text(json.dumps(request))
    with pytest.raises((ValueError, RuntimeError), match="backend|cpu or gpu"):
        load_request(path)


@pytest.mark.parametrize("backend", ["fortran", "jax"])
def test_job_selects_managed_runtime(tmp_path, monkeypatch, backend):
    commands = []

    def popen(command, **kwargs):
        commands.append(command)
        return SimpleNamespace(pid=123)

    monkeypatch.setattr("tuner.job_runtime.subprocess.Popen", popen)
    job = TunerJob.create(
        {"cases": ["bomex"], "backend": backend, "jax_options": "gpu"},
        job_dir=tmp_path / "job",
    )
    job.start()
    command = commands[0]
    if backend == "jax":
        assert Path(command[1]).name == "run_jax.py"
        assert "-options=gpu" in command
        assert "-module=tuner.tune_clubb" in command
    else:
        assert command[1:3] == ["-m", "tuner.tune_clubb"]


def test_worker_payload_retains_backend(tmp_path, monkeypatch):
    request = cli.build_request(request_args(["-jax"]))
    path = tmp_path / "request.json"
    path.write_text(json.dumps(request))
    scheduler = TuningScheduler(
        request=load_request(path),
        job_dir=tmp_path,
        control_path=tmp_path / "control.json",
        status_path=tmp_path / "status.json",
        results_path=tmp_path / "results.json",
    )
    payloads = []

    class Process:
        pid = 123

        def __init__(self, *, target, args):
            payloads.append(args[1])

        def start(self):
            pass

    conn = SimpleNamespace(close=lambda: None)
    scheduler.ctx = SimpleNamespace(Pipe=lambda: (conn, conn), Process=Process)
    scheduler._start_worker("bomex")
    assert payloads[0]["backend"] == "jax"


@pytest.mark.parametrize('backend', [[], ['-jax=gpu']])
def test_top_result_reruns_keep_backend_and_physics_override(tmp_path, monkeypatch, backend):
    args = request_args([
        *backend, "-output_run_dir", str(tmp_path),
        "-override", "l_diag_Lscale_from_tau=.true.",
    ])
    request = cli.build_request(args)
    result_path = tmp_path / "results.json"
    result_path.write_text(json.dumps({
        "request": request,
        "best_results": [{"selected_params": {"C8": 0.4}}],
    }))
    monkeypatch.setattr(
        cli, "write_top_params_file", lambda *args, **kwargs: tmp_path / "top.in",
    )
    commands = []

    def run_command(command, **kwargs):
        commands.append(command)
        return 0

    monkeypatch.setattr(cli, "run_command", run_command)
    assert cli.run_top_results(
        args, tmp_path / "request.json", result_path, "both",
    ) == 0
    assert len(commands) == 2
    for command in commands:
        assert ("-jax=gpu" in command) == bool(backend)
        assert command[command.index('-override') + 1] == request['override']


def test_invalid_jax_options_fail_before_launch_with_terminal_results(tmp_path, monkeypatch):
    def reject_launch(*args, **kwargs):
        raise AssertionError('Invalid options must be checked before starting a process')

    monkeypatch.setattr('tuner.job_runtime.subprocess.Popen', reject_launch)
    job = TunerJob.create(
        {'cases': ['bomex'], 'backend': 'jax', 'jax_options': 'other'},
        job_dir=tmp_path / 'job',
    )
    with pytest.raises(ValueError, match='cpu or gpu'):
        job.start()
    for payload in (job.status(), job.results()):
        assert payload['state'] == 'error'
        assert 'cpu or gpu' in payload['error_message']


def test_managed_launcher_records_errors_without_controller_polling(tmp_path):
    job = TunerJob.create({'cases': ['bomex'], 'backend': 'jax'}, job_dir=tmp_path / 'job')
    launcher = Path(cli.__file__).resolve().parents[1] / 'clubb_jax' / 'run_jax.py'
    result = subprocess.run(
        [sys.executable, str(launcher), '-options=other', '-module=tuner.tune_clubb',
         '-job_dir', str(job.job_dir)],
        capture_output=True, text=True, timeout=30,
    )
    assert result.returncode == 1
    for payload in (job.status(), job.results()):
        assert payload['state'] == 'error'
        assert 'got: other' in payload['error_message']


@pytest.mark.parametrize('state', ['created', 'running', 'finished', 'stopped', 'error'])
def test_process_exit_preserves_terminal_state_and_results(tmp_path, state):
    job = TunerJob.create(
        {'cases': ['bomex'], 'backend': 'jax'}, job_dir=tmp_path / 'job', initial_state=state,
    )
    result = job.results()
    result['best_results'] = [{'selected_params': {'C8': .4}, 'total_loss': .2}]
    job.results_path.write_text(json.dumps(result))
    job.proc = SimpleNamespace(poll=lambda: 1)
    assert job.poll() == 1
    expected = state if state in {'finished', 'stopped', 'error'} else 'error'
    assert job.status()['state'] == job.results()['state'] == expected
    assert job.results()['best_results'] == result['best_results']
