"""Exercise persistent Dash worker transport and ownership without model runs.

The existing clubb_dash pytest stage owns these provisional checks. Calls use
standard-library functions so persistence, task errors and shutdown are checked
independently of JAX compilation or CLUBB case setup.
"""
import operator
import os
from pathlib import Path
import subprocess
import sys
import time

import pytest

from dash_app.workers import RuntimeExecutor


def test_worker_retains_process_after_calls_and_reports_task_errors():
    """Repeated calls and an intervening exception reuse one owned process."""
    worker = RuntimeExecutor()
    try:
        pid = worker.submit(os.getpid).result(timeout=10)
        assert pid != os.getpid()
        assert worker.submit(os.getpid).result(timeout=10) == pid
        with pytest.raises(RuntimeError, match='ZeroDivisionError'):
            worker.submit(operator.truediv, 1, 0).result(timeout=10)
        assert worker.submit(os.getpid).result(timeout=10) == pid
    finally:
        worker.terminate_workers()
    with pytest.raises(ProcessLookupError):
        os.kill(pid, 0)
    with pytest.raises(RuntimeError, match='stopping'):
        worker.submit(os.getpid)


def test_shutdown_stops_an_active_worker_and_resolves_its_future():
    """Terminate an active call, reap its process and settle the pending result."""
    worker = RuntimeExecutor()
    try:
        pid = worker.submit(os.getpid).result(timeout=10)
        future = worker.submit(time.sleep, 300)
        deadline = time.monotonic() + 5
        while not future.running() and time.monotonic() < deadline:
            time.sleep(.01)
        assert future.running()
        worker.terminate_workers()
        assert future.done()
        with pytest.raises(RuntimeError, match='worker exited'):
            future.result()
        with pytest.raises(ProcessLookupError):
            os.kill(pid, 0)
    finally:
        worker.terminate_workers()


def test_maps_and_sweeps_retain_separate_workers_with_shared_admission():
    """Map/sweep sessions stay distinct while only one scientific task runs."""
    from concurrent.futures import Future
    from dash_app.plot_tab.tasks import PlotTasks

    pools = []

    class Pool:
        def __init__(self):
            """Collect requests submitted to this simulated runtime pool."""
            self.calls = []
            pools.append(self)

        def submit(self, fn, request):
            """Retain a pending Future until the test completes the request."""
            future = Future()
            self.calls.append((request, future))
            return future

        def shutdown(self, **kwargs):
            """Allow broker cleanup without starting a process."""
            pass

    tasks = PlotTasks(Pool)
    try:
        mapping = dict(scope='page:map', revision=1, task='gradient_matrix', args=[])
        sweep = dict(scope='page:sweep', revision=1, task='gradient_matrix_response', args=[])
        tasks.submit(mapping)
        tasks.submit(sweep)
        assert len(pools) == 1
        assert tasks.poll('page:sweep', 1)['state'] == 'queued'
        pools[0].calls[0][1].set_result('{}')
        assert len(pools) == 2
        tasks.submit(dict(mapping, revision=2))
        assert len(pools[0].calls) == 1
        pools[1].calls[0][1].set_result('{}')
        assert len(pools) == 2
        assert len(pools[0].calls) == 2
        assert tasks.poll('page:map', 2)['state'] == 'running'
        pools[0].calls[1][1].set_result('{}')
        tasks.submit(dict(sweep, scope='page:points', task='gradient_matrix_response_gradients'))
        assert len(pools) == 2 and len(pools[1].calls) == 2
        pools[1].calls[1][1].set_result('{}')
    finally:
        tasks.close()


def test_backend_edits_invalidate_runtime_but_tests_do_not(tmp_path):
    """Retained workers refresh on JAX source edits while test edits stay harmless."""
    from dash_app.shared.provenance import runtime_source_fingerprint

    backend = tmp_path / 'clubb_jax'
    backend.mkdir()
    source = backend / 'gradient_response.py'
    source.write_text('version = 1\n')
    first = runtime_source_fingerprint(tmp_path)
    source.write_text('version = 2\n')
    second = runtime_source_fingerprint(tmp_path)
    assert first != second
    tests = backend / 'auto_llm_generated_tests'
    tests.mkdir()
    (tests / 'test_explorer.py').write_text('version = 3\n')
    assert runtime_source_fingerprint(tmp_path) == second


def test_explorer_services_import_without_dashboard_or_model_dependencies():
    """Pickled explorer tasks must import in the numerical worker environment."""
    script = """
import sys
from dash_app.services.gradient_matrix import calculate_matrix
from dash_app.services.gradient_response import execute_gradient_request
assert not any(name in sys.modules for name in ('dash', 'plotly', 'pydantic', 'jax', 'numpy'))
try:
    execute_gradient_request({'task': 'invalid'})
except ValueError as exc:
    assert str(exc) == 'Unknown gradient task'
else:
    raise AssertionError('Unknown task was accepted')
"""
    result = subprocess.run([sys.executable, '-S', '-c', script],
                            cwd=Path(__file__).resolve().parents[4],
                            capture_output=True, text=True, timeout=10)
    assert result.returncode == 0, result.stdout + result.stderr
