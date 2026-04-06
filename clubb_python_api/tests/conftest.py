"""Pytest helpers for the F2PY test suite."""

from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys
import tempfile

import pytest



@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item, call):
    """Store per-phase test outcomes on the item for fixture teardown."""
    outcome = yield
    report = outcome.get_result()
    setattr(item, f"rep_{report.when}", report)


@pytest.fixture(autouse=True)
def silence_fortran_stdio(request):
    """Suppress low-level stdio noise unless the current test fails.

    Several tests intentionally trigger CLUBB fatal-error paths and then
    validate the returned error code. The underlying Fortran writes directly
    to stdout/stderr, which produces noisy suite output even for passing tests.
    Redirect file descriptors 1 and 2 for the duration of each test and replay
    the buffered output only when that test fails.
    """
    sys.stdout.flush()
    sys.stderr.flush()
    original_stdout_fd = os.dup(1)
    original_stderr_fd = os.dup(2)
    with tempfile.TemporaryFile(mode="w+b") as capture_out, tempfile.TemporaryFile(mode="w+b") as capture_err:
        os.dup2(capture_out.fileno(), 1)
        os.dup2(capture_err.fileno(), 2)
        try:
            yield
        finally:
            sys.stdout.flush()
            sys.stderr.flush()
            os.dup2(original_stdout_fd, 1)
            os.dup2(original_stderr_fd, 2)
            os.close(original_stdout_fd)
            os.close(original_stderr_fd)
            report = getattr(request.node, "rep_call", None)
            if report is not None and report.failed:
                capture_out.seek(0)
                stdout_data = capture_out.read()
                if stdout_data:
                    sys.stdout.write(stdout_data.decode("utf-8", errors="replace"))
                    sys.stdout.flush()
                capture_err.seek(0)
                stderr_data = capture_err.read()
                if stderr_data:
                    sys.stderr.write(stderr_data.decode("utf-8", errors="replace"))
                    sys.stderr.flush()


@pytest.fixture
def run_quiet_python():
    """Run a short Python snippet in a subprocess and capture its stdio.

    This is useful for tests that intentionally trigger Fortran fatal-error
    diagnostics but only need to assert on the returned Python-visible state.
    """
    repo_root = Path(__file__).resolve().parents[2]
    api_root = Path(__file__).resolve().parents[1]

    def _run(script: str) -> subprocess.CompletedProcess[str]:
        env = os.environ.copy()
        env["PYTHONPATH"] = f"{repo_root}:{api_root}{os.pathsep}{env['PYTHONPATH']}" if "PYTHONPATH" in env else f"{repo_root}:{api_root}"
        return subprocess.run(
            [sys.executable, "-c", script],
            cwd=api_root,
            env=env,
            capture_output=True,
            text=True,
            check=False,
        )

    return _run
