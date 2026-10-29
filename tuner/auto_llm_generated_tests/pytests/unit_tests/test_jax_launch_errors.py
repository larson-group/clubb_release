"""Direct job-directory launches must persist invalid-request failures.

The managed launch path records terminal errors even without controller polling.
Exercise the equivalent direct tuner entry point before runtime/model startup.
"""

import json
import subprocess
import sys

import pytest

from tuner.job_runtime import TunerJob
from tuner.paths import REPO_ROOT


@pytest.mark.parametrize("raw_request", [
    {"backend": "jax", "jax_options": None},
    {"backend": "jax", "jax_options": ["cpu"]},
    None,
    [],
])
def test_direct_launch_records_invalid_request_without_controller(tmp_path, raw_request):
    job = TunerJob.create({"cases": ["bomex"]}, job_dir=tmp_path / "job")
    job.request_path.write_text(json.dumps(raw_request))
    result = subprocess.run(
        [sys.executable, "-m", "tuner.tune_clubb", "-job_dir", str(job.job_dir)],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode != 0
    for payload in (job.status(), job.results()):
        assert payload["state"] == "error", result.stderr
        assert payload["error_message"]
