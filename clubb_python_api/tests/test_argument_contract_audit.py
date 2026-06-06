"""Argument-contract audit tests."""

from __future__ import annotations

from pathlib import Path
import subprocess
import sys


REPO_ROOT = Path(__file__).resolve().parents[2]
ENFORCER_DIR = REPO_ROOT / "clubb_python_api" / "tests" / "argument_list_enforcer"
AUDIT_SCRIPT = ENFORCER_DIR / "argument_contract_audit.py"


def test_argument_contract_audit_cli_passes():
    result = subprocess.run(
        [sys.executable, str(AUDIT_SCRIPT), "--repo-root", str(REPO_ROOT)],
        check=False,
        text=True,
        capture_output=True,
    )
    assert result.returncode == 0, (
        result.stdout
        + ("\nSTDERR:\n" + result.stderr if result.stderr else "")
    )
