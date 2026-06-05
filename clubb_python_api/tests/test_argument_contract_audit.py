"""Argument-contract audit tests."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[2]
ENFORCER_DIR = REPO_ROOT / "clubb_python_api" / "tests" / "argument_list_enforcer"
AUDIT_SCRIPT = ENFORCER_DIR / "argument_contract_audit.py"


def _load_audit_module():
    spec = importlib.util.spec_from_file_location("argument_contract_audit", AUDIT_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError("Failed to load argument_contract_audit module")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_argument_contract_audit_reports_current_drift():
    audit = _load_audit_module()
    issues = audit.find_contract_issues(REPO_ROOT)
    assert issues

    public_mismatches = [
        issue for issue in issues
        if issue.issue.startswith("PUBLIC_API_SIGNATURE_MISMATCH")
    ]
    assert public_mismatches
    assert any(
        issue.pyf_symbol == "f2py_stats_accumulate"
        and issue.python_api_fn.endswith("stats_clubb_utilities.py::stats_accumulate")
        for issue in public_mismatches
    )

    leaked_internal_args = [
        arg
        for issue in issues
        if issue.issue == "PUBLIC_API_FORBIDDEN_INTERNAL_ARG"
        for arg in issue.python_order
    ]
    assert "sclr_dim_transport" not in leaked_internal_args
    assert "edsclr_dim_transport" not in leaked_internal_args
