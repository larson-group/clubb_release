"""Argument-order checker tests."""

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


def test_argument_order_checker_finds_no_issues():
    audit = _load_audit_module()
    issues = audit.find_order_issues(REPO_ROOT)
    assert not issues, "\n".join(
        f"{issue.issue}: {issue.python_api_fn} -> {issue.pyf_symbol}"
        for issue in issues
    )
