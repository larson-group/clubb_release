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

    blocking_issues = [
        issue for issue in issues
        if issue.issue.split("::", 1)[0] != "GENERIC_DISPATCH_WRAPPER"
    ]
    generic_dispatch = [
        issue for issue in issues
        if issue.issue.startswith("GENERIC_DISPATCH_WRAPPER")
    ]
    assert generic_dispatch
    assert any(
        issue.issue.startswith("GENERIC_DISPATCH_WRAPPER")
        and issue.python_api_fn.endswith("matrix_solver_wrapper.py::tridiag_solve")
        for issue in issues
    )
    assert not any(
        issue.pyf_symbol == "f2py_stats_finalize"
        and issue.python_api_fn.endswith("stats_netcdf.py::finalize_stats")
        and issue.issue.startswith("PUBLIC_API_SIGNATURE_")
        for issue in issues
    )
    assert not any(
        issue.pyf_symbol.startswith("f2py_lapack_")
        for issue in issues
    )
    assert not any(
        issue.issue == "MISSING_FORTRAN_SOURCE_MAPPING"
        and issue.pyf_symbol in {
            "get_err_code",
            "f2py_get_param_names",
            "f2py_get_stats_config",
            "f2py_get_stats_var_meta",
            "f2py_get_stats_var_data",
            "f2py_set_simplified_radiation_params",
            "f2py_reset_err_code",
            "f2py_initialize_error_headers",
        }
        for issue in issues
    )

    leaked_internal_args = [
        arg
        for issue in issues
        if issue.issue == "PUBLIC_API_FORBIDDEN_INTERNAL_ARG"
        for arg in issue.python_order
    ]
    assert "sclr_dim_transport" not in leaked_internal_args
    assert "edsclr_dim_transport" not in leaked_internal_args
    assert not blocking_issues, "\n".join(
        f"{issue.issue}: {issue.python_api_fn} -> {issue.pyf_symbol}"
        for issue in blocking_issues
    )
