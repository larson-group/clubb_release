#!/usr/bin/env python3
"""Verify the standalone CLUBB settings validator against compiled F2PY.

Run after ``./compile.py -python``.  This is intentionally a plain executable
test (rather than a Dash or pytest test) so the existing Jenkins Python job can
run it directly with the same compiled API it tests elsewhere.
"""

from __future__ import annotations

import math
import os
from pathlib import Path
import sys
from contextlib import contextmanager


ROOT = Path(__file__).resolve().parents[1]
_api_candidates = [
    ROOT / "install" / "latest" / "python",
    *ROOT.glob("build/*_PYTHON/clubb_python_api/f2py_runtime"),
]
_fresh_api_dir = max(
    (path for path in _api_candidates if any(path.glob("clubb_f2py*.so"))),
    key=lambda path: max(item.stat().st_mtime for item in path.glob("clubb_f2py*.so")),
    default=None,
)
for path in (ROOT, ROOT / "clubb_python_api", _fresh_api_dir):
    if path is None:
        continue
    if path.is_dir() and str(path) not in sys.path:
        sys.path.insert(0, str(path))

from utilities.clubb_settings_validation import (  # noqa: E402
    PARAMETER_NAMES,
    parameter_hard_bounds,
    validation_errors,
)


def _require_compiled_api():
    try:
        from clubb_python import clubb_api
        from clubb_python.derived_types.err_info import ErrInfo
    except ImportError as exc:  # pragma: no cover - a Jenkins setup failure
        raise RuntimeError(
            "This test requires ./compile.py -python first; no compiled CLUBB Python API was found."
        ) from exc
    return clubb_api, ErrInfo


def _assert(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def _assert_bounds_match(clubb_api) -> None:
    python_bounds = parameter_hard_bounds()
    nparams = len(clubb_api.get_param_names())
    fortran_bounds = clubb_api.get_parameter_hard_bounds(nparams)
    _assert([item["name"] for item in fortran_bounds] == list(PARAMETER_NAMES), "parameter ordering differs from Fortran")
    for item in fortran_bounds:
        expected = python_bounds[item["name"]]
        for endpoint in ("min", "max"):
            actual = item[endpoint]
            desired = expected[endpoint]
            if actual is None or desired is None:
                _assert(actual is desired, f"{item['name']} {endpoint}: {actual!r} != {desired!r}")
            else:
                _assert(math.isclose(actual, desired, rel_tol=0.0, abs_tol=0.0), f"{item['name']} {endpoint}: {actual!r} != {desired!r}")


def _baseline(clubb_api):
    params = clubb_api.init_clubb_params(1, iunit=10, filename="")
    flags = clubb_api.get_default_config_flags()
    names = clubb_api.get_param_names()
    values = {name: float(params[0, index]) for index, name in enumerate(names)}
    return params, flags, names, values


@contextmanager
def _silence_fortran_stdio():
    """Keep expected rejecting cases from obscuring the Jenkins result."""
    saved_stdout, saved_stderr = os.dup(1), os.dup(2)
    with open(os.devnull, "w", encoding="utf-8") as null:
        try:
            os.dup2(null.fileno(), 1)
            os.dup2(null.fileno(), 2)
            yield
        finally:
            os.dup2(saved_stdout, 1)
            os.dup2(saved_stderr, 2)
            os.close(saved_stdout)
            os.close(saved_stderr)


def _fortran_rejects(clubb_api, ErrInfo, params, flags, *, l_input_fields: bool = False) -> bool:
    with _silence_fortran_stdio():
        result = clubb_api.check_clubb_settings(
            ngrdcol=1,
            params=params,
            l_implemented=False,
            l_input_fields=l_input_fields,
            clubb_config_flags=flags,
            err_info=ErrInfo(ngrdcol=1),
        )
    return any(int(code) != 0 for code in result.err_code)


def _mutated_case(clubb_api, param_updates=None, flag_updates=None):
    params, flags, names, values = _baseline(clubb_api)
    indices = {name: index for index, name in enumerate(names)}
    for name, value in (param_updates or {}).items():
        params[0, indices[name]] = value
        values[name] = value
    flags = flags._replace(**(flag_updates or {}))
    return params, flags, values


def _assert_case_matches(clubb_api, ErrInfo, label: str, *, param_updates=None, flag_updates=None, l_input_fields=False) -> None:
    params, flags, values = _mutated_case(clubb_api, param_updates, flag_updates)
    python_rejects = bool(validation_errors(values, flags._asdict(), lmin=1.0, l_input_fields=l_input_fields))
    fortran_rejects = _fortran_rejects(clubb_api, ErrInfo, params, flags, l_input_fields=l_input_fields)
    _assert(python_rejects == fortran_rejects, f"{label}: Python rejects={python_rejects}, Fortran rejects={fortran_rejects}")


def run_test() -> None:
    clubb_api, ErrInfo = _require_compiled_api()
    _assert_bounds_match(clubb_api)

    # One representative case for each rejecting branch reachable through the
    # public F2PY settings API.  The static bounds comparison above covers all
    # 114 scalar endpoints; these cover cross-setting logic.
    _assert_case_matches(clubb_api, ErrInfo, "valid baseline")
    _assert_case_matches(clubb_api, ErrInfo, "hard upper bound", param_updates={"C_uu_shr": 1.01})
    _assert_case_matches(clubb_api, ErrInfo, "hard lower bound", param_updates={"C_uu_shr": -0.01})
    _assert_case_matches(clubb_api, ErrInfo, "linked parameter equality", param_updates={"C6rt": 2.0, "C6thl": 1.0})
    _assert_case_matches(clubb_api, ErrInfo, "EM damping relation", param_updates={"C1": 2.0, "C14": 1.0}, flag_updates={"l_damp_wp2_using_em": True})
    _assert_case_matches(clubb_api, ErrInfo, "EM damping stability flag", flag_updates={"l_damp_wp2_using_em": True, "l_stability_correct_tau_zm": True})
    _assert_case_matches(clubb_api, ErrInfo, "PDF type range", flag_updates={"iiPDF_type": 11})
    _assert_case_matches(clubb_api, ErrInfo, "PDF needs input fields", flag_updates={"iiPDF_type": 2})
    _assert_case_matches(clubb_api, ErrInfo, "input field PDF accepted", flag_updates={"iiPDF_type": 2}, l_input_fields=True)
    _assert_case_matches(clubb_api, ErrInfo, "saturation formula range", flag_updates={"saturation_formula": 5})
    _assert_case_matches(clubb_api, ErrInfo, "PDF placement range", flag_updates={"ipdf_call_placement": 3})
    _assert_case_matches(clubb_api, ErrInfo, "predicted wind PDF relation", flag_updates={"l_predict_upwp_vpwp": True, "iiPDF_type": 4})
    _assert_case_matches(clubb_api, ErrInfo, "xp2 clipping relation", flag_updates={"l_min_xp2_from_corr_wx": True, "l_enable_relaxed_clipping": True})
    print("PASS CLUBB Python/Fortran settings-validation parity")


if __name__ == "__main__":
    run_test()
