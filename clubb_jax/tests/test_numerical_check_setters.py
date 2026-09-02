#!/usr/bin/env python3
"""test_numerical_check_setters.py — pin the check_nan / check_negative err_code setters.

`check_nan(arr, name, location, err_code)` and `check_negative(arr, name, location, err_code)` (numerical_check.py ↔
numerical_check.F90) are the low-level validity guards: they set `err_code` to CLUBB_FATAL_ERROR (in place) iff `arr`
contains a NaN/Inf (resp. a negative value), else leave it untouched. They are exercised on the live path only when
something is already wrong (so passing cases never trip them), and were the last untested numerical_check routines
(sfc_varnce_check/length_check/pdf_closure_check/rad_check are in test_validation_checks.py). This pins both
directions — clean input leaves err_code unchanged; a bad value sets CLUBB_FATAL_ERROR — so a regression that silently
disables the guard is caught. Pure-Python; never SKIPs. (iter 573)
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import numpy as np

from clubb_jax.src.CLUBB_core.numerical_check import check_nan, check_negative, CLUBB_FATAL_ERROR


def _ec():
    return np.array([0], dtype=np.int64)


def test_check_negative():
    ec = _ec()
    check_negative(np.array([0.0, 1.0, 5.0, 1e-9]), "x", "loc", ec)
    assert ec[0] == 0, "all-nonnegative wrongly flagged"
    ec = _ec()
    check_negative(np.array([0.0, 1.0, -1e-12, 3.0]), "x", "loc", ec)
    assert ec[0] == CLUBB_FATAL_ERROR, "negative value not flagged"
    print("  check_negative: clean→unchanged, negative→CLUBB_FATAL_ERROR  PASS")


def test_check_nan():
    ec = _ec()
    check_nan(np.array([0.0, -3.0, 5.0, 1e10]), "y", "loc", ec)
    assert ec[0] == 0, "all-finite wrongly flagged"
    for bad in (np.nan, np.inf, -np.inf):
        ec = _ec()
        check_nan(np.array([1.0, 2.0, bad]), "y", "loc", ec)
        assert ec[0] == CLUBB_FATAL_ERROR, f"non-finite {bad} not flagged"
    print("  check_nan: finite→unchanged, NaN/±Inf→CLUBB_FATAL_ERROR  PASS")


def main():
    print("test_numerical_check_setters:")
    test_check_negative()
    test_check_nan()
    print("All check_nan/check_negative checks PASSED")


if __name__ == "__main__":
    main()
