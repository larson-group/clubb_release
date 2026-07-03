#!/usr/bin/env python3
"""Regression checks for Python loss metric helpers.

The helper functions in utilities.loss_metrics intentionally mirror
src/clubb_loss_driver.F90.  These small checks cover the Taylor-metric edge
cases that matter for flat or nearly flat benchmark profiles, where ordinary
correlation and standard-deviation ratios are not mathematically defined but
the Fortran loss driver defines neutral Taylor components.
"""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


TESTS_DIR = Path(__file__).resolve().parent
CLUBB_ROOT = TESTS_DIR.parent
if str(CLUBB_ROOT) not in sys.path:
    sys.path.insert(0, str(CLUBB_ROOT))

from utilities.loss_metrics import calculate_taylor_metrics  # noqa: E402


def assert_close(label: str, actual: float, expected: float) -> None:
    if not np.isclose(actual, expected, rtol=0.0, atol=1.0e-12):
        raise AssertionError(f"{label}: expected {expected}, got {actual}")


def run_test() -> None:
    metrics = calculate_taylor_metrics(
        np.array([2.0, 2.0, 2.0]),
        np.array([2.0, 2.0, 2.0]),
    )
    assert_close("flat correlation", metrics["correlation"], 1.0)
    assert_close("flat std_ratio", metrics["std_ratio"], 1.0)
    assert_close("flat centered_rmse_norm", metrics["centered_rmse_norm"], 0.0)
    assert_close("flat bias_norm", metrics["bias_norm"], 0.0)

    metrics = calculate_taylor_metrics(
        np.array([1.0, 2.0, 3.0]),
        np.array([2.0, 2.0, 2.0]),
    )
    assert_close("flat benchmark correlation", metrics["correlation"], 1.0)
    assert_close("flat benchmark std_ratio", metrics["std_ratio"], 1.0)
    assert_close(
        "flat benchmark centered_rmse_norm",
        metrics["centered_rmse_norm"],
        np.sqrt(2.0 / 3.0),
    )
    assert_close("flat benchmark bias_norm", metrics["bias_norm"], 0.0)

    print("PASS loss metrics regression test")


if __name__ == "__main__":
    run_test()
