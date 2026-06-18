#!/usr/bin/env python3
"""Python/f2py regression harness for the in-memory CLUBB loss driver.

Usage
-----
Normally run through run_scm_loss.py:

    python run_scripts/run_scm_loss.py -python -driver_test -cases bomex

run_scm_loss.py builds the loss-driver namelist, sets PYTHONPATH so the repo and
clubb_python_api package are visible, and passes the generated input file to
this module. For direct debugging:

    python -m tuner.clubb_loss_driver_test <namelist_path>

Like the normal standalone drivers, exit code 6 means success.

Understanding
-------------
This is the Python counterpart to src/clubb_loss_driver_test.F90. The Fortran
test exercises the native loss-driver API directly; this module runs the same
basic contract through the f2py layer used by the tuner.

The test starts with clubb_get_loss(), which is the one-shot path. That result
is the baseline.

It then uses the reusable path: init_clubb_loss(), clubb_get_loss_for_params(),
and finalize_clubb_loss(). With unchanged parameters, the reusable path must
match the one-shot baseline BFB. Repeated reusable calls must also stay BFB,
changing C8 must change the loss, and finalizing/reinitializing must return to
the same baseline.

Those checks are intentionally practical. They catch wrapper lifetime mistakes,
stale Fortran state, ignored parameter updates, array shape/order mismatches,
and Python/Fortran conversion bugs before they show up as confusing tuner
behavior.
"""

from __future__ import annotations

import sys

import numpy as np

from clubb_python import clubb_api
from tuner.clubb_loss_driver import (
    clubb_get_loss,
    clubb_get_loss_for_params,
    finalize_clubb_loss,
    init_clubb_loss,
    print_loss_matrix,
)
from tuner.taylor_metrics import LOSS_METRIC_NAMES


def assert_loss_outputs_match(
    test_name: str,
    test_var_names,
    test_metrics: tuple[np.ndarray, ...],
    reference_var_names,
    reference_metrics: tuple[np.ndarray, ...],
):
    """Require two loss-driver outputs to match bit-for-bit."""
    if list(test_var_names) != list(reference_var_names):
        raise RuntimeError(f"{test_name} failed: variable-name ordering changed")
    if len(test_metrics) != len(LOSS_METRIC_NAMES) or len(reference_metrics) != len(LOSS_METRIC_NAMES):
        raise RuntimeError(f"{test_name} failed: unexpected metric tuple length")
    for metric_name, test_metric, reference_metric in zip(LOSS_METRIC_NAMES, test_metrics, reference_metrics):
        if test_metric.shape != reference_metric.shape:
            raise RuntimeError(f"{test_name} failed: {metric_name} shapes changed")
        if not np.array_equal(test_metric, reference_metric):
            raise RuntimeError(f"{test_name} failed: {metric_name} changed")


def run_manual_loss_bfb_test(
    namelist_filename: str,
    reference_var_names,
    reference_metrics: tuple[np.ndarray, ...],
):
    """Repeat the reusable loss path twice with the same parameters and require BFB output."""
    manual_var_names, clubb_params_all = init_clubb_loss(
        namelist_filename,
        return_default_params=True,
    )
    try:
        manual_metrics_1 = clubb_get_loss_for_params(clubb_params_all)
        manual_metrics_2 = clubb_get_loss_for_params(clubb_params_all)
    finally:
        finalize_clubb_loss()

    assert_loss_outputs_match(
        "Manual loss BFB check (first rerun)",
        manual_var_names,
        manual_metrics_1,
        reference_var_names,
        reference_metrics,
    )
    assert_loss_outputs_match(
        "Manual loss BFB check (second rerun)",
        manual_var_names,
        manual_metrics_2,
        reference_var_names,
        reference_metrics,
    )
    for metric_name, manual_metric_1, manual_metric_2 in zip(LOSS_METRIC_NAMES, manual_metrics_1, manual_metrics_2):
        if not np.array_equal(manual_metric_1, manual_metric_2):
            raise RuntimeError(
                f"Manual loss BFB check failed: repeated manual {metric_name} calls are not identical"
            )

    print("Manual loss BFB check passed")


def run_tweaked_parameter_test(
    namelist_filename: str,
    reference_var_names,
    reference_scaled_rmse: np.ndarray,
):
    """Verify that perturbing C8 changes the loss."""
    tweaked_var_names, clubb_params_all = init_clubb_loss(
        namelist_filename,
        return_default_params=True,
    )
    try:
        if list(tweaked_var_names) != list(reference_var_names):
            raise RuntimeError("Tweaked-parameter test failed: variable-name ordering changed")

        param_names = clubb_api.get_param_names()
        try:
            c8_idx = param_names.index("C8")
        except ValueError as exc:
            raise RuntimeError("Tweaked-parameter test failed: could not locate C8 in param_names") from exc

        clubb_params_all = np.array(clubb_params_all, copy=True)
        clubb_params_all[0, c8_idx] = clubb_params_all[0, c8_idx] + 0.5
        tweaked_scaled_rmse, _, _, _, _ = clubb_get_loss_for_params(clubb_params_all)
    finally:
        finalize_clubb_loss()

    if np.array_equal(tweaked_scaled_rmse, reference_scaled_rmse):
        raise RuntimeError("Tweaked-parameter test failed: changing C8 did not change scaled_rmse")

    print("Tweaked-parameter test passed")


def run_reinitialize_loss_bfb_test(
    namelist_filename: str,
    reference_var_names,
    reference_metrics: tuple[np.ndarray, ...],
):
    """Reinitialize the reusable loss path and require it to reproduce the baseline."""
    reinit_var_names_1, clubb_params_all_1 = init_clubb_loss(
        namelist_filename,
        return_default_params=True,
    )
    try:
        reinit_metrics_1 = clubb_get_loss_for_params(clubb_params_all_1)
    finally:
        finalize_clubb_loss()

    reinit_var_names_2, clubb_params_all_2 = init_clubb_loss(
        namelist_filename,
        return_default_params=True,
    )
    try:
        reinit_metrics_2 = clubb_get_loss_for_params(clubb_params_all_2)
    finally:
        finalize_clubb_loss()

    assert_loss_outputs_match(
        "Reinitialize loss BFB check (first init)",
        reinit_var_names_1,
        reinit_metrics_1,
        reference_var_names,
        reference_metrics,
    )
    assert_loss_outputs_match(
        "Reinitialize loss BFB check (second init)",
        reinit_var_names_2,
        reinit_metrics_2,
        reference_var_names,
        reference_metrics,
    )

    if list(reinit_var_names_1) != list(reinit_var_names_2):
        raise RuntimeError("Reinitialize loss BFB check failed: variable names changed across init/finalize")
    if not np.array_equal(clubb_params_all_1, clubb_params_all_2):
        raise RuntimeError("Reinitialize loss BFB check failed: default parameter matrix changed across init/finalize")

    print("Reinitialize loss BFB check passed")


def main():
    namelist_filename = "clubb.in"
    if len(sys.argv) >= 2 and sys.argv[1] not in ("-h", "--help"):
        namelist_filename = sys.argv[1]
    elif len(sys.argv) >= 2:
        print("Usage: python -m tuner.clubb_loss_driver_test <namelist_path>")
        sys.exit(0)

    clubb_var_names, scaled_rmse, correlation, std_ratio, centered_rmse_norm, bias_norm = clubb_get_loss(namelist_filename)
    reference_metrics = (scaled_rmse, correlation, std_ratio, centered_rmse_norm, bias_norm)
    run_manual_loss_bfb_test(namelist_filename, clubb_var_names, reference_metrics)
    run_tweaked_parameter_test(namelist_filename, clubb_var_names, scaled_rmse)
    run_reinitialize_loss_bfb_test(namelist_filename, clubb_var_names, reference_metrics)
    print_loss_matrix(clubb_var_names, scaled_rmse, correlation, std_ratio, centered_rmse_norm, bias_norm)
    sys.exit(6)


if __name__ == "__main__":
    main()
