#!/usr/bin/env python3
"""Python front-end for the reusable CLUBB loss driver."""

from __future__ import annotations

import sys

import numpy as np

from clubb_python import clubb_api


def init_clubb_loss(namelist_path: str, return_default_params: bool = False):
    """Initialize the reusable Fortran loss driver."""
    return clubb_api.init_clubb_loss(namelist_path, return_default_params=return_default_params)


def clubb_get_loss_for_params(clubb_params_all: np.ndarray):
    """Score one full parameter matrix with the reusable Fortran loss driver."""
    return clubb_api.clubb_get_loss_for_params(clubb_params_all)


def finalize_clubb_loss():
    """Release the active reusable Fortran loss session."""
    clubb_api.finalize_clubb_loss()


def clubb_get_loss(namelist_path: str):
    """Run the one-shot loss path using the default parameter matrix."""
    return clubb_api.clubb_get_loss(namelist_path)


def print_loss_matrix(
    clubb_var_names,
    scaled_rmse: np.ndarray,
    correlation: np.ndarray,
    std_ratio: np.ndarray,
    centered_rmse_norm: np.ndarray,
    bias_norm: np.ndarray,
):
    """Print the loss table in the same row-oriented format as the Fortran tools."""
    header = ["variable"]
    if scaled_rmse.shape[0] > 1:
        header.append("window")
    header.extend(f"scaled_rmse_col{i + 1}" for i in range(scaled_rmse.shape[2]))
    header.extend(f"corr_col{i + 1}" for i in range(scaled_rmse.shape[2]))
    header.extend(f"std_ratio_col{i + 1}" for i in range(scaled_rmse.shape[2]))
    header.extend(f"crmse_norm_col{i + 1}" for i in range(scaled_rmse.shape[2]))
    header.extend(f"bias_norm_col{i + 1}" for i in range(scaled_rmse.shape[2]))
    header.append("best_col")
    print(" ".join(header))

    for window_idx in range(scaled_rmse.shape[0]):
        for row_idx, var_name in enumerate(clubb_var_names):
            best_col_idx = int(np.argmin(scaled_rmse[window_idx, row_idx, :])) + 1
            row = [str(var_name).strip()]
            if scaled_rmse.shape[0] > 1:
                row.append(f"window{window_idx + 1}")
            row.extend(f"{value:16.8E}" for value in scaled_rmse[window_idx, row_idx, :])
            row.extend(f"{value:16.8E}" for value in correlation[window_idx, row_idx, :])
            row.extend(f"{value:16.8E}" for value in std_ratio[window_idx, row_idx, :])
            row.extend(f"{value:16.8E}" for value in centered_rmse_norm[window_idx, row_idx, :])
            row.extend(f"{value:16.8E}" for value in bias_norm[window_idx, row_idx, :])
            row.append(f"col{best_col_idx}")
            print(" ".join(row))


def main():
    if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help"):
        print("Usage: python -m tuner.clubb_loss_driver <namelist_path>")
        sys.exit(0 if "--help" in sys.argv else 1)

    namelist_path = sys.argv[1]
    clubb_var_names, scaled_rmse, correlation, std_ratio, centered_rmse_norm, bias_norm = clubb_get_loss(namelist_path)
    print_loss_matrix(clubb_var_names, scaled_rmse, correlation, std_ratio, centered_rmse_norm, bias_norm)
    sys.exit(6)


if __name__ == "__main__":
    main()
