#!/usr/bin/env python3

import argparse
import subprocess
import sys
from pathlib import Path

# --------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------

SCRIPT_DIR = Path(__file__).resolve().parent
EXECUTABLE = SCRIPT_DIR / ".." / "install" / "latest" / "G_unit_tests"

# --------------------------------------------------------------------------
# Hard-coded default test settings (from your original G_unit_tests.in)
# --------------------------------------------------------------------------

DEFAULT_FLAGS = {
    "l_w_up_in_cloud_test": True,
    "l_KK_unit_tests": True,
    "l_corr_cholesky_mtx_tests": True,
    "l_hole_filling_tests": True,
    "l_Nc_Ncn_test": True,
    "l_read_corr_mtx_test": True,
    "l_silhs_category_test": True,
    "l_mu_sigma_hm_tests": True,
    "l_pdf_parameter_tests": True,
    "l_spurious_source_test": True,
    "l_tuner_tests": True,
    "l_smooth_heaviside_test": True,
    "l_smooth_min_max_test": True,
    "l_rev_direction_grid_test": False,
    "l_fill_holes_test": True,
}

# Mapping from CLI flag → namelist variable
TEST_OPTIONS = {
    "KK_unit_tests": "l_KK_unit_tests",
    "corr_cholesky_mtx_tests": "l_corr_cholesky_mtx_tests",
    "hole_filling_tests": "l_hole_filling_tests",
    "Nc_Ncn_test": "l_Nc_Ncn_test",
    "read_corr_mtx_test": "l_read_corr_mtx_test",
    "silhs_category_test": "l_silhs_category_test",
    "mu_sigma_hm_tests": "l_mu_sigma_hm_tests",
    "pdf_parameter_tests": "l_pdf_parameter_tests",
    "spurious_source_test": "l_spurious_source_test",
    "tuner_tests": "l_tuner_tests",
    "smooth_heaviside_test": "l_smooth_heaviside_test",
    "smooth_min_max_test": "l_smooth_min_max_test",
    "rev_direction_grid_test": "l_rev_direction_grid_test",
    "fill_holes_test": "l_fill_holes_test",
    "w_up_in_cloud_test": "l_w_up_in_cloud_test",
}

# --------------------------------------------------------------------------
# Write temporary namelist
# --------------------------------------------------------------------------

def write_namelist(flags_dict):
    lines = ["&G_unit_namelist"]
    for key, value in flags_dict.items():
        lines.append(f"  {key} = {'.true.' if value else '.false.'}")
    lines.append("/")
    return "\n".join(lines) + "\n"

# --------------------------------------------------------------------------
# CLI setup
# --------------------------------------------------------------------------

def build_arg_parser():
    parser = argparse.ArgumentParser(
        description="Run CLUBB G_unit_tests with selectable subsets of tests."
    )

    parser.add_argument("--all", action="store_true", help="Enable all tests")

    # Add individual test flags
    for flag in TEST_OPTIONS:
        parser.add_argument(
            f"--{flag}",
            action="store_true",
            help=f"Enable only {flag.replace('_',' ')}"
        )

    return parser

# --------------------------------------------------------------------------
# Main logic
# --------------------------------------------------------------------------

def main():
    parser = build_arg_parser()
    args = parser.parse_args()

    # Determine whether any specific test flags were used
    specific_flags_used = any(
        getattr(args, flag) for flag in TEST_OPTIONS
    )

    if args.all:
        # --all → force everything on
        final_flags = {k: True for k in DEFAULT_FLAGS}

    elif specific_flags_used:
        # Specific tests were requested → ONLY those are TRUE
        final_flags = {k: False for k in DEFAULT_FLAGS}
        for flag, namelist_key in TEST_OPTIONS.items():
            if getattr(args, flag):
                final_flags[namelist_key] = True

    else:
        # No flags specified → use built-in defaults
        final_flags = DEFAULT_FLAGS.copy()

    # Write local namelist
    tmp_namelist = SCRIPT_DIR / "G_unit_tests.in"
    tmp_namelist.write_text(write_namelist(final_flags))

    # Run the executable
    print(f"Running {EXECUTABLE} ...")
    result = subprocess.run([str(EXECUTABLE)], cwd=str(SCRIPT_DIR))

    # Cleanup
    try:
        tmp_namelist.unlink()
    except FileNotFoundError:
        pass

    # Handle return code
    if result.returncode != 0:
        print("Error: A G unit test has failed.")
        sys.exit(1)
    else:
        print("All tests have succeeded.")
        sys.exit(0)


# --------------------------------------------------------------------------

if __name__ == "__main__":
    main()
