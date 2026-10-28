#!/usr/bin/env python3
"""Compare reused-process case output with fresh-process controls.

Exercise the full initialization/advance/cleanup path through clubb_case_runner,
using BOMEX (two columns, five steps), RICO_SILHS (one column, three steps), then
BOMEX again. The Python SCM uses BOMEX for the middle request because it does
not implement SILHS. Default physics, statistics and timesteps are retained;
only column counts, iteration caps and debug 0 are supplied explicitly.

Usage:
    python tests/case_runner_lifecycle_test.py -output_dir output/runner_check
    python tests/case_runner_lifecycle_test.py -jax=cpu -output_dir output/runner_jax
    - main(): prepare inputs, run reused/fresh processes and print comparisons.
    - compare_stats(actual, expected): return the number of numeric variables
      checked, raising on a variable-name, numeric shape, mask or value mismatch.

Understanding:
    Require one PID throughout the reused sequence, plus exact equality of all
    numeric stats variables against independent same-backend processes. Compare
    unmasked values and masks; NetCDF attributes and file bytes are not compared.
    Output directories must be new, and logs/inputs/stats are retained there.
    This model workflow runs in the first stage of the clubb_driver Jenkins
    pipeline for native Fortran, Python and CPU JAX; it is separate from pytests.

    Extends src/clubb_driver_test.F90's same-case reinitialization check to
    different cases and column counts. Returning to BOMEX checks that
    microphysics settings, notably Nc0_in_cloud after RICO_SILHS, are reset
    rather than inherited from the preceding case.
"""

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

if __name__ == "__main__":
    from utilities.setup_python_venv import ensure_python_venv

    ensure_python_venv()

import netCDF4
import numpy as np

from utilities.clubb_case_runner import (
    add_runner_arguments,
    clubb_case_runner,
    extract_jax_options,
    parse_case_arguments,
    validate_runtime_options,
)
from input.scripts.create_case_namelist import create_case_namelist_from_args


def compare_stats(actual_path, expected_path):
    """Return the numeric variable count after comparing two NetCDF paths exactly.

    Require matching variable names and numeric shapes/masks/unmasked values.
    Attributes and byte-level file representation are outside this comparison.
    """
    with netCDF4.Dataset(actual_path) as actual, netCDF4.Dataset(expected_path) as expected:
        assert set(actual.variables) == set(expected.variables), "Stats variable schemas differ"
        count = 0
        for name, variable in expected.variables.items():
            if variable.dtype.kind not in "fiu":
                continue
            observed = actual.variables[name][:]
            reference = variable[:]
            assert observed.shape == reference.shape, f"{name}: stats shapes differ"
            np.testing.assert_array_equal(
                np.ma.getmaskarray(observed), np.ma.getmaskarray(reference), err_msg=name,
            )
            np.testing.assert_array_equal(observed.compressed(), reference.compressed(), err_msg=name)
            count += 1
        assert count, "No numeric stats variables compared"
        return count


def main():
    """Run the selected backend's case sequence and compare every stats output with fresh controls."""
    parser = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    add_runner_arguments(parser)
    parser.add_argument("-output_dir", type=Path, default=ROOT / "output/tests/case_runner_lifecycle",
                        help="New directory for input, logs and stats output.")
    normalized, jax_options, occurrences = extract_jax_options(sys.argv[1:])
    options = parser.parse_args(normalized)
    if occurrences > 1:
        parser.error("-jax may be specified only once.")
    options.jax_options = jax_options
    try:
        validate_runtime_options(options)
    except ValueError as exc:
        parser.error(str(exc))
    if options.gdb or options.driver_test:
        parser.error("This workflow requires the standalone, without -gdb or -driver_test.")
    output = options.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    middle_case = "bomex" if options.python else "rico_silhs"
    scenarios = (("bomex", 2, 5), (middle_case, 1, 3), ("bomex", 2, 5))

    # Keep one worker across the sequence, with a fresh process for each control.
    comparisons = []
    with clubb_case_runner(
        options, persistent=True,
        log_file=output / "startup.log", timeout=180,
    ) as runner:
        runner.init()
        pid = runner.status()["pid"]
        for index, (case, columns, steps) in enumerate(scenarios):
            paths = {}
            for label in ("persistent", "control"):
                directory = output / f"{index}_{case}" / label
                directory.mkdir(parents=True)
                args = parse_case_arguments([
                    case, "-multicol", str(columns), "-max_iters", str(steps), "-debug", "0",
                ])
                with (directory / "prepare.log").open("w") as log:
                    namelist = create_case_namelist_from_args(args, directory, log_file=log)
                if label == "persistent":
                    result = runner.run_and_wait(namelist, log_file=directory / "run.log")
                    assert result["state"] == "completed", result
                    assert result["pid"] == pid, "Reusable process was replaced"
                else:
                    with clubb_case_runner(options, persistent=False, timeout=180) as control:
                        result = control.run_and_wait(namelist, log_file=directory / "run.log")
                assert result["state"] == "completed", result
                paths[label] = directory / f"{case}_stats.nc"
            comparisons.append((index, case, columns, steps, paths))
    # Compare only closed files, after every reused/control run has completed.
    for index, case, columns, steps, paths in comparisons:
        count = compare_stats(paths["persistent"], paths["control"])
        print(f"{index + 1}: {case}, {columns} columns, {steps} steps: "
              f"{count} numeric variables exact; PID {pid}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
