#!/usr/bin/python3

# =============================================================================
# Test column independence by running generated parameters in ABC and CBA order
# and comparing matching columns. See tests/README for usage examples.
#
# Setup: reconfigure and incrementally build gfortran CPU/double Release in
# build/column_mirror, enabling SILHS_MULTI_COL_RAND_DUPLICATE through a toolchain
# wrapper. Requires CMake, gfortran and normal CLUBB dependencies; leaves normal
# toolchains/install links unchanged. -exe or another explicit runtime skips
# building and must provide equivalent SILHS sampling. Both runs force
# l_lh_straight_mc=.true. Extra model options are forwarded; -- is optional.
#
# Runs: generate parameters once, then reverse their arrays exactly to avoid
# resampling roundoff. Retain parameters, build/params/forward/reverse logs and
# forward/reverse output in a fresh run_* under output/column_mirror_test (or
# -out_dir). Stop on build/run failure and report the log. Two positional output
# directories select comparison only, without building or running the model.
#
# Comparison: require matching, nonempty *_stats.nc sets with >=2 columns.
# Compare numeric time/col fields after reversal; require matching names,
# dimensions, shapes and masks. Reject unmasked NaN/Inf and files without valid
# data; skip fully masked fields. Each column pair's mean absolute difference
# must be <= -t/--tol (default 0); diagnostics also report maximum differences.
# Print stages, cases, parameter order and per-case PASS/FAIL; -v lists fields.
# Failures return nonzero. Long parameter lists are abbreviated on screen.
#
# Debugging: failures save differences/comparison.log and aligned copies of
# failed CBA files in differences/reverse_aligned (temporary storage when
# comparing saved output). Copies reverse col-dependent data and parameters,
# preserving raw values/masks and leaving the col coordinate and originals
# unchanged. The printed run_bindiff_all.py -v 2 command uses these copies;
# comparison.log also covers invalid/missing data that bindiff cannot diagnose.
# Diagnostic indices refer to concatenated ABC CBA: col0 vs col5 compares A to A.
# Check incorrect column indexing or state carried between columns; bugs in
# code paths the selected configuration never executes will not be detected.
# =============================================================================

import argparse
import io
import os
import platform
import shutil
from contextlib import ExitStack, redirect_stdout
import sys
import shlex
import subprocess
import tempfile
import textwrap
from pathlib import Path

import netCDF4
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
BUILD_DIR = REPO_ROOT / "build" / "column_mirror"
sys.path.insert(0, str(REPO_ROOT))

from utilities.create_case_namelist import resolve_tunable_config_dir
from utilities.create_multi_col_params import parse_hypergrid_range_spec
from utilities.output_paths import resolve_output_dir
from run_scripts.run_scm_all import STANDARD_CASES


def _is_numeric_netcdf_var(var):
    return np.issubdtype(var.dtype, np.number)


def _compare_mirrored_columns_vectorized(data_with_col_last):
    """Compare mirrored columns in one vectorized pass.

    data_with_col_last is expected to have the column axis at -1.
    Returns average and max absolute differences for each mirrored pair.
    """
    ncol = data_with_col_last.shape[-1]
    npairs = ncol // 2
    if npairs == 0:
        return np.array([]), np.array([])

    # Compare first half against reversed second half:
    # col 0 vs col n-1, col 1 vs col n-2, ...
    lhs = data_with_col_last[..., :npairs]
    rhs = np.flip(data_with_col_last, axis=-1)[..., :npairs]
    abs_diff = np.ma.abs(lhs - rhs)

    # Reduce over every axis except mirrored-pair axis (last axis).
    if abs_diff.ndim == 1:
        avg_by_pair = abs_diff
        max_by_pair = abs_diff
    else:
        # After moving col to the last axis, all leading axes are sample axes
        # (e.g., time/zt/zm/...), so reduce them and keep one value per col-pair.
        reduce_axes = tuple(range(abs_diff.ndim - 1))
        avg_by_pair = np.ma.mean(abs_diff, axis=reduce_axes)
        max_by_pair = np.ma.max(abs_diff, axis=reduce_axes)

    # Fill masked entries to keep output stable.
    return np.ma.filled(avg_by_pair, 0.0), np.ma.filled(max_by_pair, 0.0)


def check_file(file_path, reverse_file_path, tolerance, verbose):
    with ExitStack() as stack:
        try:
            dset = stack.enter_context(netCDF4.Dataset(file_path))
            reverse_dset = stack.enter_context(netCDF4.Dataset(reverse_file_path))
        except OSError as exc:
            print(f"Error opening comparison files: {exc}")
            return True
        for dataset in (dset, reverse_dset):
            if "col" not in dataset.dimensions or len(dataset.dimensions["col"]) < 2:
                print(f"{dataset.filepath()} must contain at least two columns.")
                return True
        ngrdcol = len(dset.dimensions["col"])
        if len(reverse_dset.dimensions["col"]) != ngrdcol:
            print(f"Column-count mismatch between {file_path} and {reverse_file_path}.")
            return True

        def fields(dataset):
            return {name for name, var in dataset.variables.items()
                    if "col" in var.dimensions and "time" in var.dimensions
                    and _is_numeric_netcdf_var(var)}

        names = fields(dset)
        if not names or names != fields(reverse_dset):
            print(f"Missing or mismatched time-dependent column fields in {file_path} and {reverse_file_path}.")
            return True
        differing_var_count = 0
        valid_var_count = 0
        for var_name in sorted(names):
            var = dset.variables[var_name]
            reverse_var = reverse_dset.variables[var_name]
            if var.dimensions != reverse_var.dimensions or var.shape != reverse_var.shape:
                print(f" - {var_name}: dimensions or shape differ between runs")
                differing_var_count += 1
                continue
            data = np.moveaxis(np.ma.array(var[:]), var.dimensions.index("col"), -1)
            reverse_data = np.moveaxis(np.ma.array(reverse_var[:]), reverse_var.dimensions.index("col"), -1)
            mirrored = np.flip(reverse_data, axis=-1)
            if not np.array_equal(np.ma.getmaskarray(data), np.ma.getmaskarray(mirrored)):
                print(f" - {var_name}: missing-value masks differ between matching columns")
                differing_var_count += 1
                continue
            if not np.isfinite(data.compressed()).all() or not np.isfinite(mirrored.compressed()).all():
                print(f" - {var_name}: non-finite values in output")
                differing_var_count += 1
                continue
            if not data.count():
                continue
            valid_var_count += 1
            if verbose:
                print(f"checking {var_name} dims={var.dimensions}")
            avg_by_pair, max_by_pair = _compare_mirrored_columns_vectorized(
                np.ma.concatenate((data, reverse_data), axis=-1)
            )
            variable_has_difference = False
            for col, (avg_abs_diff, max_abs_diff) in enumerate(zip(avg_by_pair, max_by_pair)):
                if avg_abs_diff > tolerance:
                    print(
                        f" - {var_name} col{col} vs col{2 * ngrdcol - 1 - col}: "
                        f"avg_abs_diff = {avg_abs_diff}, max_abs_diff = {max_abs_diff}"
                    )
                    variable_has_difference = True
            differing_var_count += int(variable_has_difference)
        if not valid_var_count:
            print("No valid column data were compared.")
        return differing_var_count > 0 or valid_var_count == 0


def report_differences(forward_directory, reverse_directory, failed_names, details, report_dir):
    """Retain diagnostics and align failed CBA files for the ordinary bindiff tool."""
    report_dir = Path(report_dir) if report_dir else Path(tempfile.mkdtemp(prefix="column_mirror_diff_"))
    report_dir.mkdir(parents=True, exist_ok=True)
    log_path = report_dir / "comparison.log"
    log_path.write_text(details)
    print(f"  Details: {log_path}", flush=True)
    aligned = report_dir / "reverse_aligned"
    aligned.mkdir()
    for name in failed_names:
        if not (forward_directory / name).is_file() or not (reverse_directory / name).is_file():
            continue
        target = aligned / name
        try:
            # Copy before editing, and reverse raw data to avoid packing/rounding changes.
            shutil.copyfile(reverse_directory / name, target)
            with netCDF4.Dataset(target, "r+") as dataset:
                dataset.set_auto_maskandscale(False)
                for var_name, var in dataset.variables.items():
                    if "col" in var.dimensions and var_name != "col":
                        var[:] = np.flip(var[:], axis=var.dimensions.index("col"))
        except (OSError, RuntimeError, ValueError) as exc:
            target.unlink(missing_ok=True)
            with log_path.open("a") as log:
                log.write(f"Could not align {target.name}: {exc}\n")
    if any(aligned.glob("*_stats.nc")):
        print("  Inspect differences (CBA copies aligned to ABC):", flush=True)
        print("  " + shlex.join([
            sys.executable, str(REPO_ROOT / "run_scripts/run_bindiff_all.py"),
            str(forward_directory.resolve()), str(aligned.resolve()), "-v", "2",
        ]), flush=True)


def compare_directories(forward_directory, reverse_directory, tolerance=0.0, verbose=False,
                        report_dir=None):
    forward_directory, reverse_directory = Path(forward_directory), Path(reverse_directory)
    for directory in (forward_directory, reverse_directory):
        if not Path(directory).is_dir():
            print(f"Error: {directory} is not a valid directory.")
            return 1

    forward_files = {p.name for p in forward_directory.glob("*_stats.nc")}
    reverse_files = {p.name for p in reverse_directory.glob("*_stats.nc")}
    if not forward_files and not reverse_files:
        print("FAIL: no *_stats.nc files to compare.")
        return 1

    failed_names, details = [], []
    for name in sorted(forward_files | reverse_files):
        with redirect_stdout(io.StringIO()) as output:
            if name not in forward_files or name not in reverse_files:
                missing = "ABC" if name not in forward_files else "CBA"
                print(f"Missing {missing} output: {name}")
                failed = True
            else:
                failed = check_file(forward_directory / name, reverse_directory / name,
                                    tolerance, verbose)
        print(f"  {'FAIL' if failed else 'PASS'} {name.removesuffix('_stats.nc')}", flush=True)
        if verbose:
            print(output.getvalue(), end="", flush=True)
        if failed:
            failed_names.append(name)
            details.append(f"{name}\n{output.getvalue()}\n")
    if failed_names:
        report_differences(forward_directory, reverse_directory, failed_names,
                           "".join(details), report_dir)
    return int(bool(failed_names))


def generate_parameter_files(run_dir, spec, base_params):
    """Generate once, then reverse complete columns without resampling the ranges."""
    forward = run_dir / "params_forward.in"
    reverse = run_dir / "params_reverse.in"
    log_path = run_dir / "params.log"
    with log_path.open("w") as log:
        result = subprocess.run([
            sys.executable, str(REPO_ROOT / "utilities/create_multi_col_params.py"),
            "-hr", spec, "-param_file", str(base_params), "-out_file", str(forward),
        ], stdout=log, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f"Parameter generation failed; see {log_path}")
    lines = []
    in_params = False
    for line in forward.read_text().splitlines():
        if line.strip().lower() == "&clubb_params_nl":
            in_params = True
        elif line.strip() == "/":
            in_params = False
        elif in_params and "=" in line:
            name, values = line.split("=", 1)
            # The generator writes one comma-separated parameter array per line.
            line = name + "= " + ", ".join(reversed([v.strip() for v in values.split(",")]))
        lines.append(line)
    reverse.write_text("\n".join(lines) + "\n")
    return forward, reverse


def print_run_parameters(param_file, spec):
    """Display sampled parameters from the actual input, including hypergrid order."""
    targets = {target for item in parse_hypergrid_range_spec(spec) for target in item["targets"]}
    for line in param_file.read_text().splitlines():
        if "=" not in line:
            continue
        name, values = line.split("=", 1)
        if name.strip() not in targets:
            continue
        values = [f"{float(value):.6g}" for value in values.split(",")]
        count = len(values)
        if count > 12:
            values = values[:4] + ["..."] + values[-4:]
        suffix = f" ({count} columns; full list in {param_file.name})" if count > 12 else ""
        print(f"  {name.strip()} = [{', '.join(values)}]{suffix}", flush=True)


def parse_args(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    forwarded = []
    if "--" in argv:
        split = argv.index("--")
        argv, forwarded = argv[:split], argv[split + 1:]
    parser = argparse.ArgumentParser(
        description="Build and run a column-order independence test, or compare two saved output directories.",
        epilog="The default gfortran build is refreshed incrementally. Extra run_scm.py options are forwarded (-- is optional); -exe PATH skips compilation.",
        allow_abbrev=False,
    )
    parser.add_argument("directories", nargs="*", metavar="DIRECTORY",
                        help="Two saved output directories to compare without running the model")
    cases = parser.add_mutually_exclusive_group()
    cases.add_argument("-case", help="Run one case (default: run_scm_all.py standard case set)")
    cases.add_argument("-cases", help="Comma-separated case list, passed to run_scm_all.py")
    columns = parser.add_mutually_exclusive_group()
    columns.add_argument("-n", type=int, help="Number of columns with C8 evenly spaced from 0.2 to 0.8 (default: 3)")
    columns.add_argument("-hr", help="Parameter range specification accepted by create_multi_col_params.py")
    parser.add_argument("-params", help="Base single-column parameter file (default: from -config)")
    parser.add_argument("-config", help="Tunable configuration name or directory (default: default)")
    parser.add_argument("-out_dir", default="column_mirror_test",
                        help="Output parent; a fresh run directory is retained inside it")
    parser.add_argument("-nproc", type=int, default=2, help="Concurrent cases for a case list (default: 2)")
    parser.add_argument("-max_iters", type=int, default=200, help="Timesteps per run (default: 200)")
    parser.add_argument("-override", help="Namelist overrides for both runs; straight SILHS sampling is always enabled")
    parser.add_argument("-t", "--tol", type=float, default=0.0,
                        help="Average absolute difference tolerance (default: 0)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Print every checked variable")
    # Partition before parsing: parse_known_args would mistake -tout for -t,
    # -nzmax for -n, and forwarded option values for positional directories.
    test_args, extra = [], []
    index = 0
    while index < len(argv):
        token = argv[index]
        option = token.split("=", 1)[0]
        action = parser._option_string_actions.get(option)
        # Preserve argparse's compact numeric forms, e.g. -n3 and -t1e-8.
        compact = False
        if action is None and token[:2] in {"-n", "-t"} and len(token) > 2:
            candidate = parser._option_string_actions[token[:2]]
            try:
                candidate.type(token[2:])
            except ValueError:
                pass
            else:
                action, compact = candidate, True
        if action is not None:
            test_args.append(token)
            if action.nargs != 0 and "=" not in token and not compact and index + 1 < len(argv):
                index += 1
                test_args.append(argv[index])
        elif token.startswith("-"):
            extra.append(token)
            while index + 1 < len(argv) and not argv[index + 1].startswith("-"):
                index += 1
                extra.append(argv[index])
        else:
            test_args.append(token)
        index += 1
    forwarded = extra + forwarded
    args = parser.parse_args(test_args)
    if not np.isfinite(args.tol) or args.tol < 0:
        parser.error("tolerance must be finite and nonnegative")
    if args.directories:
        if len(args.directories) != 2:
            parser.error("provide two saved directories, or use -case to run a case")
        comparison_options = {"-t", "--tol", "-v", "--verbose"}
        if forwarded or any(token.startswith("-") and token.split("=", 1)[0] not in comparison_options
                            for token in argv):
            parser.error("saved-directory comparison accepts only -t and -v")
        return args, forwarded
    if args.n is not None and args.n < 2:
        parser.error("-n must be at least 2")
    if args.nproc < 1 or args.max_iters < 1:
        parser.error("-nproc and -max_iters must be positive")
    args.spec = args.hr or f"C8/0.2:0.8/{args.n or 3}"
    try:
        specs = parse_hypergrid_range_spec(args.spec)
    except ValueError as exc:
        parser.error(str(exc))
    if not all(np.isfinite(s[key]) for s in specs for key in ("min", "max")):
        parser.error("parameter range endpoints must be finite")
    if not any(s["npoints"] > 1 and s["min"] != s["max"] for s in specs):
        parser.error("the test needs at least two distinct parameter columns")
    reserved = {"-params", "-config", "-out_dir", "-multicol", "-batch_size", "-override",
                "-max_iters", "-nproc", "-cases", "-all", "-min_cases", "-short_cases", "-priority_cases"}
    if any(token.split("=", 1)[0] in reserved for token in forwarded):
        parser.error("parameter generation, output paths, overrides, and case selection must use test-runner options before --")
    return args, forwarded


def build_test_executable(log_path):
    """Refresh an isolated test build; let CMake reuse up-to-date objects."""
    cmake = shutil.which("cmake")
    compiler = shutil.which("gfortran")
    if not cmake or not compiler:
        raise RuntimeError("Automatic test builds require CMake and gfortran on PATH; "
                           "load/install them, or supply -exe PATH.")
    platform_name = f"{platform.system().lower()}_{platform.machine().lower()}"
    base_toolchain = REPO_ROOT / "cmake" / "toolchains" / f"{platform_name}_gcc.cmake"
    if not base_toolchain.is_file():
        raise RuntimeError(f"No gfortran toolchain for {platform_name}; supply -exe PATH.")
    BUILD_DIR.mkdir(parents=True, exist_ok=True)
    wrapper = BUILD_DIR / "toolchain.cmake"
    content = (f'include([=[{base_toolchain.as_posix()}]=])\n'
               'add_compile_definitions(SILHS_MULTI_COL_RAND_DUPLICATE)\n')
    if not wrapper.exists() or wrapper.read_text() != content:
        wrapper.write_text(content)
    commands = [
        [cmake, "-S", str(REPO_ROOT), "-B", str(BUILD_DIR),
         f"-DCMAKE_TOOLCHAIN_FILE={wrapper}", f"-DCMAKE_Fortran_COMPILER={compiler}",
         "-DCMAKE_BUILD_TYPE=Release", "-DGPU=none", "-DPRECISION=double",
         "-DENABLE_TESTS=OFF", "-DENABLE_F2PY=OFF", "-DENABLE_OMP=OFF",
         "-DUSE_GPTL=OFF", "-DTUNING=OFF"],
        [cmake, "--build", str(BUILD_DIR), "--target", "clubb_standalone",
         "--parallel", str(min(8, os.cpu_count() or 1))],
    ]
    with log_path.open("w") as log:
        for phase, command in zip(("configuration", "build"), commands):
            log.write(shlex.join(command) + "\n")
            log.flush()
            result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
            if result.returncode:
                log.flush()
                tail = "\n".join(log_path.read_text(errors="replace").splitlines()[-20:])
                raise RuntimeError(f"CMake {phase} failed (exit {result.returncode}); "
                                   f"see {log_path}\n{tail}")
    executable = BUILD_DIR / "src" / "clubb_standalone"
    if not executable.is_file():
        raise RuntimeError(f"Build did not produce {executable}; see {log_path}")
    return executable


def prepare_runtime(forwarded, log_path):
    selectors = {"-exe", "-install_dir", "-driver_test", "-python", "-jax"}
    if any(token.split("=", 1)[0] in selectors for token in forwarded):
        print("[1/4] Configuring and compiling: skipped (supplied runtime)", flush=True)
        return list(forwarded)
    print("[1/4] Configuring and compiling (gfortran, CPU)", flush=True)
    return [*forwarded, "-exe", str(build_test_executable(log_path))]


def run_test(args, forwarded):
    config_dir = Path(resolve_tunable_config_dir(args.config))
    base_params = Path(args.params).resolve() if args.params else config_dir / "tunable_parameters.in"
    output_parent = resolve_output_dir(args.out_dir)
    output_parent.mkdir(parents=True, exist_ok=True)
    run_dir = Path(tempfile.mkdtemp(prefix="run_", dir=output_parent))
    print(f"Logs and output: {run_dir}", flush=True)
    params = generate_parameter_files(run_dir, args.spec, base_params)
    forwarded = prepare_runtime(forwarded, run_dir / "build.log")
    override = ",".join(filter(None, [args.override, "l_lh_straight_mc=.true."]))
    script = "run_scm.py" if args.case else "run_scm_all.py"
    cases = [args.case] if args.case else (args.cases.split(",") if args.cases else STANDARD_CASES)
    for stage, order, label, param_file in zip((2, 3), ("forward", "reverse"), ("ABC", "CBA"), params):
        print(f"[{stage}/4] {label} run", flush=True)
        print(textwrap.fill("Cases: " + ", ".join(case.strip() for case in cases),
                            width=100, initial_indent="  ", subsequent_indent="    "), flush=True)
        print(f"  Grid: -hr {args.spec}" + (" (columns reversed)" if order == "reverse" else ""), flush=True)
        print_run_parameters(param_file, args.spec)
        command = [sys.executable, str(REPO_ROOT / "run_scripts" / script), *forwarded,
                   "-config", str(config_dir), "-params", str(param_file),
                   "-out_dir", str(run_dir / order), "-max_iters", str(args.max_iters),
                   "-override", override]
        if args.case:
            command.append(args.case)
        else:
            command += ["-nproc", str(args.nproc)]
            if args.cases:
                command += ["-cases", args.cases]
        log_path = run_dir / f"{order}.log"
        with log_path.open("w") as log:
            log.write(shlex.join(command) + "\n")
            log.flush()
            result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
        if result.returncode:
            print(f"  FAIL {label} run (exit {result.returncode}); see {log_path}", flush=True)
            return 1
    print("[4/4] Comparison", flush=True)
    return compare_directories(run_dir / "forward", run_dir / "reverse", args.tol, args.verbose,
                               run_dir / "differences")


def main(argv=None):
    args, forwarded = parse_args(argv)
    if args.directories:
        print("Comparison (saved output)", flush=True)
        return compare_directories(*args.directories, args.tol, args.verbose)
    try:
        return run_test(args, forwarded)
    except (OSError, RuntimeError, ValueError, subprocess.CalledProcessError) as exc:
        print(f"Column mirror test failed: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
