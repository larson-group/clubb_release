#!/usr/bin/env python3
"""Require the real JAX/Fortran comparison to detect injected model errors.

Run a passing control, then the same comparison with a small JAX-only
error in heating, parameter loading, column indexing, or stats output. Each
mutation runs independently against the control, using private source and output trees.
The comparison script, Fortran executable, tolerances, and inputs are unchanged.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from clubb_jax.run_jax import ensure_environment  # noqa: E402

HARNESS = Path("tests/run_jax_vs_fortran_cases.py")
RESULTS = Path("output/tests/jax_driver_test_results")
CASE = "bomex"
TIMESTEPS_TO_RUN = 4


@dataclass(frozen=True)
class Mutation:
    name: str
    description: str
    target: Path
    original: str
    mutated: str
    field: str
    later_columns_only: bool = False
    expected_c8: tuple[float, ...] = (0.2, 0.4, 0.6, 0.8)
    kind: str = "numerical"


# 1. Wrong forcing: add 1e-6 K/s of heating instead of zero in every column.
HEATING = Mutation(
    name="heating",
    description="Add 1e-6 K/s of spurious heating to every column.",
    target=Path(f"clubb_jax/src/Benchmark_cases/{CASE}.py"),
    original="thlm_forcing = jnp.zeros((ngrdcol, gr.nzt), dtype=rtm.dtype)",
    mutated="thlm_forcing = jnp.full((ngrdcol, gr.nzt), 1.0e-6, dtype=rtm.dtype)",
    field="thlm",
)

# 2. Wrong parameter loading: halve the C8 spread while keeping column 1.
# The harness supplies [0.2, 0.4, 0.6, 0.8]; JAX instead uses [0.2, 0.3, 0.4, 0.5].
# Require a wp3 failure, not just a difference in the saved parameter metadata.
PARAMETER_HANDLING = Mutation(
    name="parameter_handling",
    description="Load C8 as [0.2, 0.3, 0.4, 0.5] instead of [0.2, 0.4, 0.6, 0.8].",
    target=Path("clubb_jax/src/CLUBB_core/parameters_tunable.py"),
    original="values[column, _NAME_TO_IDX[key]] = float(value)",
    mutated=(
        'values[column, _NAME_TO_IDX[key]] = '
        'values[0, _NAME_TO_IDX[key]] + 0.5 * (float(value) - values[0, _NAME_TO_IDX[key]]) '
        'if key == "c8" and column > 0 else float(value)'
    ),
    field="wp3",
    later_columns_only=True,
    expected_c8=(0.2, 0.3, 0.4, 0.5),
)

# 3. Wrong column in a physical calculation: broadcast old wp3(1,k) into the
# time-tendency equation for every column instead of using each column's wp3(i,k).
# Keep the singleton axis (0:1) so this is a numerical bug, not a shape error.
FIRST_COLUMN = Mutation(
    name="first_column",
    description="Use column 1's previous wp3 in every column's timestep equation.",
    target=Path("clubb_jax/src/CLUBB_core/advance_wp2_wp3_module.py"),
    original="rhs = rhs.at[:, 3:-2:2].add(invrs_dt * wp3[:, 1:-1])",
    mutated="rhs = rhs.at[:, 3:-2:2].add(invrs_dt * wp3[0:1, 1:-1])",
    field="wp3",
    later_columns_only=True,
)

# 4. JAX computes wprtp but fails to create its NetCDF output variable.
# The writer already skips uncreated variables, so both models still finish.
MISSING_STAT = Mutation(
    name="missing_stat",
    description="Omit wprtp from JAX statistics output while keeping its calculation.",
    target=Path("clubb_jax/src/CLUBB_core/stats_netcdf.py"),
    original="if dims is None:",
    mutated='if dims is None or name == "wprtp":',
    field="wprtp",
    kind="missing_stat",
)

MUTATIONS = (HEATING, PARAMETER_HANDLING, FIRST_COLUMN, MISSING_STAT)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def interrupt(signum, _frame):
    # Unwind run_comparison so a parent-only TERM also stops the real harness.
    raise KeyboardInterrupt(f"Interrupted by signal {signum}")


def inject_error(source: str, mutation: Mutation) -> str:
    lines = source.splitlines(keepends=True)
    sites = [i for i, line in enumerate(lines) if line.strip() == mutation.original]
    require(len(sites) == 1,
            f"Mutation anchor changed in {mutation.target}; review the test's injection site.")
    # Match a complete line: an altered expression may still start with ORIGINAL.
    i = sites[0]
    lines[i] = lines[i].replace(mutation.original, mutation.mutated, 1)
    return "".join(lines)


def make_snapshot(destination: Path) -> None:
    """Copy current source (including local edits), never hard-link mutable files."""
    destination.mkdir()
    ignore = shutil.ignore_patterns("__pycache__", "*.pyc", ".pytest_cache")
    # The comparison and batch runner share worker defaults from tuner.
    for name in ("clubb_jax", "run_scripts", "utilities", "tuner"):
        shutil.copytree(ROOT / name, destination / name, ignore=ignore)
    # Master now bootstraps shared Python tools as well as JAX. Resolve that
    # helper in the real checkout so it finds its requirements and reuses one
    # environment; model sources and comparison code remain private copies.
    setup = Path("utilities/setup_python_venv.py")
    (destination / setup).unlink()
    (destination / setup).symlink_to(ROOT / setup)
    (destination / "tests").mkdir()
    shutil.copy2(ROOT / HARNESS, destination / HARNESS)
    # Shared inputs and installed binaries are only read by these runs.
    for name in ("input", "install"):
        (destination / name).symlink_to(ROOT / name, target_is_directory=True)


def run_comparison(snapshot: Path, timeout: float) -> int:
    command = [sys.executable, str(snapshot / HARNESS),
               "-cases", CASE, "-workers", "1", "-max_iters", str(TIMESTEPS_TO_RUN),
               "-dt_main", "60", "-stats", "input/stats/multi_col_stats.in"]
    (snapshot / "command.json").write_text(json.dumps(command, indent=2) + "\n")
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    # Use the copied package, even if the invoking shell has a repo PYTHONPATH.
    env.pop("PYTHONPATH", None)
    with (snapshot / "comparison.log").open("w") as log:
        process = subprocess.Popen(command, cwd=snapshot, env=env,
                                   stdout=log, stderr=subprocess.STDOUT)
        try:
            return process.wait(timeout=timeout)
        finally:
            if process.poll() is None:
                # The real harness handles TERM and reaps its model groups.
                process.send_signal(signal.SIGTERM)
                try:
                    process.wait(timeout=15)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()


def validate_summary(summary: dict, returncode: int, *, mutation: Mutation | None) -> None:
    mutated = mutation is not None
    expected_rc = 1 if mutated else 0
    expected_status = "diff" if mutated else "match"
    require(returncode == expected_rc,
            f"Expected comparison exit {expected_rc}, got {returncode}.")
    require(summary.get("total_cases") == 1 and len(summary.get("cases", [])) == 1,
            "Expected exactly one compared case.")
    require(summary.get("match") == int(not mutated)
            and summary.get("diff") == int(mutated), "Unexpected comparison counts.")
    require(all(summary.get(key) == 0 for key in
                ("jax_failed", "fortran_failed", "both_failed")),
            "A model/launcher failure does not count as detecting the mutation.")
    case = summary["cases"][0]
    require(case.get("case") == CASE and case.get("status") == expected_status,
            f"Expected {CASE} status {expected_status}.")
    require(case.get("jax_rc") == case.get("fortran_rc") == 0,
            "Both models must exit successfully.")
    require(case.get("jax_timesteps") == case.get("fortran_timesteps") == TIMESTEPS_TO_RUN,
            f"Both models must complete exactly {TIMESTEPS_TO_RUN} timesteps.")
    require(case.get("bindiff_rc") == summary.get("final_bindiff_rc") == expected_rc,
            "Per-case and final comparisons must agree with the expected verdict.")
    for key in ("bindiff_threshold", "bindiff_percent_threshold"):
        value = case.get(key)
        require(isinstance(value, (int, float)) and math.isfinite(value) and value >= 0,
                f"Missing or invalid comparison threshold: {key}.")
    if mutated:
        first_failure = case.get("first_failing_timestep")
        if mutation.kind == "missing_stat":
            require(first_failure is None,
                    "A missing-stat failure should have no numerical failure timestep.")
        else:
            require(first_failure in range(1, TIMESTEPS_TO_RUN + 1),
                    "The numerical failure must be located within the completed run.")


def read_fields(snapshot: Path, driver: str) -> dict:
    import netCDF4
    import numpy as np

    path = snapshot / RESULTS / f"{driver}_output/default/{CASE}_stats.nc"
    with netCDF4.Dataset(path) as dataset:
        dataset.set_auto_mask(False)  # A zero fill value is also valid model data.
        fields = {}
        for name in ("thlm", "rtm", "wp3"):
            require(name in dataset.variables, f"Missing {name}: {path}")
            var = dataset[name]
            values = np.asarray(var[:])
            require(var.dimensions == ("time", "zt", "col")
                    and values.shape[0] == TIMESTEPS_TO_RUN and values.shape[-1] == 4,
                    f"Expected {TIMESTEPS_TO_RUN} saved records and four columns for {name}: {path}")
            require(values.size > 0 and np.isfinite(values).all(),
                    f"Empty/nonfinite {name}: {path}")
            fields[name] = values
        return fields


def check_parameters(snapshot: Path, driver: str, expected: tuple) -> None:
    import netCDF4
    import numpy as np

    path = snapshot / RESULTS / f"{driver}_output/default/{CASE}_stats.nc"
    with netCDF4.Dataset(path) as dataset:
        # Decode fixed-width character rows directly; chartostring assigns to
        # ndarray.shape internally, which is deprecated in NumPy 2.5.
        dataset.set_auto_mask(False)
        dataset.set_auto_chartostring(False)
        names = [row.tobytes().decode("utf-8").strip("\x00 ")
                 for row in dataset["param_name"][:]]
        values = np.asarray(dataset["clubb_params"][names.index("C8"), :])
        require(values.shape == (4,) and np.allclose(values, expected, rtol=0, atol=1e-14),
                f"Unexpected {driver} C8 values: {values}; expected {expected}.")


def check_phase(snapshot: Path, returncode: int, *,
                mutation: Mutation | None = None) -> tuple[tuple[dict, dict], dict]:
    summary = json.loads((snapshot / RESULTS / "case_compare_summary.json").read_text())
    validate_summary(summary, returncode, mutation=mutation)
    jax_fields, fortran_fields = (read_fields(snapshot, driver) for driver in ("jax", "fortran"))
    check_parameters(snapshot, "fortran", HEATING.expected_c8)
    check_parameters(snapshot, "jax", mutation.expected_c8 if mutation else HEATING.expected_c8)
    for name in jax_fields:
        require(jax_fields[name].shape == fortran_fields[name].shape,
                f"JAX/Fortran shapes differ for {name}.")
    if mutation is not None and mutation.kind == "numerical":
        log = (snapshot / RESULTS / f"logs/default/{CASE}_bindiff.log").read_text()
        require(re.search(rf"^\s*{mutation.field}\s+\d", log, re.MULTILINE) is not None,
                f"The comparator must report {mutation.field} as exceeding its thresholds.")
    return (jax_fields, fortran_fields), summary["cases"][0]


def check_missing_stat(control: Path, changed: Path, mutation: Mutation) -> dict:
    """Require a one-variable output omission, with every other value unchanged."""
    import netCDF4
    import numpy as np

    report_path = changed / RESULTS / f"logs/default/{CASE}_bindiff.json"
    report = json.loads(report_path.read_text())
    require(report.get("strict") is True, "The structural mutation requires strict bindiff.")
    case = report.get("cases", {}).get(CASE, {})
    files = [file for file in case.get("files", []) if file.get("name") == f"{CASE}_stats.nc"]
    require(len(files) == 1 and files[0].get("status") == "diff",
            "The missing-stat file must be reported as different.")
    file = files[0]
    categories = file.get("variables", {})
    require(categories.get("only_right") == [mutation.field],
            f"Strict bindiff must report only {mutation.field} as missing from JAX output.")
    require(not any(categories.get(kind) for kind in ("only_left", "shape_mismatch", "different"))
            and not file.get("issues"),
            "An unrelated variable or file difference also caused the comparison to fail.")

    for driver in ("jax", "fortran"):
        control_path = control / RESULTS / f"{driver}_output/default/{CASE}_stats.nc"
        changed_path = changed / RESULTS / f"{driver}_output/default/{CASE}_stats.nc"
        with netCDF4.Dataset(control_path) as baseline, netCDF4.Dataset(changed_path) as mutated:
            for dataset in (baseline, mutated):
                dataset.set_auto_mask(False)
                dataset.set_auto_chartostring(False)
            before = set(baseline.variables)
            after = set(mutated.variables)
            expected = before - {mutation.field} if driver == "jax" else before
            require(mutation.field in before and after == expected,
                    f"{driver} output did not have exactly the expected variable set.")
            for name in after:
                require(np.array_equal(baseline[name][:], mutated[name][:]),
                        f"{driver} variable {name} changed alongside the missing stat.")
    return {"missing_stat": mutation.field, "other_stats_unchanged": True}


def error_metrics(a, b, axis=None) -> tuple:
    import numpy as np

    absolute = np.mean(np.abs(a - b), axis=axis)
    # Match bindiff's clipping convention, including for signed wp3 values.
    a_clip, b_clip = (np.clip(x, 1e-7, 9999999.0) for x in (a, b))
    percent = np.mean(200 * np.abs(a_clip - b_clip) / (a_clip + b_clip), axis=axis)
    return absolute, percent


def check_signal(control: tuple[dict, dict], changed: tuple[dict, dict],
                 mutation: Mutation, comparison: dict) -> dict:
    """Independently establish a finite physical signal, not just a failing exit."""
    import numpy as np

    # Use the effective defaults recorded by the real comparison harness.
    abs_threshold = comparison["bindiff_threshold"]
    percent_threshold = comparison["bindiff_percent_threshold"]
    for name, values in control[1].items():
        require(np.array_equal(values, changed[1][name]),
                f"Fortran {name} changed between control and mutation runs.")
    a, b = changed[0][mutation.field], changed[1][mutation.field]
    abs_diff, percent_diff = error_metrics(a, b)
    require(abs_diff > abs_threshold and percent_diff > percent_threshold,
            f"Injected {mutation.field} error did not exceed both comparison thresholds.")
    col_abs, col_percent = error_metrics(a, b, axis=(0, 1))
    col_fails = (col_abs > abs_threshold) & (col_percent > percent_threshold)
    if mutation.later_columns_only:
        require(not col_fails[0] and np.all(col_fails[1:]),
                "Column 1 must match and every later column must fail numerically.")
        for name in control[0]:
            # Guard against an unrelated change to column 1, independent of Fortran.
            first_abs, first_pct = error_metrics(control[0][name][..., 0], changed[0][name][..., 0])
            require(first_abs <= abs_threshold or first_pct <= percent_threshold,
                    f"Mutation unexpectedly changed column 1 of {name}.")
    result = {"field": mutation.field, "mean_abs_diff": float(abs_diff),
              "mean_abs_percent_diff": float(percent_diff),
              "column_mean_abs_diff": col_abs.tolist(),
              "column_mean_abs_percent_diff": col_percent.tolist(),
              "column_fails": col_fails.tolist()}
    if mutation.name == "heating":
        warming = float(np.mean(a - control[0]["thlm"]))
        require(warming > 0, "The JAX-only positive heating mutation had no warming effect.")
        result["jax_mean_warming_K"] = warming
    return result


def print_summary(report: dict, output: Path) -> None:
    labels = {"passed": "PASS", "failed": "FAIL", "not_run": "NOT RUN",
              "cancelled": "CANCELLED"}
    rule = "=" * 78
    print(f"\n{rule}\nExpected failure test: {labels[report['status']]}\n{rule}")
    print(f"Case: {report['case']} | {report['timesteps']} timesteps x 60 s | 4 columns | CPU/double")
    print("Statistics: input/stats/multi_col_stats.in")
    print("Tolerances: JAX-versus-Fortran harness defaults", end="")
    if "thresholds" in report:
        limits = report["thresholds"]
        print(f" (absolute {limits['bindiff_threshold']:g}; "
              f"percentage {limits['bindiff_percent_threshold']:g}%)")
    else:
        print()
    control_status = report["control_status"]
    print(f"\n[{labels[control_status]}] Unmodified control"
          + (" — JAX and Fortran matched." if control_status == "passed" else ""))
    for result in report["mutations"]:
        print(f"\n[{labels[result['status']]}] {result['name']}")
        print(f"  {result['description']}")
        if result["status"] == "passed":
            if result["kind"] == "missing_stat":
                print(f"  {result['field']} is missing only from JAX output; both models "
                      f"completed {report['timesteps']} timesteps.")
                print("  Strict bindiff detected the missing variable; other stats were unchanged.")
            else:
                columns = ", ".join(str(i + 1) for i, failed in
                                    enumerate(result["column_fails"]) if failed)
                print(f"  Numerical difference detected; both models completed "
                      f"{report['timesteps']} timesteps.")
                print(f"  {result['field']}: mean absolute difference {result['mean_abs_diff']:.6g}; "
                      f"mean percentage difference {result['mean_abs_percent_diff']:.6g}%")
                print(f"  Failing columns: {columns}; first comparison failure: "
                      f"timestep {result['first_failing_timestep']}.")
        elif result["status"] == "not_run":
            print("  Not reached because the test stopped earlier.")
        else:
            print("  Did not complete all expected-failure checks; see error below.")
    if "error" in report:
        print(f"\nError: {report['error']}")
    passed = sum(result["status"] == "passed" for result in report["mutations"])
    print(f"\nResult: {labels[report['status']]} — {passed}/{len(report['mutations'])} "
          f"mutations verified | {report['elapsed_seconds']:.1f} s")
    print(f"Logs and outputs: {output}\n{rule}", flush=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, add_help=False)
    parser.add_argument("-h", "-help", action="help", help="Show this help message and exit.")
    parser.add_argument('-output_dir', dest='output_dir', type=Path,
                        help="New directory for source copies, logs and outputs (must not exist).", metavar='DIR')
    parser.add_argument("-timeout", type=float, default=300,
                        help="Maximum seconds per comparison, including compilation (default: 300).")
    parser.add_argument("-mutations", nargs="+", choices=[m.name for m in MUTATIONS],
                        help="Run selected mutations (default: all four, sharing one control).")
    args = parser.parse_args()
    if args.timeout <= 0:
        parser.error("-timeout must be positive")
    # Pin this small regression test to CPU/double; retain any custom CPU venv.
    os.environ.update(CLUBB_JAX_ACCELERATOR="cpu", CLUBB_JAX_PRECISION="double")
    os.environ["CLUBB_JAX_VENV"] = str(Path(os.environ.get(
        "CLUBB_JAX_VENV", ROOT / ".venv-jax")).resolve())
    os.environ["CLUBB_JAX_TOOLS_DIR"] = str(ROOT / ".clubb-jax-tools")
    ensure_environment()
    if args.output_dir is None:
        parent = ROOT / "output/tests"
        parent.mkdir(parents=True, exist_ok=True)
        output = Path(tempfile.mkdtemp(prefix="jax_comparison_mutation_", dir=parent))
    else:
        output = args.output_dir.resolve()
        output.mkdir(parents=True, exist_ok=False)
    print(f"Mutation test artifacts: {output}", flush=True)
    selected = [m for m in MUTATIONS if args.mutations is None or m.name in args.mutations]
    report = {"status": "failed", "case": CASE, "timesteps": TIMESTEPS_TO_RUN,
              "stats": "input/stats/multi_col_stats.in", "control_status": "not_run",
              "mutations": [
                  {"name": m.name, "description": m.description, "kind": m.kind,
                   "field": m.field, "status": "not_run", "target": str(m.target),
                   "original": m.original, "mutated": m.mutated}
                  for m in selected]}
    started = time.monotonic()
    previous_handlers = {sig: signal.signal(sig, interrupt)
                         for sig in (signal.SIGINT, signal.SIGTERM)}
    try:
        # All copies are taken before any run; each contains only its own mutation.
        control = output / "control"
        make_snapshot(control)
        for mutation, result in zip(selected, report["mutations"]):
            snapshot = output / mutation.name
            make_snapshot(snapshot)
            source = (snapshot / mutation.target).read_text()
            (snapshot / mutation.target).write_text(inject_error(source, mutation))
            result["source_sha256"] = hashlib.sha256(source.encode()).hexdigest()
        print(f"Control: running unmodified {CASE} comparison...", flush=True)
        report["control_status"] = "failed"
        clean_rc = run_comparison(control, args.timeout)
        report["control_rc"] = clean_rc
        clean_fields, clean_comparison = check_phase(control, clean_rc)
        report["control_status"] = "passed"
        report["thresholds"] = {key: clean_comparison[key] for key in
                                ("bindiff_threshold", "bindiff_percent_threshold")}
        print("Control passed.", flush=True)
        for mutation, result in zip(selected, report["mutations"]):
            print(f"Mutation {mutation.name}: running comparison...", flush=True)
            result["status"] = "failed"
            snapshot = output / mutation.name
            result["comparison_rc"] = run_comparison(snapshot, args.timeout)
            fields, comparison = check_phase(snapshot, result["comparison_rc"], mutation=mutation)
            require(all(comparison[key] == value for key, value in report["thresholds"].items()),
                    "Control and mutation must use the same comparison thresholds.")
            result["first_failing_timestep"] = comparison["first_failing_timestep"]
            if mutation.kind == "missing_stat":
                result.update(check_missing_stat(control, snapshot, mutation))
                print(f"PASS {mutation.name}: {mutation.field} absent only from JAX output; "
                      "strict comparison rejected it.", flush=True)
            else:
                result.update(check_signal(clean_fields, fields, mutation, comparison))
                print(f"PASS {mutation.name}: {mutation.field} mean absolute difference "
                      f"{result['mean_abs_diff']:.6g}; failing columns "
                      f"{[i + 1 for i, failed in enumerate(result['column_fails']) if failed]}.", flush=True)
            result["status"] = "passed"
        report["status"] = "passed"
    except (RuntimeError, OSError, ValueError, subprocess.TimeoutExpired) as exc:
        report["error"] = str(exc)
        print(f"FAIL: {exc}\nSee comparison.log and model logs under {output}", flush=True)
        return 1
    except KeyboardInterrupt as exc:
        report["status"] = "cancelled"
        report["error"] = str(exc)
        print(f"CANCELLED: {exc}", flush=True)
        return 130
    finally:
        for sig, handler in previous_handlers.items():
            signal.signal(sig, handler)
        report["elapsed_seconds"] = time.monotonic() - started
        (output / "mutation_test_summary.json").write_text(json.dumps(report, indent=2) + "\n")
        print_summary(report, output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
