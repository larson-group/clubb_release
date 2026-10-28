#!/usr/bin/env python3
"""Prepare comparison inputs and run the CLUBB loss frontend once per case.

A loss compares selected model fields with configured benchmark profiles.
This command prepares the comparison namelist and chooses the native, Python
or JAX frontend through the shared case runner, then prints its loss table.

Usage:
    ./run_scripts/run_scm_loss.py bomex
    ./run_scripts/run_scm_loss.py -jax=cpu -fields thlm,wp2 bomex
    - parse_csv_values(values, label): split nonempty CLI field/case lists.
    - case_default_overrides_from_args(args, case): collect explicit comparison
      window/height overrides for the configured case defaults.
    - build_loss_case_namelist(args, output_dir, case, fields): prepare the
      aggregate namelist and return its path and selected fields.
    - run_loss_case(args, case, fields, output_dir): wait for one frontend and
      return its success/failure status.
    - main(): parse the CLI and run every selected case.

Understanding:
    Benchmark locations come from tuner/case_defaults.json. Time ranges are
    absolute seconds and altitude ranges are meters. Each frontend invocation
    initializes, evaluates and cleans up once; this command does not retain a
    loss session. The tuner independently reuses allocations through the loss
    APIs. -disable_stats_storage keeps in-memory scoring and omits NetCDF output.
"""

import argparse
import json
import os
import sys
from pathlib import Path

CLUBB_ROOT = str(Path(__file__).resolve().parents[1])
if CLUBB_ROOT not in sys.path:
    sys.path.insert(0, CLUBB_ROOT)

from input.scripts.create_case_namelist import (
    create_loss_case_namelist,
    validate_multicol,
)
from utilities.output_paths import resolve_output_dir
from utilities.clubb_case_runner import (
    clubb_case_runner, extract_jax_options, validate_runtime_options,
)
from tuner.case_defaults import DEFAULT_LOSS_FIELDS, read_case_defaults


def parse_csv_values(values, label):
    """Return stripped, nonempty items from CLI strings, allowing commas in each.

    label names the option in the SystemExit message if no items remain.
    """
    parsed = []
    for value in values or []:
        for item in str(value).split(","):
            item = item.strip()
            if item:
                parsed.append(item)
    if not parsed:
        sys.exit(f"{label} must include at least one value")
    return parsed


def case_default_overrides_from_args(args, case_name):
    """Return case_name's comparison-window overrides from parsed args.

    Altitudes are meters and time ranges are absolute model seconds. Per-case
    config values override general CLI ranges; num_time_windows fills a missing
    case value. Benchmark-file overrides exit because case_defaults owns them.
    """
    overrides = {}
    if args.les_stats_file:
        sys.exit("LES benchmark files are configured in tuner/case_defaults.json and cannot be overridden here.")
    if args.altitude_range:
        overrides["altitude_comparison_range"] = [
            float(args.altitude_range[0]),
            float(args.altitude_range[1]),
        ]
    if args.time_range:
        overrides["time_average_range"] = [int(args.time_range[0]), int(args.time_range[1])]
    case_config = (args.case_configs or {}).get(case_name, {})
    for key in ("altitude_comparison_range", "time_average_range", "average_time_seconds", "num_time_windows"):
        if key in case_config:
            overrides[key] = case_config[key]
    if "les_stats_file" in case_config:
        sys.exit("case config files must not contain les_stats_file; edit tuner/case_defaults.json instead.")
    if args.num_time_windows is not None and "num_time_windows" not in overrides:
        overrides["num_time_windows"] = int(args.num_time_windows)
    return overrides


def build_loss_case_namelist(args, output_dir, case_name, fields):
    """Return (aggregate namelist path, selected fields) for case_name.

    Apply args' model options and case-default loss windows under output_dir.
    Empty fields selects DEFAULT_LOSS_FIELDS; driver tests disable disk stats
    while keeping in-memory loss data. Preparation errors exit with a diagnostic.
    """
    try:
        case_defaults = read_case_defaults(case_name, overrides=case_default_overrides_from_args(args, case_name))
        requested_fields = fields or list(DEFAULT_LOSS_FIELDS)
        clubb_input_namelist, requested_fields, _ = create_loss_case_namelist(
            case_name,
            output_dir,
            requested_fields,
            case_defaults=case_defaults,
            num_time_windows=case_defaults.get("num_time_windows", args.num_time_windows or 1),
            config=args.config,
            params=args.params,
            flags=args.flags,
            stats=getattr(args, "stats", None) or "all_tuning_stats",
            silhs_params=getattr(args, "silhs_params", None),
            multicol=args.multicol,
            batch_size=args.batch_size,
            disable_stats_storage=args.disable_stats_storage or args.driver_test,
            override=args.override,
            verbose=True,
        )
    except Exception as exc:
        sys.exit(str(exc))

    return clubb_input_namelist, requested_fields


def run_loss_case(args, case_name, fields, output_dir):
    """Generate and run case_name's loss namelist; return 0 for success, 1 for failure.

    args selects runtime/model settings, fields chooses comparison variables,
    and output_dir receives inputs, the console log and optional NetCDF output.
    The frontend prints its loss table through echo=True; the runner closes
    after this evaluation. Model/setup failures are reported to stderr.
    """
    clubb_input_namelist, _ = build_loss_case_namelist(
        args,
        output_dir,
        case_name,
        fields,
    )

    # Loss frontends own one complete evaluation; reusable tuning stays in the API.
    from input.scripts.create_case_namelist import case_run_name
    case_name = case_run_name(case_name)
    print(f"=================== Running {case_name} loss ===================")
    disable_stats_storage = args.disable_stats_storage or args.driver_test
    try:
        with clubb_case_runner(args, persistent=False, loss=True, echo=True) as runner:
            result = runner.run_and_wait(
                clubb_input_namelist, log_file=Path(output_dir) / f"{case_name}_log",
            )
    except (OSError, RuntimeError, ValueError) as exc:
        print(str(exc), file=sys.stderr)
        return 1
    if result['state'] != 'completed':
        print(result['error'], file=sys.stderr)
        return 1
    stats_file = Path(output_dir) / f"{case_name}_stats.nc"
    if not disable_stats_storage and not stats_file.is_file():
        print(f"WARNING: stats output not found: {stats_file}")
    return 0


def main():
    """Parse loss inputs, run the selected cases, and exit nonzero if any frontend fails."""
    parser = argparse.ArgumentParser(description="Run the in-memory CLUBB loss driver.", add_help=False, allow_abbrev=False)
    parser.add_argument("-h", "-help", action="help", help="Show this help and exit.")
    from input.scripts.create_case_namelist import add_input_arguments
    add_input_arguments(parser, cases=False)
    parser.add_argument('-output_dir', dest='out_dir', metavar="[DIR]",
        help=(
            "Output directory for generated namelists, logs, and stats. Bare names are "
            "rooted under output/; defaults to output."
        ))
    parser.add_argument("-multicol", metavar="[NUM|SPEC]", type=validate_multicol,
        help=("Generate a multi-column parameter file. "
              "Use an integer for dup_tweak mode, e.g. -multicol 4, or a parameter-grid spec like "
              "-multicol C8/0.2:0.8/4"))
    parser.add_argument("-batch_size", metavar="[NUM]", type=int,
        help=(
            "Runtime batch size written to &multicol_def. Requires -multicol. "
            "Values larger than ngrdcol are clipped to ngrdcol."
        ))
    parser.add_argument(
        "-disable_stats_storage",
        action="store_true",
        help="Disable stats NetCDF file output while keeping stats enabled in memory.",
    )
    parser.add_argument(
        "-driver_test",
        action="store_true",
        help="Run the CLUBB loss-driver test executable and force stats NetCDF output off.",
    )
    parser.add_argument(
        "-python",
        action="store_true",
        help="Run the Python loss-driver front-end instead of the compiled executable.",
    )
    parser.add_argument(
        "-override", action="append", metavar="ASSIGNMENTS|JSON",
        help=(
            "Assignments, inline JSON, or a JSON file, with all/case settings. "
            "Repeat to apply in order; later values win."
        ),
    )
    parser.add_argument(
        "-cases",
        nargs="+",
        help="Case names separated by spaces or commas, e.g. -cases atex arm.",
    )
    parser.add_argument(
        "-fields",
        help="Comma-separated CLUBB fields to compare, e.g. -fields cloud_frac,rcm",
    )
    parser.add_argument(
        "-les_stats_file",
        help="Deprecated; LES benchmark files must be configured in tuner/case_defaults.json.",
    )
    parser.add_argument(
        "-altitude_range",
        nargs=2,
        type=float,
        metavar=("Z_MIN", "Z_MAX"),
        help="Override the case-default altitude comparison range.",
    )
    parser.add_argument(
        "-time_range",
        nargs=2,
        type=int,
        metavar=("START_SEC", "END_SEC"),
        help="Override the case-default absolute time averaging range.",
    )
    parser.add_argument(
        "-time_window_mode",
        choices=("single_average", "split_average"),
        default=None,
        help="Deprecated compatibility option. Use -num_time_windows; 1 window replaces single_average.",
    )
    parser.add_argument(
        "-num_time_windows",
        type=int,
        default=None,
        help="Number of equal time subwindows. Use 1 for the former single-average behavior.",
    )
    parser.add_argument(
        "-case_config_file",
        help="JSON file containing per-case time/height/average-window overrides.",
    )
    parser.add_argument("case_name", nargs="?", help="Name of the case to run")
    parser.add_argument(
        "-jax", action="store_true",
        help="Run the JAX loss driver; accepts -jax=cpu or -jax=gpu.",
    )
    normalized, options, occurrences = extract_jax_options(sys.argv[1:])
    if occurrences > 1:
        parser.error("-jax may be specified only once")
    args = parser.parse_args(normalized)
    args.jax_options = options
    try:
        validate_runtime_options(args, loss=True)
    except ValueError as exc:
        parser.error(str(exc))

    if args.batch_size is not None and args.multicol is None:
        sys.exit("-batch_size requires -multicol")
    if args.num_time_windows is not None and args.num_time_windows < 1:
        sys.exit("-num_time_windows must be >= 1")
    if args.time_window_mode == "single_average":
        args.num_time_windows = 1
    args.case_configs = {}
    if args.case_config_file:
        try:
            raw_case_config_data = json.loads(Path(args.case_config_file).read_text(encoding="utf-8"))
        except Exception as exc:
            sys.exit(f"Could not read case config file: {exc}")
        raw_case_configs = raw_case_config_data.get("case_configs", raw_case_config_data)
        if not isinstance(raw_case_configs, list):
            sys.exit("-case_config_file must contain a list or an object with case_configs")
        for raw_config in raw_case_configs:
            case_name = str((raw_config or {}).get("case_name", "")).strip()
            if not case_name:
                sys.exit("Each case config requires case_name")
            args.case_configs[case_name] = dict(raw_config)

    try:
        output_dir = str(resolve_output_dir(args.out_dir))
    except ValueError as exc:
        parser.error(str(exc))
    os.makedirs(output_dir, exist_ok=True)
    print(f"Output directory: {output_dir}")

    if args.cases:
        case_names = parse_csv_values(args.cases, "-cases")
    elif args.case_name:
        case_names = [args.case_name]
    else:
        sys.exit("Provide either a case_name or -cases")
    fields = parse_csv_values([args.fields], "-fields") if args.fields else []

    exit_code = 0
    for case_name in case_names:
        result = run_loss_case(args, case_name, fields, output_dir)
        if result != 0:
            exit_code = result

    sys.exit(exit_code)


if __name__ == "__main__":
    from utilities.setup_python_venv import ensure_python_venv

    # Loss input preparation needs managed Python dependencies for every backend.
    ensure_python_venv()
    main()
