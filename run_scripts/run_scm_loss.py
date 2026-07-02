#!/usr/bin/env python3
import argparse
import json
import os
import sys
from pathlib import Path

if os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")) not in sys.path:
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")))

from utilities.create_case_namelist import (
    create_loss_case_namelist,
    validate_multicol,
)
from run_scm import choose_install_dir, python_runtime_dir_from_install, run_case
from tuner.case_defaults import DEFAULT_LOSS_FIELDS, read_case_defaults

RUN_SCRIPTS = os.path.dirname(os.path.abspath(__file__))
CLUBB_ROOT = os.path.join(RUN_SCRIPTS, "..")
DEFAULT_OUTPUT_DIR = os.path.join(CLUBB_ROOT, "output")


def choose_run_command(args):
    if args.python:
        module_name = "tuner.clubb_loss_driver_test" if args.driver_test else "tuner.clubb_loss_driver"
        python_driver = os.path.join(CLUBB_ROOT, "tuner", f"{module_name.split('.')[-1]}.py")
        if not os.path.isfile(python_driver):
            sys.exit(f"Python loss driver not found: {python_driver}")

        executable = f"{sys.executable} -m {module_name}"
        run_cmd = [sys.executable, "-m", module_name]
        run_env = os.environ.copy()
        existing_pythonpath = run_env.get("PYTHONPATH", "")
        install_dir, install_source = choose_install_dir()
        f2py_runtime_dir = python_runtime_dir_from_install(install_dir)
        pythonpath_entries = [f2py_runtime_dir, CLUBB_ROOT]
        if existing_pythonpath:
            pythonpath_entries.append(existing_pythonpath)
        run_env["PYTHONPATH"] = os.pathsep.join(pythonpath_entries)
        print(f" - using Python runtime dir ({install_source}): {os.path.realpath(f2py_runtime_dir)}")
        print(f" - using executable: {executable}")
        return run_cmd, RUN_SCRIPTS, run_env

    executable_name = "clubb_loss_driver_test" if args.driver_test else "clubb_standalone_loss"
    executable = os.path.join(CLUBB_ROOT, f"install/latest/{executable_name}")

    if not os.path.isfile(executable):
        sys.exit(f"{executable} not found (did you re-compile?)")

    print(f" - using executable: {executable}")
    return [executable], RUN_SCRIPTS, None


def parse_csv_values(values, label):
    """Parse one or more comma-separated CLI values."""
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
    """Return explicit CLI overrides for case comparison defaults."""
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
    """Build the loss-driver aggregate namelist using the shared namelist builder."""
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
            multicol=args.multicol,
            batch_size=args.batch_size,
            disable_stats_storage=args.disable_stats_storage or args.driver_test,
            override=args.override,
            verbose=True,
        )
    except Exception as exc:
        sys.exit(str(exc))

    return clubb_input_namelist, requested_fields


def run_loss_case(args, case_name, fields, output_dir, run_cmd, run_cwd, run_env):
    model_file = os.path.join(CLUBB_ROOT, f"input/case_setups/{case_name}_model.in")
    if not os.path.isfile(model_file):
        sys.exit(f"{model_file} does not exist")

    clubb_input_namelist, _ = build_loss_case_namelist(
        args,
        output_dir,
        case_name,
        fields,
    )

    print(f"=================== Running {case_name} loss ===================")
    disable_stats_storage = args.disable_stats_storage or args.driver_test
    return run_case(
        run_cmd,
        run_cwd,
        case_name,
        str(clubb_input_namelist),
        output_dir,
        run_env,
        expect_stats_output=not disable_stats_storage,
    )


def main():
    parser = argparse.ArgumentParser(description="Run the in-memory CLUBB loss driver.")
    parser.add_argument("-config", metavar="[DIR]",
        help=("Directory containing tunable_parameters.in, configurable_model_flags.in, "
              "and silhs_parameters.in. Defaults to input/tunable_parameters."))
    parser.add_argument("-params", metavar="[FILE]",
        help="Define the tunable parameters. Used to override params file defined by --config")
    parser.add_argument("-flags", metavar="[FILE]",
        help="Model flags file. Used to override flags file defined by --config")
    parser.add_argument("-out_dir", metavar="[DIR]",
        help="Output directory for generated namelists, logs, and stats. Defaults to output.")
    parser.add_argument("-multicol", metavar="[NUM|SPEC]", type=validate_multicol,
        help=("Generate a multi-column parameter file. "
              "Use an integer for dup_tweak mode, e.g. -multicol 4, or an hr spec like "
              "-multicol C8/0.2:0.8/4"))
    parser.add_argument("-batch_size", metavar="[NUM]", type=int,
        help="Runtime batch size written to &multicol_def. Requires -multicol.")
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
        "-override",
        help=(
            "Comma-separated key=value overrides. Values may be column lists when used "
            "with -multicol, e.g. -override C8=0.8,0.7,C11=1.0,1.1"
        ),
    )
    parser.add_argument(
        "-cases",
        nargs="+",
        help="Comma-separated case names to run, e.g. -cases atex,arm",
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
    args = parser.parse_args()

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

    output_dir = os.path.abspath(args.out_dir) if args.out_dir else os.path.abspath(DEFAULT_OUTPUT_DIR)
    os.makedirs(output_dir, exist_ok=True)

    run_cmd, run_cwd, run_env = choose_run_command(args)

    if args.cases:
        case_names = parse_csv_values(args.cases, "-cases")
    elif args.case_name:
        case_names = [args.case_name]
    else:
        sys.exit("Provide either a case_name or -cases")
    fields = parse_csv_values([args.fields], "-fields") if args.fields else []

    exit_code = 0
    for case_name in case_names:
        result = run_loss_case(args, case_name, fields, output_dir, run_cmd, run_cwd, run_env)
        if result != 0:
            exit_code = result

    sys.exit(exit_code)


if __name__ == "__main__":
    main()
