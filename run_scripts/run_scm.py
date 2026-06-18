#!/usr/bin/env python3
import argparse
import os
import subprocess
import sys

from create_case_namelist import validate_multicol

# Directory where this script lives, assumes clubb/run_scripts, which is important
# since this is used to find CLUBB_ROOT
RUN_SCRIPTS = os.path.dirname(os.path.abspath(__file__))
CLUBB_ROOT = os.path.join(RUN_SCRIPTS, "..")
DEFAULT_OUTPUT_DIR = os.path.join(CLUBB_ROOT, "output")
CREATE_CASE_NAMELIST = os.path.join(RUN_SCRIPTS, "create_case_namelist.py")


def run_case(run_cmd, run_cwd, case_name, namelist_file, output_dir, run_env=None, expect_stats_output=True):

    if not run_cmd:
        print("No run command was provided.")
        return 1

    os.makedirs(output_dir, exist_ok=True)

    log_path = os.path.join(output_dir, f"{case_name}_log")
    with open(log_path, "w", encoding="utf-8") as log_file:
        process = subprocess.Popen(
            run_cmd + [namelist_file],
            cwd=run_cwd,
            env=run_env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            errors="replace",
            bufsize=1,
        )

        if process.stdout is not None:
            for line in process.stdout:
                print(line, end="", flush=True)
                log_file.write(line)
            process.stdout.close()

        process.wait()

    if process.returncode not in (0, 6):
        return 1

    stats_file = os.path.join(output_dir, f"{case_name}_stats.nc")
    if expect_stats_output and not os.path.isfile(stats_file):
        print(f"WARNING: stats output not found: {stats_file}")

    return 0


def choose_run_command(args):
    run_cwd = RUN_SCRIPTS
    run_env = None

    if args.exe:
        if args.legacy:
            print(f"-legacy overriden by -exe entry: {args.exe}")
        executable = os.path.abspath(args.exe)
        if not os.path.isfile(executable):
            sys.exit(f"{executable} not found (did you re-compile?)")
        run_cmd = [executable]
    elif args.legacy:
        executable = os.path.join(CLUBB_ROOT, "bin/clubb_standalone")
        if not os.path.isfile(executable):
            sys.exit(f"{executable} not found (did you re-compile?)")
        run_cmd = [executable]
    elif args.python:
        python_driver = os.path.join(CLUBB_ROOT, "clubb_python_driver", "clubb_standalone.py")
        if not os.path.isfile(python_driver):
            sys.exit(f"Python standalone driver not found: {python_driver}")
        clubb_python_api_dir = os.path.join(CLUBB_ROOT, "clubb_python_api")
        if not os.path.isdir(clubb_python_api_dir):
            sys.exit(f"Python API directory not found: {clubb_python_api_dir}")
        executable = f"{sys.executable} -m clubb_python_driver.clubb_standalone"
        run_cmd = [sys.executable, "-m", "clubb_python_driver.clubb_standalone"]
        run_env = os.environ.copy()
        existing_pythonpath = run_env.get("PYTHONPATH", "")
        pythonpath_entries = [CLUBB_ROOT, clubb_python_api_dir]
        if existing_pythonpath:
            pythonpath_entries.append(existing_pythonpath)
        run_env["PYTHONPATH"] = os.pathsep.join(pythonpath_entries)
    elif args.jax:
        jax_driver = os.path.join(CLUBB_ROOT, "clubb_jax", "clubb_standalone.py")
        if not os.path.isfile(jax_driver):
            sys.exit(f"JAX standalone driver not found: {jax_driver}")
        clubb_python_api_dir = os.path.join(CLUBB_ROOT, "clubb_python_api")
        if not os.path.isdir(clubb_python_api_dir):
            sys.exit(f"Python API directory not found: {clubb_python_api_dir}")
        executable = f"{sys.executable} -m clubb_jax.clubb_standalone"
        run_cmd = [sys.executable, "-m", "clubb_jax.clubb_standalone"]
        run_env = os.environ.copy()
        existing_pythonpath = run_env.get("PYTHONPATH", "")
        pythonpath_entries = [CLUBB_ROOT, clubb_python_api_dir]
        if existing_pythonpath:
            pythonpath_entries.append(existing_pythonpath)
        run_env["PYTHONPATH"] = os.pathsep.join(pythonpath_entries)
    else:
        if args.driver_test:
            executable = os.path.join(CLUBB_ROOT, "install/latest/clubb_driver_test")
        else:
            executable = os.path.join(CLUBB_ROOT, "install/latest/clubb_standalone")
        if not os.path.isfile(executable):
            sys.exit(f"{executable} not found (did you re-compile?)")
        run_cmd = [executable]

    print(f" - using executable: {executable}")
    return run_cmd, run_cwd, run_env


def create_case_namelist(args, output_dir):
    if not os.path.isfile(CREATE_CASE_NAMELIST):
        sys.exit(f"{CREATE_CASE_NAMELIST} not found")

    cmd = [sys.executable, CREATE_CASE_NAMELIST, "-out_dir", output_dir]

    forwarded_opts = (
        ("-config", args.config),
        ("-params", args.params),
        ("-flags", args.flags),
        ("-silhs_params", args.silhs_params),
        ("-stats", args.stats),
        ("-multicol", args.multicol),
        ("-batch_size", args.batch_size),
        ("-zt_grid", args.zt_grid),
        ("-zm_grid", args.zm_grid),
        ("-nzmax", args.nzmax),
        ("-debug", args.debug),
        ("-max_iters", args.max_iters),
        ("-dt_main", args.dt_main),
        ("-dt_rad", args.dt_rad),
        ("-tout", args.tout),
        ("-stats_tstart", args.stats_tstart),
        ("-stats_tend", args.stats_tend),
        ("-override", args.override),
    )
    for opt, value in forwarded_opts:
        if value is not None:
            cmd.extend([opt, str(value)])

    cmd.append(args.case_name)

    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as exc:
        sys.exit(f"create_case_namelist failed (exit {exc.returncode}). Command was:\n  {' '.join(cmd)}")

    return os.path.join(output_dir, f"{args.case_name}.in")


def main():

    parser = argparse.ArgumentParser(description="Run the standalone CLUBB model")

    run_group = parser.add_argument_group("Run options handled by run_scm.py")
    namelist_group = parser.add_argument_group("Namelist options passed to create_case_namelist.py")

    namelist_group.add_argument("-config", metavar="[DIR]",
        help=("Directory containing all three tunable files:\n"
              "  tunable_parameters.in\n"
              "  configurable_model_flags.in\n"
              "  silhs_parameters.in\n"
              "Defaults to input/tunable_parameters if not given."))

    namelist_group.add_argument("-params", metavar="[FILE]",
        help=("Define the tunable parameters.\n"
              "Used to override params file defined by --config"))
    namelist_group.add_argument("-flags", metavar="[FILE]",
        help=("Model flags file.\n"
              "Used to override flags file defined by --config"))
    namelist_group.add_argument("-silhs_params", metavar="[FILE]",
        help=("SILHS parameters file.\n"
              "Used to override silhs_params file defined by --config"))

    namelist_group.add_argument("-zt_grid", metavar="[FILE]",
        help="Specify a zt grid file from input/grid.\nDefault: unused")
    namelist_group.add_argument("-zm_grid", metavar="[FILE]",
        help="Specify a zm grid file from input/grid.\nDefault: unused")
    namelist_group.add_argument("-nzmax", metavar="[NUM]", type=int,
        help="Max number of levels (required if specifying a zt/zm grid)")

    namelist_group.add_argument("-stats", metavar="[FILE]",
        help=("Stats file defining fields to output.\n"
              "Default: input/stats/standard_stats.in.\n"
              "Use 'none' to disable stats output."))

    run_group.add_argument("-exe", metavar="[EXECUTABLE]",
        help="CLUBB executable to use.\nDefault: install/clubb_standalone")

    run_group.add_argument("-driver_test", action="store_true",
        help="Runs the clubb_driver_test executable instead of clubb_standalone")

    run_group.add_argument("-python", action="store_true",
        help="Run the Python standalone driver (python -m clubb_python_driver.clubb_standalone)")

    run_group.add_argument("-jax", action="store_true",
        help="Run the JAX standalone driver (python -m clubb_jax.clubb_standalone)")

    run_group.add_argument(
        "-legacy",
        action="store_true",
        help="Runs the legacy compiled version of clubb_standalone (with compile.bash)"
    )

    run_group.add_argument("-out_dir", metavar="[DIR]",
        help="Output directory for results.\nDefault: output")

    namelist_group.add_argument("-debug", metavar="[NUM]",
        help="Debug level (0–3) that controls CLUBB's runtime checks (0 is no checks).\nDefault specified in model file.")
    namelist_group.add_argument("-max_iters", metavar="[NUM]", type=int,
        help="Maximum number of iterations")
    namelist_group.add_argument("-dt_main", metavar="[SECONDS]", type=int,
        help="Main timestep (s).\nDefault from model file.")
    namelist_group.add_argument("-dt_rad", metavar="[SECONDS]", type=int,
        help="Radiation timestep (s).\nDefault from model file.")
    namelist_group.add_argument("-tout", metavar="[SECONDS]", type=int,
        help="Stats output interval (s). Use 0 to disable.\nDefault from model file.")
    namelist_group.add_argument("-stats_tstart", metavar="[SECONDS]", type=float,
        help="Stats output window start time (s). Default from model file or driver.")
    namelist_group.add_argument("-stats_tend", metavar="[SECONDS]", type=float,
        help="Stats output window end time (s). Default from model file or driver.")

    namelist_group.add_argument("-multicol", metavar="[NUM|SPEC]", type=validate_multicol,
        help=("Generate a multi-column parameter file. "
              "Use an integer for dup_tweak mode, e.g. -multicol 4, or an hr spec like "
              "-multicol C8/0.2:0.8/4"))
    namelist_group.add_argument("-batch_size", metavar="[NUM]", type=int,
        help="Runtime batch size written to &multicol_def. Requires -multicol.")

    namelist_group.add_argument(
        "-override",
        help="Comma-separated key=value pairs, e.g. -override FLAG1=true,C2=2.0,...",
    )

    parser.add_argument("case_name", help="Name of the case to run")
    args = parser.parse_args()

    ndefined = sum(bool(x) for x in [args.exe, args.legacy, args.driver_test, args.python, args.jax])
    if ndefined > 1:
        parser.error("Only one of -exe, -legacy, -driver_test, -python, or -jax may be specified.")
    if args.zt_grid and args.zm_grid:
        sys.exit("\n\033[91mERROR: Cannot specify both a ZT grid and a ZM grid\033[0m")
    if args.nzmax and not (args.zt_grid or args.zm_grid):
        print("\n\033[93mWARNING: Specifying --nzmax will have no effect without specifying a --zm_grid or --zt_grid\033[0m")
    if args.batch_size is not None and args.multicol is None:
        parser.error("-batch_size requires -multicol.")
    if (args.stats_tstart is None) != (args.stats_tend is None):
        parser.error("-stats_tstart and -stats_tend must be provided together.")

    output_dir = os.path.abspath(args.out_dir) if args.out_dir else os.path.abspath(DEFAULT_OUTPUT_DIR)
    os.makedirs(output_dir, exist_ok=True)

    clubb_input_namelist = create_case_namelist(args, output_dir)
    run_cmd, run_cwd, run_env = choose_run_command(args)

    print(f"=================== Running {args.case_name} ===================")
    expect_stats_output = ((args.stats or "").strip().lower() != "none") and args.tout != 0
    result = run_case(run_cmd, run_cwd, args.case_name, clubb_input_namelist, output_dir, run_env,
                      expect_stats_output=expect_stats_output)

    sys.exit(result)


if __name__ == "__main__":
    main()
