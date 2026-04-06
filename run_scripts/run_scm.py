#!/usr/bin/env python3
import argparse
import os
import re
import subprocess
import sys

# Directory where this script lives, assumes clubb/run_scripts, which is important
# since this is used to find CLUBB_ROOT
RUN_SCRIPTS = os.path.dirname(os.path.abspath(__file__))
CLUBB_ROOT  = os.path.join(RUN_SCRIPTS, "..") 
DEFAULT_OUTPUT_DIR  = os.path.join(CLUBB_ROOT, "output")

multi_col_params_script = os.path.join(RUN_SCRIPTS, "create_multi_col_params.py") 

HR_SPEC_RE = re.compile(
    r"^[A-Za-z_]\w*/[-+]?\d*\.?\d+(?:[eEdD][-+]?\d+)?:[-+]?\d*\.?\d+(?:[eEdD][-+]?\d+)?/\d+"
    r"(?:,[A-Za-z_]\w*/[-+]?\d*\.?\d+(?:[eEdD][-+]?\d+)?:[-+]?\d*\.?\d+(?:[eEdD][-+]?\d+)?/\d+)*$"
)


def strip_comments_and_remove_keys(content: str, keys_to_remove=None) -> str:
    """Remove Fortran namelist comments (!) and specified keys."""
    if keys_to_remove is None:
        keys_to_remove = []
    lines = []
    for line in content.splitlines():
        line = re.sub(r"!.*", "", line)  # strip comments
        if any(key in line for key in keys_to_remove):
            continue
        lines.append(line)
    return "\n".join(lines)

def run_case(run_cmd, run_cwd, case_name, namelist_file, output_dir, run_env=None):

    if not run_cmd:
        print("No run command was provided.")
        return 1

    # clubb requires the output directory to exist prior to running
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

        # Stream model output to terminal and file at the same time.
        if process.stdout is not None:
            for line in process.stdout:
                print(line, end="", flush=True)
                log_file.write(line)
            process.stdout.close()

        process.wait()

    if process.returncode not in (0, 6):
        return 1

    stats_file = os.path.join(output_dir, f"{case_name}_stats.nc")
    if not os.path.isfile(stats_file):
        print(f"WARNING: stats output not found: {stats_file}")

    return 0


def read_model_times(model_file):
    """Read time_initial, time_final, dt_main from a model file if present."""
    values = {}
    line_re = re.compile(r'^\s*([a-zA-Z_]\w*)\s*=\s*([-+0-9.eEdD]+)')
    with open(model_file) as f:
        for line in f:
            m = line_re.match(line)
            if m:
                key, val = m.groups()
                val = val.replace("D", "E").replace("d", "e")
                try:
                    values[key.lower()] = float(val)
                except ValueError:
                    pass
    return values


def parse_multicol_arg(value: str):
    """Parse -multicol as either an integer ngrdcol or an hr-style spec string."""
    stripped = value.strip()
    if not stripped:
        raise argparse.ArgumentTypeError("must provide a non-empty -multicol value")

    try:
        ngrdcol = int(stripped)
    except ValueError:
        if HR_SPEC_RE.fullmatch(stripped):
            return {"ngrdcol": None, "hr_spec": stripped}
        raise argparse.ArgumentTypeError(
            "must be either an integer column count or an hr spec like C8/0.2:0.8/4"
        )

    if ngrdcol < 1:
        raise argparse.ArgumentTypeError("integer -multicol value must be >= 1")

    return {"ngrdcol": ngrdcol, "hr_spec": None}


def convert_to_multi_col(
    params_file: str,
    case_name: str,
    output_dir: str,
    ngrdcol: int | None,
    hr_spec: str | None = None,
) -> str:
    """
    Create a temporary multi-column params file and return its path.

    If ``hr_spec`` is provided, it is forwarded directly to
    ``create_multi_col_params.py -hr``. Otherwise ``ngrdcol`` uses the legacy
    dup_tweak path.
    """

    if not os.path.isfile(multi_col_params_script):
        sys.exit(f"Missing helper script: {multi_col_params_script}")

    out_file = os.path.join(output_dir, f"{case_name}_multicol_params.in")

    cmd = [sys.executable, multi_col_params_script, "-param_file", params_file, "-out_file", out_file]
    if hr_spec:
        cmd.extend(["-hr", hr_spec])
    else:
        cmd.extend(["-mode", "dup_tweak", "-n", str(ngrdcol)])

    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        sys.exit(f"create_multi_col_params failed (exit {e.returncode}). Command was:\n  {' '.join(cmd)}")

    if not os.path.isfile(out_file):
        sys.exit(f"Expected output file not created: {out_file}")

    return out_file


def override_value(override_string, clubb_in_text):
    """
    Apply overrides from -override KEY1=val1,KEY2=val2,... to the aggregate text.
    """
    for pair in override_string.split(","):
        if "=" not in pair:
            continue
        key, val = pair.split("=", 1)
        key, val = key.strip(), val.strip()
        replacement = f"{key} = {val}"

        assignment_re = re.compile(
            rf"(?im)^(\s*){re.escape(key)}\s*=.*$"
        )

        clubb_in_text, replacements = assignment_re.subn(
            lambda match: f"{match.group(1)}{key} = {val}",
            clubb_in_text,
        )

        if replacements == 0:
            clubb_in_text += f"\n{replacement}\n"
    return clubb_in_text


def setup_files_and_aggregate(args, output_dir):
    """Resolve file paths, validate, and create the aggregate namelist."""

    # Model file
    model_file = os.path.join(CLUBB_ROOT, f"input/case_setups/{args.case_name}_model.in")
    if not os.path.isfile(model_file):
        sys.exit(f"{model_file} does not exist")

    # Config dir, default is input/tunable_parameters
    config_dir = os.path.abspath(args.config) if args.config else os.path.join(CLUBB_ROOT, "input/tunable_parameters")

    if not os.path.isdir(config_dir):
        sys.exit(f"--config directory does not exist: {config_dir}")

    # Files (respect overrides)
    params_file       = args.params       or os.path.join(config_dir, "tunable_parameters.in")
    flags_file        = args.flags        or os.path.join(config_dir, "configurable_model_flags.in")
    silhs_params_file = args.silhs_params or os.path.join(config_dir, "silhs_parameters.in")
    stats_arg = (args.stats or "").strip()
    disable_stats = stats_arg.lower() == "none"
    stats_file = None if disable_stats else (args.stats or os.path.join(CLUBB_ROOT, "input/stats/standard_stats.in"))

    run_cwd = RUN_SCRIPTS
    run_env = None
    run_cmd = None

    if args.exe:
        # Use the users input from -exe to determine which executable to use
        if args.legacy:
            print(f"-legacy overriden by -exe entry: {args.exe}")
        executable  = os.path.abspath(args.exe)
        if not os.path.isfile(executable):
            sys.exit(f"{executable} not found (did you re-compile?)")
        run_cmd = [executable]
    elif args.legacy:
        # The legacy install location is /bin/clubb_standalone
        executable  = os.path.join(CLUBB_ROOT, f"bin/clubb_standalone")
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
        run_cwd = RUN_SCRIPTS
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
        run_cwd = RUN_SCRIPTS
        executable = f"{sys.executable} -m clubb_jax.clubb_standalone"
        run_cmd = [sys.executable, "-m", "clubb_jax.clubb_standalone"]
        run_env = os.environ.copy()
        existing_pythonpath = run_env.get("PYTHONPATH", "")
        pythonpath_entries = [CLUBB_ROOT, clubb_python_api_dir]
        if existing_pythonpath:
            pythonpath_entries.append(existing_pythonpath)
        run_env["PYTHONPATH"] = os.pathsep.join(pythonpath_entries)
    else:
        # Find the executable based on compiler by default (install/COMPILER/clubb_standalone), otherwise
        # default to (install/lastest/clubb_standalone) which points to the last one compiled
        #compiler    = os.path.basename(os.environ["FC"]) if "FC" in os.environ else "latest"

        if args.driver_test:
            executable  = os.path.join(CLUBB_ROOT, f"install/latest/clubb_driver_test")
        else:
            executable  = os.path.join(CLUBB_ROOT, f"install/latest/clubb_standalone")
        if not os.path.isfile(executable):
            sys.exit(f"{executable} not found (did you re-compile?)")
        run_cmd = [executable]

    print(f" - using executable: {executable}")

    # Expand multi-column params if requested
    if args.multicol is not None:
        params_file = convert_to_multi_col(
            params_file,
            args.case_name,
            output_dir,
            args.multicol["ngrdcol"],
            args.multicol["hr_spec"],
        )

    # Validate files
    files_to_validate = [
        ("-params", params_file),
        ("-flags", flags_file),
        ("-silhs_params", silhs_params_file),
    ]
    if not (args.python or args.jax):
        files_to_validate.append(("-exe", executable))
    if not disable_stats:
        files_to_validate.append(("-stats", stats_file))

    for opt, f in files_to_validate:
        if not os.path.isfile(f):
            sys.exit(f"Required file for {opt} not found: {f}")

    # Aggregate into <output_dir>/CASE.in
    os.makedirs(output_dir, exist_ok=True)
    clubb_input_namelist = os.path.join(output_dir, f"{args.case_name}.in")
    with open(clubb_input_namelist, "w") as out:
        files_to_aggregate = [params_file, silhs_params_file, flags_file, model_file]
        if not disable_stats:
            files_to_aggregate.append(stats_file)
        for f in files_to_aggregate:
            with open(f) as src:
                out.write(strip_comments_and_remove_keys(src.read()))
                out.write("\n")

    return clubb_input_namelist, model_file, run_cmd, run_cwd, run_env


def _set_stats_output_dir(clubb_in: str, output_dir: str) -> str:
    """Ensure &stats_setting contains output_dir."""
    out_norm = output_dir.replace("\\", "/")
    stats_match = re.search(r"(?is)(&\s*stats_setting\b)(.*?)(/)", clubb_in)
    if not stats_match:
        clubb_in += (
            "\n&stats_setting\n"
            f"output_dir = '{out_norm}',\n"
            "/\n"
        )
        return clubb_in

    header = stats_match.group(1)
    body = stats_match.group(2)
    end = stats_match.group(3)
    if re.search(r"(?im)^\s*output_dir\s*=", body):
        body = re.sub(
            r"(?im)^\s*output_dir\s*=.*$",
            f"output_dir = '{out_norm}',",
            body,
        )
    else:
        body = body.rstrip() + f"\noutput_dir = '{out_norm}',\n"

    return clubb_in[:stats_match.start()] + header + body + end + clubb_in[stats_match.end():]


def edit_namelist(args, clubb_input_namelist, model_file, output_dir):
    """Apply modifications to the aggregate namelist."""

    with open(clubb_input_namelist) as f:
        clubb_in = f.read()

    # Timestep modifications
    if args.dt_main is not None:
        clubb_in = re.sub(r"dt_main\s*=.*", f"dt_main = {args.dt_main}.0", clubb_in)

    if args.dt_rad is not None:
        clubb_in = re.sub(r"dt_rad\s*=.*", f"dt_rad = {args.dt_rad}.0", clubb_in)

    # Explicit grids
    if args.zt_grid is not None:
        clubb_in += f"\nnzmax = {args.nzmax}\nzt_grid_fname = '{args.zt_grid}'\ngrid_type = 2\n"

    if args.zm_grid is not None:
        clubb_in += f"\nnzmax = {args.nzmax}\nzm_grid_fname = '{args.zm_grid}'\ngrid_type = 3\n"

    # Stats output control, args.tout defines output frequency, and setting to 0 disables output.
    disable_stats = ((args.stats or "").strip().lower() == "none")
    if args.tout is not None:
        if args.tout == 0:
            disable_stats = True
        else:
            clubb_in = re.sub(r"stats_tout\s*=.*", f"stats_tout = {args.tout}.0", clubb_in)

    if disable_stats:
        if re.search(r"l_stats\s*=.*", clubb_in):
            clubb_in = re.sub(r"l_stats\s*=.*", "l_stats = .false.", clubb_in)
        else:
            clubb_in += "\nl_stats = .false.\n"

    # Debug level
    if args.debug is not None:
        clubb_in = re.sub(r"debug_level\s*=.*", f"debug_level = {args.debug}", clubb_in)

    # Iteration control
    if args.max_iters is not None:
        vals = read_model_times(model_file)
        time_initial = vals["time_initial"]
        dt_main_val = args.dt_main if args.dt_main else vals.get("dt_main")
        new_time_final = time_initial + dt_main_val * args.max_iters

        # only update time_final if it's less than the current one\
        if ( vals["time_final"] >= new_time_final ):
            clubb_in = re.sub(r"time_final\s*=.*", f"time_final = {new_time_final}", clubb_in)

    # Overrides
    if args.override:
        clubb_in = override_value(args.override, clubb_in)

    # Route all CLUBB output files into the selected directory.
    clubb_in = _set_stats_output_dir(clubb_in, output_dir)

    # Save back
    with open(clubb_input_namelist, "w") as f:
        f.write(clubb_in)


def main():

    parser = argparse.ArgumentParser( description="Run the standalone CLUBB model")

    # The "config" folder is assumed to contain the tunable parameters, silhs parameters, and config flags.
    # This is meant to be a quick way to input all these files without having to specify them individually.
    # Individually specificed files will overwrite the ones found from here
    parser.add_argument("-config", metavar="[DIR]",
        help=("Directory containing all three tunable files:\n"
            "  tunable_parameters.in\n"
            "  configurable_model_flags.in\n"
            "  silhs_parameters.in\n"
            "Defaults to input/tunable_parameters if not given."))

    # per-file overrides (params/flags already existed; add silhs_params)
    parser.add_argument("-params", metavar="[FILE]",
        help=("Define the tunable parameters.\n"
            "Used to override params file defined by --config"))
    parser.add_argument("-flags", metavar="[FILE]",
        help=("Model flags file.\n"
            "Used to override flags file defined by --config"))
    parser.add_argument("-silhs_params", metavar="[FILE]",
        help=("SILHS parameters file.\n"
            "Used to override silhs_params file defined by --config"))

    # Options to input grid files, and nzmax defines the maximum number of vertical levels when inputting grid files like this
    parser.add_argument("-zt_grid", metavar="[FILE]",
        help="Specify a zt grid file from input/grid.\nDefault: unused")
    parser.add_argument("-zm_grid", metavar="[FILE]",
        help="Specify a zm grid file from input/grid.\nDefault: unused")
    parser.add_argument("-nzmax", metavar="[NUM]", type=int,
        help="Max number of levels (required if specifying a zt/zm grid)")

    # Stats fiules define which fields to output, can be overriden here
    parser.add_argument("-stats", metavar="[FILE]",
        help=("Stats file defining fields to output.\n"
              "Default: input/stats/standard_stats.in.\n"
              "Use 'none' to disable stats output."))

    # This script will try to figure out the right executable to use based on the compiler in the environment
    # but inputting a specific executable will override that with the specified one
    parser.add_argument("-exe", metavar="[EXECUTABLE]",
        help="CLUBB executable to use.\nDefault: install/clubb_standalone")
        
    parser.add_argument("-driver_test", action="store_true",
        help="Runs the clubb_driver_test executable instead of clubb_standalone"
    )

    parser.add_argument("-python", action="store_true",
        help="Run the Python standalone driver (python -m clubb_python_driver.clubb_standalone)")

    parser.add_argument("-jax", action="store_true",
        help="Run the JAX standalone driver (python -m clubb_jax.clubb_standalone)")

    # The old method of compile clubb resulted in the executable "clubb/bin/clubb_standalone"
    # this option causes that to be the prefered executable, unless -exe is specified
    parser.add_argument(
        "-legacy",
        action="store_true",
        help="Runs the legacy compiled version of clubb_standalone (with compile.bash)"
    )

    # Allow a custom output directory to be used for all generated files.
    parser.add_argument("-out_dir", metavar="[DIR]",
        help="Output directory for results.\nDefault: output")

    # Runtime options
    parser.add_argument("-debug", metavar="[NUM]",
        help="Debug level (0–3) that controls CLUBB's runtime checks (0 is no checks).\nDefault specified in model file.")
    parser.add_argument("-max_iters", metavar="[NUM]", type=int,
        help="Maximum number of iterations")
    parser.add_argument("-dt_main", metavar="[SECONDS]", type=int,
        help="Main timestep (s).\nDefault from model file.")
    parser.add_argument("-dt_rad", metavar="[SECONDS]", type=int,
        help="Radiation timestep (s).\nDefault from model file.")
    parser.add_argument("-tout", metavar="[SECONDS]", type=int,
        help="Stats output interval (s). Use 0 to disable.\nDefault from model file.")

    # Setting -multicol will call create_multi_col_params.py to generate a multi-column parameter file.
    # Integer input uses the legacy dup_tweak path. A string matching the hr syntax is forwarded to
    # create_multi_col_params.py -hr for hypergrid generation.
    parser.add_argument("-multicol", metavar="[NUM|SPEC]", type=parse_multicol_arg,
        help=("Generate a multi-column parameter file. "
              "Use an integer for dup_tweak mode, e.g. -multicol 4, or an hr spec like "
              "-multicol C8/0.2:0.8/4"))

    # This can be used to override pretty much any settings in the aggregate namelist
    parser.add_argument(
        "-override",
        help="Comma-separated key=value pairs, e.g. -override FLAG1=true,C2=2.0,...",
    )

    parser.add_argument("case_name", help="Name of the case to run")
    args = parser.parse_args()

    # Error check
    ndefined = sum(bool(x) for x in [args.exe, args.legacy, args.driver_test, args.python, args.jax])
    if ndefined > 1:
        parser.error("Only one of -exe, -legacy, -driver_test, -python, or -jax may be specified.")
    # Validate grid options
    if args.zt_grid and args.zm_grid:
        sys.exit(f"\n\033[91mERROR: Cannot specify both a ZT grid and a ZM grid\033[0m")
    if args.nzmax and not (args.zt_grid or args.zm_grid):
        print("\n\033[93mWARNING: Specifying --nzmax will have no effect without specifying a --zm_grid or --zt_grid\033[0m")

    output_dir = os.path.abspath(args.out_dir) if args.out_dir else os.path.abspath(DEFAULT_OUTPUT_DIR)
    os.makedirs(output_dir, exist_ok=True)

    # Step 1: setup and aggregate namelist files into <output_dir>/CASE.in
    clubb_input_namelist, model_file, run_cmd, run_cwd, run_env = setup_files_and_aggregate(args, output_dir)

    # Step 2: edit clubb_input_namelist based on input specifications
    edit_namelist(args, clubb_input_namelist, model_file, output_dir)

    # Step 3: run model
    print(f"=================== Running {args.case_name} ===================")
    result = run_case(run_cmd, run_cwd, args.case_name, clubb_input_namelist, output_dir, run_env)

    sys.exit(result)


if __name__ == "__main__":
    main()
