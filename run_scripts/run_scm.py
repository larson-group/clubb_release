#!/usr/bin/env python3
import argparse
import glob
import os
import shutil
import subprocess
import sys
import time

# Directory where this script lives, assumes clubb/run_scripts, which is important
# since this is used to find CLUBB_ROOT
RUN_SCRIPTS = os.path.dirname(os.path.abspath(__file__))
CLUBB_ROOT = os.path.join(RUN_SCRIPTS, "..")
if CLUBB_ROOT not in sys.path:
    sys.path.insert(0, CLUBB_ROOT)

from utilities.create_case_namelist import add_namelist_arguments  # noqa: E402
from utilities.output_paths import resolve_output_dir  # noqa: E402

CREATE_CASE_NAMELIST = os.path.join(CLUBB_ROOT, "utilities", "create_case_namelist.py")
INSTALL_DIR = os.path.join(CLUBB_ROOT, "install")
SELECTED_INSTALL = os.path.join(INSTALL_DIR, "selected")
LATEST_INSTALL = os.path.join(INSTALL_DIR, "latest")


def extract_jax_options(argv):
    # argparse's nargs="?" would consume CASE in "-jax CASE". Extract only
    # attached values (-jax=VALUE); the JAX launcher interprets their contents.
    normalized = []
    value = None
    occurrences = 0
    for token in argv:
        option, separator, attached = token.partition("=")
        if option == "-jax":
            occurrences += 1
            value = attached if separator else None
            normalized.append("-jax")
        else:
            normalized.append(token)
    return normalized, value, occurrences


def run_case(
    run_cmd,
    run_cwd,
    case_name,
    namelist_file,
    output_dir,
    run_env=None,
    expect_stats_output=True,
    passthrough_output=False,
):

    if not run_cmd:
        print("No run command was provided.")
        return 1

    os.makedirs(output_dir, exist_ok=True)

    if passthrough_output:
        result = subprocess.run(
            run_cmd + [namelist_file],
            cwd=run_cwd,
            env=run_env,
            check=False,
        )
        return 0 if result.returncode in (0, 6) else 1

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


def choose_install_dir(args=None):
    """Choose the install tree used by run commands."""
    install_dir = getattr(args, "install_dir", None)
    if install_dir:
        return os.path.abspath(install_dir), "explicit"
    if os.path.lexists(SELECTED_INSTALL):
        return SELECTED_INSTALL, "selected"
    return LATEST_INSTALL, "latest"


def python_runtime_dir_from_install(install_dir):
    """Return the self-contained Python runtime directory for an install tree."""
    runtime_dir = os.path.join(install_dir, "python")
    missing_reason = None
    backend_candidates = [
        os.path.join(runtime_dir, "libclubb_f2py_backend.so"),
        os.path.join(runtime_dir, "libclubb_f2py_backend.dylib"),
        os.path.join(runtime_dir, "clubb_f2py_backend.dll"),
    ]

    if not os.path.isdir(runtime_dir):
        missing_reason = "directory is missing"
    elif not glob.glob(os.path.join(runtime_dir, "clubb_f2py*.so")):
        missing_reason = "clubb_f2py extension is missing"
    elif not any(os.path.isfile(candidate) for candidate in backend_candidates):
        missing_reason = "libclubb_f2py_backend shared library is missing"
    elif not os.path.isdir(os.path.join(runtime_dir, "clubb_python")):
        missing_reason = "clubb_python package is missing"

    if missing_reason:
        sys.exit(
            f"{runtime_dir} is not a complete Python runtime ({missing_reason}). "
            "Rebuild that install with ./compile.py -python or choose a Python-enabled install."
        )

    return runtime_dir


def choose_run_command(args):
    run_cwd = RUN_SCRIPTS
    run_env = None
    install_dir = None
    install_source = None
    show_install_dir = False
    f2py_runtime_dir = None

    if args.exe:
        executable = os.path.abspath(args.exe)
        if not os.path.isfile(executable):
            sys.exit(f"{executable} not found (did you re-compile?)")
        run_cmd = [executable]
    elif args.python:
        python_driver = os.path.join(CLUBB_ROOT, "clubb_python_driver", "clubb_standalone.py")
        if not os.path.isfile(python_driver):
            sys.exit(f"Python standalone driver not found: {python_driver}")
        install_dir, install_source = choose_install_dir(args)
        executable = f"{sys.executable} -m clubb_python_driver.clubb_standalone"
        run_cmd = [sys.executable, "-m", "clubb_python_driver.clubb_standalone"]
        run_env = os.environ.copy()
        existing_pythonpath = run_env.get("PYTHONPATH", "")
        f2py_runtime_dir = python_runtime_dir_from_install(install_dir)
        pythonpath_entries = [f2py_runtime_dir, CLUBB_ROOT]
        if existing_pythonpath:
            pythonpath_entries.append(existing_pythonpath)
        run_env["PYTHONPATH"] = os.pathsep.join(pythonpath_entries)
    elif args.jax:
        jax_launcher = os.path.join(CLUBB_ROOT, "clubb_jax", "run_jax.py")
        if not os.path.isfile(jax_launcher):
            sys.exit(f"JAX launcher not found: {jax_launcher}")
        executable = jax_launcher
        run_cmd = [jax_launcher]
        if args.jax_options is not None:
            # Keep profile/modifier parsing in the launcher so CLI and Dash
            # share its validation and environment setup rules.
            run_cmd.append(f"-options={args.jax_options}")
    else:
        install_dir, install_source = choose_install_dir(args)
        show_install_dir = True
        if os.path.islink(install_dir) and not os.path.exists(install_dir):
            sys.exit(f"{install_dir} points to a missing install directory.")
        if not os.path.isdir(install_dir):
            sys.exit(f"{install_dir} is not an install directory (did you re-compile?)")
        executable_name = "clubb_driver_test" if args.driver_test else "clubb_standalone"
        executable = os.path.join(install_dir, executable_name)
        if not os.path.isfile(executable):
            sys.exit(f"{executable} not found (did you re-compile?)")
        run_cmd = [executable]

    if args.gdb:
        gdb_path = shutil.which("gdb")
        if gdb_path is None:
            sys.exit("gdb not found on PATH")
        run_cmd = [gdb_path, "--args", *run_cmd]

    if install_dir and (show_install_dir or install_source == "explicit"):
        print(f" - using install dir ({install_source}): {os.path.realpath(install_dir)}")
    if f2py_runtime_dir:
        print(f" - using Python runtime dir ({install_source}): {os.path.realpath(f2py_runtime_dir)}")
    print(f" - using executable: {executable}")
    if args.gdb:
        print(f" - launching with debugger: {gdb_path}")
    return run_cmd, run_cwd, run_env


def create_case_namelist(args, output_dir):
    if not os.path.isfile(CREATE_CASE_NAMELIST):
        sys.exit(f"{CREATE_CASE_NAMELIST} not found")

    cmd = [sys.executable, CREATE_CASE_NAMELIST, "-output_dir", output_dir]

    forwarded_opts = (
        ("-config", args.config),
        ("-params_file", args.params),
        ("-flags", args.flags),
        ("-silhs_params_file", args.silhs_params),
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
    )
    for opt, value in forwarded_opts:
        if value is not None:
            cmd.extend([opt, str(value)])

    # Preserve repeated overrides for the namelist generator to resolve in order.
    overrides = args.override if isinstance(args.override, list) else [args.override]
    for value in overrides:
        if value is not None:
            cmd.extend(["-override", value])
    cmd.append(args.case_name)

    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as exc:
        sys.exit(f"create_case_namelist failed (exit {exc.returncode}). Command was:\n  {' '.join(cmd)}")

    return os.path.join(output_dir, f"{args.case_name}.in")


def main():

    parser = argparse.ArgumentParser(description="Run the standalone CLUBB model", add_help=False, allow_abbrev=False)
    parser.add_argument("-h", "-help", action="help", help="Show this help and exit.")

    run_group = parser.add_argument_group("Run options handled by run_scm.py")
    namelist_group = parser.add_argument_group("Model settings passed to create_case_namelist.py")
    add_namelist_arguments(namelist_group)

    run_group.add_argument("-exe", metavar="[EXECUTABLE]",
        help="CLUBB executable to use. Overrides -install_dir and selected/latest install dirs.")

    run_group.add_argument("-install_dir", metavar="[DIR]",
        help="Install directory containing CLUBB executables.\nDefault: install/selected if present, otherwise install/latest.")

    run_group.add_argument("-driver_test", action="store_true",
        help="Runs the clubb_driver_test executable instead of clubb_standalone")

    run_group.add_argument("-python", action="store_true",
        help="Run the Python standalone driver (python -m clubb_python_driver.clubb_standalone)")

    run_group.add_argument("-jax", action="store_true",
        help=("Run through the JAX wrapper. An optional attached -jax=VALUE is "
              "forwarded unchanged; see clubb_jax/run_jax.py -launcher_help."))

    run_group.add_argument(
        "-gdb",
        action="store_true",
        help="Launch the selected compiled executable with gdb."
    )

    parser.add_argument("case_name", help="Name of the case to run")
    normalized_argv, jax_options, jax_occurrences = extract_jax_options(sys.argv[1:])
    args = parser.parse_args(normalized_argv)
    if jax_occurrences > 1:
        parser.error("-jax may be specified only once.")
    args.jax_options = jax_options

    ndefined = sum(bool(x) for x in [args.exe, args.driver_test, args.python, args.jax])
    if ndefined > 1:
        parser.error("Only one of -exe, -driver_test, -python, or -jax may be specified.")
    if args.gdb and (args.python or args.jax):
        parser.error("-gdb can only be used with compiled executables, not -python or -jax.")
    if args.zt_grid and args.zm_grid:
        sys.exit("\n\033[91mERROR: Cannot specify both a ZT grid and a ZM grid\033[0m")
    if args.nzmax and not (args.zt_grid or args.zm_grid):
        print("\n\033[93mWARNING: Specifying -nzmax will have no effect without specifying a -zm_grid or -zt_grid\033[0m")
    if args.batch_size is not None and args.multicol is None:
        parser.error("-batch_size requires -multicol.")
    if (args.stats_tstart is None) != (args.stats_tend is None):
        parser.error("-stats_tstart and -stats_tend must be provided together.")

    try:
        output_dir = str(resolve_output_dir(args.out_dir))
    except ValueError as exc:
        parser.error(str(exc))
    os.makedirs(output_dir, exist_ok=True)
    print(f"Output directory: {output_dir}")

    clubb_input_namelist = create_case_namelist(args, output_dir)
    run_cmd, run_cwd, run_env = choose_run_command(args)

    print(f"=================== Running {args.case_name} ===================")
    expect_stats_output = ((args.stats or "").strip().lower() != "none") and args.tout != 0
    result = run_case(run_cmd, run_cwd, args.case_name, clubb_input_namelist, output_dir, run_env,
                      expect_stats_output=expect_stats_output, passthrough_output=args.gdb)

    return result


if __name__ == "__main__":
    if "-python" in sys.argv[1:]:
        from utilities.setup_python_venv import ensure_python_venv

        ensure_python_venv()
    start_time = time.perf_counter()
    try:
        exit_code = main()
    finally:
        elapsed_time = time.perf_counter() - start_time
        print("-" * 50)
        print(f"run_scm.py total runtime: {elapsed_time:.2f} s")
        print("-" * 50)
    sys.exit(exit_code)
