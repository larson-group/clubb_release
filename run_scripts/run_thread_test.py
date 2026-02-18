#!/usr/bin/env python3
import os
import sys
import shutil
import subprocess
import argparse
from pathlib import Path

# -----------------------------------------------------------------------------
# Description:
#   Python script to run the test whether CLUBB is still threadsafe.
#   Ported from the original Bash script.
# -----------------------------------------------------------------------------

THREAD_COUNT = "8"

RUN_CASES = ["fire", "bomex", "rico_silhs", "gabls3"]
NAMELISTS = ["clubb_1.in", "clubb_2.in", "clubb_3.in", "clubb_4.in"]

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
INPUT_DIR = REPO_ROOT / "input"
OUTPUT_DIR = REPO_ROOT / "output"

FLAGS_FILE = INPUT_DIR / "tunable_parameters" / "configurable_model_flags.in"
PARAMS_FILE = INPUT_DIR / "tunable_parameters" / "tunable_parameters.in"
SILHS_PARAMS_FILE = INPUT_DIR / "tunable_parameters" / "silhs_parameters.in"
STATS_FILE = INPUT_DIR / "stats" / "standard_stats.in"
SERIAL = OUTPUT_DIR / "serial"
PARALLEL = OUTPUT_DIR / "parallel"
AGGREGATE_FILES = [PARAMS_FILE, SILHS_PARAMS_FILE, FLAGS_FILE]

EXEC_NAME_STANDALONE = "clubb_standalone"
EXEC_NAME_THREAD_TEST = "clubb_thread_test"
INSTALL_LATEST_DIR = REPO_ROOT / "install" / "latest"
BINDIFF_SCRIPT = SCRIPT_DIR / "run_bindiff_all.py"


def strip_comments(filename: str):
    """Remove everything after '!' on each line (Fortran namelist comment)."""
    path = Path(filename)
    lines = path.read_text().splitlines()
    path.write_text("\n".join(line.split("!")[0].rstrip() for line in lines) + "\n")


def find_executable_pair() -> tuple[Path, Path]:
    """
    Find a matching pair of executables for standalone and thread tests.
    Use install/latest only.
    """
    latest_dir = INSTALL_LATEST_DIR.resolve()
    standalone = latest_dir / EXEC_NAME_STANDALONE
    thread_test = latest_dir / EXEC_NAME_THREAD_TEST
    if standalone.exists() and thread_test.exists():
        return standalone, thread_test
    raise FileNotFoundError(
        f"Could not find both executables in {latest_dir}"
    )


def copy_case_outputs(destination: Path):
    destination.mkdir(parents=True, exist_ok=True)
    for case in RUN_CASES:
        for src in OUTPUT_DIR.glob(f"{case}*"):
            dst = destination / src.name
            if dst.exists():
                if dst.is_dir():
                    shutil.rmtree(dst)
                else:
                    dst.unlink()
            if src.is_dir():
                shutil.copytree(src, dst)
            else:
                shutil.copy2(src, dst)


def clear_old_outputs():
    for directory in [SERIAL, PARALLEL]:
        if directory.exists():
            shutil.rmtree(directory)
    for case in RUN_CASES:
        for file_path in OUTPUT_DIR.glob(f"{case}*"):
            if file_path.is_dir():
                shutil.rmtree(file_path)
            else:
                file_path.unlink()


def write_aggregate_namelist(output_path: Path, input_files: list[Path]):
    """
    Concatenate namelist fragments while forcing a newline between files.
    This avoids accidental token merges when an input file lacks trailing newline.
    """
    with open(output_path, "w") as out:
        for src_path in input_files:
            text = src_path.read_text()
            out.write(text)
            if not text.endswith("\n"):
                out.write("\n")


def run_serial(clubb_standalone: Path, run_env):
    print("Running CLUBB in serial... ", end="", flush=True)
    l_failed = False
    for case in RUN_CASES:
        model_file = INPUT_DIR / "case_setups" / f"{case}_model.in"
        write_aggregate_namelist(
            SCRIPT_DIR / "clubb.in",
            [*AGGREGATE_FILES, model_file, STATS_FILE],
        )

        strip_comments(SCRIPT_DIR / "clubb.in")
        result = subprocess.run([str(clubb_standalone)], cwd=SCRIPT_DIR, env=run_env)
        if result.returncode not in (0, 6):
            print(f"\nSerial run failed for case {case} with exit code {result.returncode}")
            l_failed = True

        (SCRIPT_DIR / "clubb.in").unlink()
    print("Done!")
    return not l_failed


def run_parallel(clubb_thread_test: Path, run_env):
    print("Running CLUBB in parallel... ", end="", flush=True)
    # Some code paths append to case setup files with status='old'.
    # Seed them so parallel startup is robust even if setup-file initialization
    # is skipped in a thread.
    for case in RUN_CASES:
        (OUTPUT_DIR / f"{case}_setup.txt").touch(exist_ok=True)

    for case, namelist in zip(RUN_CASES, NAMELISTS):
        model_file = INPUT_DIR / "case_setups" / f"{case}_model.in"
        namelist_path = SCRIPT_DIR / namelist
        write_aggregate_namelist(
            namelist_path,
            [*AGGREGATE_FILES, model_file, STATS_FILE],
        )

        print(f"Configuring file: {model_file.name}")
        strip_comments(namelist_path)

    print("Checking existence of clubb_thread_test")
    result = subprocess.run([str(clubb_thread_test)], cwd=SCRIPT_DIR, env=run_env)
    l_failed = result.returncode not in (0, 6)
    if l_failed:
        print(f"\nParallel thread-test failed with exit code {result.returncode}")
    print("clubb_thread_test executable exists")

    for namelist in NAMELISTS:
        print(f"Removing file(s): {namelist}")
        (SCRIPT_DIR / namelist).unlink()
    print("Done!")
    return not l_failed


def parse_args():
    parser = argparse.ArgumentParser(description="Run CLUBB thread-safety regression test.")
    parser.add_argument(
        "--threads",
        type=int,
        default=int(THREAD_COUNT),
        help="OMP thread count to use (default: 8).",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    run_env = os.environ.copy()
    run_env["OMP_NUM_THREADS"] = str(args.threads)
    # Helps on systems where OpenMP shared-memory segments are restricted.
    run_env["KMP_USE_SHM"] = run_env.get("KMP_USE_SHM", "FALSE")

    try:
        clubb_standalone, clubb_thread_test = find_executable_pair()
    except FileNotFoundError as err:
        print(err)
        raise SystemExit(1)

    print(f"Using standalone executable: {clubb_standalone}")
    print(f"Using thread-test executable: {clubb_thread_test}")
    print(f"OMP_NUM_THREADS={run_env['OMP_NUM_THREADS']}")

    clear_old_outputs()

    serial_ok = run_serial(clubb_standalone, run_env)
    copy_case_outputs(SERIAL)

    parallel_ok = run_parallel(clubb_thread_test, run_env)
    copy_case_outputs(PARALLEL)

    diff_path = SCRIPT_DIR / "diff.txt"
    print("Diffing the output... ", end="", flush=True)
    # Use numeric NetCDF comparison via bindiff to avoid false failures from
    # byte-level metadata/layout differences in netCDF files.
    bindiff_cmd = [
        sys.executable,
        str(BINDIFF_SCRIPT),
        "-v",
        "0",
        str(SERIAL),
        str(PARALLEL),
    ]
    with open(diff_path, "w") as diff_file:
        bindiff_result = subprocess.run(
            bindiff_cmd,
            cwd=SCRIPT_DIR,
            stdout=diff_file,
            stderr=subprocess.STDOUT,
            check=False,
        )
    print("Done!")

    result = 0
    if not serial_ok:
        print("One or more serial executions returned non-zero status.")
        result = 1
    if not parallel_ok:
        print("Parallel thread-test execution returned non-zero status.")
        result = 1
    if bindiff_result.returncode != 0:
        print("run_bindiff_all.py reported numerical differences.")
        print(diff_path.read_text())
        result = 1
    elif diff_path.stat().st_size > 0:
        # Keep summary output from run_bindiff_all.py for visibility.
        print(diff_path.read_text())
    else:
        print("No differences found")

    diff_path.unlink(missing_ok=True)
    shutil.rmtree(SERIAL, ignore_errors=True)
    shutil.rmtree(PARALLEL, ignore_errors=True)

    raise SystemExit(result)


if __name__ == "__main__":
    main()
