#!/usr/bin/env python3
import os
import shutil
import subprocess
import fileinput
from pathlib import Path

# -----------------------------------------------------------------------------
# Description:
#   Python script to run the test whether CLUBB is still threadsafe.
#   Ported from the original Bash script.
# -----------------------------------------------------------------------------

# Set the number of threads
os.environ["OMP_NUM_THREADS"] = "8"

RUN_CASES = ["fire", "bomex", "rico_silhs", "gabls3"]
NAMELISTS = ["clubb_1.in", "clubb_2.in", "clubb_3.in", "clubb_4.in"]

FLAGS_FILE = "../input/tunable_parameters/configurable_model_flags.in"
PARAMS_FILE = "../input/tunable_parameters/tunable_parameters.in"
SILHS_PARAMS_FILE = "../input/tunable_parameters/silhs_parameters.in"
STATS_FILE = "../input/stats/standard_stats.in"
SERIAL = Path("../output_serial")
PARALLEL = Path("../output_parallel")

# these "latest" paths will use the most recent cmake compiled version
clubb_thread_test = "../install/latest/clubb_thread_test"
clubb_standalone = "../install/latest/clubb_standalone"


def strip_comments(filename: str):
    """Remove everything after '!' on each line (Fortran namelist comment)."""
    for line in fileinput.input(filename, inplace=True):
        print(line.split("!")[0].rstrip())


def run_serial():
    print("Running CLUBB in serial... ", end="", flush=True)
    for case in RUN_CASES:
        model_file = f"../input/case_setups/{case}_model.in"
        with open("clubb.in", "w") as out:
            for f in [model_file, FLAGS_FILE, PARAMS_FILE, SILHS_PARAMS_FILE, STATS_FILE]:
                with open(f) as src:
                    out.write(src.read())

        strip_comments("clubb.in")

        if Path(clubb_standalone).exists():
            subprocess.run([clubb_standalone])
        else:
            print("clubb_standalone not found (did you re-compile?)")
            exit(1)

        os.remove("clubb.in")
    print("Done!")


def run_parallel():
    print("Running CLUBB in parallel... ", end="", flush=True)
    for case, namelist in zip(RUN_CASES, NAMELISTS):
        model_file = f"../input/case_setups/{case}_model.in"
        with open(namelist, "w") as out:
            for f in [model_file, FLAGS_FILE, PARAMS_FILE, SILHS_PARAMS_FILE, STATS_FILE]:
                with open(f) as src:
                    out.write(src.read())

        print(f"Configuring file: {model_file}")
        strip_comments(namelist)

    print("Checking existence of clubb_thread_test")
    if Path(clubb_thread_test).exists():
        subprocess.run([clubb_thread_test], check=True)
        print("clubb_thread_test executable exists")
    else:
        print("clubb_thread_test not found (did you re-compile?)")
        exit(1)

    for namelist in NAMELISTS:
        print(f"Removing file(s): {namelist}")
        os.remove(namelist)
    print("Done!")


def main():
    # Save working dir and change to script dir
    restore_dir = Path.cwd()
    os.chdir(Path(__file__).parent)

    run_serial()
    SERIAL.mkdir(parents=True, exist_ok=True)
    for f in Path("../output").glob("*.*"):
        shutil.move(str(f), SERIAL)

    run_parallel()
    PARALLEL.mkdir(parents=True, exist_ok=True)
    for f in Path("../output").glob("*.*"):
        shutil.move(str(f), PARALLEL)

    print("Diffing the output... ", end="", flush=True)
    with open("diff.txt", "w") as diff_file:
        for case in RUN_CASES:
            for suffix in ["_zt.nc", "_zm.nc", "_sfc.nc"]:
                f1 = SERIAL / f"{case}{suffix}"
                f2 = PARALLEL / f"{case}{suffix}"
                subprocess.run(["diff", str(f1), str(f2)], stdout=diff_file, stderr=subprocess.STDOUT)
    print("Done!")

    result = 0
    if os.path.getsize("diff.txt") > 0:
        print(Path("diff.txt").read_text())
        result = 1
    else:
        print("No differences found")

    setup_files = list(SERIAL.glob("*_setup.txt"))
    if len(setup_files) != len(RUN_CASES):
        print("One or more simulations failed to run in serial")
        result = 1

    setup_files = list(PARALLEL.glob("*_setup.txt"))
    if len(setup_files) != len(RUN_CASES):
        print("One or more simulations failed to run in parallel")
        result = 1

    os.remove("diff.txt")
    shutil.rmtree(SERIAL)
    shutil.rmtree(PARALLEL)

    os.chdir(restore_dir)
    exit(result)


if __name__ == "__main__":
    main()
