#!/usr/bin/env python3
import argparse
import os
import subprocess
import sys
import shutil
import glob

# Directory this script lives in, which is required to be clubb root
CLUBB_ROOT = os.path.dirname(os.path.abspath(__file__))

threadprivate_script  = os.path.join(CLUBB_ROOT, "utilities/check_for_missing_threadprivate.py")
clubbstandards_script = os.path.join(CLUBB_ROOT, "utilities/CLUBBStandardsCheck.py") # perl version is utilities/CLUBBStandardsCheck.pl

clubb_src_dir         = os.path.join(CLUBB_ROOT, "src/")
CLUBB_core_dir        = os.path.join(CLUBB_ROOT, "src/CLUBB_core/")
SILHS_dir             = os.path.join(CLUBB_ROOT, "src/SILHS/")
Benchmark_cases_dir   = os.path.join(CLUBB_ROOT, "src/Benchmark_cases/")
KK_microphys_dir      = os.path.join(CLUBB_ROOT, "src/KK_microphys/")
G_unit_test_types_dir = os.path.join(CLUBB_ROOT, "src/G_unit_test_types/")

def run_and_log(cmd, logfile, term_out=True):
    """Run a subprocess command, log output, optionally print to terminal, return exit code."""

    with open(logfile, "a") as f:
        f.write("\n================= Running: " + " ".join(cmd) + " =================\n")
        f.flush()

        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True
        )

        for line in process.stdout:
            if term_out:          # only print if enabled
                print(line, end="")
            f.write(line)         # always log
            f.flush()

        process.wait()
        f.write("\n")
        f.flush()

        return process.returncode

def to_on_off(flag: bool) -> str:
    """Convert Python bool to CMake ON/OFF string."""
    return "ON" if flag else "OFF"

def configure_cmake(args, compiler, inst_dir, build_type):

    """Set up toolchain file and run cmake configure step."""

    # Set toolchain file
    if args.toolchain:
        toolchain_file = args.toolchain
    else:
        kernel = os.uname().sysname.lower()
        arch = os.uname().machine
        toolchain_file = os.path.join(CLUBB_ROOT, f"cmake/toolchains/{kernel}_{arch}_{compiler}.cmake") 
        print(f"Using inferred toolchain file: {toolchain_file}")

    print(f"about to cmnake {os.getcwd()}")

    # Configure CMake
    cmake_cmd = [
        "cmake",
        f"-S{CLUBB_ROOT}",
        f"-DUSE_NetCDF={to_on_off(not args.disable_netcdf)}",  # default ON
        f"-DSILHS={to_on_off(not args.disable_silhs)}",        # default ON
        f"-DENABLE_OMP={to_on_off(args.openmp)}",              # default OFF
        f"-DTUNING={to_on_off(args.tuning)}",                  # default OFF
        f"-DUSE_GPTL={to_on_off(args.gptl)}",                      # default OFF
        f"-DCMAKE_TOOLCHAIN_FILE={toolchain_file}",
        f"-DCMAKE_INSTALL_PREFIX={inst_dir}",
        f"-DPRECISION={args.precision}",
        f"-DCMAKE_BUILD_TYPE={build_type}",
        f"-DGPU={args.gpu}",
        f"-DPython_EXECUTABLE={shutil.which('python')}",
    ]

    if shutil.which("ninja"):
        cmake_cmd += ["-G", "Ninja"]

    cmake_cmd += args.extra_args

    print("Running CMake configure...")

    try:
        subprocess.run(cmake_cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"\n\033[91mERROR: cmake configure failed.\033[0m")
        sys.exit(e.returncode)
        

def run_cmake_build(logfile):

    nproc = os.cpu_count() or 4  # fallback to 4 if detection fails

    # This script command allows us to log output to "logfile" but also get nice looking
    # output from cmake
    cmd = f'script -qfec "cmake --build . --target install" -- -j{nproc} /dev/stdout | tee -a {logfile}'

    # this one might work better, need to figure out
    #cmd = f'script -q -c "cmake --build . --target install -j{nproc}" /dev/stdout | tee -a {logfile}' 

    # Run the command
    subprocess.call(cmd, shell=True)

    # Parse the logfile, and save the reported COMMAND_EXIT_CODE into retcode
    retcode = 0
    with open(logfile, "r") as f:
        for line in f:
            if "COMMAND_EXIT_CODE=" in line:
                retcode = int(line.split("=")[1].strip('"]\n'))

    if retcode != 0:
        print(f"\n\033[91mBuild failed. See {logfile} for details.\033[0m")
        sys.exit(retcode)
    

def run_threadprivate_check(logfile):

    # Run check_for_missing_threadprivate check on CLUBB_core and SILHS. 
    # This is important because these are the libraries used within host models
    retcode = run_and_log(
        ["python", threadprivate_script, CLUBB_core_dir, SILHS_dir],
        logfile
    )

    return retcode


def run_clubb_standards_check(logfile):

    retcode = 0

    # Run CLUBBStandardsCheck on all the source code we create.
    # This is important to keep Vince happy
    for directory in [
        clubb_src_dir,
        CLUBB_core_dir,
        SILHS_dir,
        Benchmark_cases_dir,
        KK_microphys_dir,
        G_unit_test_types_dir
    ]:
        files = glob.glob(f"{directory}/*.F90")

        if not files:
            print(f"No matches for {directory}")
            continue

        retcode = run_and_log(
            ["python", clubbstandards_script] + files,
            logfile
        )

        # Save an error
        if retcode != 0: 
          retcode = 1

    return retcode


def run_ctests(logfile):
    """
    Run ctest in the current working directory (assumed to be build dir).
    """

    print(f"\nRunning CMake ctests")
    retcode = run_and_log( ["ctest"], logfile, term_out=True)

    if retcode != 0:
        print("\033[91mSome tests failed. See ctest output above.\033[0m")
    else:
        print("\033[92mAll tests passed.\033[0m")

    return retcode


def main():

    parser = argparse.ArgumentParser(description="Compile CLUBB with CMake")
    
    # Core build options
    parser.add_argument("-install", metavar="DIR", help="Install directory for CLUBB")
    parser.add_argument("-toolchain", metavar="FILE", help="Path to CMake toolchain file")
    parser.add_argument("-gpu", choices=["none", "openacc", "openmp"], default="none", help="GPU option for build")
    parser.add_argument("-precision", choices=["single", "double", "quad"], default="double", help="Floating-point precision")
    parser.add_argument("-debug", action="store_true", help="Compile in debug mode")
    parser.add_argument("-run_tests", action="store_true", help="Run ctests after compilation")

    # Feature toggles
    parser.add_argument("-disable_netcdf", action="store_true", help="Disable NetCDF output support (default: enabled)")
    parser.add_argument("-disable_silhs", action="store_true", help="Disable SILHS (default: enabled)")
    parser.add_argument("-openmp", action="store_true", help="Enable OpenMP threading (default: disabled)")
    parser.add_argument("-tuning", action="store_true", help="Enable TUNING mode with extra runtime checks (default: disabled)")
    parser.add_argument("-gptl", action="store_true", help="Enable GPTL timing library (default: disabled)")

    # Passthrough
    parser.add_argument("extra_args", nargs=argparse.REMAINDER, help="Extra arguments passed to CMake")

    args = parser.parse_args()

    # Our CMake files distinguish between "Debug" and "Release" for CMAKE_BUILD_TYPE
    build_type = "Debug" if args.debug else "Release"


    # Determine compiler using "FC" (fortran compiler) environment variable.
    # This is set on larson-group computers through the use of lmod (module system).
    # FC is not required for cmake
    if "FC" in os.environ:
      fc = os.environ.get("FC")
      compiler = os.path.basename(shutil.which(fc))
    else:

      if args.toolchain:
        print(f"WARNING: No Fortran compiler (FC) set in user environment. Setting compiler name to 'unknown'")
        print(f"         Relying on specified toolchain for cmake {args.toolchain}")
      else:
        # If we can't find a compiler, and no toolchain is specified, this is an error
        print(f"ERROR: No Fortran compiler (FC) detected and no entry specified for -toolchain")
        sys.exit(1)

    subdir_suffix =  ""
    subdir_suffix += f"_DEBUG" if args.debug else ""
    subdir_suffix += f"_GPU{args.gpu}" if args.gpu != "none" else ""
    subdir_suffix += f"_PREC{args.precision}"

    # Create build directory and cd into it
    build_dir = os.path.join(CLUBB_ROOT, f"build/{compiler}{subdir_suffix}") 
    os.makedirs(build_dir, exist_ok=True)
    os.chdir(build_dir)
    print(f"Setting CLUBB installation dir: {build_dir}")

    # Reset build log
    build_log = os.path.join(build_dir, f"cmake_build_output.txt") 
    open(build_log, "w").close()

    inst_dir = args.install if args.install else os.path.join(CLUBB_ROOT, f"install/{compiler}{subdir_suffix}") 
    print(f"Setting CLUBB installation dir: {inst_dir}")

    # Run configure step and save installation directory
    configure_cmake(args, compiler, inst_dir, build_type)

    # Compile with cmake command
    run_cmake_build(build_log)

    # If build finished, set symlink "latest" to new install directory
    link_path = os.path.join(CLUBB_ROOT, "install/latest") 
    if os.path.lexists(link_path): os.remove(link_path)
    os.symlink(inst_dir, link_path)

    # Run the CLUBB threadprivate check, which ensures that all module variables are 
    # declared "threadprivate" 
    if run_threadprivate_check(build_log) != 0:
        print(f"\033[91m===============================================================")
        print(f"\033[91mTHREADPRIVATE CHECK FAILED")
        print(f"\033[91m  THIS IS PRINTED IN ALL RED, CAPITAL LETTERS, AND USES")
        print(f"\033[91m  AN EXCLAMATION MARK TO ENSURE THE DEVELOPERS FEEL SHAME!")
        print(f"\033[91m  IF YOU ARE ONE OF THESE \"DEVELOPERS\" CHECK THE")
        print(f"\033[91m  LOG FILE FOR DETAILS: {build_log}")
        print(f"\033[91m===============================================================\033[0m")
        source_check_failed = True

    # Run the CLUBB standards check, which looks for various code issues
    if run_clubb_standards_check(build_log) != 0:
        print(f"\033[91m===============================================================")
        print(f"\033[91mCLUBBStandardsCheck FAILED")
        print(f"\033[91m  THIS IS PRINTED IN ALL RED, CAPITAL LETTERS, AND USES")
        print(f"\033[91m  AN EXCLAMATION MARK TO ENSURE THE DEVELOPERS FEEL SHAME!")
        print(f"\033[91m  IF YOU ARE ONE OF THESE \"DEVELOPERS\" CHECK THE")
        print(f"\033[91m  LOG FILE FOR DETAILS: {build_log}")
        print(f"\033[91m===============================================================\033[0m")
        source_check_failed = True

    if source_check_failed:
        # Build passed, but one of the checks failed
        print("\n\033[93mBuild completed successfully, but some source code checks have failed.\033[0m")
    else:
        # Successful build and no failed checks 
        print("\n\033[92mBuild completed successfully, and all source code checks passed.\033[0m")

    if args.run_tests:
        # If we run the ctests, consider this the pass/fail criteria
        return run_ctests(build_log)
    else:
        # Otherwise just return 0, since getting here means the build has completed
        return 0


if __name__ == "__main__":
    sys.exit(main())
