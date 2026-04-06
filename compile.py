#!/usr/bin/env python3
import argparse
import os
import subprocess
import sys
import shutil
import glob

# Directory this script lives in, which is required to be clubb root
CLUBB_ROOT = os.path.dirname(os.path.abspath(__file__))

clubbstandards_script = os.path.join(CLUBB_ROOT, "utilities/CLUBBStandardsCheck.py") # perl version is utilities/CLUBBStandardsCheck.pl

clubb_driver_src      = os.path.join(CLUBB_ROOT, "src/clubb_driver.F90")
CLUBB_core_dir        = os.path.join(CLUBB_ROOT, "src/CLUBB_core/")
SILHS_dir             = os.path.join(CLUBB_ROOT, "src/SILHS/")
Benchmark_cases_dir   = os.path.join(CLUBB_ROOT, "src/Benchmark_cases/")
Radiation_dir         = os.path.join(CLUBB_ROOT, "src/Radiation/")
Microphys_dir         = os.path.join(CLUBB_ROOT, "src/Microphys/")
KK_microphys_dir      = os.path.join(CLUBB_ROOT, "src/Microphys/KK_microphys/")
G_unit_test_types_dir = os.path.join(CLUBB_ROOT, "src/G_unit_test_types/")

# Maps compiler/family names (from FC or LMOD_FAMILY_COMPILER) to canonical toolchain names.
# This lets multiple names (e.g. ifx, ifort, intel-oneapi) resolve to the same toolchain file
COMPILER_NAME_MAP = {
    "intel":          "intel",
    "intel-oneapi":   "intel",
    "intel-classic":  "intel",
    "ifx":            "intel",
    "ifort":          "intel",
    "nvhpc":          "nvhpc",
    "nvfortran":      "nvhpc",
    "gcc":            "gcc",
    "gnu":            "gcc",
    "gfortran":       "gcc",
    "ftn":            "cce", # this is uncomfortably generic, but that's what cce sets $FC to
    "cce":            "cce",
    "crayftn":        "cce",
}

def canonical_compiler(name):
    """Map a compiler or family name to its canonical toolchain name."""
    return COMPILER_NAME_MAP.get(name.lower(), name.lower())

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

def resolve_compiler_and_toolchain(args):
    """Identify the canonical compiler name and toolchain file to use.

    If -toolchain is given by the user, it is used directly and no file lookup occurs.
    Otherwise, LMOD_FAMILY_COMPILER is tried first, then FC, using the canonical name
    from COMPILER_NAME_MAP to locate a matching toolchain file.
    Exits with an error if no toolchain file can be found.

    Returns (compiler_name, toolchain_file).
    """
    kernel = os.uname().sysname.lower()
    arch   = os.uname().machine

    lmod_family = os.environ.get("LMOD_FAMILY_COMPILER")
    fc          = os.environ.get("FC")

    # Derive canonical compiler name for build/install directory naming
    if lmod_family:
        compiler = canonical_compiler(lmod_family)
    elif fc and shutil.which(fc):
        compiler = canonical_compiler(os.path.basename(shutil.which(fc)))
    else:
        compiler = "unknown"

    # User-specified toolchain: use it directly, skip file lookup
    if args.toolchain:
        return compiler, args.toolchain

    # Build a list of canonical names to try, in preference order
    candidates = []
    if lmod_family:
        candidates.append(canonical_compiler(lmod_family))
    if fc and shutil.which(fc):
        fc_canonical = canonical_compiler(os.path.basename(shutil.which(fc)))
        if fc_canonical not in candidates:
            candidates.append(fc_canonical)

    if not candidates:
        print(f"ERROR: No compiler detected (FC or LMOD_FAMILY_COMPILER) and no -toolchain specified")
        sys.exit(1)

    for name in candidates:
        toolchain_file = os.path.join(CLUBB_ROOT, f"cmake/toolchains/{kernel}_{arch}_{name}.cmake")
        if os.path.isfile(toolchain_file):
            print(f"Using inferred toolchain file: {toolchain_file}")
            return compiler, toolchain_file

    print(f"ERROR: No toolchain file found for detected compilers/kernel/arch combination:")
    for name in candidates:
        print(f"  cmake/toolchains/{kernel}_{arch}_{name}.cmake")
    print(f"  Use -toolchain to specify one explicitly.")
    sys.exit(1)


def configure_cmake(args, toolchain_file, inst_dir, build_type):

    """Run the cmake configure step."""

    print(f"about to cmnake {os.getcwd()}")

    # Configure CMake
    cmake_cmd = [
        "cmake",
        f"-S{CLUBB_ROOT}",
        f"-DUSE_NetCDF={to_on_off(not args.disable_netcdf)}",  # default ON
        f"-DSILHS={to_on_off(not args.disable_silhs)}",        # default ON
        f"-DENABLE_F2PY={to_on_off(args.python)}",             # default OFF
        f"-DENABLE_OMP={to_on_off(args.openmp)}",              # default OFF
        f"-DTUNING={to_on_off(args.tuning)}",                  # default OFF
        f"-DUSE_GPTL={to_on_off(args.gptl)}",                  # default OFF
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

    nproc = os.cpu_count() or 4

    cmd = ["cmake", "--build", ".", "--target", "install", f"-j{nproc}"]

    with open(logfile, "a") as log:
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True
        )

        # Stream output live to console and logfile
        for line in process.stdout:
            print(line, end="")
            log.write(line)

    retcode = process.wait()

    if retcode != 0:
        print(f"\n\033[91mBuild failed. See {logfile} for details.\033[0m")
        sys.exit(retcode)

def run_clubb_standards_check(logfile):

    retcode = 0

    # Run CLUBBStandardsCheck on all the source code we create.
    # This is important to keep Vince happy
    for path in [
        clubb_driver_src,
        CLUBB_core_dir,
        SILHS_dir,
        Benchmark_cases_dir,
        Radiation_dir,
        Microphys_dir,
        KK_microphys_dir,
        G_unit_test_types_dir
    ]:
        if os.path.isdir(path):
            files = glob.glob(f"{path}/*.F90")
        else:
            files = [path] if os.path.isfile(path) else []

        if not files:
            print(f"No matches for {path}")
            continue

        # Just sum all error codes together
        retcode += run_and_log(
            ["python", clubbstandards_script] + files,
            logfile
        )

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
    parser.add_argument("-python", action="store_true", help="Enable F2PY Python extension build")

    # Feature toggles
    parser.add_argument("-disable_netcdf", action="store_true", help="Disable NetCDF output support (default: enabled)")
    parser.add_argument("-disable_silhs", action="store_true", help="Disable SILHS (default: enabled)")
    parser.add_argument("-openmp", action="store_true", help="Enable OpenMP threading (default: disabled)")
    parser.add_argument("-tuning", action="store_true", help="Enable TUNING mode with extra runtime checks (default: disabled)")
    parser.add_argument("-gptl", action="store_true", help="Enable GPTL timing library (default: disabled)")

    # Passthrough
    parser.add_argument("extra_args", nargs=argparse.REMAINDER, help="Extra arguments passed to CMake")

    args = parser.parse_args()

    if args.python and args.disable_netcdf:
        print("ERROR: -python requires NetCDF. Remove -disable_netcdf.")
        sys.exit(1)

    # Our CMake files distinguish between "Debug" and "Release" for CMAKE_BUILD_TYPE
    build_type = "Debug" if args.debug else "Release"

    # Determine compiler and toolchain file to use based on environment and arguments
    compiler, toolchain_file = resolve_compiler_and_toolchain(args)

    subdir_suffix =  ""
    subdir_suffix += f"_DEBUG" if args.debug else ""
    subdir_suffix += f"_GPU{args.gpu}" if args.gpu != "none" else ""
    subdir_suffix += f"_PREC{args.precision}"
    subdir_suffix += "_PYTHON" if args.python else ""

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
    configure_cmake(args, toolchain_file, inst_dir, build_type)

    # Compile with cmake command
    run_cmake_build(build_log)

    # If build finished, set symlink "latest" to new install directory
    link_path = os.path.join(CLUBB_ROOT, "install/latest") 
    if os.path.lexists(link_path): os.remove(link_path)
    os.symlink(inst_dir, link_path)

    source_check_failed = False

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
