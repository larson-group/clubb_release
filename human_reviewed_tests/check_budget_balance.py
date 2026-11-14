#!/usr/bin/env python3
"""Run SCM cases and check that their reported budgets balance.

Run from the repository root with a double-precision CLUBB build available:
    python human_reviewed_tests/check_budget_balance.py
    python human_reviewed_tests/check_budget_balance.py /path/to/clubb -workers 4
    python human_reviewed_tests/check_budget_balance.py -check_only -output_dir /path/to/stats

Usage:
    main(): run the fixed CASES suite, then check the output directory; return
        zero on success or one when budgets fail or no testable files exist.
    check_budgets(output_dir, iteration=0): check sorted *_stats.nc files,
        excluding SKIP_LIST cases. Zero selects all iterations; positive values
        select one iteration, numbered from one.
    checkNetcdfBudgets(fileName, iteration=0): check one NetCDF stats file and
        return whether all selected budgets and completeness checks pass.
    inferBudgetUnits(): use budget units or infer them from component metadata.
    calcTolerance(): choose the absolute tolerance from units and the term name.
    calcPercentError(): calculate the scaled discrepancy used by dispError().
    dispError(): report offending levels and update the failure count.
    run_cmd(): run the SCM command and stop if it fails.

Checks and limits:
    Cases run at a 300-second model/radiation timestep and output interval, with
    stats sampled every 300 seconds and l_lmm_stepping disabled. The case runner
    owns the worker default (half the available logical CPUs). By default this
    script checks the source tree's output/; -check_only checks existing output
    without running cases. -output_dir selects the run/check destination.

    Each *_bt tendency must equal the sum of components whose long_name contains
    "budget:". Terms are discovered by name and metadata, so they need not be
    contiguous in the file. Completeness also compares rtm/thlm tendencies with
    state changes between outputs, starting at the second sample. Time units
    must be seconds or minutes, with at least two output times.

    Checks retain the existing first-column and vertical-boundary exclusions.
    Tolerances depend on units, term-specific limits and per-level magnitude;
    budget closure uses a 60-second scale, completeness the output interval.
    calcTolerance() records the limits (velocity 0.02 m/s, temperature 0.01 K,
    mixing ratio 1e-8 kg/kg, number concentration 1e-5 per kg) and special cases.
    The existing discrepancy is (left-right)/(100*tolerance), with failure at
    magnitude >= TEST_LENIENCY (100). Non-finite checked values raise an error;
    MAX_FAILURES (1000) stops excessive diagnostics. A pass establishes these
    checks for the selected files, times, levels and first column.

Budget calculations originated in Cavyn VonDeylen's September 2010 checker
(Larson Group, UWM), updated to Python 3 by Tyler Cernik. NetCDF input uses
netCDF4 and NumPy.
"""

import argparse
import os
from pathlib import Path
import re
import subprocess
import sys
from time import strftime

TESTS_DIR = Path(__file__).resolve().parent
DEFAULT_CLUBB_ROOT = TESTS_DIR.parent

if __name__ == "__main__":
    sys.path.insert(0, str(DEFAULT_CLUBB_ROOT))
    from utilities.setup_python_venv import ensure_python_venv

    ensure_python_venv()

import numpy as np
from netCDF4 import Dataset


CASES = [
    "arm",
    "atex",
    "bomex",
    "dycoms2_rf01",
    "dycoms2_rf02_do",
    "dycoms2_rf02_ds",
    "dycoms2_rf02_nd",
    "dycoms2_rf02_so",
    "fire",
    "gabls2",
    "gabls3",
    "mpace_a",
    "mpace_b",
    "rico",
    "wangara",
]


# Set this to false to skip the completeness tests
COMPLETENESS_TEST = True

# Scale for calculating budget balance tolerance since we cannot easily access the
# model timestep (dt_main). Completeness tests use the frequency of statistical
# output which is obtained from the output files
TIME_SCALE_DENOMINATOR = 60 # Seconds

# How strict the test should be. Anything above this number will be considered
# a failure.
TEST_LENIENCY = 100 # Percent Error

# Any case name containing these strings will be skipped. Usually done due to time constraints
SKIP_LIST = ["nov11_altocu", \
             "twp_ice", \
             "cloud_feedback", \
             "arm_97", \
             "gabls3_night"]

# If more failures than this are reported, the test will exit to avoid a huge failure log.
MAX_FAILURES = 1000
numFails = 0 # Global variable indicating the current number of errors found.


def checkNetcdfBudgets(fileName, iteration=0):
    """
    Checks the balance of CLUBB budgets. This is done by assuming
    the budget term (henceforth: leftHandValue)
    is the sum of the budget term components(RightHandValue).
    For example, vm_bt = vm_ma + vm_gf + vm_cf + vm_ta + vm_f + vm_sdmp.

    Reads NetCDF stats through netCDF4 and compares the first column.

    Input: fileName: Path to a NetCDF stats file.
           iteration: The iteration to look at when balancing the budgets. 0 for all iterations
    """

    testSuccess = True

    # Open File
    try:
        ncFile = Dataset(fileName, 'r')
    except IOError:
        sys.stderr.write("\nCannot find file " + str(fileName) + "\n\n")
        sys.exit(1)

    # Find number of z levels
    try:
        if "zt" in ncFile.dimensions:
            numLevels = len(ncFile.dimensions["zt"])
        elif "zm" in ncFile.dimensions:
            numLevels = len(ncFile.dimensions["zm"])
        else:
            raise KeyError("Expected zt or zm dimension")
    except Exception:
        sys.stderr.write("Error parsing number of z levels\n")
        sys.exit(1)

    # Find the number of iterations
    try:
        numIterations = len(ncFile.dimensions['time'])
    except Exception:
        sys.stderr.write("Error parsing number of iterations\n")
        sys.exit(1)

    # Find the timestep
    try:
        if "seconds" in ncFile.variables['time'].units.lower():
            timestep = (ncFile.variables['time'][numIterations-1:numIterations] - \
                        ncFile.variables['time'][numIterations-2:numIterations-1])
        elif "minutes" in ncFile.variables['time'].units.lower():
           timestep = 60 * (ncFile.variables['time'][numIterations-1:numIterations] - \
                            ncFile.variables['time'][numIterations-2:numIterations-1])
        else:
            sys.stderr.write("Invalid units for timestep\n")
            sys.exit(1)
    except Exception:
        sys.stderr.write("Error parsing timestep\n")
        sys.exit(1)

    # Check that user entered a valid iteration
    if iteration > numIterations:
        sys.stderr.write(str(fileName) + " only indicates " + str(numIterations) + " iterations\n")
        sys.exit(1)

    # Discover the budgets once, then read only their data and state variables.
    varList = sorted(ncFile.variables)
    budgets = []
    required = set()
    for varName in varList:
        match = re.match(r"\w+_bt", varName)
        if not match:
            continue
        budgetName = match.group()
        budgetVar = ncFile.variables[varName]
        components = []
        for name in varList:
            prefix = re.match(r"\w+?_", name)
            if prefix and prefix.group() == budgetName[:-2] and name[-2:] != "bt":
                var = ncFile.variables[name]
                if "budget:" in getattr(var, "long_name", ""):
                    components.append(name)
        units = inferBudgetUnits(ncFile, varList, budgetName, budgetVar.units)
        state = budgetName[:-3] if budgetName[:-3] in ("rtm", "thlm") else None
        if state not in ncFile.variables:
            state = None
        budgets.append((varName, budgetName, units, "zt" in budgetVar.dimensions, components, state))
        required.add(varName)
        required.update(components)
        if state:
            required.add(state)

    profiles = {}
    for name in required:
        variable = ncFile.variables[name]
        if len(variable.shape) == 3:
            profiles[name] = np.array(variable[:, :, 0])
        elif len(variable.shape) == 2:
            profiles[name] = np.array(variable[:, :])
        else:
            raise ValueError("Expected 2D or 3D stats variable, got shape {}".format(variable.shape))
    ncFile.close()

    iterations = range(1, numIterations + 1) if iteration <= 0 else (iteration,)
    for t in iterations:
        for varName, budgetName, units, on_zt, components, state in budgets:
            left = profiles[varName][t-1]
            if t != 1 and COMPLETENESS_TEST and state:
                right = (profiles[state][t-1] - profiles[state][t-2]) / timestep
                tolerance = calcTolerance(units, timestep, state)
                testSuccess = dispError(left, right, left - right, tolerance, t, 0,
                                        state + " completeness test", units, testSuccess)

            right = None
            for name in components:
                if right is None:
                    right = np.zeros_like(left)
                right = right + profiles[name][t-1]
            difference = left - right
            tolerance = calcTolerance(units, TIME_SCALE_DENOMINATOR, budgetName[:-3])
            zLevel = 1 if on_zt else 0
            if on_zt:
                difference = difference[1:]
            testSuccess = dispError(left, right, difference, tolerance, t, zLevel,
                                    budgetName[:-3], units, testSuccess)
    return testSuccess

#--------------------------------------------------------------------------------------------------
def inferBudgetUnits(ncFile, varList, budgetVarName, budgetUnits):
    """
    Resolve units for budget variables.
    New stats output may leave *_bt units blank, but component terms usually
    retain units. In that case, infer units from the first matching component.
    """
    if budgetUnits is not None and str(budgetUnits).strip() != "":
        return budgetUnits

    prefix = budgetVarName[:-3]
    for variableName in varList:
        if variableName.startswith(prefix + "_") and not variableName.endswith("_bt"):
            componentUnits = getattr(ncFile.variables[variableName], "units", "")
            if componentUnits is not None and str(componentUnits).strip() != "":
                return componentUnits

    return budgetUnits

#--------------------------------------------------------------------------------------------------
def dispError(leftHandValue, rightHandValue, errorDifference, allowedTolerance, iteration, zLevel, termName, termUnits, testSuccess):
    """
    Find and display any budget failures to stdout
    Input: leftHandValue: list of floats representing the left side of equation being examined
           rightHandValue: list of floats representing the right side of equation being examined
           errorDifference: list of differences between two lists being compared.
           allowedTolerance: Tolerance the errorDifference may be off before considered a fail
           zLevel: Starting z level. For zm grid this should be 0, for zt it should be 1.
           termName: Name of variable that is being tested
           termUnits: Units of term being tested
           testSuccess: Whether the test is succeeding or failing
    """
    global numFails

    # Each level uses its own magnitude rather than inheriting earlier tolerances.
    start = max(zLevel, 1)
    left = leftHandValue[start:]
    right = rightHandValue[start:]
    if not np.all(np.isfinite(left)) or not np.all(np.isfinite(right)):
        raise ValueError("Non-finite value in {} at iteration {}".format(termName, iteration))
    magnitude = np.minimum(np.where(left == 0, 1e-40, np.abs(left)),
                        np.where(right == 0, 1e-40, np.abs(right)))
    with np.errstate(divide="ignore", under="ignore"):
        precision = np.power(10.0, np.floor(np.log10(magnitude)))
    tolerance = np.maximum(np.abs(allowedTolerance), np.abs(precision))
    percentError = calcPercentError(left, right, tolerance)
    for index in np.flatnonzero(np.abs(percentError) >= TEST_LENIENCY):
        level = start + index + 1
        value = errorDifference[level - zLevel - 1]
        testSuccess = False
        numFails += 1
        print(" ".join([termName, "fails at t=", str(iteration), "and z=", str(level),
                        "with a difference of", "%e" % value, termUnits, "and error",
                        "%.9f" % percentError[index], "%"]))
        if numFails >= MAX_FAILURES:
            print("Too many failures: exiting test. (Change MAX_FAILURES variable to view more)")
            sys.exit(1)

    return testSuccess

#--------------------------------------------------------------------------------------------------
def calcPercentError(experimental, accepted, denominator):
    """
    Finds the percent error using the given inputs
    Input: experimental: Your results in list form
           accepted: The correct values in list form
           denominator: Integer same as accepted unless accepted is extremely small (i.e. 0)
    """
    return (experimental - accepted) / (denominator * 100)

#--------------------------------------------------------------------------------------------------
def calcTolerance(termUnits, timestep, termName):
    """
    Finds the tolerance above which we consider the budgets do not balance. Which one is chosen
    is dependant on a terms units.
    Input: termUnits: String containting the units of the variable in question.
           timestep: Integer of real time in seconds between model iterations
           termName: Name of variable being looked at. Used for finding tolerance for special cases
    """

    # These were copied from constants_clubb.F90
    w_tol = 2e-2 # m/s
    thl_tol = 1e-2 # K
    rt_tol = 1e-8 # kg/kg
    Nr_tol = 1e-5 # num/kg

    if termUnits is None:
        termUnits = ""
    termUnits = termUnits.strip()

    # Special Cases for tolerance
    if termName == "Ncm":
        tol = 2e2
    elif termName == "Nrm":
        tol = 5e3
    # New stats files may omit units; fall back to known budget term names.
    elif termUnits == "":
        if termName in ("Nsm", "Nim", "Ngm"):
            tol = Nr_tol
        elif termName in ("rtm", "rrm", "rsm", "rim", "rgm"):
            tol = rt_tol
        elif termName == "rtp2":
            tol = rt_tol * rt_tol * 100
        elif termName == "rtp3":
            tol = rt_tol * rt_tol * rt_tol
        elif termName == "thlm":
            tol = thl_tol
        elif termName == "thlp2":
            tol = thl_tol * thl_tol
        elif termName == "thlp3":
            tol = thl_tol * thl_tol * thl_tol
        elif termName in ("um", "vm"):
            tol = w_tol
        elif termName in ("wp2", "up2", "vp2", "upwp", "vpwp"):
            tol = w_tol * w_tol
        elif termName == "wp3":
            tol = w_tol * w_tol * w_tol
        elif termName == "wprtp":
            tol = w_tol * rt_tol
        elif termName == "wpthlp":
            tol = thl_tol * w_tol
        elif termName == "rtpthlp":
            tol = thl_tol * rt_tol
        else:
            sys.stderr.write("Error parsing units: " + termUnits + "\nCheck this script's calcTolerance method")
            sys.exit(1)
    # Get error tolerance based on variable units
    elif termUnits == "(num/kg)/s":                                    # num/kg
        tol = Nr_tol
    elif termUnits == "(kg/kg)/s" or termUnits == "kg kg^{-1} s^{-1}": # kg/kg
        tol = rt_tol
    elif termUnits == "(kg^2)/(kg^2 s)": 			       # kg^2/kg^2
        tol = rt_tol * rt_tol * 100 # Multiply by 100 because otherwise it's too small for tests to pass
    elif termUnits == "m s^{-2}":                                      # m/s
        tol = w_tol
    elif termUnits == "m^2/s^3":                                       # m^2/s^2
        tol = w_tol * w_tol
    elif termUnits == "m^{3} s^{-4}":                                  # m^3/s^3
        tol = w_tol * w_tol * w_tol
    elif termUnits == "K s^{-1}":                                      # K
        tol = thl_tol
    elif termUnits == "(K^2)/s":                                       # K^2
        tol = thl_tol * thl_tol
    elif termUnits == "(m kg)/(s^2 kg)":                               # m/s kg/kg
        tol = w_tol * rt_tol
    elif termUnits == "(kg K)/(kg s)":                                 # K kg/kg
        tol = thl_tol * rt_tol
    elif termUnits == "(m K)/s^2":                                     # K m/s
        tol = thl_tol * w_tol
    elif termUnits == "kg^{3} kg^{-3} s^{-1}": 		               # kg^3/kg^3
        tol = rt_tol * rt_tol * rt_tol
    elif termUnits == "K^{3} s^{-1}":                                  # K^3
        tol = thl_tol * thl_tol * thl_tol
    else:
        sys.stderr.write("Error parsing units: " + termUnits + "\nCheck this script's calcTolerance method")
        sys.exit(1)

    return tol / timestep


def check_budgets(output_dir, iteration=0):
    """Check eligible stats files in output_dir; return zero only when all pass."""
    global numFails
    numFails = 0
    testable_files = sorted(
        path for path in Path(output_dir).iterdir()
        if path.is_file() and path.name.endswith("_stats.nc")
        and not any(case in path.name for case in SKIP_LIST)
    )
    test_success = bool(testable_files)
    if not testable_files:
        print("Unable to find testable data")
    for path in testable_files:
        print("\n" + strftime("%H:%M:%S") + " - Testing " + path.name)
        if not checkNetcdfBudgets(path, iteration):
            test_success = False
    print("Budgets successfully balance!" if test_success else "Budgets fail to balance")
    return 0 if test_success else 1


def run_cmd(cmd, *, cwd=None, env=None):
    """Run a case command, raising on failure before checking output."""
    print(f"+ {' '.join(str(part) for part in cmd)}")
    return subprocess.run(cmd, cwd=cwd, env=env, check=True)


def main() -> int:
    """Run the configured SCM cases unless checking existing output, then validate budgets."""
    parser = argparse.ArgumentParser(
        description="Run SCM cases and check budget closure, or check existing stats output.",
        add_help=False, allow_abbrev=False
    )
    parser.add_argument("-h", "-help", action="help", help="Show this help and exit.")
    parser.add_argument(
        "clubb_source",
        nargs="?",
        default=None,
        help="Path to CLUBB source root (default: derived from this script location).",
    )
    parser.add_argument(
        "-workers", metavar="N",
        help="Maximum concurrent cases (default: half the available logical CPUs).",
    )
    parser.add_argument(
        "-check_only", action="store_true",
        help="Check existing stats files without running SCM cases.",
    )
    parser.add_argument(
        "-output_dir", type=Path,
        help="Stats directory to write/check (default: <clubb_source>/output).",
    )
    parser.add_argument(
        "-iteration", type=int, default=0,
        help="Output iteration to check, numbered from one (default: 0 checks all).",
    )
    args = parser.parse_args()
    if args.iteration < 0:
        parser.error("-iteration must be zero or a positive integer")

    clubb_source = Path(args.clubb_source).resolve() if args.clubb_source else DEFAULT_CLUBB_ROOT.resolve()
    env = os.environ.copy()
    env["LAPACK"] = "/usr/lib64"

    run_scripts_dir = clubb_source / "run_scripts"
    run_scm_py = run_scripts_dir / "run_scm.py"
    output_dir = args.output_dir.resolve() if args.output_dir else clubb_source / "output"
    if args.check_only:
        return check_budgets(output_dir, args.iteration)
    if not run_scm_py.exists():
        raise SystemExit(f"Missing run_scm.py at: {run_scm_py}")

    run_cmd(
        [
            sys.executable,
            str(run_scm_py),
            "-output_dir", str(output_dir),
            *(["-workers", args.workers] if args.workers is not None else []),
            "-dt_main",
            "300",
            "-dt_rad",
            "300",
            "-tout",
            "300",
            # Keep l_lmm_stepping disabled and match sampling to output cadence.
            "-override",
            "l_lmm_stepping=.false.,stats_tsamp=300.0",
            "-cases",
            ",".join(CASES),
        ],
        cwd=str(run_scripts_dir),
        env=env,
    )

    return check_budgets(output_dir, args.iteration)


if __name__ == "__main__":
    sys.exit(main())
