#!/usr/bin/env python3

######################################################################################
# Description:
#   This script tests the differences between the CPU and GPU results of
#   monotonic_turbulent_flux_limit, located in mono_flux_limiter.F90. The usual
#   CPU vs GPU diffness test isn't sufficient to test this function because only
#   the most unstable cases make use of monotonic_turbulent_flux_limit(), so it's
#   impossible to tell the differences between an error or accumulated bit changes
#   after enough timesteps to trigger the use of this function.
#
#   This test works by making use of NVIDIA's PCAST compare feature, where we use
#   the flag "-gpu=redundant" to compile the code so that it's run on both
#   CPU and GPU at the same time, then use "!$acc compare(var)" commands to test
#   the differences between the CPU and GPU results. The PCAST feature keeps the
#   CPU and GPU results in sync from timestep to timstep, so the differences
#   we measure here will not be mixed with accumulating bit changes.
#
#   The PCAST feature relies on an environment variable "PCAST_COMPARE", that for
#   this test we set to "abs=10", which means that only differences above "10^(-10)"
#   will be reported. This was chosen by some small experiments. Usually correct
#   changes are more in the range of 10^-20 or smaller.
#
# Note:
#   This test does NOT confirm that changes made in monotonic_turbulent_flux_limit
#   are correct. This only confirms that there are no errors present that only
#   affect GPU results. Bit changing changes should be tested separately.
#
# References:
#   https://developer.nvidia.com/blog/detecting-divergence-using-pcast-to-compare-gpu-to-cpu-results/
#   https://github.com/larson-group/cam/issues/175#issuecomment-2259291053
#
# Author:
#   Gunther Huebler
######################################################################################

import os
import signal
import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent

# File paths (cmake workflow version)
TOOLCHAIN_FILE = REPO_ROOT / "cmake" / "toolchains" / "linux_x86_64_nvfortran.cmake"
SOURCE_FILE = REPO_ROOT / "src" / "CLUBB_core" / "mono_flux_limiter.F90"
PARAMS_FILE = REPO_ROOT / "input" / "tunable_parameters" / "tunable_parameters.in"
MULTICOL_PARAMS = SCRIPT_DIR / "clubb_params_multi_col.in"

# Script paths
COMPILE_SCRIPT = REPO_ROOT / "compile.py"
RUN_SCM_SCRIPT = SCRIPT_DIR / "run_scm.py"
MULTICOL_SCRIPT = SCRIPT_DIR / "create_multi_col_params.py"

# Case setup
RUN_CASE = "mc3e"      # unstable/noisy case is required
CASE_MAX_ITERS = "60"  # run exactly 60 timesteps

# Preserve original file contents for robust cleanup
ORIGINAL_TEXT: dict[Path, str | None] = {}

ERR_CODE = 1  # default to error


def save_original(path: Path) -> None:
    if path not in ORIGINAL_TEXT:
        if path.exists():
            ORIGINAL_TEXT[path] = path.read_text(encoding="utf-8")
        else:
            ORIGINAL_TEXT[path] = None


def restore_backups() -> None:
    # This should be called to restore the backups
    for path, text in ORIGINAL_TEXT.items():
        if text is None:
            if path.exists():
                path.unlink()
        else:
            path.write_text(text, encoding="utf-8")


# This is called when kill signals are given.
def signal_handler(sig, frame):
    print("\n======== Script killed. Restoring and exiting. ========\n")
    restore_backups()
    raise SystemExit(0)


# Setup restore_backups to be called if Ctrl+C is hit
signal.signal(signal.SIGINT, signal_handler)   # Handle Ctrl+C (SIGINT)
signal.signal(signal.SIGTERM, signal_handler)  # Handle termination (SIGTERM)


def run_command(cmd: list[str], cwd: Path | None = None, env: dict[str, str] | None = None) -> tuple[int, str]:
    process = subprocess.Popen(
        cmd,
        cwd=str(cwd) if cwd else None,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        errors="replace",
        bufsize=1,
    )
    output_lines: list[str] = []
    assert process.stdout is not None
    for line in process.stdout:
        print(line, end="", flush=True)
        output_lines.append(line)
    process.wait()
    return process.returncode, "".join(output_lines)


def patch_monoflux_source(text: str) -> str:
    # Remove the specific test comment
    return text.replace("! MONOFLUX TEST COMMENT DO NOT REMOVE ", "")


def patch_toolchain_for_pcast(text: str) -> str:
    # Enable Nvidia PCAST features with (-gpu=redundant)
    # Also include math_uniform and nofma to prevent CPU vs GPU divergence as much as possible
    #
    # SILHS uses a CUDA call that doesn't seem to work with -gpu=redundant,
    # so CUDA is removed temporarily for this test.
    replacements = [
        (
            "set(GPU_DEFINITIONS CLUBB_GPU CUDA)",
            "set(GPU_DEFINITIONS CLUBB_GPU)",
        ),
        (
            'set(GPU_COMPILE_FLAGS   "${GPU_COMPILE_FLAGS}" "-acc")',
            'set(GPU_COMPILE_FLAGS   "${GPU_COMPILE_FLAGS}" "-acc=gpu" "-gpu=redundant,math_uniform,nofma" "-nofma")',
        ),
        (
            'set(GPU_LINK_FLAGS      "-acc")',
            'set(GPU_LINK_FLAGS      "-acc=gpu" "-gpu=redundant,math_uniform,nofma" "-nofma")',
        ),
    ]

    updated = text
    for old, new in replacements:
        if old in updated:
            updated = updated.replace(old, new)
            continue
        if new in updated:
            continue
        raise RuntimeError(f"Failed to patch toolchain (missing expected snippet): {old}")

    return updated


try:
    ##################################################################
    #               Backup files to be changed
    ##################################################################
    save_original(SOURCE_FILE)
    save_original(TOOLCHAIN_FILE)
    save_original(MULTICOL_PARAMS)

    ##################################################################
    #         Modify the source file (mono_flux_limiter.F90)
    #         and toolchain file for this test
    ##################################################################

    SOURCE_FILE.write_text(
        patch_monoflux_source(ORIGINAL_TEXT[SOURCE_FILE] or ""),
        encoding="utf-8",
    )

    TOOLCHAIN_FILE.write_text(
        patch_toolchain_for_pcast(ORIGINAL_TEXT[TOOLCHAIN_FILE] or ""),
        encoding="utf-8",
    )

    ##################################################################
    #                           Compile
    ##################################################################

    # Compile with cmake/compile.py in OpenACC mode.
    # Use debug mode so optimization is low (similar intent to old -O0 path).
    compile_cmd = [sys.executable, str(COMPILE_SCRIPT), "-gpu", "openacc", "-debug"]
    rc, _ = run_command(compile_cmd, cwd=REPO_ROOT)
    if rc != 0:
        raise RuntimeError("Compilation failed.")

    ##################################################################
    #                           Run Case
    ##################################################################

    # Create a multicol param
    create_multicol_cmd = [
        sys.executable,
        str(MULTICOL_SCRIPT),
        "-n",
        "16",                    # use 16 columns
        "-mode",
        "dup_tweak",             # duplicate and tweak initial param values
        "-param_file",
        str(PARAMS_FILE),        # location of param file
        "-out_file",
        str(MULTICOL_PARAMS),    # define output file name
    ]
    rc, _ = run_command(create_multicol_cmd, cwd=SCRIPT_DIR)
    if rc != 0:
        raise RuntimeError("Failed to create multi-column parameter file.")

    # Setting "PCAST_COMPARE=abs=n" causes the "acc compare" clauses used in this test
    # to only report if the GPU results differ from CPU results by more than 10^-n
    # see https://developer.nvidia.com/blog/detecting-divergence-using-pcast-to-compare-gpu-to-cpu-results/
    run_env = os.environ.copy()
    run_env["PCAST_COMPARE"] = "abs=10"

    ##################################################################
    #                       Edit Input Files
    ##################################################################
    # Old script edited *_model.in and configurable_model_flags.in directly.
    # New script passes equivalent settings via run_scm.py options to avoid file edits:
    #   - max_iters to control run length
    #   - override to set penta/tridiag/l_lh_straight_mc
    #   - stats none because we only care about the output/log text

    ##################################################################
    #                           Run Case
    ##################################################################

    # Run the case with stats disabled, we only care about the output text.
    run_cmd = [
        sys.executable,
        str(RUN_SCM_SCRIPT),
        "-exe",
        str(REPO_ROOT / "install" / "latest" / "clubb_standalone"),
        "-debug",
        "0",
        "-max_iters",
        CASE_MAX_ITERS,
        "-params",
        str(MULTICOL_PARAMS),
        "-stats",
        "none",
        "-override",
        "penta_solve_method=2,tridiag_solve_method=2,l_lh_straight_mc=true",
        RUN_CASE,
    ]
    rc, run_output = run_command(run_cmd, cwd=SCRIPT_DIR, env=run_env)
    if rc != 0:
        raise RuntimeError("Case run failed.")

    ##################################################################
    #                   Check Output and Report Result
    ##################################################################

    # Prefer the case log from run_scm.py if it exists, otherwise use captured output.
    run_log = REPO_ROOT / "output" / f"{RUN_CASE}_log"
    if run_log.exists():
        output_text = run_log.read_text(encoding="utf-8", errors="replace")
    else:
        output_text = run_output

    # Check for relevant lines in the output
    monoflux_wpxp_found = "MONOFLUX: wpxp adjusted" in output_text
    monoflux_xm_found = "MONOFLUX: xm adjusted" in output_text
    fail_abs_found = "FAIL ABS" in output_text

    print("\n==================================== RESULT ====================================")

    # If "FAIL ABS" appears, then we fail no matter what
    # If there are no "adjusting" messages then monoflux_limiter was never really tested, we consider this a fail
    # If "adjusting" messages are found and "FAIL ABS" isnt' then we pass
    if fail_abs_found:
        print("\n\tTEST FAILED: 'FAIL ABS' found in the output.")
        print("\tGPU results differ too significantly from CPU results.\n")
    elif not (monoflux_wpxp_found and monoflux_xm_found):
        print("\n\tTEST FAILED: Neither 'MONOFLUX: wpxp adjusted' nor 'MONOFLUX: xm adjusted'")
        print("\twere found in the output. This means the flux limiter wasn't tested.\n")
    else:
        ERR_CODE = 0
        print("\nTEST PASSED: mono_flux_limiter did modify fields, and CPU results match GPU results.\n")

except Exception as e:
    print("\n==================================== Test Incomplete ====================================")
    print(f"ERROR: {e}")

finally:
    restore_backups()

raise SystemExit(ERR_CODE)
