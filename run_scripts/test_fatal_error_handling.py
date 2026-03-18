#!/usr/bin/env python3

###############################################################################
# $Id$
#
#   Description:
#       This script injects errors into clubb then runs a case and captures the
#       error output. The injected errors need to be labeled in the fortran code
#       with a name, then added to the TESTS list in this script.
#
#   Example:
#
#       In fortran code:
#                           ! advance_wp2_wp3_bad_wp2
#
#       In this script:
#                           TestCase("advance_wp2_wp3_bad_wp2", "wp2",
#                                    "Error calling advance_wp2_wp3"),
#
#       This results in "! advance_wp2_wp3_bad_wp2" in the code being replaced
#       with a line that sets "wp2" to NaN, then the error messages that clubb
#       produces are saved in a file named "advance_wp2_wp3_bad_wp2" in run_scripts
#
#   Notes:
#           THIS SCRIPT NO LONGER NEEDS TO BE RUN FROM WITHIN 'run_scripts'
#           BECAUSE IT RESOLVES PATHS RELATIVE TO ITS OWN LOCATION
#
#       To add error tests, simply append a TestCase to the TESTS list with a
#       unique name matching a comment tag in the Fortran source.
#
###############################################################################

from __future__ import annotations

import shutil
import signal
import subprocess
import sys
from pathlib import Path
from typing import NamedTuple


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent

SOURCE_FILE = REPO_ROOT / "src" / "CLUBB_core" / "advance_clubb_core_module.F90"
TOOLCHAIN_FILE = REPO_ROOT / "cmake" / "toolchains" / "linux_x86_64_gcc.cmake"
COMPILE_SCRIPT = REPO_ROOT / "compile.py"
RUN_SCM_SCRIPT = SCRIPT_DIR / "run_scm.py"
CASE_LOG = REPO_ROOT / "output" / "arm_97_log"

class TestCase(NamedTuple):
    name: str            # matches the "! <name>" comment tag in the Fortran source
    variable: str        # the Fortran variable that will be set to NaN
    expected_error: str  # substring expected in CLUBB's error output


########### TEST DEFINITIONS  ##################
TESTS: list[TestCase] = [
    TestCase("advance_wp2_wp3_bad_wp2", "wp2", "Error calling advance_wp2_wp3"),
    TestCase("advance_xm_wpxp_bad_wp2", "wp2", "Error calling advance_xm_wpxp"),
]

ORIGINAL_TEXT: dict[Path, str | None] = {}


def save_original(path: Path) -> None:
    """Snapshot a file's contents before we mutate it, so restore_backups()
    can put it back even if the script crashes or is interrupted mid-edit."""
    if path not in ORIGINAL_TEXT:
        # Store None for missing files so restore_backups() knows to delete them
        ORIGINAL_TEXT[path] = path.read_text(encoding="utf-8") if path.exists() else None


def restore_backups() -> None:
    """Undo every mutation made during the run.  Called both on normal exit
    (finally block) and on SIGINT/SIGTERM so the source tree is never left
    in a half-patched state that would break subsequent builds."""
    for path, text in ORIGINAL_TEXT.items():
        if text is None:
            if path.exists():
                path.unlink()
        else:
            path.write_text(text, encoding="utf-8")


def cleanup_handler(sig: int, frame: object) -> None:
    restore_backups()
    raise SystemExit(1)


def run_command(cmd: list[str], cwd: Path) -> tuple[int, str]:
    """Stream output line-by-line so long compiles/runs show progress in CI,
    while also capturing everything for later log-file analysis."""
    process = subprocess.Popen(
        cmd,
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,  # merge stderr into stdout so we capture everything
        text=True,
        errors="replace",  # don't crash on non-UTF-8 bytes from Fortran runtime output
        bufsize=1,  # line-buffered so we can stream in real time
    )

    output_lines: list[str] = []
    assert process.stdout is not None  # guaranteed by PIPE, but keeps mypy happy
    for line in process.stdout:
        print(line, end="", flush=True)  # mirror to terminal for CI visibility
        output_lines.append(line)

    process.wait()  # stdout is already drained; just collect the exit code
    return process.returncode, "".join(output_lines)


def replace_once(path: Path, old: str, new: str) -> None:
    """Single-occurrence substitute that hard-fails if the marker is missing,
    so we notice immediately when the Fortran source moves the comment tag."""
    text = path.read_text(encoding="utf-8")

    if old not in text:
        raise RuntimeError(f"Expected text not found in {path}: {old}")

    path.write_text(text.replace(old, new, 1), encoding="utf-8")  # count=1: only first match


def patch_toolchain_for_nan_test(path: Path) -> None:
    """GCC's default flags trap on NaN/Inf and refuse out-of-range bit patterns.
    We need both disabled: -fno-range-check lets the transfer() NaN literal
    compile, and removing -ffpe-trap lets the NaN propagate at runtime so
    CLUBB's own error-checking (not the FPE signal) is what catches it."""
    text = path.read_text(encoding="utf-8")

    # Inject -fno-range-check so transfer(2143289344, 1.0) compiles without error
    common_flags = 'set(CMAKE_Fortran_FLAGS         "${COMMON_FLAGS} ${WARNINGS}")'
    common_flags_patched = 'set(CMAKE_Fortran_FLAGS         "-fno-range-check ${COMMON_FLAGS} ${WARNINGS}")'
    if common_flags in text:
        text = text.replace(common_flags, common_flags_patched, 1)
    elif common_flags_patched not in text:  # already patched is fine; anything else is wrong
        raise RuntimeError(f"Expected toolchain line not found in {path}")

    # Strip the FPE trap flag — whitespace varies by position so check both forms
    debug_trap = "-ffpe-trap=invalid,zero,overflow "
    debug_trap_alt = " -ffpe-trap=invalid,zero,overflow"
    if debug_trap in text:
        text = text.replace(debug_trap, "", 1)
    elif debug_trap_alt in text:
        text = text.replace(debug_trap_alt, "", 1)
    # If neither is present it was already stripped — that's fine

    path.write_text(text, encoding="utf-8")


def write_test_log(test_name: str, captured_output: str) -> None:
    """Prefer the on-disk case log (written by CLUBB itself) over captured
    stdout, because CLUBB may flush diagnostics to the log that never reach
    stdout when it aborts early."""
    log_path = SCRIPT_DIR / test_name

    if CASE_LOG.exists():
        shutil.copy2(CASE_LOG, log_path)
    else:
        log_path.write_text(captured_output, encoding="utf-8")


def expected_error_detected(test_name: str, expected_error: str) -> bool:
    """Check the saved log for the expected error string. Uses the log file
    (not captured output) so it stays consistent with write_test_log's
    preference for CLUBB's own log."""
    log_path = SCRIPT_DIR / test_name

    if not log_path.exists():
        return False

    return expected_error in log_path.read_text(encoding="utf-8", errors="replace")


def main() -> int:

    # Make signals (like Ctrl-C) trigger cleanup of the Fortran source
    signal.signal(signal.SIGINT, cleanup_handler)
    signal.signal(signal.SIGTERM, cleanup_handler)

    results: list[tuple[str, str, bool]] = []

    try:

        # Save original versions of files we will mutate
        save_original(SOURCE_FILE)
        save_original(TOOLCHAIN_FILE)

        # Save compiler flags so they can be restored after the test
        patch_toolchain_for_nan_test(TOOLCHAIN_FILE)

        for test in TESTS:

            # Find test name and replace with line that sets test variable to NaN
            replace_once(
                SOURCE_FILE,
                f"! {test.name}",
                f"{test.variable} = transfer( 2143289344, 1.0 )",  # 2143289344 = 0x7FC00000 = quiet NaN
            )

            #compile
            compile_cmd = [sys.executable, str(COMPILE_SCRIPT)]
            compile_rc, _ = run_command(compile_cmd, REPO_ROOT)

            if compile_rc != 0:
                raise RuntimeError(f"Compilation failed for {test.name}")

            # Run a case (any non error causing case will work) and write error output to a file
            run_cmd = [sys.executable, str(RUN_SCM_SCRIPT), "arm_97"]
            run_rc, run_output = run_command(run_cmd, REPO_ROOT)
            write_test_log(test.name, run_output)

            # Nonzero exit is expected (CLUBB should abort on the NaN).
            # Only fail if we got nothing at all — no log file AND no stdout.
            if run_rc != 0 and not CASE_LOG.exists() and not run_output.strip():
                raise RuntimeError(f"run_scm.py failed for {test.name}")

            detected = expected_error_detected(test.name, test.expected_error)
            results.append((test.name, test.expected_error, detected))

            # Revert the NaN injection so the next test starts from clean source
            replace_once(
                SOURCE_FILE,
                f"{test.variable} = transfer( 2143289344, 1.0 )",
                f"! {test.name}",
            )

        print("\nError detection summary:")
        all_detected = True

        for test_name, expected_error, detected in results:
            if detected:
                print(f"PASS: {test_name}: found '{expected_error}'")
            else:
                print(f"FAIL: {test_name}: missing '{expected_error}'")
                all_detected = False

        return 0 if all_detected else 1

    finally:
        # Restore compiler flags
        restore_backups()


if __name__ == "__main__":
    raise SystemExit(main())  # raise instead of sys.exit so finally blocks run
