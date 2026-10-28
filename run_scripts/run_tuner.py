#!/usr/bin/env python3
"""Run CLUBB tuner using Python workflows and install/latest executables."""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path


RUN_SCRIPTS = Path(__file__).resolve().parent
CLUBB_ROOT = RUN_SCRIPTS.parent
if str(CLUBB_ROOT) not in sys.path:
    sys.path.insert(0, str(CLUBB_ROOT))
from input.scripts.create_case_namelist import (
    DEFAULT_OUTPUT_DIR, case_run_name, create_case_settings, resolve_case_model_input,
    write_case_namelist,
)
from input.scripts.stats_json_to_namelist import select_stats_entries
INPUT_DIR = CLUBB_ROOT / "input"
OUTPUT_DIR = CLUBB_ROOT / "output"
TUNABLE_DIR = INPUT_DIR
TUNER_MISC_DIR = CLUBB_ROOT / "input_misc" / "tuner"

RUN_SCM = RUN_SCRIPTS / "run_scm.py"
TUNER_EXE = CLUBB_ROOT / "install" / "latest" / "clubb_tuner"


def require_file(path: Path) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"{path} does not exist")


def native_controls_path(selection):
    """Resolve a legacy control name or explicit namelist path without reading JSON."""
    path = Path(selection).expanduser()
    return path.resolve() if path.suffix == ".in" or path.is_absolute() else TUNER_MISC_DIR / f"error_{selection}.in"


def create_native_tuner_inputs(
    case_names, controls, output_dir, *, stats="tuning_stats", override=None, params=None,
):
    """Write final native inputs from legacy controls and modern model settings.

    controls selects error_NAME.in under input_misc/tuner or an explicit .in
    path. Model cases and parameters stay in Python until their final files are
    written. Return error.in, the case input paths and the native tuning type.
    """
    from clubb_python_driver.io.namelist import _read_namelist_groups
    from input.scripts.namelist import atomic_text
    from input.scripts.parameters_json_to_namelist import parameter_settings, parameters_to_namelist

    control_file = native_controls_path(controls)
    control_text = control_file.read_text(encoding="utf-8")
    # Match the ESA default declared in src/error.F90 when a custom file omits it.
    tune_type = int(_read_namelist_groups(control_file).get("stats", {}).get("tune_type", 2))
    parameter_text = parameters_to_namelist(parameter_settings(params))
    destination = Path(output_dir)
    cases = [(case_run_name(case), create_case_settings(case, DEFAULT_OUTPUT_DIR, stats=stats, override=override))
             for case in case_names]
    paths = [write_case_namelist(settings, destination / f"{name}_hoc.in") for name, settings in cases]
    error_path = destination / "error.in"
    atomic_text(error_path, control_text.rstrip() + "\n\n" + parameter_text)
    return error_path, paths, tune_type


def newest_matching(pattern: str, directory: Path) -> Path | None:
    matches = sorted(directory.glob(pattern), key=lambda p: p.stat().st_mtime, reverse=True)
    return matches[0] if matches else None


def run_and_tee(
    cmd: list[str],
    cwd: Path,
    log_path: Path | None = None,
    stdin_text: str | None = None,
) -> int:
    print(f"+ {' '.join(cmd)}")
    cwd.mkdir(parents=True, exist_ok=True)
    if log_path is not None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
    run_env = os.environ.copy()
    run_env["PWD"] = str(cwd.resolve())

    with (log_path.open("a", encoding="utf-8") if log_path is not None else open("/dev/null", "w", encoding="utf-8")) as log:
        proc = subprocess.Popen(
            cmd,
            cwd=str(cwd),
            env=run_env,
            stdin=subprocess.PIPE if stdin_text is not None else None,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            errors="replace",
        )
        if stdin_text is not None and proc.stdin is not None:
            proc.stdin.write(stdin_text)
            proc.stdin.close()
        assert proc.stdout is not None
        for line in proc.stdout:
            print(line, end="")
            if log_path is not None:
                log.write(line)
        proc.wait()
        return proc.returncode


def run_case_with_run_scm(
    case_name: str,
    params_file: Path | str | None,
    flags_file: Path | str | None,
    silhs_params_file: Path | str | None,
    stats_file: str,
    log_path: Path,
    override=None,
) -> int:
    """Run one SCM case, forwarding explicit selections and inheriting omitted defaults."""
    cmd = [
        sys.executable,
        str(RUN_SCM),
        "-stats",
        str(stats_file),
    ]
    for option, selection in (("-params", params_file), ("-flags", flags_file),
                              ("-silhs_params_file", silhs_params_file)):
        if selection is not None:
            cmd.extend([option, str(selection)])
    cmd.append(case_name)
    if override:
        cmd += ["-override", override]
    return run_and_tee(cmd, cwd=CLUBB_ROOT, log_path=log_path)


def move_case_outputs(case_name: str, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    patterns = [
        f"{case_name}_stats*",
        f"{case_name}_zt*",
        f"{case_name}_zm*",
        f"{case_name}_sfc*",
        f"{case_name}_lh_zt*",
        f"{case_name}_lh_sfc*",
        f"{case_name}_rad*",
        f"{case_name}_setup*",
    ]
    for pattern in patterns:
        for path in OUTPUT_DIR.glob(pattern):
            shutil.move(str(path), str(destination / path.name))


def copy_nightly_artifacts(log_path: Path) -> None:
    nightly_dir = Path.home() / "tuner_output"
    nightly_dir.mkdir(parents=True, exist_ok=True)

    for pattern in ("tuning_run_results*", "error_*"):
        for path in INPUT_DIR.glob(pattern):
            shutil.copy2(path, nightly_dir / path.name)

    for path in TUNABLE_DIR.glob("tunable_parameters_*"):
        shutil.copy2(path, nightly_dir / path.name)

    if log_path.is_file():
        shutil.copy2(log_path, nightly_dir / log_path.name)


def parse_args() -> argparse.Namespace:
    """Read native execution/preparation options and validate the output-directory scope."""
    parser = argparse.ArgumentParser(description="Run CLUBB tuner with Python tooling.", add_help=False, allow_abbrev=False)
    parser.add_argument("-h", "-help", action="help", help="Show this help and exit.")
    parser.add_argument('-initial_output', dest='initial_output', action="store_true", help="Save pre-tuning outputs.")
    parser.add_argument('-nightly', dest='nightly', action="store_true", help="Enable nightly mode behavior.")
    parser.add_argument(
        '-run_type', dest='run_type',
        choices=("single", "multiple"),
        default="single",
        help="Tuning run type.",
    )
    parser.add_argument('-run_case', dest='run_case', default="fire", help="Case name used for tuner input selection.")
    parser.add_argument(
        '-model_mult', dest='model_mult',
        nargs="+",
        default=["fire", "atex"],
        help="Case names for multiple-case tuning.",
    )
    parser.add_argument(
        '-stats_tune', dest='stats_tune',
        default="tuning_stats",
        help="Stats file used in *_hoc.in for tuner.",
    )
    parser.add_argument(
        '-stats_opt', dest='stats_opt',
        default="standard",
        help="Stats file used for standalone output runs.",
    )
    parser.add_argument("-override", help="Model namelist edits applied to tuning and output runs.")
    parser.add_argument("-controls", metavar="NAME|IN", help="Legacy tuner controls; defaults to the selected case's benchmark name.")
    parser.add_argument("-params", metavar="NAME|JSON", help="Initial parameter preset or JSON overrides.")
    parser.add_argument("-prepare_only", action="store_true", help="Write error.in and CASE_hoc.in without launching the tuner.")
    parser.add_argument("-output_dir", type=Path, help="Prepared-input directory; requires -prepare_only (default: run_scripts).")
    args = parser.parse_args()
    if args.output_dir is not None and not args.prepare_only:
        parser.error("-output_dir requires -prepare_only")
    return args


def main() -> int:
    """Prepare native inputs, run the tuner and pass optimized JSON back to SCM."""
    args = parse_args()

    stats_opt_file = args.stats_opt
    select_stats_entries(args.stats_tune)
    select_stats_entries(stats_opt_file)
    comparison_case = args.controls or resolve_case_model_input(args.run_case).get("benchmark_case") or args.run_case
    controls_file = native_controls_path(comparison_case)
    rand_seed_file = TUNER_MISC_DIR / "rand_seed.dat"
    cases_for_runs = [args.run_case] if args.run_type == "single" else args.model_mult
    if args.prepare_only:
        error_path, _, _ = create_native_tuner_inputs(
            cases_for_runs, comparison_case, args.output_dir or RUN_SCRIPTS,
            stats=args.stats_tune, override=args.override, params=args.params)
        print(error_path)
        return 0
    try:
        for required in (RUN_SCM, TUNER_EXE, controls_file, rand_seed_file):
            require_file(required)
    except FileNotFoundError as exc:
        print(str(exc))
        return 1

    if args.run_type == "multiple" and not args.model_mult:
        print("No cases specified for -model_mult")
        return 1

    hoc_files: list[Path] = []
    error_tmp = RUN_SCRIPTS / "error.in"
    rand_seed_tmp = RUN_SCRIPTS / "rand_seed.dat"
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    date_str = datetime.now().strftime("%Y-%m-%d")
    log_name = f"tuner_{args.run_case}_{date_str}.log" if args.nightly else "tuner_tmp.log"
    log_path = OUTPUT_DIR / log_name
    tuner_exit = 1
    overall_fail = False

    try:
        error_tmp, hoc_files, tune_type = create_native_tuner_inputs(
            cases_for_runs, comparison_case, RUN_SCRIPTS, stats=args.stats_tune,
            override=args.override, params=args.params)
        shutil.copy2(rand_seed_file, rand_seed_tmp)

        if args.initial_output and not args.nightly:
            initial_output_dir = OUTPUT_DIR / "tuner" / "initial_output"
            initial_output_dir.mkdir(parents=True, exist_ok=True)
            for case in cases_for_runs:
                status = run_case_with_run_scm(
                    case,
                    params_file=args.params,
                    flags_file=None,
                    silhs_params_file=None,
                    stats_file=stats_opt_file,
                    log_path=log_path,
                    override=args.override,
                )
                if status != 0:
                    overall_fail = True
                move_case_outputs(case, initial_output_dir)

        print(f"Tuning {args.run_case}")
        # Non-interactive mode: answer the tuner "continue?" prompt with "n".
        tuner_exit = run_and_tee(
            [str(TUNER_EXE)],
            cwd=RUN_SCRIPTS,
            log_path=log_path,
            stdin_text="n\n",
        )
        if tuner_exit != 0:
            overall_fail = True

        if args.nightly:
            copy_nightly_artifacts(log_path)

        newest_tuner_log = newest_matching("tuning_run_results_*", INPUT_DIR)
        if newest_tuner_log is not None:
            with newest_tuner_log.open("a", encoding="utf-8") as f:
                f.write("\n********************************\n")
                f.write(f"contents of {controls_file}\n")
                f.write(controls_file.read_text(encoding="utf-8"))
        else:
            print("Warning: No tuning_run_results_* file found in input/")

        print("Running with the optimal parameter set")
        tuned_params = newest_matching("tunable_parameters_*.json", TUNABLE_DIR) if tune_type != 3 else args.params
        tuned_flags = newest_matching("configurable_model_flags_*.json", TUNABLE_DIR) if tune_type == 3 else None
        if tuner_exit != 0 or (tune_type != 3 and tuned_params is None) or (tune_type == 3 and tuned_flags is None):
            print("Native tuning failed or did not produce the expected parameter/flag JSON")
            overall_fail = True
        else:
            for case in cases_for_runs:
                status = run_case_with_run_scm(
                    case,
                    params_file=tuned_params,
                    flags_file=tuned_flags,
                    silhs_params_file=None,
                    stats_file=stats_opt_file,
                    log_path=log_path,
                    override=args.override,
                )
                if status != 0:
                    overall_fail = True

        if log_path.is_file():
            log_text = log_path.read_text(encoding="utf-8", errors="replace")
            if "Fatal error" in log_text:
                overall_fail = True

        return 1 if overall_fail else tuner_exit
    finally:
        if error_tmp.exists():
            error_tmp.unlink()
        if rand_seed_tmp.exists():
            rand_seed_tmp.unlink()
        for hoc_file in hoc_files:
            if hoc_file.exists():
                hoc_file.unlink()


if __name__ == "__main__":
    sys.exit(main())
