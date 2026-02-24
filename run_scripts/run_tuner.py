#!/usr/bin/env python3
"""Run CLUBB tuner using Python workflows and install/latest executables."""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path


RUN_SCRIPTS = Path(__file__).resolve().parent
CLUBB_ROOT = RUN_SCRIPTS.parent
INPUT_DIR = CLUBB_ROOT / "input"
OUTPUT_DIR = CLUBB_ROOT / "output"
TUNABLE_DIR = INPUT_DIR / "tunable_parameters"
CASE_SETUP_DIR = INPUT_DIR / "case_setups"
STATS_DIR = INPUT_DIR / "stats"
TUNER_MISC_DIR = CLUBB_ROOT / "input_misc" / "tuner"

RUN_SCM = RUN_SCRIPTS / "run_scm.py"
TUNER_EXE = CLUBB_ROOT / "install" / "latest" / "clubb_tuner"


def strip_fortran_comments(text: str) -> str:
    lines = []
    for line in text.splitlines():
        lines.append(re.sub(r"!.*", "", line))
    return "\n".join(lines) + "\n"


def concat_stripped_files(sources: list[Path], destination: Path) -> None:
    text = "".join(strip_fortran_comments(path.read_text(encoding="utf-8")) for path in sources)
    destination.write_text(text, encoding="utf-8")


def require_file(path: Path) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"{path} does not exist")


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

    with (log_path.open("a", encoding="utf-8") if log_path is not None else open("/dev/null", "w", encoding="utf-8")) as log:
        proc = subprocess.Popen(
            cmd,
            cwd=str(cwd),
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
    params_file: Path,
    flags_file: Path,
    silhs_params_file: Path,
    stats_file: Path,
    log_path: Path,
) -> int:
    cmd = [
        sys.executable,
        str(RUN_SCM),
        "-params",
        str(params_file),
        "-flags",
        str(flags_file),
        "-silhs_params",
        str(silhs_params_file),
        "-stats",
        str(stats_file),
        case_name,
    ]
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
    parser = argparse.ArgumentParser(description="Run CLUBB tuner with Python tooling.")
    parser.add_argument("-i", "--initial-output", action="store_true", help="Save pre-tuning outputs.")
    parser.add_argument("-n", "--nightly", action="store_true", help="Enable nightly mode behavior.")
    parser.add_argument(
        "--run-type",
        choices=("single", "multiple"),
        default="single",
        help="Tuning run type.",
    )
    parser.add_argument("--run-case", default="fire", help="Case name used for tuner input selection.")
    parser.add_argument(
        "--model-mult",
        nargs="+",
        default=["fire", "atex"],
        help="Case names for multiple-case tuning.",
    )
    parser.add_argument(
        "--stats-tune",
        default=str(STATS_DIR / "tuning_stats.in"),
        help="Stats file used in *_hoc.in for tuner.",
    )
    parser.add_argument(
        "--stats-opt",
        default=str(STATS_DIR / "standard_stats.in"),
        help="Stats file used for standalone output runs.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    params_file = TUNABLE_DIR / "tunable_parameters.in"
    silhs_params_file = TUNABLE_DIR / "silhs_parameters.in"
    flags_file = TUNABLE_DIR / "configurable_model_flags.in"
    stats_tune_file = Path(args.stats_tune).resolve()
    stats_opt_file = Path(args.stats_opt).resolve()
    error_in_file = TUNER_MISC_DIR / f"error_{args.run_case}.in"
    rand_seed_file = TUNER_MISC_DIR / "rand_seed.dat"

    try:
        for required in (
            RUN_SCM,
            TUNER_EXE,
            params_file,
            silhs_params_file,
            flags_file,
            stats_tune_file,
            stats_opt_file,
            error_in_file,
            rand_seed_file,
        ):
            require_file(required)
    except FileNotFoundError as exc:
        print(str(exc))
        return 1

    cases_for_runs = [args.run_case] if args.run_type == "single" else args.model_mult
    if args.run_type == "multiple" and not args.model_mult:
        print("No cases specified for --model-mult")
        return 1

    model_files: dict[str, Path] = {}
    try:
        for case in cases_for_runs:
            model_file = CASE_SETUP_DIR / f"{case}_model.in"
            require_file(model_file)
            model_files[case] = model_file
    except FileNotFoundError as exc:
        print(str(exc))
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
        for case in cases_for_runs:
            hoc_path = RUN_SCRIPTS / f"{case}_hoc.in"
            concat_stripped_files(
                [model_files[case], stats_tune_file, flags_file, silhs_params_file],
                hoc_path,
            )
            hoc_files.append(hoc_path)

        error_text = strip_fortran_comments(error_in_file.read_text(encoding="utf-8"))
        params_text = strip_fortran_comments(params_file.read_text(encoding="utf-8"))
        error_tmp.write_text(error_text + params_text, encoding="utf-8")
        shutil.copy2(rand_seed_file, rand_seed_tmp)

        if args.initial_output and not args.nightly:
            initial_output_dir = CLUBB_ROOT / "initial_output"
            initial_output_dir.mkdir(parents=True, exist_ok=True)
            for case in cases_for_runs:
                status = run_case_with_run_scm(
                    case,
                    params_file=params_file,
                    flags_file=flags_file,
                    silhs_params_file=silhs_params_file,
                    stats_file=stats_opt_file,
                    log_path=log_path,
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
                f.write(f"contents of {error_in_file}\n")
                f.write(error_in_file.read_text(encoding="utf-8"))
        else:
            print("Warning: No tuning_run_results_* file found in input/")

        print("Running with the optimal parameter set")
        tuned_params = newest_matching("tunable_parameters_*", TUNABLE_DIR)
        tuned_flags = newest_matching("configurable_model_flags*", TUNABLE_DIR)
        if tuned_params is None or tuned_flags is None:
            print("Could not find tuned parameter/flag files in input/tunable_parameters")
            overall_fail = True
        else:
            for case in cases_for_runs:
                status = run_case_with_run_scm(
                    case,
                    params_file=tuned_params,
                    flags_file=tuned_flags,
                    silhs_params_file=silhs_params_file,
                    stats_file=stats_opt_file,
                    log_path=log_path,
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
