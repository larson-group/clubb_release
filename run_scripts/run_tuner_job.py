#!/usr/bin/env python3
"""Command-line front door for the CLUBB tuner.

Usage
-----
The smallest useful tuner run names a case, one or more loss fields, and one or
more parameter ranges:

    python run_scripts/run_tuner_job.py \
      -cases bomex \
      -fields cloud_frac \
      -params C8:0.2:0.8

Cases can be plain names, fixed windows, or split windows:

    bomex
    bomex:7200:18000
    bomex:7200:18000:2700
    case:t_start:t_end:t_out

Parameters use PARAM:MIN:MAX. Strategies use MODE:SETTING, for example
random:16 or resolve:0.9. Interactive runs ask whether to rerun the best result
after tuning; Jenkins and other non-interactive callers should set -run_top to
never, complete, window, or both.

Understanding
-------------
This script is a controller, not a loss calculator.

It creates a TunerJob directory under output_tuner and writes the small files
that describe the job:

  - request.json: cases, fields, parameters, strategy, and loss settings
  - control.json: stop requests and the controller keepalive lease
  - status.json: live scheduler progress
  - results.json: completed samples and rankings

The actual work is done by the scheduler launched as:

    python -m tuner.tune_clubb --job-dir <job_dir>

While the scheduler runs, this wrapper polls status.json/results.json, prints a
compact status line, and renews the keepalive in control.json. If the wrapper or
Dash app disappears and the keepalive is not refreshed for five minutes, the
scheduler treats that as a graceful stop request.

Ctrl-C is intentionally two-stage: the first press asks the scheduler to stop
cleanly, and the second press terminates the subprocess owned by this wrapper.

After a finished or stopped run, the optional post-run step writes a temporary
tunable_parameters-style file for the best selected parameters and forwards it
to run_scm.py, run_scm_loss.py, or both.
"""

from __future__ import annotations

import argparse
from datetime import datetime
import math
import os
from pathlib import Path
import subprocess
import sys
import time


RUN_SCRIPTS = Path(__file__).resolve().parent
CLUBB_ROOT = RUN_SCRIPTS.parent
OUTPUT_TUNER_DIR = CLUBB_ROOT / "output_tuner"
CLUBB_OUTPUT_DIR = CLUBB_ROOT / "output"
DEFAULT_TUNABLE_PARAMS = CLUBB_ROOT / "input" / "tunable_parameters" / "tunable_parameters.in"
RUN_SCM = RUN_SCRIPTS / "run_scm.py"
RUN_SCM_LOSS = RUN_SCRIPTS / "run_scm_loss.py"

if str(CLUBB_ROOT) not in sys.path:
    sys.path.insert(0, str(CLUBB_ROOT))

from tuner.job_runtime import TERMINAL_STATES, TunerJob, tuner_worker_env  # noqa: E402
from tuner.status import read_json_or_default  # noqa: E402
from tuner.taylor_metrics import DEFAULT_AGGREGATION_MODE, DEFAULT_LOSS_MODE  # noqa: E402


def _float_from_text(value: str, label: str) -> float:
    try:
        return float(str(value).replace("D", "E").replace("d", "e"))
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"{label} must be numeric: {value}") from exc


def _int_from_text(value: str, label: str) -> int:
    numeric = _float_from_text(value, label)
    if int(numeric) != numeric:
        raise argparse.ArgumentTypeError(f"{label} must be an integer: {value}")
    return int(numeric)


def _split_values(values: list[str] | None) -> list[str]:
    result = []
    for value in values or []:
        for item in str(value).split(","):
            item = item.strip()
            if item:
                result.append(item)
    return result


def parse_case_spec(spec: str) -> dict:
    """Parse CASE[:T_START:T_END[:T_INTERVAL]]."""
    parts = [part.strip() for part in str(spec).split(":")]
    if len(parts) not in (1, 3, 4) or any(part == "" for part in parts):
        raise argparse.ArgumentTypeError(
            f"Invalid case spec '{spec}'. Use CASE, CASE:T_START:T_END, "
            "or CASE:T_START:T_END:T_INTERVAL."
        )

    case_config = {"case_name": parts[0]}
    if len(parts) >= 3:
        t_start = _int_from_text(parts[1], f"{parts[0]} t_start")
        t_end = _int_from_text(parts[2], f"{parts[0]} t_end")
        if t_end <= t_start:
            raise argparse.ArgumentTypeError(f"{parts[0]} requires t_end > t_start")
        case_config["time_average_range"] = [t_start, t_end]

    if len(parts) == 4:
        average_time_seconds = _int_from_text(parts[3], f"{parts[0]} t_interval")
        if average_time_seconds < 1:
            raise argparse.ArgumentTypeError(f"{parts[0]} t_interval must be >= 1")
        window_seconds = case_config["time_average_range"][1] - case_config["time_average_range"][0]
        if window_seconds % average_time_seconds != 0:
            raise argparse.ArgumentTypeError(f"{parts[0]} t_interval must divide t_end - t_start evenly")
        case_config["average_time_seconds"] = average_time_seconds
        case_config["num_time_windows"] = window_seconds // average_time_seconds

    return case_config


def parse_param_spec(spec: str) -> dict:
    """Parse PARAM:MIN:MAX."""
    parts = [part.strip() for part in str(spec).split(":")]
    if len(parts) != 3 or any(part == "" for part in parts):
        raise argparse.ArgumentTypeError(f"Invalid param spec '{spec}'. Use PARAM:MIN:MAX.")
    min_value = _float_from_text(parts[1], f"{parts[0]} min")
    max_value = _float_from_text(parts[2], f"{parts[0]} max")
    if min_value > max_value:
        raise argparse.ArgumentTypeError(f"{parts[0]} requires min <= max")
    return {"name": parts[0], "min": min_value, "max": max_value}


def parse_strategy_spec(spec: str) -> dict:
    """Parse random:MAX_SAMPLES or resolve:SPACING."""
    parts = [part.strip() for part in str(spec).split(":")]
    name = parts[0].lower() if parts and parts[0] else ""
    if name == "random":
        if len(parts) != 2 or not parts[1]:
            raise argparse.ArgumentTypeError("random strategy requires random:MAX_SAMPLES")
        max_samples = _int_from_text(parts[1], "random max_samples")
        if max_samples < 1:
            raise argparse.ArgumentTypeError("random max_samples must be >= 1")
        return {"name": "random", "options": {"max_samples": max_samples}}

    if name == "resolve":
        if len(parts) != 2 or not parts[1]:
            raise argparse.ArgumentTypeError("resolve strategy requires resolve:SPACING")
        spacing = _float_from_text(parts[1], "resolve spacing")
        if spacing <= 0.0:
            raise argparse.ArgumentTypeError("resolve spacing must be > 0")
        return {"name": "resolve", "options": {"spacing": spacing}}

    raise argparse.ArgumentTypeError("strategy must be random:MAX_SAMPLES or resolve:SPACING")


def build_request(args: argparse.Namespace) -> dict:
    case_configs = [parse_case_spec(spec) for spec in _split_values(args.cases)]
    if not case_configs:
        raise argparse.ArgumentTypeError("-cases must include at least one case")

    selected_fields = _split_values(args.fields)
    if not selected_fields:
        raise argparse.ArgumentTypeError("-fields must include at least one field")

    parameter_ranges = [parse_param_spec(spec) for spec in _split_values(args.params)]
    if not parameter_ranges:
        raise argparse.ArgumentTypeError("-params must include at least one PARAM:MIN:MAX range")

    request = {
        "cases": [config["case_name"] for config in case_configs],
        "case_configs": case_configs,
        "selected_fields": selected_fields,
        "parameter_ranges": parameter_ranges,
        "batch_size": int(args.batch_size),
        "max_workers": int(args.max_workers),
        "strategy": parse_strategy_spec(args.strategy),
        "loss_mode": args.loss_mode,
        "aggregation_mode": args.aggregation_mode,
        "time_window_aggregation_mode": args.aggregation_mode,
    }
    if args.seed is not None:
        request["seed"] = int(args.seed)
    return request


def format_seconds(seconds: float | int | None) -> str:
    try:
        total = int(float(seconds or 0))
    except (TypeError, ValueError):
        total = 0
    hours, remainder = divmod(total, 3600)
    minutes, sec = divmod(remainder, 60)
    if hours:
        return f"{hours:d}:{minutes:02d}:{sec:02d}"
    return f"{minutes:02d}:{sec:02d}"


def format_params(params: dict, *, limit: int = 4) -> str:
    items = list((params or {}).items())
    pieces = [f"{name}={float(value):.6g}" for name, value in items[:limit]]
    if len(items) > limit:
        pieces.append("...")
    return ", ".join(pieces)


def format_strategy(strategy: dict) -> str:
    """Return a compact human-readable strategy description."""
    name = strategy.get("name", "unknown")
    options = dict(strategy.get("options") or {})
    if name == "random":
        return f"random, max_samples={options.get('max_samples', 'unbounded')}"
    if name == "resolve":
        return f"resolve, spacing={options.get('spacing', '--')}"
    return str(name)


def estimate_total_samples(request: dict) -> str:
    """Return the sample count implied by the CLI request, if known."""
    strategy = request.get("strategy") or {}
    options = strategy.get("options") or {}
    if strategy.get("name") == "random":
        return str(options.get("max_samples", "unknown"))
    if strategy.get("name") == "resolve":
        spacing = float(options.get("spacing"))
        total = 1
        for spec in request.get("parameter_ranges", []):
            span = float(spec["max"]) - float(spec["min"])
            intervals = 1 if span == 0.0 else max(1, int(math.ceil(span / spacing)))
            total *= intervals + 1
        return str(total)
    return "unknown"


def render_status_line(status: dict) -> str:
    state = status.get("state", "unknown")
    samples = int(status.get("samples_evaluated", 0) or 0)
    total = status.get("total_samples")
    sample_text = f"{samples}/{total}" if total is not None else f"{samples}"
    elapsed = format_seconds(status.get("elapsed_seconds", 0.0))
    active = int(status.get("active_evaluations", 0) or 0)
    idle = int(status.get("idle_workers", 0) or 0)
    initialized = int(status.get("initialized_workers", 0) or 0)
    queued = int(status.get("queued_case_jobs", 0) or 0)
    best = status.get("best_total_loss")
    if best is None:
        best_text = "best: --"
    else:
        params = {}
        top_results = status.get("top_results") or []
        if top_results:
            params = top_results[0].get("params", {})
        best_text = f"best: {float(best):.6g}"
        if params:
            best_text += f" ({format_params(params)})"
    return (
        f"state: {state} | samples: {sample_text} | elapsed: {elapsed} | "
        f"workers active/idle/init: {active}/{idle}/{initialized} | queued cases: {queued} | {best_text}"
    )


def print_table_rows(rows: list[tuple[str, str]], *, indent: str = "  ") -> None:
    """Print two-column rows with the left column aligned."""
    width = max((len(label) for label, _value in rows), default=0)
    for label, value in rows:
        print(f"{indent}{label:<{width}}  {value}")


def print_request_summary(
    request: dict,
    job_dir: Path,
    request_path: Path,
    log_path: Path,
    args: argparse.Namespace,
) -> None:
    print()
    print("Tuner Job")
    print("---------")
    print_table_rows(
        [
            ("job_dir", str(job_dir)),
            ("request", str(request_path)),
            ("log", str(log_path)),
        ]
    )

    print()
    print("Configuration")
    print("  cases")
    for config in request["case_configs"]:
        case_name = config["case_name"]
        time_range = config.get("time_average_range")
        windows = config.get("num_time_windows", "default")
        average_time = config.get("average_time_seconds")
        if time_range is None:
            detail = f"default window, windows={windows}"
        elif average_time is not None:
            detail = f"{time_range[0]}-{time_range[1]} s, avg={average_time} s, windows={windows}"
        else:
            detail = f"{time_range[0]}-{time_range[1]} s, windows={windows}"
        print(f"    {case_name:<12} {detail}")

    print("  fields")
    print("    " + ", ".join(request["selected_fields"]))

    print("  parameters")
    max_param_len = max((len(spec["name"]) for spec in request["parameter_ranges"]), default=0)
    for spec in request["parameter_ranges"]:
        print(f"    {spec['name']:<{max_param_len}}  {spec['min']:g} -> {spec['max']:g}")

    print("  run")
    print_table_rows(
        [
            ("strategy", format_strategy(request["strategy"])),
            ("estimated_samples", estimate_total_samples(request)),
            ("batch_size", str(request["batch_size"])),
            ("max_workers", str(request["max_workers"])),
            ("loss", request["loss_mode"]),
            ("aggregation", request["aggregation_mode"]),
            ("poll_interval", f"{args.poll_interval:g} s"),
            ("post_run", f"{args.run_top}, top_n={args.top_n}, out={Path(args.run_out_dir).resolve()}"),
        ],
        indent="    ",
    )
    print()


def print_output_locations(
    job_dir: Path,
    request_path: Path,
    status_path: Path,
    results_path: Path,
    log_path: Path,
    *,
    run_output_dirs: dict[str, Path] | None = None,
) -> None:
    print()
    print("Saved Output")
    print("------------")
    rows = [
        ("job_dir", str(job_dir)),
        ("results", str(results_path)),
        ("status", str(status_path)),
        ("request", str(request_path)),
        ("log", str(log_path)),
    ]
    for label, output_dir in (run_output_dirs or {}).items():
        rows.append((f"{label}_output", str(output_dir)))
    print_table_rows(rows)
    print()


def poll_tuner(job: TunerJob, poll_interval: float) -> tuple[int, str | None]:
    stop_requested = False
    last_line = None
    returncode = None
    while True:
        try:
            job.heartbeat()
            status = job.status({})
            line = render_status_line(status)
            if line != last_line:
                print(line, flush=True)
                last_line = line

            returncode = job.poll()
            if status.get("state") in TERMINAL_STATES and returncode is not None:
                return int(returncode), last_line
            if returncode is not None:
                return int(returncode), last_line
            time.sleep(max(0.1, float(poll_interval)))
        except KeyboardInterrupt:
            if not stop_requested:
                print("\nCtrl-C received: requesting graceful tuner stop.", flush=True)
                print("Waiting for active evaluations to finish. Press Ctrl-C again to terminate immediately.", flush=True)
                job.request_stop()
                stop_requested = True
                continue

            print("\nSecond Ctrl-C received: terminating tuner process.", flush=True)
            job.terminate(timeout=5)
            return 130, last_line


def brief_error_message(message: str) -> str:
    """Return the first useful line from a possibly long traceback message."""
    for line in str(message or "").splitlines():
        stripped = line.strip()
        if stripped:
            return stripped
    return ""


def read_default_tunable_params() -> list[tuple[str, float]]:
    params = []
    for raw_line in DEFAULT_TUNABLE_PARAMS.read_text(encoding="utf-8").splitlines():
        line = raw_line.split("!", 1)[0].strip()
        if not line or line.startswith("&") or line.startswith("/"):
            continue
        if "=" not in line:
            continue
        name, value_text = line.split("=", 1)
        name = name.strip()
        value_text = value_text.split(",", 1)[0].strip()
        try:
            value = float(value_text.replace("D", "E").replace("d", "e"))
        except ValueError:
            continue
        params.append((name, value))
    return params


def write_top_params_file(results: list[dict], top_n: int, output_dir: Path) -> Path:
    selected_results = list(results[:top_n])
    if not selected_results:
        raise RuntimeError("No top results are available")
    # all_params contains the full Fortran API vector; only selected_params
    # should become explicit overrides in the rerun namelist.
    param_sets = [dict(result.get("selected_params") or {}) for result in selected_results]
    defaults = read_default_tunable_params()
    default_names = {name for name, _value in defaults}
    unknown = sorted({name for param_set in param_sets for name in param_set} - default_names)
    if unknown:
        raise RuntimeError("Top result contains unknown tunable parameter(s): " + ", ".join(unknown))

    output_dir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    params_path = output_dir / f"tuner_top_{timestamp}_{os.getpid()}_params.in"
    with open(params_path, "w", encoding="utf-8") as dst:
        if len(param_sets) > 1:
            dst.write("&multicol_def\n")
            dst.write(f"ngrdcol = {len(param_sets)}\n")
            dst.write(f"batch_size = {len(param_sets)}\n")
            dst.write("/\n\n")
        dst.write("&clubb_params_nl\n")
        for name, default_value in defaults:
            values = [
                float(param_set.get(name, default_value))
                for param_set in param_sets
            ]
            dst.write(f"{name} = {', '.join(f'{value:.15g}' for value in values)}\n")
        dst.write("/\n")
    return params_path


def run_command(command: list[str], *, cwd: Path) -> int:
    print("+ " + " ".join(command), flush=True)
    process = subprocess.Popen(
        command,
        cwd=str(cwd),
        env=tuner_worker_env(),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        errors="replace",
        bufsize=1,
    )
    assert process.stdout is not None
    for line in process.stdout:
        print(line, end="", flush=True)
    process.wait()
    return int(process.returncode)


def top_run_output_dirs(args: argparse.Namespace, mode: str) -> dict[str, Path]:
    base_output_dir = Path(args.run_out_dir).resolve()
    if mode == "both":
        return {
            "window": base_output_dir / "window",
            "complete": base_output_dir / "complete",
        }
    return {mode: base_output_dir}


def run_top_results(args: argparse.Namespace, request_path: Path, results_path: Path, mode: str) -> int:
    results_payload = read_json_or_default(results_path, {})
    best_results = list(results_payload.get("best_results") or [])
    if not best_results:
        print("No top results are available to run.", flush=True)
        return 0

    top_n = min(int(args.top_n), len(best_results))
    base_output_dir = Path(args.run_out_dir).resolve()
    output_dirs = top_run_output_dirs(args, mode)
    params_path = write_top_params_file(best_results, top_n, base_output_dir)
    request = results_payload.get("request") or {}
    cases = list(request.get("cases") or args.request_cases)
    fields = list(request.get("selected_fields") or args.request_fields)

    print(f"Running top {top_n} result(s) with params: {params_path}", flush=True)
    exit_code = 0
    if mode in {"window", "both"}:
        command = [
            sys.executable,
            str(RUN_SCM_LOSS),
            "-out_dir",
            str(output_dirs["window"]),
            "-fields",
            ",".join(fields),
            "-cases",
            ",".join(cases),
            "-params",
            str(params_path),
            "-case_config_file",
            str(request_path),
        ]
        exit_code = run_command(command, cwd=CLUBB_ROOT) or exit_code

    if mode in {"complete", "both"}:
        for case_name in cases:
            command = [
                sys.executable,
                str(RUN_SCM),
                "-out_dir",
                str(output_dirs["complete"]),
                "-params",
                str(params_path),
                str(case_name),
            ]
            exit_code = run_command(command, cwd=CLUBB_ROOT) or exit_code

    return exit_code


def choose_top_run_mode(args: argparse.Namespace, results_path: Path) -> str | None:
    if args.run_top == "never":
        return None
    if args.run_top in {"complete", "window", "both"}:
        return args.run_top

    results_payload = read_json_or_default(results_path, {})
    best_results = list(results_payload.get("best_results") or [])
    if not best_results:
        print("No top results are available to run.", flush=True)
        return None
    top = best_results[0]
    print(
        "Best result: "
        f"rank 1 | total_loss {float(top.get('total_loss', 0.0)):.6g} | "
        f"{format_params(top.get('selected_params', {}), limit=8)}",
        flush=True,
    )

    if not sys.stdin.isatty():
        print("Skipping top-result run because stdin is not interactive.", flush=True)
        return None

    answer = input("Run top result now? [y/N] ").strip().lower()
    if answer not in {"y", "yes"}:
        return None
    mode = input("Run mode [complete/window/both] (complete): ").strip().lower() or "complete"
    if mode not in {"complete", "window", "both"}:
        print(f"Unknown run mode '{mode}', skipping top-result run.", flush=True)
        return None
    return mode


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the CLUBB Dash-style tuner from the command line.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python run_scripts/run_tuner_job.py -cases bomex -fields cloud_frac "
            "-params C8:0.2:0.8 -strategy random:8\n"
            "  python run_scripts/run_tuner_job.py -cases arm:10800:21600:10800 bomex:7200:18000:2700 "
            "-fields cloud_frac rcm -params C8:0.2:0.8 C11:0.1:1.0 -strategy resolve:0.1\n"
        ),
    )
    parser.add_argument(
        "-cases",
        nargs="+",
        required=True,
        help="Case specs: CASE, CASE:T_START:T_END, or CASE:T_START:T_END:T_INTERVAL.",
    )
    parser.add_argument("-fields", nargs="+", required=True, help="CLUBB-facing fields, comma-separated or space-separated.")
    parser.add_argument("-params", nargs="+", required=True, help="Parameter ranges as PARAM:MIN:MAX.")
    parser.add_argument(
        "-strategy",
        default="random:8",
        help="Strategy spec: random:MAX_SAMPLES or resolve:SPACING. Default: random:8.",
    )
    parser.add_argument("-batch_size", type=int, default=8, help="Parameter batch size. Default: 8.")
    parser.add_argument("-max_workers", type=int, default=1, help="Maximum concurrent case workers. Default: 1.")
    parser.add_argument("-seed", type=int, help="Optional random seed.")
    parser.add_argument("-loss_mode", default=DEFAULT_LOSS_MODE, help=f"Loss mode. Default: {DEFAULT_LOSS_MODE}.")
    parser.add_argument(
        "-aggregation_mode",
        default=DEFAULT_AGGREGATION_MODE,
        help=f"Aggregation mode. Default: {DEFAULT_AGGREGATION_MODE}.",
    )
    parser.add_argument("-job_dir", help="Exact job directory to create. Must be empty if it exists.")
    parser.add_argument("-output_root", default=str(OUTPUT_TUNER_DIR), help="Root for generated job directories.")
    parser.add_argument("-poll_interval", type=float, default=1.0, help="Seconds between status polls. Default: 1.")
    parser.add_argument(
        "-run_top",
        choices=("ask", "never", "complete", "window", "both"),
        default="ask",
        help="What to do with top results after the tuner exits. Default: ask.",
    )
    parser.add_argument("-top_n", type=int, default=1, help="Number of top parameter sets to run. Default: 1.")
    parser.add_argument("-run_out_dir", default=str(CLUBB_OUTPUT_DIR), help="Output directory for post-tuning runs.")
    parser.add_argument("-dry_run", action="store_true", help="Write request/control/status files and exit without launching.")
    args = parser.parse_args(argv)

    if args.batch_size < 1:
        parser.error("-batch_size must be >= 1")
    if args.max_workers < 1:
        parser.error("-max_workers must be >= 1")
    if args.top_n < 1:
        parser.error("-top_n must be >= 1")
    if args.poll_interval <= 0.0:
        parser.error("-poll_interval must be > 0")
    return args


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        request = build_request(args)
        args.request_cases = list(request["cases"])
        args.request_fields = list(request["selected_fields"])
        job = TunerJob.create(
            request,
            output_root=args.output_root,
            job_dir=args.job_dir,
            prefix="cli",
        )
        print_request_summary(request, job.job_dir, job.request_path, job.log_path, args)
        if args.dry_run:
            print("Dry run requested; tuner was not launched.")
            print_output_locations(job.job_dir, job.request_path, job.status_path, job.results_path, job.log_path)
            return 0

        job.start()
        print("Launching tuner. Press Ctrl-C once for graceful stop, twice to terminate.", flush=True)
        tuner_exit, last_status_line = poll_tuner(job, args.poll_interval)
        final_status = job.status({})
        final_line = render_status_line(final_status)
        if final_line != last_status_line:
            print(final_line, flush=True)
        if final_status.get("error_message"):
            print(f"error: {brief_error_message(final_status['error_message'])}", flush=True)
            print(f"See full log: {job.log_path}", flush=True)

        top_mode = choose_top_run_mode(args, job.results_path)
        top_exit = 0
        if top_mode is not None:
            top_exit = run_top_results(args, job.request_path, job.results_path, top_mode)
        run_output_dirs = top_run_output_dirs(args, top_mode) if top_mode is not None else None
        print_output_locations(
            job.job_dir,
            job.request_path,
            job.status_path,
            job.results_path,
            job.log_path,
            run_output_dirs=run_output_dirs,
        )
        return int(tuner_exit or top_exit)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
