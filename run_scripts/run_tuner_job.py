#!/usr/bin/env python3
"""Command-line front door for the CLUBB tuner.

Usage
-----
The smallest useful tuner run names a case, one or more loss fields, and one or
more parameter ranges:

    python run_scripts/run_tuner_job.py \
      -cases bomex \
      -fields cloud_frac \
      -param_ranges C8:0.2:0.8

Cases can be plain names, fixed windows, or split windows:

    bomex
    bomex:7200:18000
    bomex:7200:18000:2700
    case:t_start:t_end:t_out

Parameters use PARAM:MIN:MAX, or PARAM=PARAM:MIN:MAX to link equal physical
parameters. Strategies include random, resolve, SimAnn, and Adam with SPSA
gradients. Interactive runs ask whether to rerun the best result
after tuning; Jenkins and other non-interactive callers should set -run_top to
never, complete, window, or both.

Understanding
-------------
This script is a controller, not a loss calculator.

It creates a TunerJob directory under output/tuner and writes the small files
that describe the job:

  - request.json: cases, fields, parameters, strategy, and loss settings
  - control.json: stop requests and the controller keepalive lease
  - status.json: live scheduler progress
  - results.json: completed samples and rankings

The actual work is done by the scheduler launched as:

    python -m tuner.tune_clubb -job_dir <job_dir>

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
import json
from datetime import datetime
import os
from pathlib import Path
import subprocess
import sys
import time


RUN_SCRIPTS = Path(__file__).resolve().parent
CLUBB_ROOT = RUN_SCRIPTS.parent
OUTPUT_TUNER_DIR = CLUBB_ROOT / "output" / "tuner"
CLUBB_OUTPUT_DIR = CLUBB_ROOT / "output"
RUN_SCM = RUN_SCRIPTS / "run_scm.py"
RUN_SCM_LOSS = RUN_SCRIPTS / "run_scm_loss.py"

if str(CLUBB_ROOT) not in sys.path:
    sys.path.insert(0, str(CLUBB_ROOT))

if __name__ == "__main__":
    from utilities.setup_python_venv import ensure_python_venv

    # Prepare dependencies before imports; follow-up runs inherit this interpreter.
    ensure_python_venv()

from clubb_jax.run_jax import runtime_arguments  # noqa: E402
from tuner.job_runtime import TERMINAL_STATES, TunerJob, tuner_worker_env  # noqa: E402
from tuner.presets import apply_preset, list_presets  # noqa: E402
from tuner.status import read_json_or_default  # noqa: E402
from tuner.system_defaults import default_max_workers  # noqa: E402
from tuner.taylor_metrics import DEFAULT_AGGREGATION_MODE, DEFAULT_LOSS_WEIGHTS  # noqa: E402
from utilities.loss_metrics import normalize_loss_weights  # noqa: E402
from input.scripts.create_case_namelist import (  # noqa: E402
    filter_tuned_overrides,
)


DEFAULT_MAX_WORKERS = default_max_workers()


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
    """Parse PARAM:MIN:MAX or linked PARAM=PARAM:MIN:MAX."""
    parts = [part.strip() for part in str(spec).split(":")]
    if len(parts) != 3 or any(part == "" for part in parts):
        raise argparse.ArgumentTypeError(f"Invalid param spec '{spec}'. Use PARAM:MIN:MAX or PARAM=PARAM:MIN:MAX.")
    targets = [item.strip() for item in parts[0].split("=") if item.strip()]
    if not targets:
        raise argparse.ArgumentTypeError(f"Invalid parameter target list '{parts[0]}'.")
    if len(targets) != len(set(targets)):
        raise argparse.ArgumentTypeError(f"Parameter target list '{parts[0]}' repeats a name.")
    min_value = _float_from_text(parts[1], f"{parts[0]} min")
    max_value = _float_from_text(parts[2], f"{parts[0]} max")
    if min_value > max_value:
        raise argparse.ArgumentTypeError(f"{parts[0]} requires min <= max")
    return {"name": targets[0], "targets": targets, "min": min_value, "max": max_value}


def parse_strategy_spec(spec: str) -> dict:
    """Parse one compact tuning-strategy specification."""
    from tuner.algorithm_registry import parse_strategy
    try:
        return parse_strategy(spec)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc


def build_request(args: argparse.Namespace) -> dict:
    """Build a canonical strategy request, including required gradients and layout."""
    case_values = _split_values(args.cases)
    field_values = _split_values(args.fields)
    param_values = _split_values(args.params)
    if not args.preset and not case_values:
        raise argparse.ArgumentTypeError("-cases must include at least one case unless -preset is used")
    if not args.preset and not field_values:
        raise argparse.ArgumentTypeError("-fields must include at least one field unless -preset is used")
    if not args.preset and not param_values:
        raise argparse.ArgumentTypeError("-param_ranges must include at least one PARAM:MIN:MAX range unless -preset is used")

    case_configs = [parse_case_spec(spec) for spec in case_values]
    selected_fields = field_values
    parameter_ranges = [parse_param_spec(spec) for spec in param_values]

    strategy = parse_strategy_spec(args.strategy)
    batch_size = int(args.batch_size)
    max_workers = int(args.max_workers)
    request = {
        "backend": "jax" if getattr(args, "jax", False) else "fortran",
        "grad": getattr(args, "grad", False),
        "gradient_mode": getattr(args, 'gradient_mode', 'auto'),
        "batch_size": batch_size,
        "max_workers": max_workers,
        "strategy": strategy,
        "loss_weights": normalize_loss_weights(json.loads(args.loss_weights)),
        "aggregation_mode": args.aggregation_mode,
        "aggregation_weights": list(args.aggregation_weights),
    }
    if args.preset:
        request["preset"] = args.preset
    if getattr(args, "jax", False):
        request["jax_options"] = args.jax_options or "cpu"
    if case_values:
        request["cases"] = [config["case_name"] for config in case_configs]
        request["case_configs"] = case_configs
    if field_values:
        request["selected_fields"] = selected_fields
    if param_values:
        request["parameter_ranges"] = parameter_ranges
    if args.config is not None:
        request["config"] = str(args.config).strip() or "default"
    if args.override is not None:
        request["override"] = args.override
    if args.seed is not None:
        request["seed"] = int(args.seed)
    from tuner.algorithm_registry import prepare_strategy
    request = apply_preset(request)
    prepare_strategy(request)
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
    from tuner.algorithm_registry import get_algorithm
    spec = get_algorithm(strategy["name"])
    options = strategy.get("options") or {}
    return spec.label + ", " + ", ".join(f"{key}={value}" for key, value in options.items())


def estimate_total_samples(request: dict) -> str:
    """Return the sample count implied by the CLI request, if known."""
    from tuner.algorithm_registry import prepare_strategy
    candidate = dict(request)
    candidate["strategy"] = {**request["strategy"], "options": dict(request["strategy"].get("options", {}))}
    try:
        prepare_strategy(candidate)
    except ValueError:
        return "invalid layout"
    count = candidate["total_samples"]
    return str(count) if count is not None else "unknown"


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
    baseline_parts = []
    for name, label in (("clubb_default", "default"), ("override_defaults", "override")):
        baseline = (status.get("baselines") or {}).get(name) or {}
        value = baseline.get("total_loss")
        if value is not None:
            baseline_parts.append(f"{label} {float(value):.6g}")
        elif baseline.get("available") is False:
            baseline_parts.append(f"{label} unavailable")
    baseline_text = "" if not baseline_parts else " | baselines: " + ", ".join(baseline_parts)
    best = status.get("best_total_loss")
    if best is None:
        best_text = "best: --"
    else:
        params = {}
        top_results = status.get("top_results") or []
        if top_results:
            params = top_results[0].get("params", {})
        best_text = f"best: {float(best):.6g}"
        if top_results:
            for key, label in (
                ("improvement_vs_clubb_default_percent", "default"),
                ("improvement_vs_override_defaults_percent", "override"),
            ):
                value = top_results[0].get(key)
                if value is not None:
                    best_text += f", vs {label} {float(value):+.2f}%"
        if params:
            best_text += f" ({format_params(params)})"
    return (
        f"state: {state} | samples: {sample_text} | elapsed: {elapsed} | "
        f"workers active/idle/init: {active}/{idle}/{initialized} | queued cases: {queued}"
        f"{baseline_text} | {best_text}"
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
            ("estimated_candidates", estimate_total_samples(request)),
            (
                "baseline_evaluations",
                str(len(request["cases"]) * (2 if str(request.get("override") or "").strip() else 1)),
            ),
            ("batch_size", str(request["batch_size"])),
            ("backend", request.get("backend", "fortran")),
            ("max_workers", str(request["max_workers"])),
            ("config", str(request.get("config") or "default")),
            ("loss", request["loss_weights"]),
            ("aggregation", request["aggregation_mode"]),
            ("aggregation_weights", ", ".join(f"{value:g}" for value in request.get("aggregation_weights", []))),
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


def read_default_tunable_params(config: str = "default") -> list[tuple[str, float]]:
    from tuner.parameter_ranges import tunable_parameter_defaults
    return list(tunable_parameter_defaults(config).items())


def write_top_params_file(results: list[dict], top_n: int, output_dir: Path, config: str = "default") -> Path:
    selected_results = list(results[:top_n])
    if not selected_results:
        raise RuntimeError("No top results are available")
    # all_params contains the full Fortran API vector; only selected_params
    # should become explicit overrides in the rerun namelist.
    param_sets = [dict(result.get("selected_params") or {}) for result in selected_results]
    defaults = read_default_tunable_params(config)
    default_names = {name for name, _value in defaults}
    unknown = sorted({name for param_set in param_sets for name in param_set} - default_names)
    if unknown:
        raise RuntimeError("Top result contains unknown tunable parameter(s): " + ", ".join(unknown))

    output_dir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    from input.scripts.parameters_json_to_namelist import parameter_settings
    settings = parameter_settings(config=config)
    params_path = output_dir / f"tuner_top_{timestamp}_{os.getpid()}_params.json"
    settings["model"] = {
        name: [float(param_set.get(name, default_value)) for param_set in param_sets]
        if len(param_sets) > 1 else float(param_sets[0].get(name, default_value))
        for name, default_value in defaults
    }
    if len(param_sets) > 1:
        settings["multicol"] = {"ngrdcol": len(param_sets), "batch_size": len(param_sets)}
    params_path.write_text(json.dumps(settings, indent=2, allow_nan=False) + "\n", encoding="utf-8")
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
    request = results_payload.get("request") or {}
    cases = list(request.get("cases") or args.request_cases)
    fields = list(request.get("selected_fields") or args.request_fields)
    config = str(request.get("config") or args.config or "default").strip() or "default"
    params_path = write_top_params_file(best_results, top_n, base_output_dir, config=config)

    # Candidate vectors replace initialized values for tuned parameters. Keep
    # that precedence in reruns while retaining the other request overrides.
    rerun_override = filter_tuned_overrides(
        request.get("override"), cases,
        {name for result in best_results[:top_n] for name in result.get("selected_params", {})},
    )

    print(f"Running top {top_n} result(s) with params: {params_path}", flush=True)
    exit_code = 0
    if mode in {"window", "both"}:
        command = [
            sys.executable,
            str(RUN_SCM_LOSS),
            "-output_dir",
            str(output_dirs["window"]),
            "-config",
            config,
            "-fields",
            ",".join(fields),
            "-cases",
            ",".join(cases),
            "-params_file",
            str(params_path),
            "-case_config_file",
            str(request_path),
        ]
        if request.get("backend") == "jax":
            command.extend(runtime_arguments(
                request.get("jax_options", "cpu"),
                device=request.get("jax_gpu", ""),
                prealloc_gpu_mem=request.get("jax_xla_prealloc"), scm=True,
            ))
        if rerun_override:
            command.extend(["-override", rerun_override])
        exit_code = run_command(command, cwd=CLUBB_ROOT) or exit_code

    if mode in {"complete", "both"}:
        for case_name in cases:
            command = [
                sys.executable,
                str(RUN_SCM),
                "-output_dir",
                str(output_dirs["complete"]),
                "-config",
                config,
                "-params_file",
                str(params_path),
                str(case_name),
            ]
            if request.get("backend") == "jax":
                command.extend(runtime_arguments(
                    request.get("jax_options", "cpu"),
                    device=request.get("jax_gpu", ""),
                    prealloc_gpu_mem=request.get("jax_xla_prealloc"), scm=True,
                ))
            if rerun_override:
                command.extend(["-override", rerun_override])
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
    for key, label in (
        ("improvement_vs_clubb_default_percent", "Improvement vs CLUBB default"),
        ("improvement_vs_override_defaults_percent", "Improvement vs override defaults"),
    ):
        if top.get(key) is not None:
            print(f"{label}: {float(top[key]):+.2f}%", flush=True)

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
    """Parse tuner job, backend and post-tuning run options from argv or the command line."""
    parser = argparse.ArgumentParser(
        description="Run the CLUBB Dash-style tuner from the command line.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python run_scripts/run_tuner_job.py -preset wpxp -strategy random:8\n"
            "  python run_scripts/run_tuner_job.py -cases bomex -fields cloud_frac "
            "-param_ranges C8:0.2:0.8 -strategy random:8\n"
            "  python run_scripts/run_tuner_job.py -cases arm:10800:21600:10800 bomex:7200:18000:2700 "
            "-fields cloud_frac rcm -param_ranges C8:0.2:0.8 C11:0.1:1.0 -strategy resolve:0.1\n"
            "  python run_scripts/run_tuner_job.py -cases bomex -fields cloud_frac "
            "-param_ranges C8:0.2:0.8 -strategy simann:2000:1.0:1e-12\n"
            "  python run_scripts/run_tuner_job.py -cases bomex -fields cloud_frac "
            "-param_ranges C8:0.2:0.8 -strategy adam:100:0.01:0.05:2\n"
        ),
        add_help=False, allow_abbrev=False
    )
    parser.add_argument("-h", "-help", action="help", help="Show this help and exit.")
    parser.add_argument(
        "-cases",
        nargs="+",
        help="Case specs: CASE, CASE:T_START:T_END, or CASE:T_START:T_END:T_INTERVAL.",
    )
    parser.add_argument(
        "-jax", action="store_true",
        help="Use the JAX loss backend; -jax=cpu or -jax=gpu selects its runtime.",
    )
    parser.add_argument("-grad", action="store_true", help="Compute and store case-loss gradients for selected parameters (JAX only).")
    from utilities.clubb_loss_runner import GRADIENT_MODES
    parser.add_argument('-gradient_mode', choices=GRADIENT_MODES, default='auto',
                        help='Gradient mode: auto uses reverse for more than one parameter targeting the case loss, otherwise forward.')
    parser.add_argument("-fields", nargs="+", help="CLUBB-facing fields, comma-separated or space-separated.")
    parser.add_argument('-param_ranges', dest='params', nargs="+", help="Ranges as PARAM:MIN:MAX or linked PARAM=PARAM:MIN:MAX.", metavar='RANGE')
    parser.add_argument('-preset', dest='preset', help="Named tuner preset; explicit -cases/-fields/-param_ranges override its pieces.")
    parser.add_argument('-list_presets', dest='list_presets', action="store_true", help="List available Tune presets and exit.")
    parser.add_argument("-list_algorithms", action="store_true", help="List algorithm requirements and typed options as JSON, then exit.")
    parser.add_argument(
        "-config",
        default=None,
        help="Shared JSON preset name. Default: default.",
    )
    parser.add_argument("-override", default=None, help="Optional comma-separated CLUBB key=value overrides.")
    parser.add_argument(
        "-strategy",
        default="random",
        help=(
            "Algorithm name followed by colon-separated positional or NAME=VALUE options. "
            "Use -list_algorithms for types, defaults and requirements. "
            "Example: adam_gradient:max_updates=10:learning_rate=0.01 (requires -jax). Default: random."
        ),
    )
    parser.add_argument("-batch_size", type=int, default=8, help="Parameter samples per CLUBB column batch. Default: 8.")
    parser.add_argument(
        '-workers', dest='max_workers',
        type=int,
        default=DEFAULT_MAX_WORKERS,
        help=f"Maximum active case evaluations (default: half the available logical CPUs, {DEFAULT_MAX_WORKERS}).",
        metavar='N',
    )
    parser.add_argument("-seed", type=int, help="Optional random seed.")
    parser.add_argument("-loss_weights", default=json.dumps(DEFAULT_LOSS_WEIGHTS),
                        help='Positive Taylor contribution weights as JSON, e.g. {"mean":1,"std":2,"corr":1}.')
    parser.add_argument(
        "-aggregation_mode",
        default=DEFAULT_AGGREGATION_MODE,
        help=f"Aggregation mode. Default: {DEFAULT_AGGREGATION_MODE}.",
    )
    parser.add_argument(
        "-aggregation_weights",
        nargs=4,
        type=float,
        default=(0.1, 0.4, 0.4, 0.1),
        metavar=("BEST", "LOWER_MID", "UPPER_MID", "WORST"),
        help="Best-to-worst quantile weights for quantile_weighted aggregation. Default: 0.1 0.4 0.4 0.1.",
    )
    parser.add_argument('-output_job_dir', dest='job_dir', help="Exact job directory to create. Must be empty if it exists.", metavar='DIR')
    parser.add_argument("-output_root", default=str(OUTPUT_TUNER_DIR), help="Root for generated job directories.", metavar='DIR')
    parser.add_argument("-poll_interval", type=float, default=1.0, help="Seconds between status polls. Default: 1.")
    parser.add_argument(
        "-run_top",
        choices=("ask", "never", "complete", "window", "both"),
        default="ask",
        help="What to do with top results after the tuner exits. Default: ask.",
    )
    parser.add_argument("-top_n", type=int, default=1, help="Number of top parameter sets to run. Default: 1.")
    parser.add_argument('-output_run_dir', dest='run_out_dir', default=str(CLUBB_OUTPUT_DIR), help="Output directory for post-tuning runs.", metavar='DIR')
    parser.add_argument("-dry_run", action="store_true", help="Write request/control/status files and exit without launching.")
    from utilities.clubb_case_runner import extract_jax_options

    normalized, jax_options, occurrences = extract_jax_options(sys.argv[1:] if argv is None else argv)
    if occurrences > 1:
        parser.error("-jax may be specified only once")
    args = parser.parse_args(normalized)
    args.jax_options = jax_options
    if args.batch_size < 1:
        parser.error("-batch_size must be >= 1")
    if args.max_workers < 1:
        parser.error("-workers must be >= 1")
    if args.top_n < 1:
        parser.error("-top_n must be >= 1")
    if args.poll_interval <= 0.0:
        parser.error("-poll_interval must be > 0")
    return args


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if args.list_algorithms:
        from tuner.algorithm_registry import list_algorithms
        print(json.dumps(list_algorithms(), indent=2))
        return 0
    if args.list_presets:
        for preset in list_presets():
            print(f"{preset['name']}: {preset['description']}")
        return 0
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
