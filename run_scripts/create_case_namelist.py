#!/usr/bin/env python3
"""Build the single namelist file that CLUBB actually reads.

Usage
-----
CLUBB case settings live in several input files: tunable parameters,
configurable flags, SILHS settings, the case model file, and a stats list. This
script combines those fragments into one generated namelist:

    python run_scripts/create_case_namelist.py [options] CASE

The command prints the generated file path. By default that file is
output/CASE.in, but callers can choose another directory with -out_dir.

Common options choose alternate input fragments, generate a multicolumn
parameter file, change dt_main/dt_rad/time_final, set the stats output window,
select a grid file, or apply raw key=value namelist overrides.

Most runs reach this code indirectly:

  - run_scm.py forwards its namelist-related options here before launching
    clubb_standalone.
  - run_scm_loss.py and tuner workers call create_loss_case_namelist() to build
    the specialized input used by the in-memory loss driver.

Understanding
-------------
Normal case setup:
  create_case_namelist_file() is the standard path. It finds the case model
  file and the requested input fragments, optionally creates a temporary
  multicolumn parameter file, strips comments, concatenates the namelists,
  applies requested edits, writes stats_setting%output_dir, and saves CASE.in.

Loss-driver setup:
  create_loss_case_namelist() starts with the same normal CASE.in so the loss
  path does not grow a separate version of case setup. It then appends
  &tuner_loss_nl, converts the LES benchmark file to CLUBB-facing variable
  names, trims the stats list to the requested loss fields, and sets the model
  and stats windows so the in-memory loss calculation and optional NetCDF output
  describe the same samples.

The important design point is that this file owns the translation from
high-level run options into namelist text. Keeping that logic here prevents the
normal scripts, loss scripts, and tuner workers from each developing a slightly
different idea of how parameters, flags, stats, multicolumn settings, and loss
windows should be assembled.
"""
import argparse
import os
import re
import subprocess
import sys
from pathlib import Path

RUN_SCRIPTS = os.path.dirname(os.path.abspath(__file__))
CLUBB_ROOT = os.path.join(RUN_SCRIPTS, "..")
DEFAULT_OUTPUT_DIR = os.path.join(CLUBB_ROOT, "output")
DEFAULT_STANDARD_STATS = os.path.join(CLUBB_ROOT, "input/stats/standard_stats.in")
DEFAULT_TUNER_STATS = os.path.join(CLUBB_ROOT, "input/stats/tuning_stats.in")
DEFAULT_TUNABLE_PARAMS = os.path.join(CLUBB_ROOT, "input/tunable_parameters/tunable_parameters.in")
MULTI_COL_PARAMS_SCRIPT = os.path.join(RUN_SCRIPTS, "create_multi_col_params.py")

HR_SPEC_RE = re.compile(
    r"^[A-Za-z_]\w*/[-+]?\d*\.?\d+(?:[eEdD][-+]?\d+)?:[-+]?\d*\.?\d+(?:[eEdD][-+]?\d+)?/\d+"
    r"(?:,[A-Za-z_]\w*/[-+]?\d*\.?\d+(?:[eEdD][-+]?\d+)?:[-+]?\d*\.?\d+(?:[eEdD][-+]?\d+)?/\d+)*$"
)


def strip_comments_and_remove_keys(content: str, keys_to_remove=None) -> str:
    if keys_to_remove is None:
        keys_to_remove = []
    lines = []
    for line in content.splitlines():
        line = re.sub(r"!.*", "", line)
        if any(key in line for key in keys_to_remove):
            continue
        lines.append(line)
    return "\n".join(lines)


def validate_multicol(value: str) -> str:
    stripped = value.strip()
    if not stripped:
        raise argparse.ArgumentTypeError("must provide a non-empty -multicol value")

    try:
        ngrdcol = int(stripped)
    except ValueError:
        if HR_SPEC_RE.fullmatch(stripped):
            return stripped
        raise argparse.ArgumentTypeError(
            "must be either an integer column count or an hr spec like C8/0.2:0.8/4"
        )

    if ngrdcol < 1:
        raise argparse.ArgumentTypeError("integer -multicol value must be >= 1")

    return stripped


def convert_to_multi_col(
    params_file: str,
    case_name: str,
    output_dir: str,
    multicol_spec: str,
    batch_size: int | None = None,
    *,
    integer_mode: str = "dup_tweak",
    out_file: str | None = None,
) -> str:
    if not os.path.isfile(MULTI_COL_PARAMS_SCRIPT):
        raise RuntimeError(f"Missing helper script: {MULTI_COL_PARAMS_SCRIPT}")

    if out_file is None:
        out_file = os.path.join(output_dir, f"{case_name}_multicol_params.in")
    cmd = [sys.executable, MULTI_COL_PARAMS_SCRIPT, "-param_file", params_file, "-out_file", out_file]

    try:
        ngrdcol = int(multicol_spec)
    except ValueError:
        cmd.extend(["-hr", multicol_spec])
    else:
        cmd.extend(["-mode", integer_mode, "-n", str(ngrdcol)])

    if batch_size is not None:
        cmd.extend(["-batch_size", str(batch_size)])

    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            f"create_multi_col_params failed (exit {exc.returncode}). Command was:\n  {' '.join(cmd)}"
        ) from exc

    if not os.path.isfile(out_file):
        raise RuntimeError(f"Expected output file not created: {out_file}")

    return out_file


def read_model_times(model_file):
    """Read time_initial, time_final, dt_main, etc. from a model file if present."""
    values = {}
    line_re = re.compile(r"^\s*([a-zA-Z_]\w*)\s*=\s*([-+0-9.eEdD]+)")
    with open(model_file, encoding="utf-8") as f:
        for line in f:
            m = line_re.match(line)
            if m:
                key, val = m.groups()
                val = val.replace("D", "E").replace("d", "e")
                try:
                    values[key.lower()] = float(val)
                except ValueError:
                    pass
    return values


def parse_override_pairs(override_string):
    """Parse KEY=value overrides while allowing comma-separated value lists."""
    pairs = []
    current_key = None
    current_values = []

    for raw_item in (override_string or "").split(","):
        item = raw_item.strip()
        if not item:
            continue
        if "=" in item:
            if current_key is not None:
                pairs.append((current_key, ", ".join(current_values)))
            current_key, value = item.split("=", 1)
            current_key = current_key.strip()
            current_values = [value.strip()]
        elif current_key is not None:
            current_values.append(item)

    if current_key is not None:
        pairs.append((current_key, ", ".join(current_values)))
    return [(key, value) for key, value in pairs if key]


def override_value(override_string, clubb_in_text):
    """
    Apply overrides from -override KEY1=val1,KEY2=val2,... to the aggregate text.
    Values may also be comma-separated column lists, e.g. C8=0.8,0.7,C11=1.0,1.1.
    """
    for key, val in parse_override_pairs(override_string):
        replacement = f"{key} = {val}"

        assignment_re = re.compile(rf"(?im)^(\s*){re.escape(key)}\s*=.*$")

        clubb_in_text, replacements = assignment_re.subn(
            lambda match: f"{match.group(1)}{key} = {val}",
            clubb_in_text,
        )

        if replacements == 0:
            clubb_in_text += f"\n{replacement}\n"
    return clubb_in_text


def set_namelist_value(clubb_in: str, namelist: str, key: str, value: str, *, quote: bool = False) -> str:
    """Set a key inside a namelist, adding the key or namelist if needed."""
    value_text = f"'{value}'" if quote else str(value)
    namelist_match = re.search(rf"(?ims)^(\s*&\s*{re.escape(namelist)}\b)(.*?)(^\s*/\s*$)", clubb_in)
    if not namelist_match:
        clubb_in += f"\n&{namelist}\n{key} = {value_text},\n/\n"
        return clubb_in

    header = namelist_match.group(1)
    body = namelist_match.group(2)
    end = namelist_match.group(3)
    if re.search(rf"(?im)^\s*{re.escape(key)}\s*=", body):
        body = re.sub(
            rf"(?im)^(\s*){re.escape(key)}\s*=.*$",
            lambda match: f"{match.group(1)}{key} = {value_text},",
            body,
        )
    else:
        body = body.rstrip() + f"\n{key} = {value_text},\n"

    return clubb_in[:namelist_match.start()] + header + body + end + clubb_in[namelist_match.end():]


def set_stats_output_dir(clubb_in: str, output_dir: str) -> str:
    """Ensure &stats_setting contains output_dir."""
    out_norm = output_dir.replace("\\", "/")
    return set_namelist_value(clubb_in, "stats_setting", "output_dir", out_norm, quote=True)


def set_stats_string(clubb_in: str, key: str, value: str) -> str:
    return set_namelist_value(clubb_in, "stats_setting", key, value, quote=True)


def set_stats_value(clubb_in: str, key: str, value: str) -> str:
    return set_namelist_value(clubb_in, "stats_setting", key, value)


def set_model_string(clubb_in: str, key: str, value: str) -> str:
    return set_namelist_value(clubb_in, "model_setting", key, value, quote=True)


def set_model_value(clubb_in: str, key: str, value: str) -> str:
    return set_namelist_value(clubb_in, "model_setting", key, value)


def apply_namelist_overrides(args, clubb_in: str, model_file: str, output_dir: str) -> str:
    """Apply CLI namelist edits after the base files have been aggregated."""
    if args.dt_main is not None:
        clubb_in = set_model_value(clubb_in, "dt_main", f"{args.dt_main}.0")

    if args.dt_rad is not None:
        clubb_in = set_model_value(clubb_in, "dt_rad", f"{args.dt_rad}.0")

    if args.zt_grid is not None:
        clubb_in = set_model_value(clubb_in, "nzmax", str(args.nzmax))
        clubb_in = set_model_string(clubb_in, "zt_grid_fname", args.zt_grid)
        clubb_in = set_model_string(clubb_in, "zm_grid_fname", "")
        clubb_in = set_model_value(clubb_in, "grid_type", "2")

    if args.zm_grid is not None:
        clubb_in = set_model_value(clubb_in, "nzmax", str(args.nzmax))
        clubb_in = set_model_string(clubb_in, "zm_grid_fname", args.zm_grid)
        clubb_in = set_model_string(clubb_in, "zt_grid_fname", "")
        clubb_in = set_model_value(clubb_in, "grid_type", "3")

    disable_stats = ((args.stats or "").strip().lower() == "none")
    if args.tout is not None:
        if args.tout == 0:
            disable_stats = True
        else:
            clubb_in = set_stats_value(clubb_in, "stats_tout", f"{args.tout}.0")

    if args.stats_tstart is not None:
        clubb_in = set_stats_value(clubb_in, "stats_tstart", str(args.stats_tstart))

    if args.stats_tend is not None:
        clubb_in = set_stats_value(clubb_in, "stats_tend", str(args.stats_tend))

    if disable_stats:
        clubb_in = set_stats_value(clubb_in, "l_stats", ".false.")
        clubb_in = set_stats_string(clubb_in, "stats_output_filename", "")
    else:
        clubb_in = set_stats_string(clubb_in, "stats_output_filename", f"{args.case_name}_stats.nc")

    if args.debug is not None:
        clubb_in = set_model_value(clubb_in, "debug_level", str(args.debug))

    if args.max_iters is not None:
        vals = read_model_times(model_file)
        time_initial = vals.get("time_initial")
        time_final = vals.get("time_final")
        dt_main_val = args.dt_main if args.dt_main is not None else vals.get("dt_main")
        if time_initial is None or time_final is None or dt_main_val is None:
            sys.exit(f"Could not determine time_initial, time_final, and dt_main from {model_file}")
        new_time_final = time_initial + dt_main_val * args.max_iters
        if time_final >= new_time_final:
            clubb_in = set_model_value(clubb_in, "time_final", str(new_time_final))

    if args.override:
        clubb_in = override_value(args.override, clubb_in)

    return set_stats_output_dir(clubb_in, output_dir)


def _require_existing_file(opt: str, path: str | None) -> None:
    if path is None:
        return
    if not os.path.isfile(path):
        raise RuntimeError(f"Required file for {opt} not found: {path}")


def _ensure_repo_root_on_path() -> None:
    repo_root = os.path.abspath(CLUBB_ROOT)
    if repo_root not in sys.path:
        sys.path.insert(0, repo_root)


def create_case_namelist_file(
    case_name: str,
    output_dir: str | os.PathLike | None = None,
    *,
    config: str | None = None,
    params: str | None = None,
    flags: str | None = None,
    silhs_params: str | None = None,
    stats: str | None = None,
    multicol: str | None = None,
    batch_size: int | None = None,
    zt_grid: str | None = None,
    zm_grid: str | None = None,
    nzmax: int | None = None,
    debug: str | None = None,
    max_iters: int | None = None,
    dt_main: int | None = None,
    dt_rad: int | None = None,
    tout: int | None = None,
    stats_tstart: float | None = None,
    stats_tend: float | None = None,
    override: str | None = None,
) -> Path:
    """Create the aggregate case namelist used by normal CLUBB runs."""
    output_dir_abs = os.path.abspath(os.fspath(output_dir) if output_dir else DEFAULT_OUTPUT_DIR)
    os.makedirs(output_dir_abs, exist_ok=True)

    model_file = os.path.join(CLUBB_ROOT, f"input/case_setups/{case_name}_model.in")
    if not os.path.isfile(model_file):
        raise RuntimeError(f"{model_file} does not exist")

    config_dir = os.path.abspath(config) if config else os.path.join(CLUBB_ROOT, "input/tunable_parameters")
    if not os.path.isdir(config_dir):
        raise RuntimeError(f"--config directory does not exist: {config_dir}")

    params_file = params or os.path.join(config_dir, "tunable_parameters.in")
    flags_file = flags or os.path.join(config_dir, "configurable_model_flags.in")
    silhs_params_file = silhs_params or os.path.join(config_dir, "silhs_parameters.in")
    stats_arg = (stats or "").strip()
    disable_stats = stats_arg.lower() == "none"
    stats_file = None if disable_stats else (stats or DEFAULT_STANDARD_STATS)

    if batch_size is not None and multicol is None:
        raise RuntimeError("-batch_size requires -multicol")
    if zt_grid and zm_grid:
        raise RuntimeError("Cannot specify both -zt_grid and -zm_grid")
    if (zt_grid or zm_grid) and nzmax is None:
        raise RuntimeError("-nzmax is required with -zt_grid or -zm_grid")
    if nzmax and not (zt_grid or zm_grid):
        print("WARNING: Specifying -nzmax has no effect without specifying -zm_grid or -zt_grid")
    if (stats_tstart is None) != (stats_tend is None):
        raise RuntimeError("-stats_tstart and -stats_tend must be provided together")

    if multicol is not None:
        params_file = convert_to_multi_col(
            params_file,
            case_name,
            output_dir_abs,
            multicol,
            batch_size=batch_size,
        )

    _require_existing_file("-params", params_file)
    _require_existing_file("-flags", flags_file)
    _require_existing_file("-silhs_params", silhs_params_file)
    if not disable_stats:
        _require_existing_file("-stats", stats_file)

    clubb_input_namelist = os.path.join(output_dir_abs, f"{case_name}.in")
    clubb_in = ""
    files_to_aggregate = [params_file, silhs_params_file, flags_file, model_file]
    if not disable_stats:
        files_to_aggregate.append(stats_file)
    for filename in files_to_aggregate:
        with open(filename, encoding="utf-8") as src:
            clubb_in += strip_comments_and_remove_keys(src.read())
            clubb_in += "\n"

    args = argparse.Namespace(
        case_name=case_name,
        stats=stats,
        zt_grid=zt_grid,
        zm_grid=zm_grid,
        nzmax=nzmax,
        debug=debug,
        max_iters=max_iters,
        dt_main=dt_main,
        dt_rad=dt_rad,
        tout=tout,
        stats_tstart=stats_tstart,
        stats_tend=stats_tend,
        override=override,
    )
    clubb_in = apply_namelist_overrides(args, clubb_in, model_file, output_dir_abs)

    with open(clubb_input_namelist, "w", encoding="utf-8") as out:
        out.write(clubb_in)

    return Path(clubb_input_namelist)


def build_tuner_namelist(
    case_defaults: dict,
    selected_fields: list[str],
    num_time_windows: int = 1,
) -> str:
    """Build the loss-driver request namelist appended to the aggregate case file."""
    field_values = ", ".join(f'"{name}"' for name in selected_fields)
    z_min, z_max = case_defaults["altitude_comparison_range"]
    t_start, t_end = case_defaults["time_average_range"]
    return (
        "&tuner_loss_nl\n"
        f'  les_stats_file = "{case_defaults["les_stats_file"]}"\n'
        f"  clubb_var_names = {field_values}\n"
        f"  benchmark_var_name = {field_values}\n"
        f"  altitude_comparison_range = {z_min}, {z_max}\n"
        f"  time_average_range = {t_start}, {t_end}\n"
        f"  num_time_windows = {int(num_time_windows)}\n"
        "/\n"
    )


def create_normalized_benchmark_file(
    case_name: str,
    case_defaults: dict,
    selected_fields: list[str],
    job_dir: Path,
) -> Path:
    """Convert the case LES benchmark to CLUBB-facing field names."""
    _ensure_repo_root_on_path()
    from utilities.benchmark_converter import convert_benchmark_file

    job_dir.mkdir(parents=True, exist_ok=True)
    output_path = job_dir / f"{case_name}_normalized_benchmark.nc"
    status = convert_benchmark_file(
        case_defaults["les_stats_file"],
        output_path,
        source_type="auto",
        fields=selected_fields,
    )
    missing = {name: result for name, result in status.items() if result != "written"}
    if missing:
        detail = ", ".join(f"{name}: {result}" for name, result in sorted(missing.items()))
        raise RuntimeError(f"Could not normalize requested benchmark field(s) for {case_name}: {detail}")
    return output_path


def prune_clubb_stats_namelist(clubb_in: str, requested_vars: list[str]) -> str:
    """Keep only the requested CLUBB stats entries needed by the loss driver."""
    stats_match = re.search(r"(?ims)^(\s*&\s*clubb_stats_nl\b)(.*?)(^\s*/\s*$)", clubb_in)
    if not stats_match:
        raise RuntimeError("Could not find &clubb_stats_nl in aggregated namelist")

    requested_set = set(requested_vars)
    entry_matches = re.findall(r'(?im)^\s*entry\(\d+\)\s*=\s*"([^"]+)"\s*$', stats_match.group(2))
    if not entry_matches:
        raise RuntimeError("No stats entries found in &clubb_stats_nl")

    kept_entries = []
    found_names = set()
    for entry_text in entry_matches:
        field_name = entry_text.split("|", 1)[0].strip()
        if field_name in requested_set:
            kept_entries.append(entry_text)
            found_names.add(field_name)

    missing_vars = [name for name in requested_vars if name not in found_names]
    if missing_vars:
        raise RuntimeError(
            "Requested clubb_var_names missing from tuning_stats.in: " + ", ".join(missing_vars)
        )

    new_body_lines = ["\n"]
    for idx, entry_text in enumerate(kept_entries, start=1):
        new_body_lines.append(f'  entry({idx}) = "{entry_text}"\n')

    new_block = stats_match.group(1) + "".join(new_body_lines) + stats_match.group(3)
    return clubb_in[: stats_match.start()] + new_block + clubb_in[stats_match.end() :]


def _loss_timing_settings(case_name: str, model_file: Path, case_defaults: dict, num_time_windows: int) -> dict:
    vals = read_model_times(str(model_file))
    dt_main_val = vals.get("dt_main")
    dt_rad_val = vals.get("dt_rad", dt_main_val)
    time_initial_seconds = vals.get("time_initial")
    original_time_final_seconds = vals.get("time_final")
    interval_start_seconds, interval_end_seconds = case_defaults["time_average_range"]

    if num_time_windows < 1:
        raise RuntimeError("num_time_windows must be >= 1")
    if dt_main_val is None or dt_rad_val is None:
        raise RuntimeError(f"Could not determine dt_main/dt_rad from {model_file}")
    if time_initial_seconds is None or original_time_final_seconds is None:
        raise RuntimeError(f"Could not determine time_initial/time_final from {model_file}")
    if interval_start_seconds < time_initial_seconds:
        raise RuntimeError(
            f"time_average_range for {case_name} starts before "
            f"time_initial={time_initial_seconds:g} s from {model_file}"
        )
    if interval_end_seconds > original_time_final_seconds:
        raise RuntimeError(
            f"time_average_range for {case_name} ends after "
            f"time_final={original_time_final_seconds:g} s from {model_file}"
        )

    requested_window_seconds = float(interval_end_seconds - interval_start_seconds)
    stats_tsamp_seconds = float(dt_main_val)
    stats_tout_seconds = requested_window_seconds / float(num_time_windows)

    if requested_window_seconds <= 0.0:
        raise RuntimeError(f"Derived non-positive stats_tout={requested_window_seconds} s for {case_name}")
    if round(requested_window_seconds / stats_tout_seconds) * stats_tout_seconds != requested_window_seconds:
        raise RuntimeError(
            f"Requested sampling window length={requested_window_seconds} s for {case_name}, "
            f"but it is not evenly divisible by num_time_windows={num_time_windows}"
        )
    if round(requested_window_seconds / stats_tsamp_seconds) * stats_tsamp_seconds != requested_window_seconds:
        raise RuntimeError(
            f"Requested sampling window length={requested_window_seconds} s for {case_name}, "
            f"but it is not an integer multiple of stats_tsamp={stats_tsamp_seconds} s"
        )
    if round(stats_tout_seconds / stats_tsamp_seconds) * stats_tsamp_seconds != stats_tout_seconds:
        raise RuntimeError(
            f"stats_tsamp={stats_tsamp_seconds} s must evenly divide stats_tout={stats_tout_seconds} s"
        )
    if round(stats_tout_seconds / dt_main_val) * dt_main_val != stats_tout_seconds:
        raise RuntimeError(
            f"Derived stats_tout={stats_tout_seconds} s for {case_name}, "
            f"but it is not an integer multiple of dt_main={dt_main_val} s"
        )
    if not (
        round(dt_rad_val / stats_tout_seconds) * stats_tout_seconds == dt_rad_val
        or round(stats_tout_seconds / dt_rad_val) * dt_rad_val == stats_tout_seconds
    ):
        raise RuntimeError(
            f"Derived stats_tout={stats_tout_seconds} s for {case_name}, "
            f"but it is incompatible with dt_rad={dt_rad_val} s"
        )

    return {
        "interval_start_seconds": float(interval_start_seconds),
        "interval_end_seconds": float(interval_end_seconds),
        "time_initial_seconds": float(time_initial_seconds),
        "time_final_seconds": float(interval_end_seconds),
        "stats_tsamp_seconds": float(stats_tsamp_seconds),
        "stats_tout_seconds": float(stats_tout_seconds),
        "num_time_windows": int(num_time_windows),
    }


def _print_loss_timing(timing: dict) -> None:
    print(
        " - loss timing:"
        f" absolute_window=[{timing['interval_start_seconds']:g}, {timing['interval_end_seconds']:g}] s,"
        f" model_time_initial={timing['time_initial_seconds']:g} s,"
        f" model_time_final={timing['time_final_seconds']:g} s,"
        f" stats_tstart={timing['interval_start_seconds']:g} s,"
        f" stats_tend={timing['interval_end_seconds']:g} s,"
        f" stats_tsamp={timing['stats_tsamp_seconds']:g} s,"
        f" stats_tout={timing['stats_tout_seconds']:g} s,"
        f" num_time_windows={timing['num_time_windows']}"
    )


def create_loss_case_namelist(
    case_name: str,
    output_dir: str | os.PathLike,
    selected_fields: list[str],
    *,
    case_defaults: dict | None = None,
    num_time_windows: int | None = None,
    config: str | None = None,
    params: str | None = None,
    flags: str | None = None,
    silhs_params: str | None = None,
    stats: str | None = DEFAULT_TUNER_STATS,
    multicol: str | None = None,
    batch_size: int | None = None,
    duplicate_params_for_batch: bool = False,
    disable_stats_storage: bool = False,
    override: str | None = None,
    verbose: bool = False,
) -> tuple[Path, list[str], dict]:
    """Create the aggregate namelist consumed by the in-memory loss driver."""
    output_path = Path(output_dir).resolve()
    output_path.mkdir(parents=True, exist_ok=True)

    if not selected_fields:
        _ensure_repo_root_on_path()
        from tuner.case_defaults import DEFAULT_LOSS_FIELDS

        selected_fields = list(DEFAULT_LOSS_FIELDS)
    else:
        selected_fields = list(selected_fields)

    if case_defaults is None:
        _ensure_repo_root_on_path()
        from tuner.case_defaults import read_case_defaults

        case_defaults = read_case_defaults(case_name)
    else:
        case_defaults = dict(case_defaults)

    if num_time_windows is None:
        num_time_windows = int(case_defaults.get("num_time_windows", 1))
    num_time_windows = int(num_time_windows)
    case_defaults["num_time_windows"] = num_time_windows

    params_for_case = params
    multicol_for_case = multicol
    batch_size_for_case = batch_size
    if duplicate_params_for_batch:
        if batch_size is None:
            raise RuntimeError("duplicate_params_for_batch requires batch_size")
        params_for_case = convert_to_multi_col(
            params or DEFAULT_TUNABLE_PARAMS,
            case_name,
            str(output_path),
            str(batch_size),
            batch_size=batch_size,
            integer_mode="duplicate",
            out_file=str(output_path / f"{case_name}_tuning_multicol_params.in"),
        )
        multicol_for_case = None
        batch_size_for_case = None

    aggregate_path = create_case_namelist_file(
        case_name,
        output_path,
        config=config,
        params=params_for_case,
        flags=flags,
        silhs_params=silhs_params,
        stats=stats,
        multicol=multicol_for_case,
        batch_size=batch_size_for_case,
    )

    normalized_benchmark_file = create_normalized_benchmark_file(
        case_name,
        case_defaults,
        selected_fields,
        output_path,
    )
    if verbose:
        print(f" - normalized benchmark: {normalized_benchmark_file}")

    normalized_case_defaults = dict(case_defaults)
    normalized_case_defaults["les_stats_file"] = str(normalized_benchmark_file)

    clubb_in = aggregate_path.read_text(encoding="utf-8")
    clubb_in += "\n" + build_tuner_namelist(
        normalized_case_defaults,
        selected_fields,
        num_time_windows=num_time_windows,
    )

    model_file = Path(CLUBB_ROOT) / "input" / "case_setups" / f"{case_name}_model.in"
    timing = _loss_timing_settings(case_name, model_file, case_defaults, num_time_windows)
    if verbose:
        _print_loss_timing(timing)

    clubb_in = set_stats_value(clubb_in, "l_stats", ".true.")
    clubb_in = set_model_value(clubb_in, "time_final", f"{timing['time_final_seconds']:.1f}")
    clubb_in = set_stats_value(clubb_in, "stats_tstart", f"{timing['interval_start_seconds']:.1f}")
    clubb_in = set_stats_value(clubb_in, "stats_tend", f"{timing['interval_end_seconds']:.1f}")
    clubb_in = set_stats_value(clubb_in, "stats_tsamp", f"{timing['stats_tsamp_seconds']:.1f}")
    clubb_in = set_stats_value(clubb_in, "stats_tout", f"{timing['stats_tout_seconds']:.1f}")
    clubb_in = set_model_value(clubb_in, "debug_level", "-1")
    if disable_stats_storage:
        clubb_in = set_stats_string(clubb_in, "stats_output_filename", "")
    else:
        clubb_in = set_stats_string(clubb_in, "stats_output_filename", f"{case_name}_stats.nc")
    clubb_in = prune_clubb_stats_namelist(clubb_in, selected_fields)
    if override:
        clubb_in = override_value(override, clubb_in)
    clubb_in = set_stats_output_dir(clubb_in, str(output_path))

    aggregate_path.write_text(clubb_in, encoding="utf-8")
    return aggregate_path, selected_fields, case_defaults


def main():
    parser = argparse.ArgumentParser(description="Create an aggregated CLUBB case namelist.")
    parser.add_argument("-config", metavar="[DIR]",
        help=("Directory containing tunable_parameters.in, configurable_model_flags.in, "
              "and silhs_parameters.in. Defaults to input/tunable_parameters."))
    parser.add_argument("-params", metavar="[FILE]",
        help="Define the tunable parameters. Used to override params file defined by --config")
    parser.add_argument("-flags", metavar="[FILE]",
        help="Model flags file. Used to override flags file defined by --config")
    parser.add_argument("-silhs_params", metavar="[FILE]",
        help="SILHS parameters file. Used to override silhs_params file defined by --config")
    parser.add_argument("-stats", metavar="[FILE]",
        help=("Stats file defining fields to output.\n"
              "Default: input/stats/standard_stats.in.\n"
              "Use 'none' to disable stats output."))
    parser.add_argument("-out_dir", metavar="[DIR]",
        help="Output directory for the aggregated namelist.\nDefault: output")
    parser.add_argument("-multicol", metavar="[NUM|SPEC]", type=validate_multicol,
        help=("Generate a multi-column parameter file. "
              "Use an integer for dup_tweak mode, e.g. -multicol 4, or an hr spec like "
              "-multicol C8/0.2:0.8/4"))
    parser.add_argument("-batch_size", metavar="[NUM]", type=int,
        help=("Runtime batch size written to &multicol_def. "
              "Only meaningful together with -multicol."))
    parser.add_argument("-zt_grid", metavar="[FILE]",
        help="Specify a zt grid file from input/grid.\nDefault: unused")
    parser.add_argument("-zm_grid", metavar="[FILE]",
        help="Specify a zm grid file from input/grid.\nDefault: unused")
    parser.add_argument("-nzmax", metavar="[NUM]", type=int,
        help="Max number of levels. Required with -zt_grid or -zm_grid.")
    parser.add_argument("-debug", metavar="[NUM]",
        help="Debug level (0-3) that controls CLUBB's runtime checks. Default specified in model file.")
    parser.add_argument("-max_iters", metavar="[NUM]", type=int,
        help="Maximum number of iterations")
    parser.add_argument("-dt_main", metavar="[SECONDS]", type=int,
        help="Main timestep (s). Default from model file.")
    parser.add_argument("-dt_rad", metavar="[SECONDS]", type=int,
        help="Radiation timestep (s). Default from model file.")
    parser.add_argument("-tout", metavar="[SECONDS]", type=int,
        help="Stats output interval (s). Use 0 to disable. Default from model file.")
    parser.add_argument("-stats_tstart", metavar="[SECONDS]", type=float,
        help="Stats output window start time (s). Default from model file or driver.")
    parser.add_argument("-stats_tend", metavar="[SECONDS]", type=float,
        help="Stats output window end time (s). Default from model file or driver.")
    parser.add_argument(
        "-override",
        help="Comma-separated key=value pairs, e.g. -override FLAG1=true,C2=2.0,...",
    )
    parser.add_argument("case_name", help="Name of the case to aggregate")
    args = parser.parse_args()

    try:
        clubb_input_namelist = create_case_namelist_file(
            args.case_name,
            args.out_dir,
            config=args.config,
            params=args.params,
            flags=args.flags,
            silhs_params=args.silhs_params,
            stats=args.stats,
            multicol=args.multicol,
            batch_size=args.batch_size,
            zt_grid=args.zt_grid,
            zm_grid=args.zm_grid,
            nzmax=args.nzmax,
            debug=args.debug,
            max_iters=args.max_iters,
            dt_main=args.dt_main,
            dt_rad=args.dt_rad,
            tout=args.tout,
            stats_tstart=args.stats_tstart,
            stats_tend=args.stats_tend,
            override=args.override,
        )
    except RuntimeError as exc:
        sys.exit(str(exc))

    print(clubb_input_namelist)


if __name__ == "__main__":
    main()
