#!/usr/bin/env python3
"""Resolve CLUBB's JSON inputs and write the final top-level runner namelist.

Usage:
    python input/scripts/create_case_namelist.py -cases arm -config e3sm -flags custom.json
    python input/scripts/create_case_namelist.py -list parameters

    load_catalog(), normalize_input_selection(), input_patch(),
    resolve_input_settings() and available_configs()
    read defaults and presets without creating files. Case resolution retains the
    source path and benchmark identity. Override helpers apply direct CLI edits
    to typed groups; arrays replace arrays and settings require non-null values.

    create_case_settings() assembles normal-run groups in memory. Loss setup uses
    build_loss_settings(), create_normalized_benchmark_file() and timing helpers
    before writing. convert_to_multi_col() prepares JSON column arrays.

    write_case_namelist() calls the four colocated converters exactly at the
    runner boundary. create_case_namelist_file() and create_loss_case_namelist()
    provide the complete Fortran/Python/JAX run and loss inputs.

    add_input_arguments(), add_namelist_arguments(), normalize_case_selection(),
    case_namelist_command(), create_case_namelist_from_args() and
    parse_forwarded_args() keep wrapper parsers and forwarding consistent.

Understanding:
    Cases, parameters and flags resolve catalog defaults, a matching -config
    preset, then an explicit name or JSON file. A shared config need not exist
    in every section; an explicit unknown selector fails. Stats independently
    select categories, complete presets or variable-name lists. Final direct
    overrides take precedence. Fortran, Python and JAX entry points all receive a namelist
    path; preparation code reads and edits JSON throughout.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT_DIR = os.path.dirname(SCRIPT_DIR)
CLUBB_ROOT = os.path.join(INPUT_DIR, "..")
DEFAULT_OUTPUT_DIR = os.path.join(CLUBB_ROOT, "output")
DEFAULT_TUNER_STATS = "all_tuning_stats"
MULTI_COL_PARAMS_SCRIPT = os.path.join(SCRIPT_DIR, "create_multi_col_params.py")

if os.path.abspath(CLUBB_ROOT) not in sys.path:
    sys.path.insert(0, os.path.abspath(CLUBB_ROOT))
from input.scripts.case_json_to_namelist import (
    available_case_names as json_case_names, case_definition_record, case_to_namelist,
)

HR_SPEC_RE = re.compile(
    r"^[A-Za-z_]\w*/[-+]?\d*\.?\d+(?:[eEdD][-+]?\d+)?:[-+]?\d*\.?\d+(?:[eEdD][-+]?\d+)?/\d+"
    r"(?:,[A-Za-z_]\w*/[-+]?\d*\.?\d+(?:[eEdD][-+]?\d+)?:[-+]?\d*\.?\d+(?:[eEdD][-+]?\d+)?/\d+)*$"
)


INPUT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = "default"
CATALOG_FILES = {
    "cases": "case_definitions.jsonc",
    "parameters": "clubb_params.jsonc",
    "flags": "clubb_config_flags.jsonc",
}


def load_catalog(kind, json_file=None):
    """Read one authoritative defaults/presets catalog without creating files."""
    from input.scripts.namelist import load_json, require_non_null
    if kind not in CATALOG_FILES:
        raise ValueError(f"Unknown input catalog {kind!r}")
    catalog = load_json(json_file or INPUT_ROOT / CATALOG_FILES[kind])
    if not isinstance(catalog, dict) or "defaults" not in catalog or not isinstance(catalog.get("presets"), dict):
        raise ValueError(f"{kind} catalog requires defaults and presets")
    require_non_null(catalog["defaults"], f"{kind}.defaults")
    for name, patch in catalog["presets"].items():
        require_non_null({key: value for key, value in patch.items() if key != "_metadata"}
                         if isinstance(patch, dict) else patch, f"{kind}.presets.{name}")
    return catalog


def available_configs(root=None):
    """Return the built-in default choice and names defined by any catalog."""
    root = Path(root or INPUT_ROOT)
    return sorted({DEFAULT_CONFIG} | {name for kind, filename in CATALOG_FILES.items()
                   for name in load_catalog(kind, root / filename)["presets"]})


def normalize_input_selection(selection):
    """Keep preset names intact and resolve JSON paths against the caller's directory."""
    value = os.fspath(selection)
    path = Path(value).expanduser()
    if path.suffix.lower() in {".json", ".jsonc"} or path.is_file() or value.startswith(("/", "~", ".")):
        return str(path.resolve())
    return value


def input_patch(kind, selection, *, json_file=None):
    """Read an explicit named preset or JSON override with identical contents."""
    from input.scripts.namelist import load_json
    value = normalize_input_selection(selection)
    path = Path(value)
    if path.is_absolute():
        return load_json(path)
    catalog = load_catalog(kind, json_file)
    if value not in catalog["presets"]:
        raise ValueError(f"Unknown {kind} preset {value!r}. Available: {', '.join(catalog['presets'])}")
    return catalog["presets"][value]


def resolve_input_settings(kind, selection=None, *, config=None, json_file=None):
    """Apply defaults, a shared config and explicit overrides; default skips presets."""
    import copy
    from input.scripts.namelist import merge_settings, require_non_null
    catalog = load_catalog(kind, json_file)
    if config == DEFAULT_CONFIG:
        config = None
    if config and json_file is None and config not in available_configs():
        raise ValueError(f"Unknown config {config!r}")
    result = copy.deepcopy(catalog["defaults"])
    patches = []
    if config in catalog["presets"]:
        patches.append(catalog["presets"][config])
    if selection is not None:
        patches.append(input_patch(kind, selection, json_file=json_file))
    for patch in patches:
        patch = {key: value for key, value in patch.items() if key != "_metadata"} if isinstance(patch, dict) else patch
        require_non_null(patch, f"{kind}.overrides")
        result = merge_settings(result, patch)
    return result


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
        out_file = os.path.join(output_dir, f"{case_name}_multicol_params.json")
    cmd = [sys.executable, MULTI_COL_PARAMS_SCRIPT, "-params_file", params_file, "-output_file", out_file]

    try:
        ngrdcol = int(multicol_spec)
    except ValueError:
        cmd.extend(["-multicol", multicol_spec])
    else:
        cmd.extend(["-mode", integer_mode, "-multicol", str(ngrdcol)])

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


def available_case_names():
    """Discover cases solely through the JSON catalog and JSON snapshots."""
    return json_case_names()


def resolve_case_model_input(case_name, config=None):
    """Return the authoritative case JSON and provenance without creating files."""
    definition = case_definition_record(case_name, config=config)
    return {**definition, "model_setting": definition["namelists"].get("model_setting", {})}


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


def normalize_override_string(override_string):
    """Return an override value, accepting an optional leading ``-override``."""
    value = "" if override_string is None else str(override_string).strip()
    if value.lower().startswith("-override"):
        value = value[len("-override"):].lstrip(" =")
    return value


def resolve_override_pairs(override, case_name=None):
    """Resolve assignments or JSON settings for the current case."""
    # Repeated -override arguments are applied in order, including mixed formats.
    if isinstance(override, (list, tuple)):
        return [pair for value in override for pair in resolve_override_pairs(value, case_name)]
    value = normalize_override_string(override)
    if value.startswith(("{", "[")):
        settings = json.loads(value)
    elif "=" not in value and (value.endswith((".json", ".jsonc")) or os.path.isfile(value)):
        from input.scripts.namelist import load_json
        settings = load_json(value)
    else:
        return parse_override_pairs(value)

    if not isinstance(settings, dict):
        raise ValueError("Override JSON must be an object of settings or case-specific objects")
    if any(isinstance(item, dict) for item in settings.values()):
        if not all(isinstance(item, dict) for item in settings.values()):
            raise ValueError("Override JSON cannot mix global settings and case-specific objects")
        known_cases = set(available_case_names()) | ({case_name} if case_name else set())
        unknown = settings.keys() - known_cases - {"all"}
        if unknown:
            raise ValueError("Unknown case(s) in override JSON: " + ", ".join(sorted(unknown)))
        if case_name is None and settings.keys() - {"all"}:
            raise ValueError("A case name is required for case-specific override JSON")
        groups = settings.values()
    else:
        groups = [settings]
    # Validate the whole profile so a bad entry cannot silently pass on another case.
    for group in groups:
        for key, item in group.items():
            if not isinstance(key, str) or not key.strip() or "=" in key or "," in key:
                raise ValueError(f"Invalid override key: {key!r}")
            if item is None or isinstance(item, (dict, list)):
                raise ValueError(f"Override {key} must have a scalar value")
    if any(isinstance(item, dict) for item in settings.values()):
        # Case-specific settings take precedence over shared defaults.
        settings = {**settings.get("all", {}), **settings.get(case_name, {})}
    return [(key, ".true." if item is True else ".false." if item is False else str(item))
            for key, item in settings.items()]


def filter_tuned_overrides(override, case_names, tuned_names):
    """Retain case-resolved physics overrides without replacing tuned columns."""
    tuned_names = {str(name).lower() for name in tuned_names}
    resolved_overrides = {
        case_name: resolve_override_pairs(override, case_name)
        for case_name in case_names
    }
    if not any(
        name.rsplit(".", 1)[-1].lower() in tuned_names
        for pairs in resolved_overrides.values()
        for name, _ in pairs
    ):
        return override
    remaining_overrides = {
        case_name: {
            name: value
            for name, value in pairs
            if name.rsplit(".", 1)[-1].lower() not in tuned_names
        }
        for case_name, pairs in resolved_overrides.items()
    }
    return json.dumps(remaining_overrides) if any(remaining_overrides.values()) else ""


def override_value(override_string, settings, case_name=None):
    """Apply CLI assignments to JSON groups, preserving types and derived components."""
    import copy
    from input.scripts.namelist import parse_value
    result = copy.deepcopy(settings)
    for key, raw in resolve_override_pairs(override_string, case_name):
        if "." in key:
            group_name, key = key.split(".", 1)
            groups = [name for name in result if name.casefold() == group_name.casefold()]
            if not groups or "." in key:
                raise ValueError(f"override failed - namelist group {group_name!r} is not present")
        else:
            groups = [name for name, values in result.items()
                      if any(field.casefold() == key.split("%", 1)[0].casefold() for field in values)]
            if not groups:
                raise ValueError(f"override failed - key {key!r} could not be matched. To add a key specify [namelist.]key=value")
        for group in groups:
            values = result[group]
            parts = key.split("%")
            for part in parts[:-1]:
                actual = next((name for name in values if name.casefold() == part.casefold()), part)
                values = values.setdefault(actual, {})
            actual = next((name for name in values if name.casefold() == parts[-1].casefold()), parts[-1])
            values[actual] = parse_value(raw)
    return result


def set_namelist_value(settings, namelist, key, value, *, quote=False):
    """Set one JSON group value before final serialization."""
    from input.scripts.namelist import parse_value
    settings.setdefault(namelist, {})[key] = str(value) if quote else parse_value(value) if isinstance(value, str) else value
    return settings


def set_stats_output_dir(settings, output_dir):
    """Set the final output directory in the resolved JSON settings."""
    return set_namelist_value(settings, "stats_setting", "output_dir", str(output_dir).replace("\\", "/"), quote=True)


def set_stats_string(settings, key, value):
    """Set a stats string before serialization."""
    return set_namelist_value(settings, "stats_setting", key, value, quote=True)


def set_stats_value(settings, key, value):
    """Set a typed stats value before serialization."""
    return set_namelist_value(settings, "stats_setting", key, value)


def set_model_string(settings, key, value):
    """Set a model string before serialization."""
    return set_namelist_value(settings, "model_setting", key, value, quote=True)


def set_model_value(settings, key, value):
    """Set a typed model value before serialization."""
    return set_namelist_value(settings, "model_setting", key, value)


def apply_namelist_overrides(args, settings, model_settings, output_dir):
    """Apply direct CLI settings to resolved JSON, then the general override."""
    model = settings.setdefault("model_setting", {})
    stats = settings.setdefault("stats_setting", {})
    for key in ("dt_main", "dt_rad"):
        if getattr(args, key) is not None:
            model[key] = float(getattr(args, key))
    for option, grid_type in (("zt_grid", 2), ("zm_grid", 3)):
        if getattr(args, option) is not None:
            model.update(nzmax=args.nzmax, grid_type=grid_type,
                         zt_grid_fname=args.zt_grid or "", zm_grid_fname=args.zm_grid or "")
    disabled = not settings.get("clubb_stats_nl", {}).get("entry")
    if args.tout is not None:
        disabled = disabled or args.tout == 0
        if args.tout:
            stats["stats_tout"] = float(args.tout)
    for key in ("stats_tstart", "stats_tend"):
        if getattr(args, key) is not None:
            stats[key] = getattr(args, key)
    if disabled:
        stats["l_stats"] = False
    stats["stats_output_filename"] = "" if disabled else f"{args.case_name}_stats.nc"
    if args.debug is not None:
        model["debug_level"] = int(args.debug)
    if args.max_iters is not None:
        new_final = model["time_initial"] + model["dt_main"] * args.max_iters
        model["time_final"] = min(model["time_final"], new_final)
    if args.override:
        settings = override_value(args.override, settings, args.case_name)
    return set_stats_output_dir(settings, output_dir)


def _ensure_repo_root_on_path() -> None:
    repo_root = os.path.abspath(CLUBB_ROOT)
    if repo_root not in sys.path:
        sys.path.insert(0, repo_root)


def create_case_settings(
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
    override: str | list[str] | None = None,
) -> dict:
    """Assemble resolved JSON groups for a normal run without emitting a namelist."""
    output_dir_abs = os.path.abspath(os.fspath(output_dir) if output_dir else DEFAULT_OUTPUT_DIR)
    os.makedirs(output_dir_abs, exist_ok=True)

    model_input = resolve_case_model_input(case_name, config=config)
    case_name = model_input["name"]
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

    from input.scripts.parameters_json_to_namelist import parameter_settings
    from input.scripts.flags_json_to_namelist import flag_settings
    from input.scripts.namelist import atomic_text
    from input.scripts.stats_json_to_namelist import select_stats_entries, ENTRY_FIELDS

    parameters = parameter_settings(params, config=config)
    if multicol is not None:
        params_file = os.path.join(output_dir_abs, f"{case_name}_base_parameters.json")
        atomic_text(params_file, json.dumps(parameters, indent=2, allow_nan=False) + "\n")
        params_file = convert_to_multi_col(params_file, case_name, output_dir_abs, multicol, batch_size=batch_size)
        parameters = parameter_settings(params_file)
    if silhs_params:
        parameters["silhs"] = parameter_settings(silhs_params, config=config)["silhs"]
    flag_values = flag_settings(flags, config=config)
    entries = select_stats_entries(stats)
    settings = {"clubb_params_nl": parameters["model"], "silhs_params_nl": parameters["silhs"],
                "configurable_clubb_flags_nl": flag_values["model"],
                "configurable_silhs_flags_nl": flag_values["silhs"], **model_input["namelists"],
                "clubb_stats_nl": {"entry": [" | ".join(entry[key] for key in ENTRY_FIELDS) for entry in entries]}}
    if "multicol" in parameters:
        settings = {"multicol_def": parameters["multicol"], **settings}

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
    return apply_namelist_overrides(args, settings, model_input["model_setting"], output_dir_abs)


def case_run_name(selection):
    """Return a safe output stem for a named case or case override path."""
    return Path(os.fspath(selection)).name.removesuffix(".case.json").removesuffix(".json")


def write_case_namelist(settings, path):
    """Call the specialized converters once, at the top-level runner boundary."""
    from input.scripts.parameters_json_to_namelist import parameters_to_namelist
    from input.scripts.flags_json_to_namelist import flags_to_namelist
    from input.scripts.stats_json_to_namelist import stats_to_namelist, ENTRY_FIELDS
    from input.scripts.namelist import atomic_text
    parameters = {"model": settings.get("clubb_params_nl", {}), "silhs": settings.get("silhs_params_nl", {})}
    if "multicol_def" in settings:
        parameters["multicol"] = settings["multicol_def"]
    flags = {"model": settings.get("configurable_clubb_flags_nl", {}), "silhs": settings.get("configurable_silhs_flags_nl", {})}
    omitted = {"clubb_params_nl", "silhs_params_nl", "multicol_def", "configurable_clubb_flags_nl", "configurable_silhs_flags_nl", "clubb_stats_nl"}
    text = parameters_to_namelist(parameters) + flags_to_namelist(flags)
    text += case_to_namelist({key: value for key, value in settings.items() if key not in omitted})
    entries = [dict(zip(ENTRY_FIELDS, [piece.strip() for piece in entry.split(" | ", 3)]))
               for entry in settings.get("clubb_stats_nl", {}).get("entry", [])]
    atomic_text(path, text + stats_to_namelist(entries))
    return Path(path)


def create_case_namelist_file(case_name, output_dir=None, **options):
    """Resolve all JSON inputs, then write the complete runner namelist once."""
    settings = create_case_settings(case_name, output_dir, **options)
    return write_case_namelist(settings, Path(output_dir or DEFAULT_OUTPUT_DIR) / f"{case_run_name(case_name)}.in")


def build_loss_settings(
    case_defaults: dict,
    selected_fields: list[str],
    num_time_windows: int = 1,
) -> dict:
    """Build typed loss-driver request settings for the final aggregate input."""
    return {"tuner_loss_nl": {"les_stats_file": case_defaults["les_stats_file"],
            "clubb_var_names": list(selected_fields), "benchmark_var_name": list(selected_fields),
            "altitude_comparison_range": list(case_defaults["altitude_comparison_range"]),
            "time_average_range": list(case_defaults["time_average_range"]),
            "num_time_windows": int(num_time_windows)}}


def create_normalized_benchmark_file(
    case_name: str,
    case_defaults: dict,
    selected_fields: list[str],
    job_dir: Path,
) -> Path:
    """Convert the case LES benchmark to CLUBB-facing field names."""
    _ensure_repo_root_on_path()
    from input.scripts.benchmark_converter import convert_benchmark_file

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


def prune_clubb_stats_namelist(settings, requested_vars):
    """Keep requested registry entries in the JSON settings before final emission."""
    entries = settings.get("clubb_stats_nl", {}).get("entry", [])
    requested = set(requested_vars)
    kept = [entry for entry in entries if entry.split("|", 1)[0].strip() in requested]
    missing = requested - {entry.split("|", 1)[0].strip() for entry in kept}
    if missing:
        raise RuntimeError("Requested clubb_var_names missing from stats selection: " + ", ".join(sorted(missing)))
    settings["clubb_stats_nl"] = {"entry": kept}
    return settings


def _loss_timing_settings(case_name: str, model_file: Path, case_defaults: dict, num_time_windows: int) -> dict:
    vals = model_file
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
    override: str | list[str] | None = None,
    verbose: bool = False,
) -> tuple[Path, list[str], dict]:
    """Create the aggregate namelist consumed by the in-memory loss driver."""
    case_selector = case_name
    case_name = case_run_name(case_selector)
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

        case_defaults = read_case_defaults(case_selector)
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
        from input.scripts.parameters_json_to_namelist import parameter_settings
        from input.scripts.namelist import atomic_text
        base_settings = parameter_settings(params, config=config)
        base_params_file = str(output_path / f"{case_name}_base_parameters.json")
        atomic_text(base_params_file, json.dumps(base_settings, indent=2, allow_nan=False) + "\n")
        params_for_case = convert_to_multi_col(
            base_params_file,
            case_name,
            str(output_path),
            str(batch_size),
            batch_size=batch_size,
            integer_mode="duplicate",
            out_file=str(output_path / f"{case_name}_tuning_multicol_params.json"),
        )
        multicol_for_case = None
        batch_size_for_case = None

    clubb_in = create_case_settings(
        case_selector,
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

    clubb_in.update(build_loss_settings(
        normalized_case_defaults,
        selected_fields,
        num_time_windows=num_time_windows,
    ))

    model_file = clubb_in["model_setting"]
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
        clubb_in = override_value(override, clubb_in, case_name)
    clubb_in = set_stats_output_dir(clubb_in, str(output_path))

    aggregate_path = write_case_namelist(clubb_in, output_path / f"{case_run_name(case_name)}.in")
    return aggregate_path, selected_fields, case_defaults


def add_input_arguments(parser, *, cases=True):
    """Declare the shared config, input overrides and independent stats selector."""
    if cases:
        parser.add_argument("-cases", dest="case_selector", metavar="NAME|JSON", help="Case preset or JSON override file; also accepted as the positional case.")
    parser.add_argument("-config", metavar="NAME", help="Intermediate preset applied wherever that name exists in cases, parameters or flags; default uses catalog defaults without a preset.")
    parser.add_argument('-params', '-params_file', dest='params', metavar="NAME|JSON",
        help="Parameter preset name or JSON overrides, applied after -config.")
    parser.add_argument("-flags", metavar="NAME|JSON",
        help="Flag preset name or JSON overrides, applied after -config.")
    parser.add_argument('-silhs_params_file', dest='silhs_params', metavar="[FILE]",
        help="Parameter JSON supplying only the SILHS category override.")
    parser.add_argument("-stats", metavar="[GROUPS|JSON]",
        help=("Stats categories or named presets from input/stats.jsonc, or a JSON list of variable names. "
              "Join categories with ',' or '+'. Paths select one subtree; var:NAME selects a variable. "
              "Use 'all' for every stat or 'none' to disable stats. Default: standard."))


def add_namelist_arguments(parser):
    """Declare shared input selectors and direct model settings."""
    add_input_arguments(parser)
    parser.add_argument('-output_dir', dest='out_dir', metavar="[DIR]",
        help="Output directory for the generated namelist and run results. Relative paths are rooted under output/ unless they already start with output/; absolute paths are used directly. Default: output")
    parser.add_argument("-multicol", metavar="[NUM|SPEC]", type=validate_multicol,
        help=("Generate a multi-column parameter file. "
              "Use an integer for dup_tweak mode, e.g. -multicol 4, or a parameter-grid spec like "
              "-multicol C8/0.2:0.8/4"))
    parser.add_argument("-batch_size", metavar="[NUM]", type=int,
        help=("Runtime batch size written to &multicol_def. "
              "Only meaningful together with -multicol. "
              "Values larger than ngrdcol are clipped to ngrdcol."))
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
        "-override", action="append", metavar="ASSIGNMENTS|JSON",
        help=("Comma-separated key=value pairs, inline JSON, or a JSON file. JSON maps "
              'settings to values, or case names to settings; "all" supplies shared defaults. Repeat to apply in order; '
              "later values win. Use namelist.key to add a key."),
    )


def normalize_case_selection(args, parser):
    """Preserve the case selector while deriving a safe name for run artifacts."""
    selection = getattr(args, "case_selector", None) or args.case_name
    if selection is None and getattr(args, "config", None) in load_catalog("cases")["presets"]:
        selection = args.config
    if not selection:
        parser.error("Provide a case preset or JSON file as CASE or -cases")
    args.case_selector = selection
    args.case_name = case_run_name(selection)
    return args


def case_namelist_command(args, output_dir):
    """Return the generator command for args.case_name and its model options.

    args uses add_namelist_arguments field names. output_dir must be the
    caller's resolved output directory; it replaces args.out_dir. Runtime
    options are omitted and repeated overrides retain their original order.
    """
    parser = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
    add_namelist_arguments(parser)
    command = [sys.executable, str(Path(__file__).resolve()), "-output_dir", str(output_dir)]
    for action in parser._actions:
        if action.dest == "out_dir":
            continue
        value = getattr(args, action.dest, None)
        if value is None:
            continue
        values = value if isinstance(action, argparse._AppendAction) else [value]
        for item in values:
            command.extend([action.option_strings[0], str(item)])
    return [*command, args.case_name]


def create_case_namelist_from_args(args, output_dir, *, log_file=None, timeout=None):
    """Generate args' aggregate namelist in resolved output_dir and return its Path.

    args supplies the case name and parsed model options. Run the generator as
    a subprocess so concurrent callers can send console output to separate open
    log_file streams (None inherits stdout). timeout is its limit in seconds.
    Generation failure/timeout propagates from subprocess.run. This helper does
    not initialize a model or select its runtime.
    """
    subprocess.run(case_namelist_command(args, output_dir), check=True,
                   stdout=log_file, stderr=subprocess.STDOUT, timeout=timeout)
    return Path(output_dir) / f"{case_run_name(args.case_name)}.in"


def parse_forwarded_args(parser, argv=None):
    """Return (wrapper Namespace, forwarded tokens) from argv or sys.argv.

    parser owns wrapper options; shared declarations identify runtime/model
    options whose values must be kept together when forwarded. Unknown options
    are rejected. -- forwards the remaining tokens unchanged to the SCM parser.
    """
    model = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
    add_namelist_arguments(model)
    # Runtime option declarations have the same owner as actual runner selection.
    from utilities.clubb_case_runner import add_runner_arguments

    add_runner_arguments(model)
    actions = dict(model._option_string_actions)
    own, forwarded = [], []
    tokens = list(sys.argv[1:] if argv is None else argv)
    index = 0
    while index < len(tokens):
        token = tokens[index]
        if token == "--":
            forwarded.extend(tokens[index + 1:])
            break
        option = token.split("=", 1)[0]
        action = parser._option_string_actions.get(option)
        target = own
        if action is None and option in actions:
            action, target = actions[option], forwarded
        elif action is None and token.startswith("-"):
            parser.error(f"unrecognized argument: {token}")
        target.append(token)
        index += 1
        if action is None or action.nargs == 0 or "=" in token:
            continue
        if action.nargs in ("+", "*"):
            while index < len(tokens) and not tokens[index].startswith("-"):
                target.append(tokens[index])
                index += 1
        else:
            count = action.nargs if isinstance(action.nargs, int) else 1
            if index + count > len(tokens):
                parser.error(f"{option} requires a value")
            target.extend(tokens[index:index + count])
            index += count
    return parser.parse_args(own), forwarded


def main():
    parser = argparse.ArgumentParser(description="Create an aggregated CLUBB case namelist.", add_help=False, allow_abbrev=False)
    parser.add_argument("-h", "-help", action="help", help="Show this help and exit.")
    add_namelist_arguments(parser)
    parser.add_argument("-list", choices=[*CATALOG_FILES, "stats", "configs"], help="List predefined selectors without creating files.")
    parser.add_argument("case_name", nargs="?", help="Case preset or JSON override file")
    args = parser.parse_args()
    if args.list:
        if args.list == "stats":
            from input.scripts.stats_json_to_namelist import DEFAULT_STATS_JSON, load_stats_categories
            from input.scripts.namelist import load_json
            names = [path.removeprefix("all/") for path in load_stats_categories()]
            names.extend(load_json(DEFAULT_STATS_JSON).get("presets", {}))
        else:
            names = available_configs() if args.list == "configs" else load_catalog(args.list)["presets"]
        print("\n".join(names))
        return
    args = normalize_case_selection(args, parser)

    try:
        clubb_input_namelist = create_case_namelist_file(
            args.case_selector,
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
    except (RuntimeError, ValueError, OSError) as exc:
        sys.exit(str(exc))

    print(clubb_input_namelist)


if __name__ == "__main__":
    main()
