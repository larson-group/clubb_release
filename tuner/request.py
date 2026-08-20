"""Request validation and case-default parsing for tuning jobs."""

from __future__ import annotations

import json
import math
from pathlib import Path
import re
from functools import lru_cache

from tuner.case_defaults import read_case_defaults
from tuner.presets import apply_preset
from utilities.benchmark_converter import supported_fields
from tuner.taylor_metrics import (
    AGGREGATION_MODE_NAMES,
    DEFAULT_AGGREGATION_WEIGHTS,
    DEFAULT_AGGREGATION_MODE,
    DEFAULT_LOSS_MODE,
    DEFAULT_NUM_TIME_WINDOWS,
    LOSS_MODE_NAMES,
    LOSS_POLICY_CONSTANTS,
    LOSS_POLICY_VERSION,
    DEFAULT_TIME_WINDOW_AGGREGATION_SCOPE,
    TIME_WINDOW_AGGREGATION_SCOPES,
    WORST_QUANTILE_FRACTION,
    normalize_aggregation_weights,
)
from tuner.tuning_strategy import normalize_strategy_config
from utilities.create_case_namelist import normalize_override_string, parse_override_pairs, resolve_tunable_config_dir
from utilities.clubb_settings_validation import (
    build_settings_schema,
    canonical_flag_name,
    canonical_parameter_name,
    evaluate_settings,
    is_independently_tunable,
    setting_key,
)


REQUIRED_TUNABLE_CONFIG_FILES = (
    "tunable_parameters.in",
    "configurable_model_flags.in",
    "silhs_parameters.in",
)


def load_request(request_path: Path) -> dict:
    """Load and validate the tuning request."""
    request = apply_preset(json.loads(request_path.read_text(encoding="utf-8")))

    request["config"] = _normalize_config(request.get("config"))
    available_targets = _tunable_parameter_names(request["config"])
    cases, case_configs, case_defaults = _normalize_case_configs(request)

    selected_fields = [str(name).strip() for name in request.get("selected_fields", []) if str(name).strip()]
    if not selected_fields:
        raise RuntimeError("Tuning request must contain at least one selected field")
    unsupported_fields = sorted(set(selected_fields) - set(supported_fields()))
    if unsupported_fields:
        raise RuntimeError("Unsupported normalized benchmark field(s): " + ", ".join(unsupported_fields))

    raw_ranges = list(request.get("parameter_ranges", []))
    if not raw_ranges:
        raise RuntimeError("Tuning request must contain at least one parameter range")
    parameter_ranges = []
    seen_params = set()
    owned_targets = set()
    for spec in raw_ranges:
        if not isinstance(spec, dict):
            raise RuntimeError("Each parameter range must be an object")
        name = str(spec.get("name", "")).strip()
        if not name:
            raise RuntimeError("Each parameter range must include a name")
        if name in seen_params:
            raise RuntimeError(f"Duplicate tuning parameter: {name}")
        seen_params.add(name)
        raw_targets = spec.get("targets", [name])
        if isinstance(raw_targets, str):
            raw_targets = [raw_targets]
        targets = [str(target).strip() for target in raw_targets or [] if str(target).strip()]
        if not targets:
            raise RuntimeError(f"{name} must include at least one physical target")
        duplicate_targets = owned_targets.intersection(targets)
        if duplicate_targets:
            raise RuntimeError(
                f"Physical tuning parameter(s) assigned by more than one range: {', '.join(sorted(duplicate_targets))}"
            )
        if len(targets) != len(set(targets)):
            raise RuntimeError(f"{name} repeats a physical target")
        unknown_targets = sorted(set(targets) - available_targets)
        if unknown_targets:
            raise RuntimeError(
                f"Unknown tunable parameter target(s) for config {request['config']}: "
                + ", ".join(unknown_targets)
            )
        owned_targets.update(targets)
        try:
            min_value = float(spec.get("min"))
            max_value = float(spec.get("max"))
        except (TypeError, ValueError):
            raise RuntimeError(f"Invalid numeric range for {name}")
        if not math.isfinite(min_value) or not math.isfinite(max_value):
            raise RuntimeError(f"{name} requires finite min and max")
        if min_value > max_value:
            raise RuntimeError(f"{name} requires min <= max")
        parameter_ranges.append({"name": name, "targets": targets, "min": min_value, "max": max_value})

    request["case_name"] = cases[0]
    request["cases"] = cases
    request["case_configs"] = case_configs
    request["case_defaults"] = case_defaults
    request["override"] = normalize_override_string(request.get("override"))
    resolution = evaluate_tune_settings(request["config"], request["override"])
    resolution_errors = [issue for issue in resolution.get("issues", []) if issue.get("severity") == "error"]
    if resolution_errors:
        raise RuntimeError("; ".join(str(issue.get("message") or "Invalid CLUBB settings.") for issue in resolution_errors))
    inactive_targets = sorted(
        target
        for spec in parameter_ranges
        for target in spec["targets"]
        if not is_independently_tunable(resolution.get("parameter_states", {}).get(target))
    )
    if inactive_targets:
        details = "; ".join(
            f"{name}: {resolution['parameter_states'][name]['reason']}" for name in inactive_targets
        )
        raise RuntimeError("Tune request selects inactive parameter(s): " + details)
    parameter_ranges = apply_required_parameter_links(parameter_ranges, resolution)
    request["selected_fields"] = selected_fields
    request["parameter_ranges"] = parameter_ranges
    request["settings_resolution"] = resolution
    request["batch_size"] = _positive_int(request.get("batch_size"), "Tuning request batch_size")
    request["max_workers"] = _positive_int(request.get("max_workers", 1), "Tuning request max_workers")
    request["case_weights"] = _normalize_weights(request.get("case_weights"), set(cases), "case")
    request["field_weights"] = _normalize_weights(request.get("field_weights"), set(selected_fields), "field")
    request["loss_mode"] = _normalize_choice(
        request.get("loss_mode", DEFAULT_LOSS_MODE),
        set(LOSS_MODE_NAMES),
        "loss_mode",
    )
    request["aggregation_mode"] = _normalize_choice(
        request.get("aggregation_mode", DEFAULT_AGGREGATION_MODE),
        set(AGGREGATION_MODE_NAMES),
        "aggregation_mode",
    )
    request["time_window_aggregation_mode"] = _normalize_choice(
        request.get("time_window_aggregation_mode", request["aggregation_mode"]),
        set(AGGREGATION_MODE_NAMES),
        "time_window_aggregation_mode",
    )
    try:
        request["aggregation_weights"] = normalize_aggregation_weights(
            request.get("aggregation_weights", DEFAULT_AGGREGATION_WEIGHTS)
        )
    except ValueError as exc:
        raise RuntimeError(str(exc)) from exc
    request["time_window_aggregation_scope"] = _normalize_choice(
        request.get("time_window_aggregation_scope", DEFAULT_TIME_WINDOW_AGGREGATION_SCOPE),
        set(TIME_WINDOW_AGGREGATION_SCOPES),
        "time_window_aggregation_scope",
    )
    request.pop("time_window_mode", None)
    request.pop("num_time_windows", None)
    requested_loss_policy_version = request.get("loss_policy_version", LOSS_POLICY_VERSION)
    if requested_loss_policy_version != LOSS_POLICY_VERSION:
        raise RuntimeError(
            f"Unsupported loss_policy_version {requested_loss_policy_version}; expected {LOSS_POLICY_VERSION}"
        )
    request["loss_policy_version"] = LOSS_POLICY_VERSION
    request["loss_policy_constants"] = dict(LOSS_POLICY_CONSTANTS)
    request["aggregation_options"] = {
        "worst_quantile_fraction": WORST_QUANTILE_FRACTION,
        "quantile_weights": list(request["aggregation_weights"]),
        "time_window_aggregation_scope": request["time_window_aggregation_scope"],
    }

    try:
        request["strategy"] = normalize_strategy_config(request)
    except ValueError as exc:
        raise RuntimeError(str(exc)) from exc
    if request["strategy"]["name"] == "random":
        request["max_samples"] = request["strategy"]["options"].get("max_samples")
    else:
        request.pop("max_samples", None)
    if request["strategy"]["name"] == "simann":
        options = request["strategy"]["options"]
        if options.get("chain_count") is None:
            options["chain_count"] = max(1, request["max_workers"] * request["batch_size"])
        request["total_samples"] = options.get("max_iters") * options.get("chain_count")
    elif request["strategy"]["name"] == "adam":
        options = request["strategy"]["options"]
        columns_per_chain = 2 * int(options["spsa_pairs"])
        if request["batch_size"] % columns_per_chain != 0:
            raise RuntimeError(
                "Adam requires batch_size to be divisible by 2 * spsa_pairs "
                f"({request['batch_size']} is not divisible by {columns_per_chain})"
            )
        chains_per_batch = request["batch_size"] // columns_per_chain
        concurrent_batches = math.ceil(request["max_workers"] / len(cases))
        chain_count = chains_per_batch * concurrent_batches
        options.update(
            {
                "chains_per_batch": chains_per_batch,
                "concurrent_batches": concurrent_batches,
                "chain_count": chain_count,
            }
        )
        request["total_samples"] = chain_count * (
            2 + columns_per_chain * int(options["max_updates"])
        )
    else:
        request["total_samples"] = request["strategy"]["options"].get("max_samples")

    if request.get("seed") is not None:
        try:
            request["seed"] = int(request["seed"])
        except (TypeError, ValueError):
            raise RuntimeError("Tuning request seed must be an integer")
    return request


def _normalize_config(raw_config) -> str:
    """Validate and normalize the tunable config name/path in a tuning request."""
    config = str(raw_config or "default").strip() or "default"
    try:
        config_dir = Path(resolve_tunable_config_dir(config))
    except RuntimeError as exc:
        raise RuntimeError(str(exc)) from exc
    missing = [filename for filename in REQUIRED_TUNABLE_CONFIG_FILES if not (config_dir / filename).is_file()]
    if missing:
        raise RuntimeError(
            f"Tunable config '{config}' is missing required file(s): " + ", ".join(missing)
        )
    return config


def _tunable_parameter_names(config: str) -> set[str]:
    """Read the config's physical parameter names without importing Dash/F2PY.

    The scheduler repeats this against the compiled API at launch time; this
    early checked-in-config validation gives CLI/Dash users a useful error
    before a worker process is created.
    """
    params_path = Path(resolve_tunable_config_dir(config)) / "tunable_parameters.in"
    pattern = re.compile(r"^\s*([A-Za-z][A-Za-z0-9_]*)\s*=")
    names = set()
    for line in params_path.read_text(encoding="utf-8").splitlines():
        match = pattern.match(line)
        if match:
            names.add(match.group(1))
    if not names:
        raise RuntimeError(f"Tunable config '{config}' has no readable parameter names")
    return names


def _parse_namelist_value(raw: str):
    """Parse the scalar config values relevant to resolver decisions."""
    text = str(raw).strip().rstrip(",")
    if text.lower() in {".true.", "true", "t"}:
        return True
    if text.lower() in {".false.", "false", "f"}:
        return False
    try:
        return int(text)
    except ValueError:
        return text


def _read_namelist_assignments(path: Path) -> dict[str, str]:
    """Read simple scalar namelist assignments without importing Dash."""
    assignment = re.compile(r"^\s*([A-Za-z_]\w*)\s*=\s*([^!,/]+)")
    values = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        match = assignment.match(line)
        if match:
            values[match.group(1)] = match.group(2).strip()
    return values


def _is_bool_literal(value: str) -> bool:
    """Return whether a namelist scalar is a standard logical literal."""
    return str(value).strip().lower() in {".true.", "true", "t", ".false.", "false", "f"}


@lru_cache(maxsize=None)
def settings_schema_for_tune_config(config: str) -> dict:
    """Load one checked-in config into the shared UI-neutral settings schema."""
    config_dir = Path(resolve_tunable_config_dir(config))
    config_values = _read_namelist_assignments(config_dir / "configurable_model_flags.in")
    flag_defaults = {
        name: _parse_namelist_value(value)
        for name, value in config_values.items()
        if _is_bool_literal(value)
    }
    parameter_defaults = {
        "flags": {
            name: value
            for name, value in config_values.items()
            if not _is_bool_literal(value)
        },
        "tunable": _read_namelist_assignments(config_dir / "tunable_parameters.in"),
        "silhs": _read_namelist_assignments(config_dir / "silhs_parameters.in"),
    }
    metadata = [
        {"file": file_name, "name": name}
        for file_name, values in parameter_defaults.items()
        for name in values
    ]
    return build_settings_schema(flag_defaults, parameter_defaults, metadata)


def evaluate_tune_settings(config: str, override: str = "") -> dict:
    """Evaluate a Tune config and its override through the shared settings engine."""
    schema = settings_schema_for_tune_config(config)
    flag_defaults = dict(schema.get("flag_defaults") or {})
    parameter_defaults = dict(schema.get("parameter_defaults") or {})
    override_flags = {}
    override_parameters = {}
    for raw_name, value in parse_override_pairs(override):
        flag_name = canonical_flag_name(raw_name)
        parameter_name = canonical_parameter_name(raw_name)
        if flag_name in flag_defaults:
            override_flags[flag_name] = _parse_namelist_value(value)
            continue
        for file_name, defaults in parameter_defaults.items():
            if parameter_name in dict(defaults or {}):
                override_parameters[setting_key(file_name, parameter_name)] = value
                break
    return evaluate_settings(schema, flag_values=override_flags, parameter_values=override_parameters)


def resolve_tune_settings(config: str, override: str = ""):
    """Deprecated compatibility alias for :func:`evaluate_tune_settings`."""
    return evaluate_tune_settings(config, override)


def apply_required_parameter_links(parameter_ranges: list[dict], resolution: dict) -> list[dict]:
    """Expand a sampled coordinate to all model-required equal targets."""
    updated = [dict(spec, targets=list(spec["targets"])) for spec in parameter_ranges]
    owned = {target for spec in updated for target in spec["targets"]}
    for spec in updated:
        targets = set(spec["targets"])
        for relation in resolution.get("coupled_parameters", []):
            members = set(relation["members"])
            if not targets.intersection(members) or members.issubset(targets):
                continue
            missing = members - targets
            conflict = missing.intersection(owned)
            if conflict:
                raise RuntimeError(
                    "Required equal parameter link conflicts with another range: "
                    + ", ".join(sorted(conflict))
                )
            spec["targets"].extend(sorted(missing))
            owned.update(missing)
    return updated


def read_case_tuner_defaults(case_name: str, overrides: dict | None = None) -> dict:
    """Compatibility wrapper for callers that still use the old function name."""
    return read_case_defaults(case_name, overrides=overrides)


def _normalize_case_configs(request: dict) -> tuple[list[str], list[dict], dict[str, dict]]:
    """Return canonical per-case tuning configs and full JSON-backed defaults."""
    raw_case_configs = request.get("case_configs")
    raw_case_overrides = request.get("case_overrides", {}) or {}
    if not isinstance(raw_case_overrides, dict):
        raise RuntimeError("Tuning request case_overrides must be an object")

    legacy_window_count = _legacy_window_count(request)
    raw_entries = []
    if raw_case_configs is not None:
        if not isinstance(raw_case_configs, list):
            raise RuntimeError("Tuning request case_configs must be a list")
        raw_entries = list(raw_case_configs)
    else:
        raw_cases = request.get("cases")
        if raw_cases is None:
            case_name = request.get("case_name")
            if not case_name:
                raise RuntimeError("Tuning request is missing case_name, cases, or case_configs")
            raw_cases = [case_name]
        elif isinstance(raw_cases, str):
            raw_cases = [raw_cases]
        for raw_case in raw_cases or []:
            case_name = str(raw_case).strip()
            if not case_name:
                continue
            raw_entry = {"case_name": case_name}
            if case_name in raw_case_overrides:
                raw_entry.update(raw_case_overrides[case_name] or {})
            if "num_time_windows" not in raw_entry:
                raw_entry["num_time_windows"] = legacy_window_count
            raw_entries.append(raw_entry)

    cases = []
    case_configs = []
    case_defaults = {}
    seen_cases = set()
    for raw_entry in raw_entries:
        if isinstance(raw_entry, str):
            raw_entry = {"case_name": raw_entry}
        if not isinstance(raw_entry, dict):
            raise RuntimeError("Each case config must be an object")
        case_name = str(raw_entry.get("case_name", "")).strip()
        if not case_name:
            raise RuntimeError("Each case config must include case_name")
        if case_name in seen_cases:
            raise RuntimeError(f"Duplicate tuning case: {case_name}")

        override = {
            key: raw_entry[key]
            for key in (
                "altitude_comparison_range",
                "time_average_range",
                "average_time_seconds",
                "num_time_windows",
            )
            if key in raw_entry
        }
        unsupported_keys = sorted(set(raw_entry) - {"case_name", *override.keys()})
        if unsupported_keys:
            raise RuntimeError(
                f"Unsupported case config key(s) for {case_name}: " + ", ".join(unsupported_keys)
            )

        defaults = read_case_tuner_defaults(case_name, overrides=override)
        case_defaults[case_name] = defaults
        normalized_config = {
            "case_name": case_name,
            "altitude_comparison_range": list(defaults["altitude_comparison_range"]),
            "time_average_range": list(defaults["time_average_range"]),
            "num_time_windows": int(defaults.get("num_time_windows", 1)),
        }
        if defaults.get("average_time_seconds") is not None:
            normalized_config["average_time_seconds"] = int(defaults["average_time_seconds"])
        case_configs.append(normalized_config)
        cases.append(case_name)
        seen_cases.add(case_name)

    if not cases:
        raise RuntimeError("Tuning request must contain at least one case")
    return cases, case_configs, case_defaults


def _legacy_window_count(request: dict) -> int:
    raw_count = request.get("num_time_windows", DEFAULT_NUM_TIME_WINDOWS)
    raw_mode = str(request.get("time_window_mode", "") or "").strip()
    if raw_mode == "single_average":
        return 1
    if raw_mode and raw_mode != "split_average":
        raise RuntimeError(f"Unknown time_window_mode: {raw_mode}")
    return _positive_int(raw_count, "Tuning request num_time_windows")


def _positive_int(raw_value, label: str) -> int:
    try:
        value = int(raw_value)
    except (TypeError, ValueError):
        raise RuntimeError(f"{label} must be an integer")
    if value < 1:
        raise RuntimeError(f"{label} must be >= 1")
    return value


def _normalize_choice(raw_value, valid_values: set[str], label: str) -> str:
    value = str(raw_value or "").strip()
    if value not in valid_values:
        raise RuntimeError(f"Unknown {label}: {value}")
    return value


def _normalize_weights(raw_weights: dict | None, valid_names: set[str], label: str) -> dict[str, float]:
    weights = {}
    for name, raw_value in (raw_weights or {}).items():
        name = str(name).strip()
        if name not in valid_names:
            raise RuntimeError(f"Unknown {label} weight key: {name}")
        try:
            value = float(raw_value)
        except (TypeError, ValueError):
            raise RuntimeError(f"{label} weight for {name} must be numeric")
        if value < 0.0:
            raise RuntimeError(f"{label} weight for {name} must be >= 0")
        weights[name] = value
    return weights
