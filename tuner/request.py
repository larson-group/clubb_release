"""Request validation and case-default parsing for tuning jobs."""

from __future__ import annotations

import json
from pathlib import Path

from tuner.case_defaults import read_case_defaults
from utilities.benchmark_converter import supported_fields
from tuner.taylor_metrics import (
    AGGREGATION_MODE_NAMES,
    DEFAULT_AGGREGATION_MODE,
    DEFAULT_LOSS_MODE,
    DEFAULT_NUM_TIME_WINDOWS,
    LOSS_MODE_NAMES,
    LOSS_POLICY_CONSTANTS,
    LOSS_POLICY_VERSION,
    WORST_QUANTILE_FRACTION,
)
from tuner.tuning_strategy import normalize_strategy_config


def load_request(request_path: Path) -> dict:
    """Load and validate the tuning request."""
    request = json.loads(request_path.read_text(encoding="utf-8"))

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
    for spec in raw_ranges:
        name = str(spec.get("name", "")).strip()
        if not name:
            raise RuntimeError("Each parameter range must include a name")
        if name in seen_params:
            raise RuntimeError(f"Duplicate tuning parameter: {name}")
        seen_params.add(name)
        try:
            min_value = float(spec.get("min"))
            max_value = float(spec.get("max"))
        except (TypeError, ValueError):
            raise RuntimeError(f"Invalid numeric range for {name}")
        if min_value > max_value:
            raise RuntimeError(f"{name} requires min <= max")
        parameter_ranges.append({"name": name, "min": min_value, "max": max_value})

    request["case_name"] = cases[0]
    request["cases"] = cases
    request["case_configs"] = case_configs
    request["case_defaults"] = case_defaults
    request["selected_fields"] = selected_fields
    request["parameter_ranges"] = parameter_ranges
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
    }

    try:
        request["strategy"] = normalize_strategy_config(request)
    except ValueError as exc:
        raise RuntimeError(str(exc)) from exc
    if request["strategy"]["name"] == "random":
        request["max_samples"] = request["strategy"]["options"].get("max_samples")
    else:
        request.pop("max_samples", None)
    request["total_samples"] = request["strategy"]["options"].get("max_samples")

    if request.get("seed") is not None:
        try:
            request["seed"] = int(request["seed"])
        except (TypeError, ValueError):
            raise RuntimeError("Tuning request seed must be an integer")
    return request


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
