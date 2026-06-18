"""Case-specific benchmark defaults for tuner and loss-driver workflows."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Mapping


CASE_DEFAULTS_PATH = Path(__file__).with_name("case_defaults.json")
OVERRIDABLE_KEYS = {
    "altitude_comparison_range",
    "time_average_range",
    "average_time_seconds",
    "num_time_windows",
}
DEFAULT_LOSS_FIELDS = ["cloud_frac"]


def load_case_defaults() -> dict[str, dict]:
    """Load and validate all case comparison defaults."""
    try:
        raw_defaults = json.loads(CASE_DEFAULTS_PATH.read_text(encoding="utf-8"))
    except OSError as exc:
        raise RuntimeError(f"Could not read case defaults: {CASE_DEFAULTS_PATH}") from exc
    except json.JSONDecodeError as exc:
        raise RuntimeError(f"Could not parse case defaults: {CASE_DEFAULTS_PATH}") from exc

    if not isinstance(raw_defaults, dict):
        raise RuntimeError(f"{CASE_DEFAULTS_PATH} must contain a JSON object")

    defaults = {}
    for case_name, raw_value in raw_defaults.items():
        defaults[str(case_name)] = _normalize_case_defaults(str(case_name), raw_value)
    return defaults


def available_case_names() -> list[str]:
    """Return case names that have benchmark-comparison defaults."""
    return sorted(load_case_defaults())


def read_case_defaults(case_name: str, overrides: Mapping | None = None) -> dict:
    """Return one case's defaults, applying explicit caller overrides."""
    all_defaults = load_case_defaults()
    if case_name not in all_defaults:
        raise RuntimeError(f"Missing tuner defaults for case {case_name}: {CASE_DEFAULTS_PATH}")
    defaults = dict(all_defaults[case_name])
    return apply_case_overrides(case_name, defaults, overrides)


def apply_case_overrides(case_name: str, defaults: dict, overrides: Mapping | None = None) -> dict:
    """Apply supported per-case comparison overrides to a defaults record."""
    if not overrides:
        return _normalize_case_defaults(case_name, defaults)
    if not isinstance(overrides, Mapping):
        raise RuntimeError(f"case_overrides for {case_name} must be an object")

    unknown_keys = sorted(set(overrides) - OVERRIDABLE_KEYS)
    if unknown_keys:
        raise RuntimeError(
            f"Unsupported case override key(s) for {case_name}: " + ", ".join(unknown_keys)
        )

    merged = dict(defaults)
    for key in OVERRIDABLE_KEYS:
        if key in overrides:
            merged[key] = overrides[key]
    return _normalize_case_defaults(case_name, merged)


def _normalize_case_defaults(case_name: str, raw_value) -> dict:
    if not isinstance(raw_value, Mapping):
        raise RuntimeError(f"Case defaults for {case_name} must be an object")

    les_stats_file = str(raw_value.get("les_stats_file", "")).strip()
    if not les_stats_file:
        raise RuntimeError(f"Case defaults for {case_name} require les_stats_file")

    altitude_range = _normalize_float_pair(
        raw_value.get("altitude_comparison_range"),
        f"altitude_comparison_range for {case_name}",
    )
    if altitude_range[1] < altitude_range[0]:
        raise RuntimeError(f"altitude_comparison_range for {case_name} requires min <= max")

    time_range = _normalize_int_pair(
        raw_value.get("time_average_range"),
        f"time_average_range for {case_name}",
    )
    if time_range[1] <= time_range[0]:
        raise RuntimeError(f"time_average_range for {case_name} requires end > start")

    average_time_seconds = None
    if raw_value.get("average_time_seconds") is not None:
        average_time_seconds = _positive_int(
            raw_value.get("average_time_seconds"),
            f"average_time_seconds for {case_name}",
        )
        window_seconds = time_range[1] - time_range[0]
        if window_seconds % average_time_seconds != 0:
            raise RuntimeError(
                f"average_time_seconds for {case_name} must divide time_average_range evenly"
            )
        num_time_windows = window_seconds // average_time_seconds
    else:
        num_time_windows = _positive_int(
            raw_value.get("num_time_windows", 1),
            f"num_time_windows for {case_name}",
        )

    normalized = {
        "les_stats_file": les_stats_file,
        "altitude_comparison_range": altitude_range,
        "time_average_range": time_range,
        "num_time_windows": num_time_windows,
    }
    if average_time_seconds is not None:
        normalized["average_time_seconds"] = average_time_seconds
    return normalized


def _normalize_float_pair(raw_value, label: str) -> list[float]:
    if not isinstance(raw_value, (list, tuple)) or len(raw_value) != 2:
        raise RuntimeError(f"{label} must be a two-value list")
    try:
        return [float(raw_value[0]), float(raw_value[1])]
    except (TypeError, ValueError) as exc:
        raise RuntimeError(f"{label} must contain numeric values") from exc


def _normalize_int_pair(raw_value, label: str) -> list[int]:
    if not isinstance(raw_value, (list, tuple)) or len(raw_value) != 2:
        raise RuntimeError(f"{label} must be a two-value list")
    try:
        return [int(raw_value[0]), int(raw_value[1])]
    except (TypeError, ValueError) as exc:
        raise RuntimeError(f"{label} must contain integer second values") from exc


def _positive_int(raw_value, label: str) -> int:
    try:
        value = int(raw_value)
    except (TypeError, ValueError) as exc:
        raise RuntimeError(f"{label} must be an integer") from exc
    if value < 1:
        raise RuntimeError(f"{label} must be >= 1")
    return value
