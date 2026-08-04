"""Shared, validated Tune experiment presets.

The JSON file is deliberately the single source of truth for the CLI, native
Dash controls, and durable tuner requests.  A preset is an experiment
baseline; explicit caller values may subsequently replace its individual
pieces before normal request validation runs.
"""

from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path

from tuner.parameter_ranges import default_parameter_ranges, resolve_parameter_specs


_PRESETS_PATH = Path(__file__).with_name("presets.json")


def _normalise_range(raw: dict, *, preset_name: str) -> dict:
    if not isinstance(raw, dict):
        raise ValueError(f"Preset '{preset_name}' has a non-object parameter range")
    name = str(raw.get("name", "")).strip()
    targets = raw.get("targets", [name])
    if isinstance(targets, str):
        targets = [targets]
    targets = [str(target).strip() for target in targets or [] if str(target).strip()]
    if not name or not targets:
        raise ValueError(f"Preset '{preset_name}' has a range without a name or target")
    # Presets intentionally name the coordinates but do not duplicate validity
    # intervals.  ``parameter_hard_bounds`` in compiled Fortran plus the
    # selected configuration's defaults resolve those immediately before the
    # request reaches the tuner.
    return {"name": name, "targets": targets}


def _validate_preset(raw: dict) -> dict:
    if not isinstance(raw, dict):
        raise ValueError("Each tuner preset must be an object")
    name = str(raw.get("name", "")).strip()
    if not name:
        raise ValueError("Each tuner preset must have a name")
    case_configs = [dict(item) for item in raw.get("case_configs", []) if isinstance(item, dict)]
    if not case_configs or any(not str(item.get("case_name", "")).strip() for item in case_configs):
        raise ValueError(f"Preset '{name}' must provide case_configs with case_name values")
    fields = [str(field).strip() for field in raw.get("selected_fields", []) if str(field).strip()]
    if not fields:
        raise ValueError(f"Preset '{name}' must provide selected_fields")
    ranges = [_normalise_range(item, preset_name=name) for item in raw.get("parameter_ranges", [])]
    if not ranges:
        raise ValueError(f"Preset '{name}' must provide parameter_ranges")
    names = [item["name"] for item in ranges]
    if len(names) != len(set(names)):
        raise ValueError(f"Preset '{name}' duplicates a logical range name")
    targets = [target for item in ranges for target in item["targets"]]
    if len(targets) != len(set(targets)):
        raise ValueError(f"Preset '{name}' assigns a physical target more than once")
    return {
        "name": name,
        "label": str(raw.get("label") or name),
        "description": str(raw.get("description") or ""),
        "case_configs": case_configs,
        "selected_fields": fields,
        "parameter_ranges": ranges,
        **({"config": str(raw["config"]).strip()} if raw.get("config") else {}),
        **({"override": str(raw["override"]).strip()} if raw.get("override") else {}),
    }


def load_presets(path: Path | None = None) -> dict[str, dict]:
    """Load all preset definitions, validating their small closed schema."""
    raw = json.loads((path or _PRESETS_PATH).read_text(encoding="utf-8"))
    if not isinstance(raw, dict) or raw.get("schema_version") != 1:
        raise ValueError("tuner preset file must use schema_version 1")
    definitions = [_validate_preset(item) for item in raw.get("presets", [])]
    names = [item["name"] for item in definitions]
    if len(names) != len(set(names)):
        raise ValueError("tuner preset file duplicates a preset name")
    return {item["name"]: item for item in definitions}


def list_presets() -> list[dict]:
    """Return compact metadata for UI and CLI discovery."""
    return [
        {"name": preset["name"], "label": preset["label"], "description": preset["description"]}
        for preset in load_presets().values()
    ]


def get_preset(name: str) -> dict:
    """Return an isolated, full definition for one named preset."""
    key = str(name or "").strip()
    try:
        return deepcopy(load_presets()[key])
    except KeyError as exc:
        choices = ", ".join(load_presets()) or "(none)"
        raise ValueError(f"Unknown Tune preset '{key}'. Available presets: {choices}") from exc


def apply_preset(request: dict, preset_name: str | None = None) -> dict:
    """Return a request baseline populated by a preset.

    Existing, explicitly supplied keys take precedence.  This lets the CLI use
    a preset alone or modify only ``-params``/``-cases``/``-fields``.
    """
    request = dict(request or {})
    name = str(preset_name or request.get("preset") or "").strip()
    if not name:
        return request
    preset = get_preset(name)
    effective_config = str(request.get("config") or preset.get("config") or "default").strip() or "default"
    preset_ranges = (
        resolve_parameter_specs(
            preset["parameter_ranges"], default_parameter_ranges(effective_config)
        )
        if request.get("parameter_ranges") is None
        else deepcopy(request["parameter_ranges"])
    )
    merged = {
        "preset": name,
        "case_configs": deepcopy(preset["case_configs"]),
        "cases": [item["case_name"] for item in preset["case_configs"]],
        "selected_fields": list(preset["selected_fields"]),
        "parameter_ranges": preset_ranges,
    }
    for key in ("config", "override"):
        if key in preset:
            merged[key] = preset[key]
    merged.update({key: value for key, value in request.items() if value is not None})
    merged["preset"] = name
    return merged
