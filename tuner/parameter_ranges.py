"""Derive Tune intervals from checked-in Fortran bounds and config defaults.

The Fortran ``parameter_hard_bounds`` table is the runtime validity oracle.
Its checked-in Python mirror is validated against the compiled table in
Jenkins, which lets Dash construct Tune controls without loading F2PY at
startup. Open sides use a quarter-to-four-times-default interval suggestion.
"""

from __future__ import annotations

from input.scripts.parameters_json_to_namelist import parameter_settings
from utilities.clubb_settings_validation import default_range_for_value, parameter_hard_bounds


def tunable_parameter_defaults(config: str | None = None) -> dict[str, float]:
    """Read scalar default values from a selected tunable namelist."""
    return {name: float(value) for name, value in parameter_settings(config=config)["model"].items()
            if not isinstance(value, list)}


def default_parameter_ranges(config: str | None = None) -> dict[str, dict[str, float]]:
    """Return UI-ready default intervals for every scalar tunable parameter."""
    bounds = parameter_hard_bounds()
    result: dict[str, dict[str, float]] = {}
    for name, default in tunable_parameter_defaults(config).items():
        try:
            low, high = default_range_for_value(default, bounds.get(name))
        except ValueError as exc:
            raise RuntimeError(f"Could not derive Tune range for {name}: {exc}") from exc
        result[name] = {"default": default, "min": low, "max": high}
    return result


def resolve_parameter_specs(
    specs: list[dict],
    default_ranges: dict[str, dict[str, float]],
) -> list[dict]:
    """Fill omitted interval endpoints and intersect linked physical bounds."""
    resolved: list[dict] = []
    for raw in specs:
        spec = dict(raw or {})
        name = str(spec.get("name") or "").strip()
        raw_targets = spec.get("targets", [name])
        if isinstance(raw_targets, str):
            raw_targets = [raw_targets]
        targets = [str(target).strip() for target in raw_targets or [] if str(target).strip()]
        if not name or not targets:
            raise ValueError("Tune parameter range requires a name and physical target")

        if spec.get("min") is not None and spec.get("max") is not None:
            low, high = float(spec["min"]), float(spec["max"])
        else:
            target_ranges = []
            for target in targets:
                derived = default_ranges.get(target)
                if derived is None:
                    raise ValueError(f"No default/hard-bound metadata for preset target {target}")
                target_ranges.append(derived)
            low = max(float(item["min"]) for item in target_ranges)
            high = min(float(item["max"]) for item in target_ranges)
        if low > high:
            raise ValueError(f"Linked Tune range {name} has no shared admissible interval")
        resolved.append({"name": name, "targets": targets, "min": low, "max": high})
    return resolved
