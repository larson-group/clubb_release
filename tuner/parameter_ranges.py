"""Derive Tune intervals from checked-in Fortran bounds and config defaults.

The Fortran ``parameter_hard_bounds`` table is the runtime validity oracle.
Its checked-in Python mirror is validated against the compiled table in
Jenkins, which lets Dash construct Tune controls without loading F2PY at
startup. Tune's historical quarter-to-four-times-default interval is retained
only for an *open* side of that table.
"""

from __future__ import annotations

import math
from pathlib import Path
import re

from utilities.create_case_namelist import resolve_tunable_config_dir
from utilities.clubb_settings_validation import default_range_for_value, parameter_hard_bounds


_ASSIGNMENT = re.compile(r"^\s*([A-Za-z][A-Za-z0-9_]*)\s*=\s*([^! ,/]+)")


def _numeric(value: str) -> float | None:
    try:
        parsed = float(str(value).replace("D", "E").replace("d", "e"))
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def tunable_parameter_defaults(config: str | None = None) -> dict[str, float]:
    """Read scalar default values from a selected tunable namelist."""
    path = Path(resolve_tunable_config_dir(config or "default")) / "tunable_parameters.in"
    defaults: dict[str, float] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.lstrip().startswith("!"):
            continue
        match = _ASSIGNMENT.match(line)
        if not match:
            continue
        value = _numeric(match.group(2))
        if value is not None:
            defaults[match.group(1)] = value
    return defaults


def compiled_parameter_hard_bounds() -> dict[str, dict[str, float | None]]:
    """Return Fortran-validated hard bounds without importing the F2PY API.

    The public function name is retained for callers that previously used the
    compiled lookup. ``tests/test_parameter_hard_bounds_parity.py`` verifies
    this mirror agrees with the compiled Fortran table.
    """
    return parameter_hard_bounds()


def default_range_for_parameter(
    default: float,
    hard_bounds: dict[str, float | None] | None,
) -> tuple[float, float]:
    """Compatibility wrapper around the shared settings range policy."""
    return default_range_for_value(default, hard_bounds)


def default_parameter_ranges(config: str | None = None) -> dict[str, dict[str, float]]:
    """Return UI-ready default intervals for every scalar tunable parameter."""
    bounds = compiled_parameter_hard_bounds()
    result: dict[str, dict[str, float]] = {}
    for name, default in tunable_parameter_defaults(config).items():
        try:
            low, high = default_range_for_parameter(default, bounds.get(name))
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
