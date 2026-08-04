"""UI-neutral Plot selection and plot-grid state transitions.

The Dash callbacks and MCP browser-handoff adapters both use this module for
the state changes that are common to the Plot tab.  Plot-family plugins remain
responsible for their own option catalogs and figure rendering.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from dash_app.plot_tab.benchmark_overlay import resolve_enabled_sources
from dash_app.plot_tab.plot_types.registry import PLOT_TYPES
from dash_app.services import profiles as profile_service


_WINDOW_PRESETS = {"loss", "pyplotgen"}
_PLOT_OPERATIONS = {"set_view", "add_budget"}


@dataclass(frozen=True)
class PlotStateTransition:
    """The complete ordered plot-store update produced by one transition."""

    plot_order: list[int]
    plot_state: dict[str, dict[str, Any]]
    next_id: int


def resolve_benchmark_sources(
    case_data: dict[str, Any] | None,
    requested_sources=None,
    *,
    strict: bool = False,
) -> list[str]:
    """Resolve common Plot overlays through the Plot service boundary."""
    return resolve_enabled_sources(case_data or {}, requested_sources, strict=strict)


def toggle_benchmark_source(case_data: dict[str, Any] | None, current_sources, source: str) -> list[str]:
    """Apply the native benchmark-button transition used by Plot handoffs."""
    current = resolve_benchmark_sources(case_data, current_sources)
    value = str(source or "").strip().lower()
    if value in current:
        current.remove(value)
    else:
        current.append(value)
    return resolve_benchmark_sources(case_data, current, strict=True)


def _option_values(case_data: dict[str, Any] | None, plot_type: str) -> list[str]:
    module = PLOT_TYPES.get(str(plot_type or ""))
    if module is None:
        raise ValueError(f"unknown plot type: {plot_type}")
    return [
        str(option.get("value"))
        for option in module.case_data_options(case_data or {})
        if option.get("value") is not None
    ]


def available_plot_values(case_data: dict[str, Any] | None, plot_type: str) -> list[str]:
    """Return the validated option values exposed by one registered family."""
    return _option_values(case_data, plot_type)


def validate_budget_group(case_data: dict[str, Any] | None, group: str | None = None) -> str:
    """Validate one budget group against the current case's available groups."""
    available = _option_values(case_data, "budget")
    if not available:
        raise ValueError("no budget groups are available for this output")
    requested = str(group or "wp2").strip()
    if requested not in available:
        raise ValueError(
            f"budget group '{requested}' is unavailable; available groups: {', '.join(available)}"
        )
    return requested


def _validate_common_selection(case_data: dict[str, Any], request: dict[str, Any]) -> dict[str, Any]:
    """Normalize time and benchmark fields shared by Plot operations."""
    normalized: dict[str, Any] = {}
    requested_start = request.get("time_start_seconds")
    legacy_start = request.get("time_seconds")
    if requested_start is not None and legacy_start is not None:
        if abs(float(requested_start) - float(legacy_start)) > 1.0e-9:
            raise ValueError("time_seconds and time_start_seconds disagree; use time_start_seconds")
    if requested_start is None:
        requested_start = legacy_start

    preset = str(request.get("window_preset") or "").strip().lower()
    if preset and preset not in _WINDOW_PRESETS:
        raise ValueError("window_preset must be 'loss' or 'pyplotgen'")
    if preset:
        start = case_data.get(f"{preset}_time_start_seconds")
        duration = case_data.get(f"{preset}_time_duration_minutes")
        if start is None or duration is None:
            raise ValueError(f"{preset} window is not available for {case_data.get('name') or 'this case'}")
        if requested_start is not None or request.get("average_minutes") is not None:
            # Published handoffs carry the exact physical values alongside
            # the symbolic preset so the browser can display both.  Accept
            # that canonical pair, but reject a conflicting custom window.
            if (
                requested_start is None
                or request.get("average_minutes") is None
                or abs(float(requested_start) - float(start)) > 1.0e-9
                or abs(float(request.get("average_minutes")) - float(duration)) > 1.0e-9
            ):
                raise ValueError("window_preset cannot be combined with custom start or average values")
        normalized.update(
            {
                "time_start_seconds": float(start),
                "average_minutes": float(duration),
                "window_preset": preset,
            }
        )
    elif requested_start is not None or request.get("average_minutes") is not None:
        normalized.update(
            profile_service.resolve_time_window(
                case_data,
                start_seconds=float(requested_start) if requested_start is not None else None,
                average_minutes=(
                    float(request.get("average_minutes"))
                    if request.get("average_minutes") is not None
                    else None
                ),
            )
        )

    if "benchmark_sources" in request:
        normalized["benchmark_sources"] = resolve_benchmark_sources(
            case_data,
            request.get("benchmark_sources"),
            strict=True,
        )
    return normalized


def validate_plot_request(case_data: dict[str, Any] | None, request: dict[str, Any]) -> dict[str, Any]:
    """Validate a typed Plot request without referring to Dash component IDs."""
    data = dict(case_data or {})
    operation = str(request.get("operation") or "").strip()
    if operation not in _PLOT_OPERATIONS:
        raise ValueError(f"unsupported Plot operation: {operation or 'missing'}")
    requested_case = str(request.get("case") or "").strip()
    case_name = str(data.get("name") or "").strip()
    if requested_case and case_name and requested_case != case_name:
        raise ValueError(f"Plot request case '{requested_case}' does not match loaded case '{case_name}'")

    normalized = {
        "operation": operation,
        "case": requested_case or case_name,
    }
    if isinstance(request.get("_normalized_common"), dict):
        # Existing profile selection already goes through the legacy-compatible
        # profile service (including immutable-run handling).  Reuse that
        # canonical result while still applying this service's family and
        # benchmark validation.
        normalized.update(dict(request["_normalized_common"]))
        if "benchmark_sources" in request:
            normalized["benchmark_sources"] = resolve_benchmark_sources(
                data,
                request.get("benchmark_sources"),
                strict=True,
            )
    else:
        normalized.update(_validate_common_selection(data, request))
    if request.get("output_dir"):
        normalized["output_dir"] = str(request["output_dir"])
    if request.get("run_id"):
        normalized["run_id"] = str(request["run_id"])

    if operation == "set_view":
        normalized["variables"] = profile_service.validate_profile_names(
            request.get("variables") or [],
            data,
        )
    else:
        normalized["budget_group"] = validate_budget_group(data, request.get("budget_group"))
    return normalized


def next_plot_id(plot_order, plot_state, next_id) -> int:
    """Choose a plot id not already mounted in either persisted store."""
    used = {
        int(value)
        for value in list(plot_order or []) + list((plot_state or {}).keys())
        if str(value).lstrip("-").isdigit()
    }
    candidate = max(int(next_id or 0), 0)
    while candidate in used:
        candidate += 1
    return candidate


def normalize_plot_id(value) -> int:
    """Normalize one externally supplied Plot-card ID."""
    if isinstance(value, bool):
        raise ValueError("plot_id must be an integer")
    try:
        normalized = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError("plot_id must be an integer") from exc
    if normalized < 0 or str(normalized) != str(value).strip() and not isinstance(value, int):
        raise ValueError("plot_id must be a non-negative integer")
    return normalized


def list_plot_instances(plot_order, plot_state) -> list[dict[str, Any]]:
    """Describe only the currently mounted, owned Plot cards.

    The returned records intentionally contain stable card IDs and safe
    family/selection metadata, not arbitrary Dash component state.
    """
    instances: list[dict[str, Any]] = []
    seen: set[int] = set()
    state = plot_state or {}
    for raw_id in plot_order or []:
        try:
            plot_id = normalize_plot_id(raw_id)
        except ValueError:
            continue
        if plot_id in seen:
            continue
        entry = state.get(str(plot_id))
        if not isinstance(entry, dict):
            continue
        plot_type = str(entry.get("plot_type") or "").strip()
        if plot_type not in PLOT_TYPES:
            continue
        seen.add(plot_id)
        instances.append(
            {
                "id": plot_id,
                "plot_type": plot_type,
                "selection": str(entry.get("var") or ""),
                "size": str(entry.get("size") or "normal"),
            }
        )
    return instances


def remove_plot_instance(plot_id, plot_order, plot_state, next_id) -> PlotStateTransition:
    """Remove one card only when it belongs to the current Plot stores."""
    normalized_id = normalize_plot_id(plot_id)
    order = [int(value) for value in (plot_order or [])]
    state = {str(key): dict(value or {}) for key, value in (plot_state or {}).items()}
    if normalized_id not in order or str(normalized_id) not in state:
        raise ValueError(f"plot_id {normalized_id} is not an active Plot card owned by this dashboard")
    order = [value for value in order if value != normalized_id]
    state.pop(str(normalized_id), None)
    return PlotStateTransition(order, state, int(next_id or 0))


def append_plot_instance(
    case_data: dict[str, Any] | None,
    plot_order,
    plot_state,
    next_id,
    *,
    plot_type: str,
    state: dict[str, Any] | None = None,
    variable: str | None = None,
) -> PlotStateTransition:
    """Append one validated plot-family instance to the ordered grid."""
    module = PLOT_TYPES.get(str(plot_type or ""))
    if module is None:
        raise ValueError(f"unknown plot type: {plot_type}")
    options = _option_values(case_data, module.plot_type_id)
    if not options:
        raise ValueError(f"plot type '{module.plot_type_id}' is unavailable for this output")
    plot_id = next_plot_id(plot_order, plot_state, next_id)
    entry = dict(state or module.make_default_state(case_data or {}, plot_id))
    entry["plot_type"] = module.plot_type_id
    selected = str(variable if variable is not None else entry.get("var") or "").strip()
    if selected not in options:
        raise ValueError(
            f"{module.plot_type_id} selection '{selected}' is unavailable; available values: {', '.join(options)}"
        )
    entry["var"] = selected
    entry["size"] = entry.get("size") or "normal"
    updated_order = [int(value) for value in (plot_order or [])]
    updated_state = {str(key): dict(value or {}) for key, value in (plot_state or {}).items()}
    updated_order.append(plot_id)
    updated_state[str(plot_id)] = entry
    return PlotStateTransition(updated_order, updated_state, plot_id + 1)


def apply_plot_request(
    case_data: dict[str, Any] | None,
    request: dict[str, Any],
    plot_order,
    plot_state,
    next_id,
) -> PlotStateTransition:
    """Apply a validated MCP/UI request to the Plot-tab stores."""
    normalized = validate_plot_request(case_data, request)
    if normalized["operation"] == "set_view":
        order = list(range(len(normalized["variables"])))
        state = {
            str(plot_id): {"plot_type": "profile", "var": var_name, "size": "normal"}
            for plot_id, var_name in enumerate(normalized["variables"])
        }
        return PlotStateTransition(order, state, len(order))
    return append_plot_instance(
        case_data,
        plot_order,
        plot_state,
        next_id,
        plot_type="budget",
        variable=normalized["budget_group"],
    )
