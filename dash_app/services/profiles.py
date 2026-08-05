"""Shared profile selection and figure-context service.

This module is deliberately UI-neutral: it supplies the exact data selection
contract consumed by both Plot-tab callbacks and immutable agent artifacts.
Plot cards remain responsible for presentation only.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from utilities.output_paths import OUTPUT_ROOT
from dash_app.shared.output_discovery import discover_stats_directories
from dash_app.shared.netcdf import file_signature
from dash_app.plot_tab.plot_types.shared import (
    build_case_data,
    scan_output_cases,
    snap_start_time_to_step,
    time_start_max_for_duration,
)

_MAX_AVERAGE_LENGTH_MINUTES = 240.0
_MAX_OUTPUT_DIRECTORY_SCAN = 1000
_MAX_DISCOVERED_OUTPUT_DIRECTORIES = 150
_EXCLUDED_DISCOVERY_ROOTS = frozenset({"agent_artifacts"})


def _average_length_bounds(case_data: dict[str, Any]) -> tuple[float, float, float]:
    minimum = max(1.0e-6, float((case_data or {}).get("time_slider_duration_min_minutes") or 1.0))
    raw_maximum = max(minimum, float((case_data or {}).get("time_slider_duration_max_minutes") or minimum))
    maximum = max(minimum, min(raw_maximum, _MAX_AVERAGE_LENGTH_MINUTES))
    step = max(1.0e-6, float((case_data or {}).get("time_slider_duration_step_minutes") or minimum))
    return minimum, maximum, step


def _default_average_length(minimum: float, _maximum: float) -> float:
    return minimum


def scan_case_outputs(output_dirs: list[str] | None = None) -> dict[str, list[str]]:
    """Discover profile-capable output cases with the same cache as Plot."""
    return scan_output_cases(output_dirs)


def _output_label(path: Path, output_root: Path) -> str:
    try:
        relative = path.relative_to(output_root).as_posix()
    except ValueError:
        return str(path)
    return "output" if not relative else f"output/{relative}"


def discover_output_directories(
    root: Path | None = None,
    selected_dirs: list[str] | None = None,
) -> list[dict[str, Any]]:
    """Find output subdirectories that directly contain readable stats files.

    Discovery is intentionally shallow and bounded enough for an interactive
    Dash refresh: it only examines directories beneath ``output/`` and checks
    a file signature only after finding a plausible ``*_stats.nc`` filename.
    """
    output_root = Path(root or OUTPUT_ROOT).expanduser().resolve()
    records = discover_stats_directories(
        [output_root],
        max_depth=5,
        exclude_dir_names=_EXCLUDED_DISCOVERY_ROOTS,
        max_directories=_MAX_DISCOVERED_OUTPUT_DIRECTORIES,
        max_scanned_directories=_MAX_OUTPUT_DIRECTORY_SCAN,
    )
    known_paths = {str(record["path"]) for record in records}
    for raw_path in selected_dirs or []:
        selected_path = Path(raw_path).expanduser().resolve()
        selected_text = str(selected_path)
        if selected_text in known_paths:
            continue
        direct_records = discover_stats_directories(
            [selected_path],
            recursive=False,
            max_depth=0,
            max_directories=1,
            max_scanned_directories=1,
        )
        if direct_records:
            record = direct_records[0]
        else:
            record = {
                "path": selected_text,
                "key": selected_text,
                "root": selected_text,
                "relative_path": "",
                "case_names": [],
                "cases": {},
                "case_count": 0,
                "case_fingerprints": {},
                "latest_stats_modified_ns": 0,
                "modified": 0.0,
                "available": False,
            }
        records.append(record)
        known_paths.add(selected_text)

    # Preserve the established Plot-picker label contract beneath output/;
    # advanced external folders retain their full path for clarity.
    for record in records:
        path = Path(str(record["path"])).resolve()
        label = _output_label(path, output_root)
        record["relative_path"] = label
        record["label"] = label
        record["catalog_origin"] = "external" if label.startswith("/") else "output"
    records.sort(key=lambda item: (-float(item.get("modified") or 0.0), str(item["label"])))
    newest_available_assigned = False
    for index, record in enumerate(records):
        record["recency_index"] = index
        record["is_newest"] = bool(record.get("available")) and not newest_available_assigned
        newest_available_assigned = newest_available_assigned or bool(record["is_newest"])
    return records


def stats_file_fingerprints(files: list[str] | None) -> list[dict[str, Any]]:
    """Return ordered, JSON-serializable fingerprints for Plot stats files."""
    fingerprints = []
    for path in files or []:
        signature_path, modified_ns, size_bytes = file_signature(path)
        fingerprints.append(
            {
                "path": str(signature_path),
                "mtime_ns": modified_ns,
                "size_bytes": size_bytes,
            }
        )
    return fingerprints


def build_case_metadata(case_name: str, files: list[str], output_dirs: list[str] | None = None) -> dict[str, Any]:
    """Build the canonical case metadata consumed by Plot controls and export."""
    metadata = build_case_data(case_name, files, output_dirs)
    metadata["stats_fingerprints"] = stats_file_fingerprints(files)
    return metadata


def load_case_data(repo_root: Path, case_name: str, *, required: bool, output_dirs: list[str] | None = None) -> dict[str, Any] | None:
    locations = output_dirs or [str(repo_root / "output")]
    cases = scan_case_outputs(locations)
    files = cases.get(case_name)
    if not files:
        if required:
            raise ValueError(f"no readable {case_name}_stats.nc file is available in the requested output location")
        return None
    return build_case_metadata(case_name, files, locations)


def resolve_time_window(
    case_data: dict[str, Any],
    *,
    start_seconds: float | None = None,
    average_minutes: float | None = None,
) -> dict[str, float]:
    """Validate, snap, and report a physical Plot-tab averaging interval."""
    minimum, maximum, _step = _average_length_bounds(case_data)
    duration = float(average_minutes if average_minutes is not None else _default_average_length(minimum, maximum))
    if not minimum <= duration <= maximum:
        raise ValueError(f"average_minutes must be between {minimum:g} and {maximum:g}")
    start_minimum = float(case_data.get("time_slider_start_min_seconds") or 0.0)
    start_maximum = time_start_max_for_duration(case_data, duration)
    start = float(start_seconds if start_seconds is not None else case_data.get("default_time_start_seconds") or start_minimum)
    if not start_minimum <= start <= start_maximum:
        raise ValueError(f"time_start_seconds must be between {start_minimum:g} and {start_maximum:g} for a {duration:g}-minute average")
    return {"time_start_seconds": float(snap_start_time_to_step(case_data, start, duration)), "average_minutes": duration}


def validate_profile_names(variables: list[str], case_data: dict[str, Any] | None) -> list[str]:
    names = [str(value).strip() for value in (variables or []) if str(value).strip()]
    if not names:
        raise ValueError("provide at least one profile variable")
    if case_data is not None:
        available = {str(item.get("value") or "") for item in (case_data.get("profile_vars") or [])}
        unknown = sorted(set(names) - available)
        if unknown:
            raise ValueError("profile variable(s) unavailable for this output: " + ", ".join(unknown))
    return list(dict.fromkeys(names))


def figure_context(case_data: dict[str, Any], selection: dict[str, Any]) -> dict[str, Any]:
    """Return the non-interactive context used by the native Profile figure."""
    minimum, maximum, _step = _average_length_bounds(case_data)
    duration = float(selection.get("average_minutes") or _default_average_length(minimum, maximum))
    start_minimum = float(case_data.get("time_slider_start_min_seconds") or 0.0)
    start = float(selection.get("time_start_seconds") or case_data.get("default_time_start_seconds") or start_minimum)
    return {
        "case_data": case_data,
        "enabled_benchmarks": [],
        "time_range": duration,
        "time_point": start,
        "height_range": case_data.get("default_height_range"),
        "relayout_data": None,
        "use_relayout_height_range": False,
        "selected_column": 0,
        "column_mode": "single",
        "column_filter_indices": None,
        "size": "normal",
        "theme_name": "dark",
    }
