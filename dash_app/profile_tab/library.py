"""Profile-library loading, interchange, and comparison metadata."""

from __future__ import annotations

import base64
import binascii
import shutil
from pathlib import Path
from typing import Any, Iterable, Sequence

from utilities.timing_profiles import (
    ProfileFormatError,
    discover_profiles,
    export_profiles,
    import_profiles,
    load_profiles,
    read_profile_manifest,
)


REPO_ROOT = Path(__file__).resolve().parents[2]


def resolve_library_path(value: Any) -> Path:
    text = str(value or "output/timing").strip()
    path = Path(text).expanduser()
    if not path.is_absolute():
        path = REPO_ROOT / path
    return path.resolve()


def profile_option(record: dict[str, Any]) -> dict[str, str]:
    run_id = str(record.get("run_id") or "")
    label = str(record.get("label") or run_id)
    case_name = str(record.get("case_name") or "unknown case")
    return {
        "label": f"{label} · {case_name}",
        "value": run_id,
    }


def discover_profile_library(value: Any) -> list[dict[str, Any]]:
    return discover_profiles(resolve_library_path(value))


def delete_profile_library_entry(output_value: Any, run_id: str) -> Path:
    """Delete one recognized profile below the selected library directory."""
    root = resolve_library_path(output_value)
    record = next(
        (
            item
            for item in discover_profiles(root)
            if str(item.get("run_id") or "") == str(run_id or "")
        ),
        None,
    )
    if record is None:
        raise ValueError(f"Profile '{run_id}' is not in the selected library")

    candidate = Path(str(record.get("bundle_path") or ""))
    if candidate.is_symlink():
        raise ValueError(f"Refusing to delete symlinked profile directory: {candidate}")
    target = candidate.resolve()
    if target == root or root not in target.parents:
        raise ValueError("Only profile directories inside the selected library can be deleted")
    if read_profile_manifest(target) is None:
        raise ValueError(f"Refusing to delete unrecognized profile directory: {target}")

    shutil.rmtree(target)
    return target


def load_profile_library_data(
    value: Any,
    run_ids: Sequence[str],
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    return load_profiles(resolve_library_path(value), run_ids)


def decode_upload(contents: str | None) -> bytes:
    if not contents or "," not in contents:
        raise ProfileFormatError("uploaded profile is empty")
    _header, encoded = contents.split(",", 1)
    try:
        return base64.b64decode(encoded, validate=True)
    except (binascii.Error, ValueError) as exc:
        raise ProfileFormatError("uploaded profile is not valid base64 data") from exc


def import_profile_upload(
    output_value: Any,
    file_name: str | None,
    contents: str | None,
) -> list[str]:
    name = str(file_name or "profile-upload")
    return import_profiles(resolve_library_path(output_value), name, decode_upload(contents))


def export_profile_library(output_value: Any, run_ids: Sequence[str]) -> bytes:
    return export_profiles(resolve_library_path(output_value), run_ids)


def enrich_summary_rows(rows: Iterable[dict[str, Any]]) -> list[dict[str, Any]]:
    """Attach comparison-only throughput metrics without changing the CSV schema."""
    prepared = [dict(row) for row in rows]
    step_counts: dict[tuple[str, str, str, str], float] = {}
    for row in prepared:
        if str(row.get("timer_name") or "") != "advance_clubb_core_api":
            continue
        key = (
            str(row.get("profile_id") or row.get("run_id") or ""),
            str(row.get("process_count") or ""),
            str(row.get("columns_per_process") or ""),
            str(row.get("repetition") or ""),
        )
        try:
            step_counts[key] = float(row.get("calls") or 0)
        except (TypeError, ValueError):
            continue

    for row in prepared:
        key = (
            str(row.get("profile_id") or row.get("run_id") or ""),
            str(row.get("process_count") or ""),
            str(row.get("columns_per_process") or ""),
            str(row.get("repetition") or ""),
        )
        steps = step_counts.get(key)
        row["model_steps"] = steps if steps else ""
        try:
            total_columns = float(row.get("total_columns") or 0)
            vertical_levels = float(row.get("vertical_levels") or 0)
            timer_seconds = float(row.get("timer_max_seconds") or 0)
            timer_mean = float(row.get("timer_mean_seconds") or 0)
        except (TypeError, ValueError):
            continue
        row["throughput_column_steps_per_second"] = (
            total_columns * steps / timer_seconds
            if steps and timer_seconds > 0
            else ""
        )
        row["throughput_grid_box_iterations_per_second"] = (
            vertical_levels * total_columns * steps / timer_seconds
            if vertical_levels > 0 and steps and timer_seconds > 0
            else ""
        )
        row["process_imbalance_ratio"] = (
            timer_seconds / timer_mean if timer_mean > 0 else ""
        )
    return prepared
