"""Profile-library loading, interchange, and comparison metadata."""

from __future__ import annotations

import base64
import binascii
from pathlib import Path
from typing import Any, Iterable, Sequence

from utilities.timing_profiles import (
    ProfileFormatError,
    discover_profiles,
    export_profiles,
    import_profiles,
    load_profiles,
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
    """Attach comparison-only metrics without changing the portable CSV schema."""
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
            timer_seconds = float(row.get("timer_max_seconds") or 0)
            timer_mean = float(row.get("timer_mean_seconds") or 0)
        except (TypeError, ValueError):
            continue
        row["throughput_column_steps_per_second"] = (
            total_columns * steps / timer_seconds
            if steps and timer_seconds > 0
            else ""
        )
        row["process_imbalance_ratio"] = (
            timer_seconds / timer_mean if timer_mean > 0 else ""
        )
    return prepared


def comparison_warnings(
    catalog: Iterable[dict[str, Any]],
    selected_ids: Sequence[str],
) -> list[str]:
    """Identify provenance differences that can invalidate direct comparisons."""
    selected = set(str(value) for value in selected_ids)
    records = [row for row in catalog if str(row.get("run_id") or "") in selected]
    if len(records) < 2:
        return []
    checks = (
        ("case_name", "cases"),
        ("git_revision", "source revisions"),
        ("forwarded_args", "run arguments"),
        ("vertical_levels", "vertical level counts"),
        ("model_steps", "model-step counts"),
        ("omp_num_threads", "OpenMP thread counts"),
        ("backends", "timer backends"),
        ("time_bases", "timer time bases"),
        ("hostname", "hosts"),
        ("executable_sha256", "executables"),
    )
    warnings = []
    for field, label in checks:
        values = {str(record.get(field) or "unknown") for record in records}
        if len(values) > 1:
            warnings.append(f"Selected profiles use different {label}.")
    return warnings
