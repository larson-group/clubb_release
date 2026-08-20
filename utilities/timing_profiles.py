"""Compact, portable storage for CLUBB timing profiles.

Each profile is one self-contained directory.  Profile-wide provenance lives in
``profile.json`` and ``README.md``; workload metadata lives in ``batches.csv``;
and ``timings.csv`` contains only raw timer observations.  Dash derives summary
statistics when profiles are loaded instead of persisting duplicate summaries.
"""

from __future__ import annotations

import csv
import io
import json
import os
import re
import shutil
import stat
import statistics
import tempfile
import uuid
import zipfile
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from typing import Any, Iterable, Sequence


PROFILE_FORMAT = "clubb-timing-profile"
PROFILE_FORMAT_VERSION = 2
EXPORT_FORMAT = "clubb-timing-profile-export"
EXPORT_FORMAT_VERSION = 2

MANIFEST_FILE = "profile.json"
README_FILE = "README.md"
BATCHES_FILE = "batches.csv"
TIMINGS_FILE = "timings.csv"
LOGS_DIRECTORY = "logs"
EXPORT_MANIFEST_FILE = "clubb-profile-export.json"
GROUP_WALL_TIMER = "__process_group_wall__"
DEFAULT_ANALYSIS_TIMER = "advance_clubb_to_end"

BATCH_FIELDS = (
    "batch_id",
    "process_count",
    "batch_size",
    "total_batch_size",
    "status",
    "warmups_completed",
    "measurements_completed",
    "failed_runs",
    "representative_phase",
    "representative_repetition",
    "representative_process",
)

TIMING_FIELDS = (
    "batch_id",
    "phase",
    "repetition",
    "process",
    "thread",
    "timer",
    "calls",
    "inclusive_s",
    "exclusive_s",
    "status",
    "return_code",
    "message",
)

_SAFE_ID = re.compile(r"[^A-Za-z0-9_.-]+")


class ProfileFormatError(ValueError):
    """Raised when imported timing data is unsafe or unrecognized."""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def safe_profile_id(value: str) -> str:
    cleaned = _SAFE_ID.sub("_", str(value or "").strip()).strip("._")
    return cleaned[:120] or "profile"


def profile_directory(output_root: Path, run_id: str) -> Path:
    """Return the direct child directory used by one compact profile."""
    return Path(output_root) / run_id


def reserve_profile_directory(
    output_root: Path,
    label: str,
    *,
    overwrite: bool = False,
) -> tuple[str, Path]:
    """Reserve a readable profile directory, optionally replacing that exact profile."""
    root = Path(output_root)
    root.mkdir(parents=True, exist_ok=True)
    base = safe_profile_id(label)
    if overwrite:
        path = profile_directory(root, base)
        if path.is_symlink():
            raise ProfileFormatError(f"refusing to overwrite symlinked profile path: {path}")
        if path.exists():
            if not path.is_dir() or read_profile_manifest(path) is None:
                raise ProfileFormatError(
                    f"refusing to overwrite unrecognized profile directory: {path}"
                )
            shutil.rmtree(path)
        path.mkdir()
        return base, path

    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    for index in range(10_000):
        suffix = "" if index == 0 else f"_{stamp}" if index == 1 else f"_{stamp}_{index}"
        run_id = f"{base}{suffix}"
        path = profile_directory(root, run_id)
        try:
            path.mkdir()
        except FileExistsError:
            continue
        return run_id, path
    raise OSError(f"unable to reserve a profile directory under {root}")


def _atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent, text=True
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, path)
    except BaseException:
        try:
            Path(temporary_name).unlink()
        except OSError:
            pass
        raise


def _read_csv(path: Path) -> list[dict[str, str]]:
    try:
        with Path(path).open(newline="", encoding="utf-8") as handle:
            return list(csv.DictReader(handle))
    except OSError:
        return []


def _csv_text(rows: Sequence[dict[str, Any]], fields: Sequence[str]) -> str:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=list(fields), extrasaction="ignore")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue()


def write_batches(profile_dir: Path, rows: Sequence[dict[str, Any]]) -> None:
    """Atomically replace the small one-row-per-workload table."""
    _atomic_write(Path(profile_dir) / BATCHES_FILE, _csv_text(rows, BATCH_FIELDS))


def write_profile_manifest(profile_dir: Path, manifest: dict[str, Any]) -> None:
    normalized = dict(manifest)
    normalized["format"] = PROFILE_FORMAT
    normalized["format_version"] = PROFILE_FORMAT_VERSION
    _atomic_write(
        Path(profile_dir) / MANIFEST_FILE,
        json.dumps(normalized, indent=2, sort_keys=True) + "\n",
    )


def read_profile_manifest(profile_dir: Path) -> dict[str, Any] | None:
    try:
        value = json.loads((Path(profile_dir) / MANIFEST_FILE).read_text(encoding="utf-8"))
    except (OSError, ValueError, TypeError):
        return None
    if not isinstance(value, dict) or value.get("format") != PROFILE_FORMAT:
        return None
    return value


def write_profile_readme(profile_dir: Path, manifest: dict[str, Any]) -> None:
    """Write a concise human-readable companion to the machine metadata."""
    benchmark = dict(manifest.get("benchmark") or {})
    environment = dict(manifest.get("environment") or {})
    build = dict(manifest.get("build") or {})
    results = dict(manifest.get("results") or {})
    lines = [
        f"# {manifest.get('label') or manifest.get('run_id') or 'CLUBB timing profile'}",
        "",
        "Compact CLUBB timing profile generated by `utilities/time_clubb.py`.",
        "",
        "## Run",
        "",
        f"- Profile ID: `{manifest.get('run_id', '')}`",
        f"- Case: `{manifest.get('case_name', '')}`",
        f"- Status: `{manifest.get('status', 'unknown')}`",
        f"- Started: `{manifest.get('started_utc', '')}`",
        f"- Completed: `{manifest.get('completed_utc', '')}`",
        f"- Git revision: `{manifest.get('git_revision', '')}`",
        f"- Git worktree dirty: `{manifest.get('git_dirty', '')}`",
        f"- Processes: `{benchmark.get('processes', [])}`",
        f"- Per-process batch sizes: `{benchmark.get('batch_sizes', [])}`",
        f"- Warmups: `{benchmark.get('warmups', '')}`",
        f"- Measured repetitions: `{benchmark.get('repetitions', '')}`",
        f"- Forwarded arguments: `{manifest.get('forwarded_args', '')}`",
        "",
        "## Provenance",
        "",
        f"- Host: `{environment.get('hostname', '')}`",
        f"- CPU: `{environment.get('cpu_model', '')}`",
        f"- OpenMP threads per process: `{manifest.get('omp_num_threads', '')}`",
        f"- Implementation: `{build.get('implementation', 'fortran')}`",
        f"- Executable: `{build.get('resolved_executable') or build.get('requested_executable', '')}`",
        f"- Executable SHA-256: `{build.get('executable_sha256', '')}`",
        f"- Timer backend(s): `{', '.join(results.get('timer_backends', []))}`",
        f"- Time basis: `{', '.join(results.get('time_bases', []))}`",
        "",
        "## Contents",
        "",
        "- `profile.json`: machine-readable profile-wide metadata.",
        "- `batches.csv`: one row per process-count/per-process-batch-size workload.",
        "- `timings.csv`: raw observations for every process, timer, warmup, and measured repetition.",
        "- `logs/`: one representative CLUBB input, setup, log, and native timing file per workload.",
        "",
        "Warmup observations are retained in `timings.csv` with `phase=warmup`; analysis views ignore them by default.",
        "",
    ]
    _atomic_write(Path(profile_dir) / README_FILE, "\n".join(lines))


def create_profile_manifest(
    *,
    run_id: str,
    label: str,
    case_name: str,
    started_utc: str,
    git_revision: str,
    git_dirty: bool | None,
    forwarded_args: str,
    omp_num_threads: str,
    benchmark: dict[str, Any],
    environment: dict[str, Any],
    build: dict[str, Any],
) -> dict[str, Any]:
    return {
        "format": PROFILE_FORMAT,
        "format_version": PROFILE_FORMAT_VERSION,
        "run_id": run_id,
        "label": label,
        "case_name": case_name,
        "started_utc": started_utc,
        "completed_utc": "",
        "status": "running",
        "git_revision": git_revision,
        "git_dirty": git_dirty,
        "forwarded_args": forwarded_args,
        "omp_num_threads": omp_num_threads,
        "benchmark": benchmark,
        "environment": environment,
        "build": build,
        "results": {
            "timing_rows": 0,
            "batches_completed": 0,
            "failed_batches": 0,
            "observed_model_steps": [],
            "timer_backends": [],
            "time_bases": [],
        },
        "storage": "compact",
        "imported_from": "",
        "files": {
            "batches": BATCHES_FILE,
            "timings": TIMINGS_FILE,
            "logs": LOGS_DIRECTORY,
        },
    }


def update_profile_manifest(
    profile_dir: Path,
    manifest: dict[str, Any],
    *,
    status: str | None = None,
    timing_rows: int | None = None,
    batches_completed: int | None = None,
    failed_batches: int | None = None,
    completed: bool = False,
) -> dict[str, Any]:
    updated = dict(manifest)
    results = dict(updated.get("results") or {})
    for key, value in (
        ("timing_rows", timing_rows),
        ("batches_completed", batches_completed),
        ("failed_batches", failed_batches),
    ):
        if value is not None:
            results[key] = value
    updated["results"] = results
    if status is not None:
        updated["status"] = status
    if completed:
        updated["completed_utc"] = utc_now()
    write_profile_manifest(profile_dir, updated)
    write_profile_readme(profile_dir, updated)
    return updated


def _number_list(value: Any) -> str:
    if isinstance(value, (list, tuple)):
        return ",".join(str(item) for item in value)
    return str(value or "")


def _catalog_record(manifest: dict[str, Any], profile_dir: Path) -> dict[str, Any]:
    benchmark = dict(manifest.get("benchmark") or {})
    results = dict(manifest.get("results") or {})
    environment = dict(manifest.get("environment") or {})
    build = dict(manifest.get("build") or {})
    return {
        "format_version": manifest.get("format_version", PROFILE_FORMAT_VERSION),
        "run_id": manifest.get("run_id", profile_dir.name),
        "label": manifest.get("label", profile_dir.name),
        "case_name": manifest.get("case_name", ""),
        "started_utc": manifest.get("started_utc", ""),
        "completed_utc": manifest.get("completed_utc", ""),
        "status": manifest.get("status", "unknown"),
        "git_revision": manifest.get("git_revision", ""),
        "git_dirty": manifest.get("git_dirty", ""),
        "forwarded_args": manifest.get("forwarded_args", ""),
        "omp_num_threads": manifest.get("omp_num_threads", ""),
        "process_counts": _number_list(benchmark.get("processes")),
        "batch_sizes": _number_list(benchmark.get("batch_sizes")),
        # Dash retains its established generic column terminology internally.
        "columns_per_process": _number_list(benchmark.get("batch_sizes")),
        "vertical_levels": benchmark.get("vertical_levels", ""),
        "model_steps": _number_list(results.get("observed_model_steps")),
        "backends": _number_list(results.get("timer_backends")),
        "time_bases": _number_list(results.get("time_bases")),
        "hostname": environment.get("hostname", ""),
        "cpu_model": environment.get("cpu_model", ""),
        "implementation": build.get("implementation", "fortran"),
        "executable_sha256": build.get("executable_sha256", ""),
        "timing_rows": results.get("timing_rows", 0),
        "batches": results.get("batches_completed", 0),
        "storage": manifest.get("storage", "compact"),
        "bundle_path": str(profile_dir),
        "imported_from": manifest.get("imported_from", ""),
    }


def _profile_dirs(output_root: Path) -> list[Path]:
    root = Path(output_root)
    if read_profile_manifest(root) is not None:
        return [root]
    try:
        return [path for path in root.iterdir() if path.is_dir()]
    except OSError:
        return []


def discover_profiles(output_root: Path) -> list[dict[str, Any]]:
    """Discover compact profile directories without a mutable catalog file."""
    records = []
    for profile_dir in _profile_dirs(Path(output_root)):
        manifest = read_profile_manifest(profile_dir)
        if manifest is not None:
            records.append(_catalog_record(manifest, profile_dir))
    return sorted(
        records,
        key=lambda record: (str(record.get("started_utc") or ""), str(record.get("run_id") or "")),
        reverse=True,
    )


def refresh_catalog(output_root: Path) -> list[dict[str, Any]]:
    """Compatibility name: discovery is now direct and writes no catalog."""
    return discover_profiles(Path(output_root))


def _find_profile_dir(output_root: Path, run_id: str) -> Path:
    root = Path(output_root)
    root_manifest = read_profile_manifest(root)
    if root_manifest is not None and str(root_manifest.get("run_id") or root.name) == run_id:
        return root
    return profile_directory(root, run_id)


def _integer(value: Any, default: int = 0) -> int:
    try:
        return int(str(value))
    except (TypeError, ValueError):
        return default


def _floating(value: Any) -> float | None:
    try:
        return float(str(value))
    except (TypeError, ValueError):
        return None


def _derived_rows(profile_dir: Path) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    manifest = read_profile_manifest(profile_dir)
    if manifest is None:
        return [], []
    batches = {row.get("batch_id", ""): row for row in _read_csv(profile_dir / BATCHES_FILE)}
    raw = [row for row in _read_csv(profile_dir / TIMINGS_FILE) if row.get("phase") == "measured"]
    benchmark = dict(manifest.get("benchmark") or {})
    results_meta = dict(manifest.get("results") or {})
    backends = list(results_meta.get("timer_backends") or [])
    time_bases = list(results_meta.get("time_bases") or [])
    backend = str(backends[0] if backends else "")
    time_basis = str(time_bases[0] if time_bases else "")
    run_id = str(manifest.get("run_id") or profile_dir.name)
    common_profile = {
        "schema_version": str(manifest.get("format_version") or PROFILE_FORMAT_VERSION),
        "run_id": run_id,
        "started_utc": str(manifest.get("started_utc") or ""),
        "label": str(manifest.get("label") or run_id),
        "case_name": str(manifest.get("case_name") or ""),
        "git_revision": str(manifest.get("git_revision") or ""),
        "forwarded_args": str(manifest.get("forwarded_args") or ""),
        "omp_num_threads": str(manifest.get("omp_num_threads") or ""),
        "vertical_levels": benchmark.get("vertical_levels", ""),
    }

    by_execution: dict[tuple[str, int], list[dict[str, str]]] = {}
    for row in raw:
        key = (str(row.get("batch_id") or ""), _integer(row.get("repetition")))
        by_execution.setdefault(key, []).append(row)

    summaries: list[dict[str, Any]] = []
    processes: list[dict[str, Any]] = []
    for (batch_id, repetition), rows in sorted(by_execution.items()):
        batch = batches.get(batch_id, {})
        process_count = _integer(batch.get("process_count"))
        batch_size = _integer(batch.get("batch_size"))
        total_batch_size = _integer(batch.get("total_batch_size"), process_count * batch_size)
        group_rows = [row for row in rows if row.get("timer") == GROUP_WALL_TIMER]
        group_wall = _floating(group_rows[-1].get("inclusive_s")) if group_rows else None
        execution_status = str(group_rows[-1].get("status") or "") if group_rows else ""
        timer_rows = [row for row in rows if row.get("timer")]
        failure_messages = list(
            dict.fromkeys(str(row.get("message") or "") for row in rows if str(row.get("message") or ""))
        )

        common = {
            **common_profile,
            "process_count": process_count,
            "columns_per_process": batch_size,
            "batch_size": batch_size,
            "total_columns": total_batch_size,
            "total_batch_size": total_batch_size,
            "repetition": repetition,
        }
        successful_processes = {
            _integer(row.get("process"), -1)
            for row in rows
            if row.get("process") not in (None, "") and row.get("status") == "success"
        }
        for row in rows:
            if not row.get("timer") or row.get("timer") == GROUP_WALL_TIMER:
                continue
            process_index = _integer(row.get("process"))
            artifact = (
                str(batch.get("representative_phase") or "") == "measured"
                and _integer(batch.get("representative_repetition"), -1) == repetition
                and _integer(batch.get("representative_process"), -1) == process_index
            )
            case_name = str(manifest.get("case_name") or "")
            artifact_root = f"logs/{batch_id}"
            processes.append(
                {
                    **common,
                    "process_index": process_index,
                    "status": str(row.get("status") or ""),
                    "return_code": str(row.get("return_code") or ""),
                    "backend": backend,
                    "time_basis": time_basis,
                    "thread": _integer(row.get("thread")),
                    "timer_name": str(row.get("timer") or ""),
                    "calls": str(row.get("calls") or ""),
                    "inclusive_seconds": str(row.get("inclusive_s") or ""),
                    "exclusive_seconds": str(row.get("exclusive_s") or ""),
                    "timing_file": f"{artifact_root}/{case_name}.timing" if artifact else "",
                    "log_file": f"{artifact_root}/{case_name}.log" if artifact else "",
                    "message": str(row.get("message") or ""),
                }
            )

        grouped: dict[str, dict[int, list[dict[str, str]]]] = {}
        for row in timer_rows:
            grouped.setdefault(str(row["timer"]), {}).setdefault(_integer(row.get("process")), []).append(row)
        if not grouped:
            grouped[DEFAULT_ANALYSIS_TIMER] = {}
        for timer_name, by_process in sorted(
            grouped.items(), key=lambda item: (item[0] == GROUP_WALL_TIMER, item[0])
        ):
            critical = [
                max(process_timer_rows, key=lambda item: _floating(item.get("inclusive_s")) or 0.0)
                for process_timer_rows in by_process.values()
            ]
            seconds = [value for row in critical if (value := _floating(row.get("inclusive_s"))) is not None]
            calls = [str(row.get("calls") or "") for row in critical]
            maximum = max(seconds) if seconds else None
            status = execution_status or ("success" if len(successful_processes) == process_count else "failed")
            summaries.append(
                {
                    **common,
                    "status": status,
                    "successful_processes": len(successful_processes),
                    "processes_reported": len(seconds),
                    "backend": backend,
                    "time_basis": time_basis,
                    "timer_name": timer_name,
                    "calls": calls[0] if calls and len(set(calls)) == 1 else "",
                    "timer_min_seconds": min(seconds) if seconds else "",
                    "timer_mean_seconds": statistics.fmean(seconds) if seconds else "",
                    "timer_max_seconds": maximum if maximum is not None else "",
                    "timer_stdev_seconds": statistics.pstdev(seconds) if seconds else "",
                    "throughput_columns_per_second": (
                        total_batch_size / maximum if maximum is not None and maximum > 0 else ""
                    ),
                    "message": "; ".join(failure_messages),
                }
            )
    return summaries, processes


def read_profile_rows(output_root: Path, run_id: str, *, processes: bool = False) -> list[dict[str, Any]]:
    summaries, process_rows = _derived_rows(_find_profile_dir(Path(output_root), run_id))
    return process_rows if processes else summaries


def load_profiles(
    output_root: Path,
    run_ids: Sequence[str],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    catalog = {str(row.get("run_id")): row for row in discover_profiles(Path(output_root))}
    selected_ids = list(dict.fromkeys(str(value) for value in run_ids if str(value).strip()))
    base_labels = [str(catalog.get(run_id, {}).get("label") or run_id) for run_id in selected_ids]
    label_counts = {label: base_labels.count(label) for label in set(base_labels)}
    summaries: list[dict[str, Any]] = []
    processes: list[dict[str, Any]] = []
    for run_id in selected_ids:
        record = catalog.get(run_id, {})
        label = str(record.get("label") or run_id)
        if label_counts.get(label, 0) > 1:
            revision = str(record.get("git_revision") or "")[:8]
            started = str(record.get("started_utc") or "").replace("T", " ")[:16]
            discriminator = " · ".join(value for value in (revision, started) if value) or run_id[-8:]
            label = f"{label} · {discriminator}"
        summaries.extend({**row, "profile_id": run_id, "profile_label": label} for row in read_profile_rows(output_root, run_id))
        processes.extend(
            {**row, "profile_id": run_id, "profile_label": label}
            for row in read_profile_rows(output_root, run_id, processes=True)
        )
    return summaries, processes


def _safe_archive_relative(value: str) -> PurePosixPath:
    path = PurePosixPath(value)
    if path.is_absolute() or not path.parts or any(part in {"", ".", ".."} for part in path.parts):
        raise ProfileFormatError(f"unsafe archive path: {value}")
    return path


def export_profiles(output_root: Path, run_ids: Sequence[str]) -> bytes:
    selected = list(dict.fromkeys(str(value) for value in run_ids if str(value).strip()))
    if not selected:
        raise ValueError("select at least one profile to export")
    stream = io.BytesIO()
    export_records = []
    with zipfile.ZipFile(stream, "w", compression=zipfile.ZIP_DEFLATED, allowZip64=True) as archive:
        for run_id in selected:
            profile_dir = _find_profile_dir(Path(output_root), run_id)
            manifest = read_profile_manifest(profile_dir)
            if manifest is None:
                raise ProfileFormatError(f"profile not found: {run_id}")
            export_records.append(
                {
                    "run_id": run_id,
                    "label": manifest.get("label", run_id),
                    "case_name": manifest.get("case_name", ""),
                }
            )
            for path in sorted(profile_dir.rglob("*")):
                if path.is_file() and not path.name.startswith(".import_"):
                    archive.write(path, f"profiles/{run_id}/{path.relative_to(profile_dir).as_posix()}")
        archive.writestr(
            EXPORT_MANIFEST_FILE,
            json.dumps(
                {
                    "format": EXPORT_FORMAT,
                    "format_version": EXPORT_FORMAT_VERSION,
                    "created_utc": utc_now(),
                    "profiles": export_records,
                },
                indent=2,
                sort_keys=True,
            )
            + "\n",
        )
    return stream.getvalue()


def _zip_profile_prefixes(archive: zipfile.ZipFile) -> dict[str, str]:
    prefixes: dict[str, str] = {}
    for info in archive.infolist():
        stripped = info.filename.rstrip("/")
        path = _safe_archive_relative(stripped) if stripped else None
        if path is not None and len(path.parts) >= 3 and path.parts[0] == "profiles" and path.parts[-1] == MANIFEST_FILE:
            prefixes[path.parts[1]] = "/".join(path.parts[:2]) + "/"
    if not prefixes:
        raise ProfileFormatError("archive contains no CLUBB timing profiles")
    return prefixes


def _available_import_id(output_root: Path, requested: str) -> str:
    candidate = safe_profile_id(requested)
    if not profile_directory(output_root, candidate).exists():
        return candidate
    return f"{candidate}_imported_{uuid.uuid4().hex[:8]}"


def import_profiles(output_root: Path, file_name: str, data: bytes) -> list[str]:
    """Import one or more complete compact-profile ZIP bundles."""
    root = Path(output_root)
    root.mkdir(parents=True, exist_ok=True)
    if not zipfile.is_zipfile(io.BytesIO(data)):
        raise ProfileFormatError("import a complete CLUBB timing profile ZIP; timings.csv alone has no batch metadata")
    imported: list[str] = []
    with zipfile.ZipFile(io.BytesIO(data), "r") as archive:
        infos = archive.infolist()
        if len(infos) > 100_000 or sum(info.file_size for info in infos) > 8 * 1024**3:
            raise ProfileFormatError("profile archive is too large")
        prefixes = _zip_profile_prefixes(archive)
        for archived_id, prefix in prefixes.items():
            try:
                manifest = json.loads(archive.read(prefix + MANIFEST_FILE).decode("utf-8"))
            except (KeyError, UnicodeDecodeError, ValueError) as exc:
                raise ProfileFormatError(f"invalid manifest for profile {archived_id}") from exc
            if not isinstance(manifest, dict) or manifest.get("format") != PROFILE_FORMAT:
                raise ProfileFormatError(f"invalid CLUBB profile manifest for {archived_id}")
            old_run_id = str(manifest.get("run_id") or archived_id)
            new_run_id = _available_import_id(root, old_run_id)
            staging = Path(tempfile.mkdtemp(prefix=".import_", dir=root))
            try:
                for info in infos:
                    if info.is_dir() or not info.filename.startswith(prefix):
                        continue
                    relative = _safe_archive_relative(info.filename[len(prefix) :])
                    mode = (info.external_attr >> 16) & 0o170000
                    if mode == stat.S_IFLNK:
                        raise ProfileFormatError("profile archives may not contain symbolic links")
                    target = staging.joinpath(*relative.parts)
                    target.parent.mkdir(parents=True, exist_ok=True)
                    target.write_bytes(archive.read(info))
                manifest["run_id"] = new_run_id
                manifest["storage"] = "compact"
                manifest["imported_from"] = file_name
                write_profile_manifest(staging, manifest)
                destination = profile_directory(root, new_run_id)
                os.replace(staging, destination)
            except BaseException:
                shutil.rmtree(staging, ignore_errors=True)
                raise
            imported.append(new_run_id)
    return imported
