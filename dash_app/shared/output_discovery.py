"""Read-only discovery of CLUBB statistics outputs.

The Plot tab selects parent output directories while file-oriented diagnostics
select individual NetCDF files.  Keep their common filesystem traversal here
without coupling either view to Dash callbacks, brokers, or job records.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Iterable


STATS_SUFFIX = "_stats.nc"


def has_netcdf_signature(path: str | Path) -> bool:
    """Return whether *path* has a classic NetCDF or NetCDF-4 signature."""
    try:
        with Path(path).open("rb") as handle:
            header = handle.read(8)
    except OSError:
        return False
    return header.startswith(b"CDF") or header == b"\x89HDF\r\n\x1a\n"


def _roots(roots: Iterable[str | Path]) -> list[Path]:
    unique: list[Path] = []
    seen: set[Path] = set()
    for raw_root in roots:
        root = Path(raw_root).expanduser().resolve()
        if root.is_dir() and root not in seen:
            unique.append(root)
            seen.add(root)
    return unique


def _label_records(records: list[dict]) -> None:
    """Assign stable display labels, widening only ambiguous root-relative labels."""
    counts: dict[str, int] = {}
    for record in records:
        label = str(record["label"])
        counts[label] = counts.get(label, 0) + 1
    for record in records:
        if counts[str(record["label"])] > 1:
            record["label"] = str(record["path"])


def _file_fingerprint(path: Path) -> dict[str, int] | None:
    """Return JSON-serializable metadata for a readable stats file."""
    try:
        stat = path.stat()
    except OSError:
        return None
    return {
        "mtime_ns": int(stat.st_mtime_ns),
        "size_bytes": int(stat.st_size),
    }


def _discover(
    roots: Iterable[str | Path],
    *,
    recursive: bool,
    max_depth: int,
    exclude_dir_names: Iterable[str],
    max_directories: int | None,
    max_scanned_directories: int | None,
) -> tuple[list[dict], list[dict]]:
    excluded = {str(name) for name in exclude_dir_names}
    directories: list[dict] = []
    files: list[dict] = []
    scanned = 0
    for root in _roots(roots):
        for directory_name, child_names, file_names in os.walk(root, topdown=True):
            scanned += 1
            if max_scanned_directories is not None and scanned > max_scanned_directories:
                child_names[:] = []
                break
            directory = Path(directory_name)
            relative = directory.relative_to(root)
            depth = len(relative.parts)
            if relative.parts and relative.parts[0] in excluded:
                child_names[:] = []
                continue
            child_names[:] = [
                child
                for child in child_names
                if not child.startswith(".")
                and child not in excluded
                and recursive
                and depth < max_depth
            ]
            stats_files = []
            for name in sorted(file_names):
                path = directory / name
                if not name.endswith(STATS_SUFFIX) or not has_netcdf_signature(path):
                    continue
                fingerprint = _file_fingerprint(path)
                if fingerprint is not None:
                    stats_files.append((path, fingerprint))
            if not stats_files:
                continue
            relative_text = relative.as_posix() if relative.parts else ""
            base_label = root.name if not relative_text else f"{root.name}/{relative_text}"
            cases = {
                path.name[: -len(STATS_SUFFIX)]: str(path)
                for path, _fingerprint in stats_files
            }
            case_fingerprints = {
                path.name[: -len(STATS_SUFFIX)]: fingerprint
                for path, fingerprint in stats_files
            }
            latest_stats_modified_ns = max(
                fingerprint["mtime_ns"] for _path, fingerprint in stats_files
            )
            directory_record = {
                "path": str(directory),
                "key": str(directory),
                "root": str(root),
                "relative_path": relative_text,
                "label": base_label,
                "case_names": sorted(cases),
                "cases": cases,
                "case_count": len(cases),
                "case_fingerprints": case_fingerprints,
                "latest_stats_modified_ns": latest_stats_modified_ns,
                "modified": latest_stats_modified_ns / 1_000_000_000.0,
                "available": True,
            }
            directories.append(directory_record)
            for path, fingerprint in stats_files:
                case = path.name[: -len(STATS_SUFFIX)]
                files.append(
                    {
                        "path": str(path),
                        "key": str(path),
                        "directory": str(directory),
                        "root": str(root),
                        "relative_path": f"{relative_text + '/' if relative_text else ''}{path.name}",
                        "label": f"{base_label}/{path.name}",
                        "case": case,
                        "fingerprint": fingerprint,
                        "modified": fingerprint["mtime_ns"] / 1_000_000_000.0,
                    }
                )
    _label_records(directories)
    _label_records(files)
    directories.sort(key=lambda item: (-float(item["modified"]), str(item["label"])))
    for index, record in enumerate(directories):
        record["recency_index"] = index
        record["is_newest"] = index == 0
    if max_directories is not None:
        directories = directories[:max_directories]
    files.sort(key=lambda item: (str(item["label"]), str(item["path"])))
    return directories, files


def discover_stats_directories(
    roots: Iterable[str | Path],
    *,
    recursive: bool = True,
    max_depth: int = 5,
    exclude_dir_names: Iterable[str] = ("agent_artifacts",),
    max_directories: int | None = None,
    max_scanned_directories: int | None = None,
) -> list[dict]:
    """Return output directories containing directly readable stats files."""
    directories, _files = _discover(
        roots,
        recursive=recursive,
        max_depth=max_depth,
        exclude_dir_names=exclude_dir_names,
        max_directories=max_directories,
        max_scanned_directories=max_scanned_directories,
    )
    return directories


def discover_stats_files(
    roots: Iterable[str | Path],
    *,
    recursive: bool = True,
    max_depth: int = 5,
    exclude_dir_names: Iterable[str] = ("agent_artifacts",),
    max_scanned_directories: int | None = None,
) -> list[dict]:
    """Return individual readable stats files with globally unique keys/labels."""
    _directories, files = _discover(
        roots,
        recursive=recursive,
        max_depth=max_depth,
        exclude_dir_names=exclude_dir_names,
        max_directories=None,
        max_scanned_directories=max_scanned_directories,
    )
    return files
