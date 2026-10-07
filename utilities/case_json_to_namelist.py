#!/usr/bin/env python3
"""Resolve CLUBB case definitions and emit their model namelists.

The catalog owns built-in cases; custom files store complete snapshots.
This module is shared by command-line runners, Dash and the local service.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import re
import tempfile

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CASES_JSON = REPO_ROOT / "input/case_setups/case_definitions.json"
CASE_NAME = re.compile(r"[A-Za-z][A-Za-z0-9_]*\Z")
MISSING = object()


def merge_case_settings(defaults, overrides):
    """Deep merge objects, replace arrays, and remove null-marked keys."""
    result = copy.deepcopy(defaults)
    for key, value in overrides.items():
        if value is None:
            result.pop(key, None)
        elif isinstance(value, dict):
            result[key] = merge_case_settings(result.get(key, {}), value)
        else:
            result[key] = copy.deepcopy(value)
    return result


def case_settings_patch(defaults, settings):
    """Express a resolved definition relative to the catalog defaults."""
    result = {}
    for key in dict.fromkeys([*settings, *defaults]):
        before, after = defaults.get(key, MISSING), settings.get(key, MISSING)
        if after is MISSING:
            result[key] = None
        elif before is MISSING:
            result[key] = copy.deepcopy(after)
        elif isinstance(before, dict) and isinstance(after, dict):
            patch = case_settings_patch(before, after)
            if patch:
                result[key] = patch
        elif before != after:
            result[key] = copy.deepcopy(after)
    return result


def load_case_catalog(json_file=None):
    path = Path(json_file or DEFAULT_CASES_JSON)
    catalog = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(catalog.get("defaults"), dict) or not isinstance(catalog.get("cases"), dict):
        raise ValueError("Case catalog requires defaults and cases objects")
    return catalog


def custom_case_path(name, json_file=None):
    if not CASE_NAME.fullmatch(str(name)):
        raise ValueError("Case names must start with a letter and use letters, numbers or underscores")
    root = Path(json_file or DEFAULT_CASES_JSON).parent / "custom"
    path = root / f"{name}.json"
    if not path.resolve().is_relative_to(root.resolve()):
        raise ValueError("Custom case file must stay inside the custom directory")
    return path


def available_case_names(json_file=None):
    catalog = load_case_catalog(json_file)
    names = set(catalog["cases"])
    root = Path(json_file or DEFAULT_CASES_JSON).parent / "custom"
    for path in root.glob("*.json"):
        if CASE_NAME.fullmatch(path.stem) and path.is_file() and path.resolve().is_relative_to(root.resolve()):
            names.add(path.stem)
    return sorted(names)


def case_definition_record(name, json_file=None):
    """Return JSON editor settings and their actual source/provenance."""
    path = custom_case_path(name, json_file)
    if path.is_file():
        saved = json.loads(path.read_text(encoding="utf-8"))
        if saved.get("name") != name or not isinstance(saved.get("namelists"), dict):
            raise ValueError(f"Invalid saved case {name!r}")
        settings = copy.deepcopy(saved["namelists"])
        kind = "custom"
    else:
        catalog = load_case_catalog(json_file)
        if name not in catalog["cases"]:
            raise ValueError(f"Unknown case {name!r}")
        path = Path(json_file or DEFAULT_CASES_JSON)
        saved = catalog.get("case_metadata", {}).get(name, {"copied_from": name, "benchmark_case": name})
        settings = merge_case_settings(catalog["defaults"], catalog["cases"][name])
        kind = "catalog"
    return {"name": name, "namelists": settings, "kind": kind, "source_file": str(path),
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "copied_from": saved.get("copied_from"), "benchmark_case": saved.get("benchmark_case")}


def resolved_case_settings(name, json_file=None):
    return case_definition_record(name, json_file)["namelists"]


def _literal(value):
    if isinstance(value, bool):
        return ".true." if value else ".false."
    if isinstance(value, str):
        return '"' + value.replace('"', '""') + '"'
    return "" if value is None else str(value)


def namelists_to_text(settings):
    """Write typed values, derived components and indexed arrays; no dependencies."""
    lines = []

    def fields(values, prefix=""):
        starts = values.get("_start_index", {})
        for key, value in values.items():
            if key == "_start_index":
                continue
            name = prefix + key
            if isinstance(value, dict):
                fields(value, name + "%")
            elif isinstance(value, list):
                array(name, value, starts.get(key, []))
            else:
                lines.append(f"  {name} = {_literal(value)}")

    def array(name, values, starts, indices=()):
        if any(isinstance(value, list) for value in values):
            axis = len(starts) - len(indices) - 1
            start = starts[axis] if starts and axis >= 0 else None
            if start is None:
                start = 1
            for offset, row in enumerate(values):
                array(name, row if isinstance(row, list) else [row], starts, (start + offset, *indices))
        else:
            lhs = name
            if indices or starts:
                first = starts[0] if starts else None
                span = ":" if first is None else f"{first}:{first + len(values) - 1}"
                lhs += "(" + ",".join([span, *map(str, indices)]) + ")"
            lines.append(f"  {lhs} = " + ", ".join("" if item is None else _literal(item) for item in values))

    if not isinstance(settings, dict):
        raise ValueError("Case namelists must be an object")
    for group, values in settings.items():
        if not isinstance(values, dict):
            raise ValueError("Namelist groups must be objects")
        lines.append("&" + group)
        fields(values)
        lines.append("/\n")
    return "\n".join(lines)


def _atomic_json(path, data, *, exclusive=False):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=".case-", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
            stream.write(_pretty_json(data) + "\n")
        if exclusive:
            os.link(temporary, path)
        else:
            os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def _pretty_json(value, level=0):
    if not isinstance(value, dict):
        return json.dumps(value, ensure_ascii=False, allow_nan=False)
    if not value:
        return "{}"
    rows = ["  " * (level + 1) + json.dumps(key) + ": " + _pretty_json(item, level + 1)
            for key, item in value.items()]
    return "{\n" + ",\n".join(rows) + "\n" + "  " * level + "}"


def save_case_definition(name, settings, *, copied_from=None, expected_sha256=None, json_file=None, source_record=None):
    """Create a local snapshot, or explicitly replace a confirmed existing case."""
    path = custom_case_path(name, json_file)
    catalog = load_case_catalog(json_file)
    source = source_record if source_record is not None else (case_definition_record(copied_from, json_file) if copied_from else {})
    benchmark_case = source.get("benchmark_case")
    if settings.get("model_setting", {}).get("runtype") != source.get("namelists", {}).get("model_setting", {}).get("runtype"):
        benchmark_case = None
    exists = name in available_case_names(json_file)
    if exists:
        current = case_definition_record(name, json_file)
        if expected_sha256 is None:
            raise FileExistsError(f"Case {name!r} already exists; confirm replacement")
        if current["sha256"] != expected_sha256:
            raise ValueError("Case changed since confirmation; review the replacement again")
        if current["kind"] == "catalog":
            catalog["cases"][name] = case_settings_patch(catalog["defaults"], settings)
            catalog.setdefault("case_metadata", {})[name] = {"copied_from": copied_from,
                                                            "benchmark_case": benchmark_case}
            _atomic_json(Path(json_file or DEFAULT_CASES_JSON), catalog)
            return case_definition_record(name, json_file)
    settings = copy.deepcopy(settings)
    if not exists and settings.get("stats_setting", {}).get("fname_prefix", MISSING) == source.get("namelists", {}).get("stats_setting", {}).get("fname_prefix"):
        settings["stats_setting"]["fname_prefix"] = name
    _atomic_json(path, {"name": name, "copied_from": copied_from,
                        "benchmark_case": benchmark_case, "namelists": settings}, exclusive=not exists)
    return case_definition_record(name, json_file)


def main():
    parser = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    parser.add_argument("case_name", nargs="?")
    parser.add_argument("-json_file", type=Path, default=DEFAULT_CASES_JSON)
    parser.add_argument("-output_file", type=Path)
    parser.add_argument("-list", action="store_true")
    args = parser.parse_args()
    try:
        if args.list:
            print("\n".join(available_case_names(args.json_file)))
        else:
            if not args.case_name:
                parser.error("Specify a case name or -list")
            text = namelists_to_text(resolved_case_settings(args.case_name, args.json_file))
            if args.output_file:
                args.output_file.parent.mkdir(parents=True, exist_ok=True)
                args.output_file.write_text(text, encoding="utf-8")
            else:
                print(text, end="")
    except (ValueError, OSError) as exc:
        parser.exit(2, f"{parser.prog}: error: {exc}\n")


if __name__ == "__main__":
    main()
