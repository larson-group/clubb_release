#!/usr/bin/env python3
"""Resolve, edit and serialize case JSON definitions.

case_definition_record() and resolved_case_settings() resolve a preset or an
arbitrary JSON override file against catalog defaults. available_case_names()
lists presets and locally saved *.case.json or *.case.jsonc files. load_case_catalog() reads
the source. merge_case_settings() and case_settings_patch() form/apply sparse
overrides. save_case_definition() writes copyable overrides with provenance;
existing entries require a matching source checksum. custom_case_path() confines
UI-created names to input/. case_to_namelist() serializes final case groups for
create_case_namelist.py; it does not create files.
"""
from __future__ import annotations

import copy
import hashlib
import os
from pathlib import Path
import re
import tempfile


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CASES_JSON = REPO_ROOT / "input/case_definitions.jsonc"
CASE_NAME = re.compile(r"[A-Za-z][A-Za-z0-9_]*\Z")
MISSING = object()


from input.scripts.namelist import (
    load_json, merge_settings as merge_case_settings, settings_patch as case_settings_patch,
    updated_json_text,
)


def load_case_catalog(json_file=None):
    """Read the case catalog through the shared non-null input validation."""
    from input.scripts.create_case_namelist import load_catalog
    catalog = load_catalog("cases", json_file or DEFAULT_CASES_JSON)
    if not isinstance(catalog["defaults"], dict):
        raise ValueError("Case catalog requires defaults and presets objects")
    return catalog


def custom_case_path(name, json_file=None):
    if not CASE_NAME.fullmatch(str(name)):
        raise ValueError("Case names must start with a letter and use letters, numbers or underscores")
    root = Path(json_file or DEFAULT_CASES_JSON).parent
    path = root / f"{name}.case.json"
    if not path.resolve().is_relative_to(root.resolve()):
        raise ValueError("Custom case file must stay inside the custom directory")
    return path


def available_case_names(json_file=None):
    catalog = load_case_catalog(json_file)
    names = set(catalog["presets"])
    root = Path(json_file or DEFAULT_CASES_JSON).parent
    for path in [*root.glob("*.case.json"), *root.glob("*.case.jsonc")]:
        name = path.stem.removesuffix(".case")
        if CASE_NAME.fullmatch(name) and path.is_file() and path.resolve().is_relative_to(root.resolve()):
            names.add(name)
    return sorted(names)


def case_definition_record(name, json_file=None, *, config=None):
    """Resolve a preset or override file and retain its source and benchmark identity."""
    from input.scripts.create_case_namelist import resolve_input_settings
    catalog = load_case_catalog(json_file)
    value = os.fspath(name)
    candidate = Path(value).expanduser()
    explicit_file = candidate.suffix.lower() in {".json", ".jsonc"} or candidate.is_file()
    path = candidate if explicit_file else custom_case_path(name, json_file)
    if not explicit_file:
        matches = [item for item in (path, path.with_suffix(".jsonc")) if item.is_file()]
        if len(matches) > 1:
            raise ValueError(f"Multiple custom definitions for case {name!r}")
        if matches:
            path = matches[0]
    if path.is_file():
        saved = load_json(path)
        if not isinstance(saved, dict):
            raise ValueError("A case override must be an object")
        metadata = saved.get("_metadata", {})
        selection, kind = str(path), "custom"
    else:
        if explicit_file:
            raise FileNotFoundError(path)
        metadata = catalog.get("case_metadata", {}).get(name, {"copied_from": name, "benchmark_case": name})
        path = Path(json_file or DEFAULT_CASES_JSON)
        selection, kind = name, "catalog"
    settings = resolve_input_settings("cases", selection, config=config, json_file=json_file or DEFAULT_CASES_JSON)
    if explicit_file and "benchmark_case" not in metadata:
        runtype = settings.get("model_setting", {}).get("runtype")
        metadata = {**metadata, "benchmark_case": runtype if runtype in catalog["presets"] else None}
    return {"name": candidate.stem.removesuffix(".case") if explicit_file else name,
            "namelists": settings, "kind": kind, "source_file": str(path),
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "copied_from": metadata.get("copied_from"), "benchmark_case": metadata.get("benchmark_case")}


def resolved_case_settings(name, json_file=None):
    return case_definition_record(name, json_file)["namelists"]


from input.scripts.namelist import namelists_to_text


def _atomic_json(path, data, *, exclusive=False):
    """Write a case edit atomically, preserving existing catalog comments."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=".case-", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
            stream.write(updated_json_text(path, data))
        if exclusive:
            os.link(temporary, path)
        else:
            os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


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
            catalog["presets"][name] = case_settings_patch(catalog["defaults"], settings)
            catalog.setdefault("case_metadata", {})[name] = {"copied_from": copied_from,
                                                            "benchmark_case": benchmark_case}
            _atomic_json(Path(json_file or DEFAULT_CASES_JSON), catalog)
            return case_definition_record(name, json_file)
        path = Path(current["source_file"])
    settings = copy.deepcopy(settings)
    if not exists and settings.get("stats_setting", {}).get("fname_prefix", MISSING) == source.get("namelists", {}).get("stats_setting", {}).get("fname_prefix"):
        settings["stats_setting"]["fname_prefix"] = name
    _atomic_json(path, {"_metadata": {"copied_from": copied_from, "benchmark_case": benchmark_case},
                        **case_settings_patch(catalog["defaults"], settings)}, exclusive=not exists)
    return case_definition_record(name, json_file)



def case_to_namelist(settings):
    """Serialize finalized case/loss groups for the central entry-point builder."""
    return namelists_to_text(settings)
