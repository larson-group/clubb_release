"""Resolve the shared case, parameter and flag catalogs without writing files.

Usage:
    load_catalog() reads defaults and presets. available_configs() lists their
    combined names, including the built-in default. normalize_input_selection()
    distinguishes preset names from paths; input_patch() reads either form.
    resolve_input_settings() merges defaults, a matching shared config and an
    explicit section selection, in that order.

Understanding:
    Custom JSON/JSONC files have the same sparse contents as preset entries.
    Objects merge recursively, arrays replace arrays, and null settings fail.
    catalog_index() retains leaf ancestry and editable locations;
    catalog_presets() provides inherited patches. resolve_selection() expands
    lists/groups before execution. split_values() and positive_int() parse shared
    selection/worker arguments. Case names are global; flag members are local to
    their owning group. Every group name ends in _group.
    Stats selection is independent and belongs to stats_json_to_namelist.py.
"""
from __future__ import annotations

import copy
import os
from pathlib import Path

from input.scripts.namelist import load_json, merge_settings, require_non_null

INPUT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = "default"
CATALOG_FILES = {
    "cases": "case_definitions.jsonc",
    "parameters": "clubb_params.jsonc",
    "flags": "clubb_config_flags.jsonc",
}


def load_catalog(kind, json_file=None):
    """Read one authoritative defaults/presets catalog without creating files."""
    if kind not in CATALOG_FILES:
        raise ValueError(f"Unknown input catalog {kind!r}")
    catalog = load_json(json_file or INPUT_ROOT / CATALOG_FILES[kind])
    if not isinstance(catalog, dict) or "defaults" not in catalog or not isinstance(catalog.get("presets"), dict):
        raise ValueError(f"{kind} catalog requires defaults and presets")
    require_non_null(catalog["defaults"], f"{kind}.defaults")
    catalog_index(kind, catalog)
    return catalog


def available_configs(root=None):
    """Return the built-in default choice and names defined by any catalog."""
    root = Path(root or INPUT_ROOT)
    return sorted({DEFAULT_CONFIG} | {name for kind, filename in CATALOG_FILES.items()
                   for name in catalog_index(kind, load_catalog(kind, root / filename))[0]})


def normalize_input_selection(selection):
    """Keep preset names intact and resolve JSON paths against the caller's directory."""
    value = os.fspath(selection)
    path = Path(value).expanduser()
    if path.suffix.lower() in {".json", ".jsonc"} or path.is_file() or value.startswith(("/", "~", ".")):
        return str(path.resolve())
    return value


def input_patch(kind, selection, *, json_file=None):
    """Read an explicit named preset or JSON override with identical contents."""
    value = normalize_input_selection(selection)
    path = Path(value)
    if path.is_absolute():
        return load_json(path)
    catalog = load_catalog(kind, json_file)
    presets, _, _ = catalog_index(kind, catalog)
    if value not in presets:
        raise ValueError(f"Unknown {kind} preset {value!r}. Available: {', '.join(presets)}")
    return presets[value]["settings"]


def resolve_input_settings(kind, selection=None, *, config=None, json_file=None):
    """Apply defaults, a shared config and explicit overrides; default skips presets."""
    catalog = load_catalog(kind, json_file)
    if config == DEFAULT_CONFIG:
        config = None
    if config and json_file is None and config not in available_configs():
        raise ValueError(f"Unknown config {config!r}")
    result = copy.deepcopy(catalog["defaults"])
    patches = []
    presets, _, _ = catalog_index(kind, catalog)
    if config in presets:
        patches.append(presets[config]["settings"])
    if selection is not None:
        patches.append(input_patch(kind, selection, json_file=json_file))
    for patch in patches:
        patch = {key: value for key, value in patch.items() if key != "_metadata"} if isinstance(patch, dict) else patch
        require_non_null(patch, f"{kind}.overrides")
        result = merge_settings(result, patch)
    return result


def catalog_index(kind, catalog):
    """Index leaf settings, group membership and editable locations in one tree.

    Case names are globally unique. Flag members are qualified by their owning
    group, allowing each flag collection to retain its own flag1, flag2, etc.
    Each record retains its inherited base and source container for case edits.
    """
    presets, groups, locations = {}, {}, {}

    def visit(node, inherited, owner=None):
        """Collect leaves and descendants while checking names and group structure."""
        if not isinstance(node, dict) or not isinstance(node.get("presets", {}), dict):
            raise ValueError("Catalog presets must be objects")
        members = []
        for name, patch in node.get("presets", {}).items():
            if not isinstance(name, str) or not name or "/" in name or not isinstance(patch, dict):
                raise ValueError(f"Invalid {kind} preset: {name!r}")
            key = f"{owner}/{name}" if kind == "flags" and owner else name
            if key in presets or key in groups:
                raise ValueError(f"Duplicate {kind} name: {key}")
            values = {k: v for k, v in patch.items() if k != "_metadata"}
            require_non_null(values, f"{kind}.{key}")
            presets[key] = {"settings": merge_settings(inherited, values),
                            "inherited": inherited, "name": name, "group": owner}
            locations[key] = node["presets"]
            members.append(key)
        children = node.get("groups", {})
        if not isinstance(children, dict):
            raise ValueError("Catalog groups must be objects")
        for name, group in children.items():
            if not isinstance(name, str) or not name.endswith("_group") or "/" in name:
                raise ValueError("Group names must end in _group and contain no slash")
            if name in groups or name in presets:
                raise ValueError(f"Duplicate {kind} name: {name}")
            if not isinstance(group, dict) or set(group) - {"settings", "presets", "groups", "_metadata"}:
                raise ValueError(f"Invalid group definition: {name}")
            groups[name] = []
            base = merge_settings(inherited, group.get("settings", {}))
            descendants = visit(group, base, name)
            if not descendants:
                raise ValueError(f"Empty group: {name}")
            groups[name] = descendants
            members.extend(descendants)
        return members

    visit(catalog, {})
    lists = catalog.get("lists", {})
    if not isinstance(lists, dict):
        raise ValueError("Catalog lists must be objects")
    for name, members in lists.items():
        if name in presets or name in groups or not isinstance(members, list) or not members:
            raise ValueError(f"Invalid or ambiguous list: {name}")
        if any(not isinstance(member, str) or not member for member in members):
            raise ValueError(f"List {name} requires nonempty selector names")
    return presets, groups, locations


def catalog_presets(kind, catalog):
    """Return all leaf patches with ancestor settings applied, without defaults."""
    return {name: record["settings"] for name, record in catalog_index(kind, catalog)[0].items()}


def split_values(values):
    """Split CLI selections at spaces or commas while preserving their order."""
    if isinstance(values, str):
        values = [values]
    return [part.strip() for value in values or [] for part in value.split(",") if part.strip()]


def resolve_selection(kind, selection, *, json_file=None):
    """Expand named leaves, groups and run lists; deduplicate in first-use order.

    Lists may reference other lists or groups. Cycles and unknown references
    fail before execution. Case JSON files and locally saved cases are leaves.
    """
    catalog = load_catalog(kind, json_file)
    presets, groups, _ = catalog_index(kind, catalog)
    lists = catalog.get("lists", {})
    result = []

    def expand(name, ancestors=()):
        """Resolve one selector and detect list cycles with a useful path."""
        if name in ancestors:
            raise ValueError("Selection cycle: " + " -> ".join((*ancestors, name)))
        if name in lists:
            for member in lists[name]:
                expand(member, (*ancestors, name))
        elif name in groups:
            result.extend(groups[name])
        elif name in presets:
            result.append(name)
        elif kind == "cases":
            from input.scripts.case_json_to_namelist import case_definition_record
            case_definition_record(name, json_file)
            result.append(name)
        else:
            raise ValueError(f"Unknown {kind} selection: {name}")

    for name in split_values(selection):
        expand(name)
    if not result:
        raise ValueError(f"Empty {kind} selection")
    return list(dict.fromkeys(result))


def positive_int(value):
    """Parse a positive worker count for the model and comparison CLIs."""
    import argparse
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be >= 1")
    return parsed


def catalog_selectors(kind, catalog):
    """List selectable leaf names, inherited groups and ordered run lists."""
    presets, groups, _ = catalog_index(kind, catalog)
    return [*presets, *groups, *catalog.get("lists", {})]
