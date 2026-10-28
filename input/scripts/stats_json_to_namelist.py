#!/usr/bin/env python3
"""Resolve stats selections and serialize the final registry.

load_stats_categories() collects editable metadata subtrees. select_stats_entries()
accepts a complete preset, a JSON list of names, or categories/var:NAME joined
with comma or plus. An omitted selection uses the standard category. Custom
lists retain their explicit order and take metadata from the catalog; unions
use catalog order and emit each variable once. stats_to_namelist(entries) serializes
the final registry for create_case_namelist.py. stats_json_to_namelist() combines
selection and serialization in memory. Empty selections disable stats output
when the central builder assembles the runner input.
"""
from __future__ import annotations

import re

from pathlib import Path
from input.scripts.namelist import load_json


DEFAULT_STATS_JSON = Path(__file__).resolve().parents[1] / "stats.jsonc"
DEFAULT_STATS_SELECTION = "standard"
ENTRY_FIELDS = ("name", "grid", "units", "long_name")


def load_stats_categories(json_file: str | Path = DEFAULT_STATS_JSON) -> dict[str, list[dict[str, str]]]:
    """Read category paths and their descendants without fixed memberships/counts."""
    catalog = load_json(json_file)
    if not isinstance(catalog, dict) or "all" not in catalog or set(catalog) - {"all", "presets"}:
        raise ValueError("Stats JSON requires the 'all' category and optional 'presets' object")
    categories = {}

    def collect(node, path):
        if isinstance(node, dict):
            entries = []
            for name, child in node.items():
                entries.extend(collect(child, path + (name,)))
        elif isinstance(node, list):
            entries = node
            for entry in entries:
                if not isinstance(entry, dict) or any(
                    not isinstance(entry.get(key), str) or not entry[key]
                    or "\n" in entry[key] or "\r" in entry[key]
                    for key in ENTRY_FIELDS
                ):
                    raise ValueError(f"Invalid stat definition in {'/'.join(path)}")
        else:
            raise ValueError(f"Expected a category or stat list at {'/'.join(path)}")
        if path:
            categories["/".join(path)] = entries
        return entries

    collect({"all": catalog["all"]}, ())
    return categories


def select_stats_entries(selection=None, json_file=DEFAULT_STATS_JSON):
    """Select only requested groups/names, defaulting to standard; preserve metadata."""
    from input.scripts.namelist import require_non_null
    if selection is None:
        selection = DEFAULT_STATS_SELECTION
    if isinstance(selection, list) or Path(str(selection)).suffix.lower() in {".json", ".jsonc"}:
        names = selection if isinstance(selection, list) else load_json(Path(selection).expanduser())
        require_non_null(names, "stats.selection")
        if not isinstance(names, list) or any(not isinstance(name, str) for name in names):
            raise ValueError("Stats overrides must be a list of variable names")
        entries = {entry["name"]: entry for entry in load_stats_categories(json_file)["all"]}
        missing = set(names) - entries.keys()
        if missing:
            raise ValueError("Unknown stats variable(s): " + ", ".join(sorted(missing)))
        return [entries[name] for name in dict.fromkeys(names)]
    if Path(str(selection)).suffix.lower() == ".in":
        raise ValueError("Stats definitions must be JSON")
    tokens = [token.strip() for token in re.split(r"[,+]", str(selection))]
    if any(not token for token in tokens):
        raise ValueError("Stats selection must contain non-empty category names")
    if any(token.lower() == "none" for token in tokens):
        if len(tokens) != 1:
            raise ValueError("'none' cannot be combined with stats categories")
        return []

    categories = load_stats_categories(json_file)
    available_names = {entry["name"] for entry in categories["all"]}
    catalog = load_json(json_file)
    presets = catalog.get("presets", {})
    if not isinstance(presets, dict):
        raise ValueError("Stats presets must be an object")
    for name, members in presets.items():
        if not isinstance(members, list) or any(not isinstance(member, str) or member not in available_names for member in members):
            raise ValueError(f"Invalid stats preset {name!r}")
    if len(tokens) == 1 and tokens[0] in presets:
        entries = {entry["name"]: entry for entry in categories["all"]}
        return [entries[name] for name in dict.fromkeys(presets[tokens[0]])]
    selected_names = set()
    for token in tokens:
        if token in presets:
            selected_names.update(presets[token])
            continue
        if token.startswith("var:"):
            name = token[4:]
            if name not in available_names:
                raise ValueError(f"Unknown stats variable '{name}'")
            selected_names.add(name)
            continue
        if "/" in token:
            path = token if token.startswith("all/") else "all/" + token
            matches = [path] if path in categories else []
        else:
            matches = [path for path in categories if path.rsplit("/", 1)[-1] == token]
        if not matches:
            raise ValueError(
                f"Unknown stats category '{token}'. "
                "Use input/scripts/create_case_namelist.py -list stats to see categories and presets."
            )
        for path in matches:
            selected_names.update(entry["name"] for entry in categories[path])

    # Catalog order makes equivalent unions produce the same namelist.
    entries = {}
    for entry in categories["all"]:
        name = entry["name"]
        if name in entries and entry != entries[name]:
            raise ValueError(f"Conflicting stat definitions for '{name}'")
        entries[name] = entry
    return [entry for entry in entries.values() if entry["name"] in selected_names]


def stats_to_namelist(entries):
    """Serialize the final selected registry for a top-level model entry point."""
    if not entries:
        return ""
    lines = ["&clubb_stats_nl"]
    for index, entry in enumerate(entries, 1):
        value = " | ".join(entry[key] for key in ENTRY_FIELDS).replace('"', '""')
        lines.append(f'  entry({index}) = "{value}"')
    return "\n".join([*lines, "/"]) + "\n"


def stats_json_to_namelist(selection=None, json_file=DEFAULT_STATS_JSON):
    """Return the selected registry with metadata from the authoritative catalog."""
    return stats_to_namelist(select_stats_entries(selection, json_file))
