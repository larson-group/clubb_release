#!/usr/bin/env python3
"""Select stats.json categories and emit a CLUBB stats registry namelist.

Examples:
    python utilities/stats_json_to_namelist.py core,radiation
    python utilities/stats_json_to_namelist.py standard/radiation -output_file /tmp/radiation.in
    python utilities/stats_json_to_namelist.py -list

Short names select every matching category; slash-separated paths select one
subtree. 'var:NAME' selects an individual stat. Commas and '+' form unions,
with each stat emitted once in catalog order. 'none' emits no registry;
create_case_namelist also disables stats output.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

DEFAULT_STATS_JSON = Path(__file__).resolve().parents[1] / "input/stats/stats.json"
ENTRY_FIELDS = ("name", "grid", "units", "long_name")


def load_stats_categories(json_file: str | Path = DEFAULT_STATS_JSON) -> dict[str, list[dict[str, str]]]:
    """Read category paths and their descendants without fixed memberships/counts."""
    with open(json_file, encoding="utf-8") as source:
        catalog = json.load(source)
    if not isinstance(catalog, dict) or set(catalog) != {"all"}:
        raise ValueError("Stats JSON must contain only the 'all' root category")
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

    collect(catalog, ())
    return categories


def select_stats_entries(selection: str, json_file: str | Path = DEFAULT_STATS_JSON) -> list[dict[str, str]]:
    """Return selected stats once, in catalog order, with their original metadata."""
    tokens = [token.strip() for token in re.split(r"[,+]", selection)]
    if any(not token for token in tokens):
        raise ValueError("Stats selection must contain non-empty category names")
    if any(token.lower() == "none" for token in tokens):
        if len(tokens) != 1:
            raise ValueError("'none' cannot be combined with stats categories")
        return []

    categories = load_stats_categories(json_file)
    available_names = {entry["name"] for entry in categories["all"]}
    selected_names = set()
    for token in tokens:
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
                "Use utilities/stats_json_to_namelist.py -list to see categories."
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


def stats_json_to_namelist(selection: str, json_file: str | Path = DEFAULT_STATS_JSON) -> str:
    """Return the selected registry, retaining templates and original metadata."""
    selected_entries = select_stats_entries(selection, json_file)
    if selection.strip().lower() == "none":
        return ""
    lines = ["&clubb_stats_nl"]
    for index, entry in enumerate(selected_entries, 1):
        value = " | ".join(entry[key] for key in ENTRY_FIELDS).replace('"', '""')
        lines.append(f'  entry({index}) = "{value}"')
    lines.append("/")
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__, allow_abbrev=False,
                                    formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("selection", nargs="?", default="standard",
                        help="Category names, paths or var:NAME joined with ',' or '+', all, or none.")
    parser.add_argument("-json_file", type=Path, default=DEFAULT_STATS_JSON, help="Stats catalog JSON file.")
    parser.add_argument("-output_file", type=Path, help="Write a namelist file instead of stdout.")
    parser.add_argument("-list", action="store_true", help="List category paths and stat counts.")
    args = parser.parse_args()
    try:
        if args.list:
            categories = load_stats_categories(args.json_file)
            for path, entries in categories.items():
                print(f"{path}: {len({entry['name'] for entry in entries})}")
        else:
            namelist = stats_json_to_namelist(args.selection, args.json_file)
            if args.output_file:
                args.output_file.parent.mkdir(parents=True, exist_ok=True)
                args.output_file.write_text(namelist, encoding="utf-8")
            else:
                sys.stdout.write(namelist)
    except (ValueError, OSError) as exc:
        parser.exit(2, f"{parser.prog}: error: {exc}\n")


if __name__ == "__main__":
    main()
