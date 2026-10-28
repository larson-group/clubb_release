"""Category selections and portable stats lists shared by Dash and MCP."""

import os
import json
import re
from pathlib import Path

from input.scripts.stats_json_to_namelist import DEFAULT_STATS_JSON, load_stats_categories, select_stats_entries
from input.scripts.namelist import load_json

STATS_DIR = DEFAULT_STATS_JSON.parent


def list_stats_files():
    """Return named JSON presets and saved JSON snapshots relative to input."""
    root = Path(STATS_DIR).resolve()
    saved = sorted(path.relative_to(root).as_posix() for path in [*root.glob("*.stats.json"), *root.glob("*.stats.jsonc")]
                   if path.is_file() and path.resolve().is_relative_to(root))
    presets = load_json(DEFAULT_STATS_JSON).get("presets", {})
    return [*presets, *saved]


def load_stats_choices():
    """Offer quick category selections and reusable JSON presets/snapshots."""
    return ["all", "standard", "core", *list_stats_files()]


def validate_stats_selection(value):
    """Validate categories, presets or a contained JSON snapshot."""
    value = str(value).strip()
    if value in list_stats_files():
        select_stats_entries(str(Path(STATS_DIR) / value) if value.endswith((".json", ".jsonc")) else value)
        return value
    if value.endswith((".in", ".json", ".jsonc")):
        raise ValueError("Choose a saved JSON stats list inside input/")
    select_stats_entries(value)
    return value


def stats_cli_argument(value, *, absolute=False):
    """Keep named selections intact and resolve saved JSON snapshots for the CLI."""
    value = validate_stats_selection(value)
    if value.endswith((".json", ".jsonc")):
        root = STATS_DIR if absolute else "input"
        return os.path.join(root, value)
    return value


def stats_input_path(value):
    """Identify the authoritative JSON definition for run provenance."""
    if str(value).lower() == "none":
        return None
    return Path(STATS_DIR) / value if value.endswith((".json", ".jsonc")) else DEFAULT_STATS_JSON


def stats_selection_names(value):
    """Return variable names from a category, preset or saved JSON snapshot."""
    if value in list_stats_files() and value.endswith((".json", ".jsonc")):
        names = load_json(Path(STATS_DIR) / value)
        if not isinstance(names, list) or any(not isinstance(name, str) for name in names):
            raise ValueError("Stats overrides must be a list of variable names")
        return set(names)
    source = value
    return {entry["name"] for entry in select_stats_entries(source)}


def stats_selection_expression(names, categories=None):
    """Express exact variables with complete categories and individual leftovers."""
    categories = load_stats_categories() if categories is None else categories
    selected = set(names)
    available = {entry["name"] for entry in categories["all"]}
    if selected - available:
        raise ValueError("Selected variables are missing from the stats catalog.")

    def collect(path):
        members = {entry["name"] for entry in categories[path]}
        if not members & selected:
            return []
        if members <= selected:
            return [path.removeprefix("all/")]
        children = [child for child in categories if "/" in child and child.rsplit("/", 1)[0] == path]
        if children:
            return [token for child in children for token in collect(child)]
        return [f"var:{entry['name']}" for entry in categories[path] if entry["name"] in selected]

    return ",".join(dict.fromkeys(collect("all"))) or "none"


def stats_list_options():
    """Label stored custom lists separately from the old registries."""
    return [
        {"label": f"{'Saved' if name.endswith(('.stats.json', '.stats.jsonc')) else 'Preset'} · {Path(name).stem}", "value": name}
        for name in list_stats_files()
    ]


def save_stats_config(name, selection):
    """Save a new portable JSON snapshot; never overwrite an existing list."""
    name = str(name or "").strip()
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", name):
        raise ValueError("Use a name with letters, numbers, underscores or hyphens.")
    content = json.dumps([entry["name"] for entry in select_stats_entries(selection)], indent=2, allow_nan=False) + "\n"
    if not select_stats_entries(selection):
        raise ValueError("Choose at least one variable before saving.")
    path = Path(STATS_DIR) / f"{name}.stats.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        with path.open("x", encoding="utf-8") as output:
            output.write(content)
    except FileExistsError:
        raise ValueError(f"'{name}' already exists. Choose a different name.") from None
    return path.relative_to(STATS_DIR).as_posix()
