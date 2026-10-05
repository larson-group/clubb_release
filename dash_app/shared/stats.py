"""Category selections and portable stats lists shared by Dash and MCP."""

import os
import re
from pathlib import Path

from utilities.stats_json_to_namelist import DEFAULT_STATS_JSON, load_stats_categories, select_stats_entries, stats_json_to_namelist

STATS_DIR = DEFAULT_STATS_JSON.parent


def list_stats_files():
    """Return legacy and saved namelists relative to input/stats."""
    root = Path(STATS_DIR).resolve()
    return sorted(
        path.relative_to(root).as_posix()
        for path in root.rglob("*.in")
        if path.is_file() and path.resolve().is_relative_to(root)
    )


def load_stats_choices():
    """Offer the three quick selections and reusable custom/legacy lists."""
    return ["all", "standard", "core", *list_stats_files()]


def validate_stats_selection(value):
    """Validate categories or a namelist contained in the local stats directory."""
    value = str(value).strip()
    if value.lower() == "none":
        return "none"
    if value in list_stats_files():
        return value
    select_stats_entries(value)
    return value


def stats_cli_argument(value, *, absolute=False):
    """Keep category expressions intact; root local list names for the SCM CLI."""
    value = validate_stats_selection(value)
    if value in list_stats_files():
        root = STATS_DIR if absolute else os.path.join("input", "stats")
        return os.path.join(root, value)
    return value


def stats_input_path(value):
    """Identify the actual registry/catalog for run provenance."""
    if str(value).lower() == "none":
        return None
    return Path(STATS_DIR) / value if value in list_stats_files() else DEFAULT_STATS_JSON


def stats_selection_names(value):
    """Read names for a category selection or a saved/legacy registry preview."""
    if value in list_stats_files():
        text = (Path(STATS_DIR) / value).read_text(encoding="utf-8")
        return set(re.findall(r'entry\s*\(\s*\d+\s*\)\s*=\s*[\"\']\s*([^|\"\']+?)\s*\|', text, re.I))
    return {entry["name"] for entry in select_stats_entries(value)}


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
        {"label": f"{'Saved' if name.startswith('custom/') else 'Legacy'} · {Path(name).stem}", "value": name}
        for name in list_stats_files()
    ]


def save_stats_config(name, selection):
    """Save a new portable namelist; never overwrite an existing list."""
    name = str(name or "").strip()
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", name):
        raise ValueError("Use a name with letters, numbers, underscores or hyphens.")
    content = stats_json_to_namelist(selection)
    if not content:
        raise ValueError("Choose at least one variable before saving.")
    path = Path(STATS_DIR) / "custom" / f"{name}.in"
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        with path.open("x", encoding="utf-8") as output:
            output.write(content)
    except FileExistsError:
        raise ValueError(f"'{name}' already exists. Choose a different name.") from None
    return path.relative_to(STATS_DIR).as_posix()
