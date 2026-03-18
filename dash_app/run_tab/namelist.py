"""Namelist parsing and override-file helpers for the run tab."""

import os
import re
import tempfile

BOOL_TRUE = {".true.", "true", "t"}
BOOL_FALSE = {".false.", "false", "f"}


def parse_line(line):
    """Parse one namelist assignment line into name/value metadata."""
    if "=" not in line:
        return None
    if line.lstrip().startswith(("!", "&")):
        return None
    if line.strip() == "/":
        return None
    comment = ""
    idx = line.find("!")
    if idx != -1:
        comment = line[idx:].rstrip("\n")
        body = line[:idx]
    else:
        body = line
    match = re.match(r"^(\s*)([A-Za-z_]\w*)\s*=\s*(.*)$", body)
    if not match:
        return None
    indent, name, rest = match.groups()
    rest = rest.rstrip()
    comma = rest.endswith(",")
    if comma:
        rest = rest[:-1].rstrip()
    return {
        "indent": indent,
        "name": name,
        "value": rest.strip(),
        "comment": comment,
        "comma": comma,
    }


def read_namelist_entries(path):
    """Read parsed namelist entries from a file, preserving line index metadata."""
    entries = []
    if not os.path.exists(path):
        return entries
    with open(path, "r", encoding="utf-8") as handle:
        for idx, line in enumerate(handle.readlines()):
            parsed = parse_line(line)
            if not parsed:
                continue
            parsed["line_index"] = idx
            entries.append(parsed)
    return entries


def is_bool_value(value):
    """Return whether a namelist value is a recognized Fortran boolean literal."""
    return value.strip().lower() in BOOL_TRUE | BOOL_FALSE


def is_true(value):
    """Return whether a namelist value is a recognized true literal."""
    return value.strip().lower() in BOOL_TRUE


def normalize_numeric_display(value):
    """Normalize numeric strings for stable UI and dirty-state comparisons."""
    text = str(value).strip()
    if not text:
        return text
    lower = text.lower()
    if lower in BOOL_TRUE or lower in BOOL_FALSE:
        return text
    if "e" in lower:
        return text
    if "." not in text:
        return text
    trimmed = text.rstrip("0")
    if trimmed.endswith("."):
        trimmed += "0"
    return trimmed


def apply_updates_to_lines(lines, updates):
    """Apply in-memory namelist updates to a list of source lines."""
    remaining = dict(updates or {})
    for idx, line in enumerate(lines):
        parsed = parse_line(line)
        if not parsed:
            continue
        name = parsed["name"]
        if name not in remaining:
            continue
        new_value = remaining.pop(name)
        comma = "," if parsed["comma"] else ""
        comment = parsed["comment"]
        sep = " " if comment else ""
        updated = f"{parsed['indent']}{name} = {new_value}{comma}{sep}{comment}".rstrip()
        lines[idx] = updated + "\n"
        if not remaining:
            break
    for name, new_value in remaining.items():
        lines.append(f"{name} = {new_value}\n")
    return lines


def update_namelist_value(path, name, new_value):
    """Rewrite one value in place in a namelist file."""
    if not os.path.exists(path):
        return
    with open(path, "r", encoding="utf-8") as handle:
        lines = handle.readlines()
    lines = apply_updates_to_lines(lines, {name: new_value})
    with open(path, "w", encoding="utf-8") as handle:
        handle.writelines(lines)


def write_temp_namelist(base_path, updates, prefix):
    """Write a temporary namelist override file and return its path."""
    if not updates or not os.path.exists(base_path):
        return None
    with open(base_path, "r", encoding="utf-8") as handle:
        lines = handle.readlines()
    lines = apply_updates_to_lines(lines, updates)
    temp_file = tempfile.NamedTemporaryFile(delete=False, prefix=prefix, suffix=".in", dir="/tmp", mode="w", encoding="utf-8")
    temp_file.writelines(lines)
    temp_file.close()
    return temp_file.name


def cleanup_temp_files(paths):
    """Delete temporary override files, ignoring missing-file and permission issues."""
    if not paths:
        return
    for path in paths:
        if not path:
            continue
        try:
            if os.path.exists(path):
                os.remove(path)
        except Exception:
            # Temporary-file cleanup should not make the dashboard brittle.
            pass


def build_override_updates(flag_values, param_values, defaults_data, flag_names_data, param_meta):
    """Build per-file namelist overrides from current UI settings versus defaults."""
    updates = {"flags": {}, "tunable": {}, "silhs": {}}
    if not defaults_data:
        return updates
    for name, values in zip(flag_names_data or [], flag_values or []):
        current = bool(values)
        default = bool(defaults_data["flags"].get(name))
        if current != default:
            updates["flags"][name] = ".true." if current else ".false."
    for meta, value in zip(param_meta or [], param_values or []):
        file_key = meta.get("file")
        name = meta.get("name")
        if not file_key or not name:
            continue
        default_value = defaults_data["params"][file_key].get(name)
        current_value = normalize_numeric_display(value)
        if current_value != normalize_numeric_display(default_value):
            updates[file_key][name] = current_value
    return updates
