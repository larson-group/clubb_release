"""JSON setting display and override-file helpers for the run tab."""

import math
import os
import re
import tempfile

BOOL_TRUE = {".true.", "true", "t"}
BOOL_FALSE = {".false.", "false", "f"}
INT_LITERAL_RE = re.compile(r"^[+-]?\d+$")


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


def format_update_value(existing_value, new_value):
    """Format a new value using the existing namelist literal as a type hint."""
    text = str(new_value).strip()
    if INT_LITERAL_RE.fullmatch(str(existing_value).strip()):
        try:
            numeric_value = float(text.replace("D", "E").replace("d", "e"))
        except (TypeError, ValueError):
            return text
        if math.isfinite(numeric_value) and numeric_value.is_integer():
            return str(int(numeric_value))
    return text


def write_temp_settings(base_settings, updates, prefix):
    """Write edited parameter/flag JSON to a temporary snapshot and return its path."""
    if not updates:
        return None
    import json
    from input.scripts.namelist import updated_settings

    categories = ("silhs",) if "silhs" in prefix else ("model",) if "params" in prefix else ("model", "silhs")
    settings = updated_settings(base_settings, updates, categories)
    with tempfile.NamedTemporaryFile(delete=False, prefix=prefix, suffix=".json", dir="/tmp", mode="w", encoding="utf-8") as stream:
        json.dump(settings, stream, indent=2, allow_nan=False)
        stream.write("\n")
        return stream.name


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
