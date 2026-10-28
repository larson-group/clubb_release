"""Read and merge JSON settings, then serialize finalized namelist groups.

Usage:
    load_json() reads JSON/JSONC values with duplicate-key and non-finite-constant
    checks. updated_json_text() edits existing values while retaining comments.
    require_non_null() rejects null input settings. merge_settings() applies
    sparse object patches; settings_patch()
    expresses resolved objects relative to defaults. atomic_text() replaces a
    generated text/JSON file without exposing a partial write.

    updated_settings() applies typed UI edits. parse_value() converts direct
    CLI assignments into JSON values. setting_entries() supplies scalar UI
    display values from resolved JSON.

    namelists_to_text() emits finalized groups, derived components, logicals,
    quoted strings and arrays with optional _start_index metadata. Call it
    through create_case_namelist.py at the model entry-point boundary. It
    preserves omitted keys and enforces no physics policy.
"""
from __future__ import annotations

import copy
import json
import os
from pathlib import Path
import re
import tempfile

def _literal(value):
    if isinstance(value, bool):
        return ".true." if value else ".false."
    if isinstance(value, str):
        return '"' + value.replace('"', '""') + '"'
    return str(value)


def require_non_null(value, path="settings"):
    """Reject null settings at any depth, reporting the offending field or index."""
    if value is None:
        raise ValueError(f"Null is not supported at {path}; use an explicit value")
    if isinstance(value, dict):
        for key, item in value.items():
            require_non_null(item, f"{path}.{key}")
    elif isinstance(value, list):
        for index, item in enumerate(value):
            require_non_null(item, f"{path}[{index}]")


def namelists_to_text(settings):
    """Write typed values, derived components and indexed arrays; no dependencies."""
    require_non_null(settings)
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
            start = starts[axis] if starts and axis >= 0 else 1
            for offset, row in enumerate(values):
                array(name, row if isinstance(row, list) else [row], starts, (start + offset, *indices))
        else:
            lhs = name
            if indices or starts:
                span = f"{starts[0]}:{starts[0] + len(values) - 1}" if starts else ":"
                lhs += "(" + ",".join([span, *map(str, indices)]) + ")"
            lines.append(f"  {lhs} = " + ", ".join(_literal(item) for item in values))

    if not isinstance(settings, dict):
        raise ValueError("Case namelists must be an object")
    for group, values in settings.items():
        if not isinstance(values, dict):
            raise ValueError("Namelist groups must be objects")
        lines.append("&" + group)
        fields(values)
        lines.append("/\n")
    return "\n".join(lines)



def _unique_object(pairs):
    """Reject duplicate JSON keys instead of silently losing a definition."""
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Duplicate JSON key: {key}")
        result[key] = value
    return result


def load_json(path):
    """Read JSON or JSONC; reject duplicate keys and non-finite constants."""
    def invalid_constant(value):
        raise ValueError(f"Non-finite JSON value: {value}")
    path = Path(path)
    if path.suffix.lower() not in {".json", ".jsonc"}:
        raise ValueError(f"Input definitions must be JSON or JSONC: {path}")
    value = json.loads(_uncomment_json(path.read_text(encoding="utf-8")), object_pairs_hook=_unique_object,
                       parse_constant=invalid_constant)
    return value


_JSON_PART = re.compile(r'"(?:[^"\\]|\\[\s\S])*"|//[^\r\n]*|/\*[\s\S]*?(?:\*/|\Z)')


def _uncomment_json(text):
    """Replace comments with whitespace, retaining string contents and error locations."""
    def replace(match):
        token = match.group()
        if token.startswith('"'):
            return token
        if token.startswith("/*") and not token.endswith("*/"):
            raise ValueError("Unterminated JSONC block comment")
        return re.sub(r"[^\r\n]", " ", token)
    return _JSON_PART.sub(replace, text)


def updated_json_text(path, settings):
    """Render JSON edits without rewriting unchanged values or removing comments.

    Existing objects are edited by member and equal-length arrays by element.
    Comments inside a removed/replaced container remain at that location as
    standalone notes. New files use ordinary JSON. Trailing commas are rejected.
    The caller owns atomic replacement and any concurrency checks.
    """
    path = Path(path)
    if not path.exists():
        return json.dumps(settings, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    original = path.read_text(encoding="utf-8")
    old = load_json(path)
    clean = _uncomment_json(original)
    tokens = list(re.finditer(r'"(?:[^"\\]|\\.)*"|[^\s{}\[\],:]+|[{}\[\],:]', clean))
    cursor = 0
    edits = []

    def node():
        """Index one validated JSON value and its children by source offsets."""
        nonlocal cursor
        first = tokens[cursor]
        cursor += 1
        children = []
        if first.group() in {"{", "["}:
            closing = "}" if first.group() == "{" else "]"
            while tokens[cursor].group() != closing:
                key = tokens[cursor] if first.group() == "{" else None
                if key:
                    cursor += 2
                child = node()
                comma = tokens[cursor] if tokens[cursor].group() == "," else None
                if comma:
                    cursor += 1
                children.append((key, child, comma))
            end = tokens[cursor].end()
            cursor += 1
        else:
            end = first.end()
        return first.start(), end, children

    def notes(start, end):
        """Keep comment text when the value surrounding it is removed."""
        return "\n".join(match.group() for match in _JSON_PART.finditer(original[start:end])
                         if not match.group().startswith('"'))

    def render(value, indent):
        """Indent newly inserted JSON relative to its containing member."""
        return json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False).replace("\n", "\n" + indent)

    def update(span, before, after):
        """Collect replacements for changed values and object membership."""
        start, end, children = span
        if type(before) is type(after) and before == after and not isinstance(before, (dict, list)):
            return
        indent = re.match(r" *", original[original.rfind("\n", 0, start) + 1:start]).group()
        if isinstance(before, dict) and isinstance(after, dict):
            retained = [json.loads(key.group()) for key, _, _ in children if json.loads(key.group()) in after]
            additions = {key: value for key, value in after.items() if key not in before}
            for key, child, comma in children:
                name = json.loads(key.group())
                if name in after:
                    update(child, before[name], after[name])
                else:
                    comments = notes(key.start(), child[1])
                    edits.append((key.start(), child[1], comments + "\n" if comments else ""))
                keep_comma = name in retained and (name != retained[-1] or bool(additions))
                if comma and not keep_comma:
                    edits.append((comma.start(), comma.end(), ""))
                elif keep_comma and not comma:
                    edits.append((child[1], child[1], ","))
            if additions:
                child_indent = indent + "  "
                rows = [child_indent + json.dumps(key) + ": " + render(value, child_indent)
                        for key, value in additions.items()]
                edits.append((end - 1, end - 1, "\n" + ",\n".join(rows) + "\n" + indent))
        elif isinstance(before, list) and isinstance(after, list) and len(before) == len(after):
            for (_, child, _), old_value, value in zip(children, before, after):
                update(child, old_value, value)
        else:
            comments = notes(start, end)
            edits.append((start, end, (comments + "\n" + indent if comments else "") + render(after, indent)))

    update(node(), old, settings)
    insertions = {}
    replacements = []
    for start, end, replacement in edits:
        if start == end:
            insertions[start] = insertions.get(start, "") + replacement
        else:
            replacements.append((start, end, replacement))
    edits = replacements + [(offset, offset, text) for offset, text in insertions.items()]
    for start, end, replacement in sorted(edits, reverse=True):
        original = original[:start] + replacement + original[end:]
    return original


def merge_settings(defaults, patch):
    """Copy defaults and apply sparse object overrides; lists replace whole lists."""
    require_non_null(defaults, "defaults")
    require_non_null(patch, "overrides")

    def merge(base, updates):
        """Merge validated objects without changing either input."""
        if not isinstance(base, dict) or not isinstance(updates, dict):
            raise ValueError("Config defaults and patches must be objects")
        result = copy.deepcopy(base)
        for key, value in updates.items():
            result[key] = merge(result.get(key, {}), value) if isinstance(value, dict) else copy.deepcopy(value)
        return result

    return merge(defaults, patch)


def atomic_text(path, text):
    """Replace a generated text file atomically; preserve already identical files."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.is_file() and path.read_text(encoding="utf-8") == text:
        return
    descriptor, temporary = tempfile.mkstemp(dir=path.parent, prefix=".json-input-")
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
            stream.write(text)
        os.replace(temporary, path)
    finally:
        Path(temporary).unlink(missing_ok=True)


def updated_settings(settings, updates, categories=("model", "silhs")):
    """Copy a snapshot and apply named scalar edits using its existing JSON types."""
    result = copy.deepcopy(settings)
    for name, raw in (updates or {}).items():
        matches = [(category, key) for category in categories for key in result[category]
                   if key.casefold() == name.casefold()]
        if len(matches) != 1:
            raise ValueError(f"Config setting not present or ambiguous: {name}")
        category, key = matches[0]
        old = result[category][key]
        text = str(raw).strip()
        if isinstance(old, bool):
            if text.lower() not in {".true.", ".false.", "true", "false", "t", "f"}:
                raise ValueError(f"{name} requires a logical value")
            value = text.lower() in {".true.", "true", "t"}
        else:
            value = float(text.replace("D", "E").replace("d", "e"))
            if isinstance(old, int):
                if not value.is_integer():
                    raise ValueError(f"{name} requires an integer")
                value = int(value)
        result[category][key] = value
    return result


def settings_patch(defaults, settings):
    """Return sparse overrides, requiring explicit values for inherited settings."""
    require_non_null(defaults, "defaults")
    require_non_null(settings)

    def difference(base, values, path):
        """Keep changed values and reject omissions that would require deletion."""
        missing = base.keys() - values.keys()
        if missing:
            name = next(key for key in base if key in missing)
            raise ValueError(f"Cannot omit inherited setting {path}.{name}; use an explicit value")
        result = {}
        for key, after in values.items():
            if key not in base:
                result[key] = copy.deepcopy(after)
            elif isinstance(base[key], dict) and isinstance(after, dict):
                patch = difference(base[key], after, f"{path}.{key}")
                if patch:
                    result[key] = patch
            elif base[key] != after:
                result[key] = copy.deepcopy(after)
        return result

    return difference(defaults, settings, "settings")


def parse_value(raw):
    """Convert a CLI scalar or comma-separated list to typed JSON values."""
    text = str(raw).strip().rstrip(",").strip()
    quote, pieces, start = None, [], 0
    for index, char in enumerate(text):
        if char in "\"'":
            if quote == char:
                quote = None
            elif quote is None:
                quote = char
        elif char == "," and quote is None:
            pieces.append(text[start:index])
            start = index + 1
    if quote:
        raise ValueError(f"Unterminated quoted override: {raw}")
    if pieces:
        return [parse_value(piece) for piece in [*pieces, text[start:]]]
    if len(text) >= 2 and text[0] == text[-1] and text[0] in "\"'":
        return text[1:-1].replace(text[0] * 2, text[0])
    if text.lower() == "null":
        raise ValueError("Null is not supported in overrides; use an explicit value")
    if text.lower() in {"true", ".true.", "t"}:
        return True
    if text.lower() in {"false", ".false.", "f"}:
        return False
    try:
        return int(text)
    except ValueError:
        try:
            value = float(text.replace("D", "E").replace("d", "e"))
            if not __import__("math").isfinite(value):
                raise ValueError("Overrides require finite numbers")
            return value
        except ValueError:
            if text.lower() in {"nan", "inf", "infinity", "-inf", "+inf"}:
                raise ValueError("Overrides require finite numbers")
            return text


def setting_entries(settings):
    """Return scalar JSON settings as editable UI strings, preserving their order."""
    return [{"name": name, "value": json.dumps(value)}
            for name, value in settings.items()]
