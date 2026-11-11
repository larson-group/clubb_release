"""Read Native procedure declarations for the C++ namelist generators.

Declarations are authoritative. Unsupported argument types or missing intents
raise errors so generated bindings cannot silently omit an argument.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

from generate_derived_types import ROOT, split_top


def statements(text: str) -> list[str]:
    """Join free-form continuations while preserving quoted exclamation marks."""
    result: list[str] = []
    pending = ""
    for raw in text.splitlines():
        if raw.lstrip().startswith("#"):
            continue
        quote = ""
        code = ""
        for char in raw:
            if char in "'\"":
                if quote == char:
                    quote = ""
                elif not quote:
                    quote = char
            if char == "!" and not quote:
                break
            code += char
        code = code.strip()
        if not code:
            continue
        continued = code.endswith("&")
        pending += " " + code.strip("& ")
        if not continued:
            result.append(pending.strip())
            pending = ""
    if pending:
        raise ValueError("unfinished Fortran continuation")
    return result


@dataclass
class Argument:
    name: str
    kind: str
    dims: tuple[str, ...]
    intent: str
    optional: bool
    allocatable: bool
    character_length: str = ""


@dataclass
class Procedure:
    name: str
    module: str
    source: str
    arguments: list[Argument]
    result: Argument | None = None


def declaration(statement: str) -> dict[str, Argument]:
    """Read typed entities from a modern declaration with a double colon."""
    if "::" not in statement:
        return {}
    attributes, variables = statement.split("::", 1)
    attrs = split_top(attributes)
    spec = re.sub(r"\s+", "", attrs[0]).lower()
    length = ""
    if spec.startswith("type("):
        kind = spec[5:-1]
    elif spec.startswith("integer"):
        kind = "integer"
    elif spec.startswith("logical"):
        kind = "logical"
    elif spec.startswith("real"):
        kind = (
            "double"
            if any(k in spec for k in ("time_precision", "stat_rknd", "dp", "c_double"))
            else "real"
        )
    elif spec.startswith("character"):
        kind = "character"
        match = re.search(r"len=([^,)]+)", spec)
        length = match.group(1) if match else "1"
    else:
        return {}
    dimension = next((a for a in attrs if re.match(r"dimension\s*\(", a, re.I)), "")
    dims = (
        tuple(split_top(dimension[dimension.index("(") + 1 : -1])) if dimension else ()
    )
    intent_match = re.search(r"intent\s*\(\s*(inout|out|in)\s*\)", attributes, re.I)
    intent = intent_match.group(1).lower() if intent_match else ""
    result = {}
    for variable in split_top(variables):
        match = re.match(r"(\w+)\s*(\(.*?\))?\s*(?:=|$)", variable)
        if match is None:
            raise ValueError(f"unsupported variable declaration: {statement}")
        name = match.group(1)
        var_dims = tuple(split_top(match.group(2)[1:-1])) if match.group(2) else dims
        result[name.lower()] = Argument(
            name,
            kind,
            var_dims,
            intent,
            bool(re.search(r"\boptional\b", attributes, re.I)),
            bool(re.search(r"\ballocatable\b", attributes, re.I)),
            length,
        )
    return result


def read_procedure(path: Path, name: str) -> Procedure:
    """Return one routine's argument inventory in its native signature order."""
    import subprocess

    source = subprocess.run(
        [
            "cpp",
            "-P",
            "-traditional-cpp",
            "-DCLUBB",
            "-DCLUBB_REAL_TYPE=8",
            "-DCOAMPS_MICRO",
            "-DUWM_MISC",
            "-Dnooverlap",
            "-Dradoffline",
            str(path),
        ],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    module_match = re.search(r"(?im)^\s*module\s+(?!procedure\b)(\w+)", source)
    if module_match is None:
        raise ValueError(f"missing module in {path}")
    code = statements(source)
    for start, line in enumerate(code):
        # Fortran permits a bare zero-argument subroutine header. Normalize only
        # its declaration form, so END SUBROUTINE cannot become a false header.
        if re.fullmatch(
            rf"(?:(?:pure|impure|elemental|recursive|non_recursive|module)\s+)*"
            rf"subroutine\s+{re.escape(name)}",
            line,
            re.I,
        ):
            line += "()"
        match = re.search(
            rf"\b(subroutine|function)\s+{name}\s*\((.*?)\)(?:\s*result\s*\(\s*(\w+)\s*\))?",
            line,
            re.I,
        )
        if match:
            break
    else:
        raise ValueError(f"missing procedure {name} in {path}")
    names = [arg for arg in split_top(match.group(2)) if arg]
    result_name = (
        (match.group(3) or name).lower()
        if match.group(1).lower() == "function"
        else None
    )
    decls = {}
    for line in code[start + 1 :]:
        if re.match(rf"end\s+(subroutine|function)\s+{name}\b", line, re.I):
            break
        decls.update(declaration(line))
    # These legacy declarations omit INTENT. Their native bodies only read
    # the listed inputs and overwrite every listed output before reading it.
    # All other missing contracts still fail generation.
    legacy_intents = {
        "lin_interp_between_grids": {
            "interpolate_altitudes": "in",
            "current_altitudes": "in",
            "current_values": "in",
        },
        "calc_sfc_varnce": {"dt": "in", "up2_sfc_coef": "in", "a_const": "in"},
        "calc_wp4_pdf": {"wp4": "out"},
        "calc_wpxpyp_pdf": {"wpxpyp": "out"},
    }
    for key, intent in legacy_intents.get(name.lower(), {}).items():
        if not decls[key].intent:
            decls[key].intent = intent
    missing = [
        arg
        for arg in names
        if arg.lower() not in decls or not decls[arg.lower()].intent
    ]
    if missing:
        raise ValueError(f"missing argument declaration/intent in {name}: {missing}")
    result = None
    if result_name:
        result = decls.get(result_name)
        if result is None:
            header_type = code[start].split("function", 1)[0].strip()
            parsed = declaration(f"{header_type} :: {result_name}")
            result = parsed.get(result_name)
        if result is None:
            raise ValueError(f"missing result type in {name}")
        result.intent = "out"
    return Procedure(
        name,
        module_match.group(1),
        str(path.relative_to(ROOT)),
        [decls[arg.lower()] for arg in names],
        result,
    )


@lru_cache(maxsize=1)
def source_files() -> tuple[tuple[Path, str], ...]:
    """Read the native source inventory once during a generator invocation."""
    return tuple((path, path.read_text()) for path in (ROOT / "src").rglob("*.F90"))


@lru_cache(maxsize=1)
def procedure_sources() -> dict[str, tuple[Path, ...]]:
    """Index definition names once so generating many bridges remains cheap."""
    pattern = re.compile(
        r"(?im)^[ \t]*(?:(?:pure|impure|elemental|recursive|module|logical|integer|real|complex|character|double[ \t]+precision)"
        r"(?:[ \t]*\([^\n)]*\))?[ \t]+)*(?:subroutine|function)[ \t]+(\w+)"
    )
    result: dict[str, set[Path]] = {}
    for path, text in source_files():
        for match in pattern.finditer(text):
            result.setdefault(match.group(1).lower(), set()).add(path)
    return {name: tuple(paths) for name, paths in result.items()}


def find_procedure(name: str) -> Procedure:
    """Locate a native routine by its definition rather than its public alias."""
    matches = procedure_sources().get(name.lower(), ())
    if len(matches) != 1:
        raise ValueError(f"expected one definition of {name}, got {matches}")
    return read_procedure(matches[0], name)
