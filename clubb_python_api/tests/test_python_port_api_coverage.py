"""Checks that src/CLUBB_core cross-module public routine usage has Python API coverage."""

import ast
import re
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
PY_API_FILE = REPO_ROOT / "clubb_python_api" / "clubb_python" / "clubb_api.py"
CLUBB_CORE_DIR = REPO_ROOT / "src" / "CLUBB_core"

def _collect_python_api_defs() -> set[str]:
    """Collect exported function names from clubb_python/clubb_api.py."""
    tree = ast.parse(PY_API_FILE.read_text())
    defs: set[str] = set()
    for node in tree.body:
        if isinstance(node, ast.FunctionDef):
            defs.add(node.name.lower())
        elif isinstance(node, ast.ImportFrom):
            for alias in node.names:
                if alias.name == "*":
                    continue
                defs.add((alias.asname or alias.name).lower())
    return defs


def _strip_fortran_comments(line: str) -> str:
    """Remove trailing Fortran comments from a line."""
    return line.split("!")[0].rstrip()


def _strip_gfdl_blocks(lines: list[str]) -> list[str]:
    """Drop GFDL-only preprocessor blocks; the Python driver does not support them."""
    out: list[str] = []
    depth = 0
    for line in lines:
        stripped = line.strip()
        if re.match(r"#\s*ifdef\s+GFDL\b", stripped, re.IGNORECASE):
            depth += 1
            continue
        if depth and re.match(r"#\s*if", stripped, re.IGNORECASE):
            depth += 1
            continue
        if depth and re.match(r"#\s*endif\b", stripped, re.IGNORECASE):
            depth -= 1
            continue
        if depth:
            continue
        out.append(line)
    return out


def _collect_use_only_imports(fortran_file: Path) -> dict[str, set[str]]:
    """Collect `use module, only: ...` imports from a Fortran file."""
    imports: dict[str, set[str]] = {}
    lines = _strip_gfdl_blocks(fortran_file.read_text().splitlines())

    statements: list[str] = []
    buf = ""
    in_use_stmt = False
    for raw in lines:
        line = _strip_fortran_comments(raw)
        if not line:
            continue
        if in_use_stmt:
            buf += " " + line.strip()
            if not line.strip().endswith("&"):
                statements.append(buf)
                buf = ""
                in_use_stmt = False
        elif re.match(r"\s*use\s+", line, re.IGNORECASE):
            buf = line.strip()
            if line.strip().endswith("&"):
                in_use_stmt = True
            else:
                statements.append(buf)
                buf = ""

    for stmt in statements:
        match = re.match(r"use\s+([A-Za-z0-9_]+)\s*,\s*only\s*:\s*(.*)", stmt, re.IGNORECASE)
        if not match:
            continue
        module_name = match.group(1).lower()
        names = []
        for raw_name in match.group(2).replace("&", " ").split(","):
            name = raw_name.strip()
            if not name:
                continue
            if "=>" in name:
                _, name = [part.strip() for part in name.split("=>", 1)]
            names.append(name)
        imports.setdefault(module_name, set()).update(names)

    return imports


def _collect_clubb_core_module_names() -> set[str]:
    """Collect module names defined within src/CLUBB_core."""
    module_names: set[str] = set()
    for path in sorted(CLUBB_CORE_DIR.glob("*.F90")):
        text = "\n".join(_strip_gfdl_blocks(path.read_text().splitlines()))
        for match in re.finditer(r"^\s*module\s+([A-Za-z0-9_]+)\b", text, re.IGNORECASE | re.MULTILINE):
            name = match.group(1).lower()
            if name == "procedure":
                continue
            module_names.add(name)
    return module_names


def _collect_fortran_cross_module_routine_usage(fortran_file: Path, ignored_modules: set[str]) -> list[tuple[str, str]]:
    """Return (module, routine) pairs used as cross-module routine calls/functions."""
    imports = _collect_use_only_imports(fortran_file)
    lines = _strip_gfdl_blocks(fortran_file.read_text().splitlines())
    body = "\n".join(_strip_fortran_comments(line) for line in lines)

    used: list[tuple[str, str]] = []
    seen: set[tuple[str, str]] = set()
    for module_name, names in sorted(imports.items()):
        if module_name in ignored_modules:
            continue
        for routine in sorted(names):
            if re.search(rf"\bcall\s+{re.escape(routine)}\b", body, re.IGNORECASE) or re.search(
                rf"\b{re.escape(routine)}\s*\(",
                body,
                re.IGNORECASE,
            ):
                pair = (module_name, routine)
                if pair not in seen:
                    used.append(pair)
                    seen.add(pair)
    return used


def _collect_all_clubb_core_cross_module_api_gaps() -> set[tuple[str, str, str, tuple[str, ...]]]:
    """Audit all src/CLUBB_core files for cross-module CLUBB_core routine usage without Python API exports."""
    api_defs = _collect_python_api_defs()
    clubb_core_modules = _collect_clubb_core_module_names()

    gaps: set[tuple[str, str, str, tuple[str, ...]]] = set()
    for fortran_file in sorted(CLUBB_CORE_DIR.glob("*.F90")):
        for module_name, routine in _collect_fortran_cross_module_routine_usage(
            fortran_file,
            ignored_modules={"advance_sclrm_nd_module"},
        ):
            if module_name not in clubb_core_modules:
                continue
            if routine.lower() == "clubb_at_least_debug_level_api":
                continue

            expected_names = _expected_api_names(routine)
            if any(name in api_defs for name in expected_names):
                continue

            gaps.add((fortran_file.name, module_name, routine, expected_names))

    return gaps


def _format_gap_report(gaps: list[tuple[str, str, str, tuple[str, ...]]]) -> str:
    """Render a readable failure/report message for missing Python API coverage."""
    lines = [
        "Missing Python API coverage for cross-module CLUBB_core public routine usage.",
        f"Gap count: {len(gaps)}",
        "",
    ]
    for fortran_file, module_name, routine, expected_names in gaps:
        expected = ", ".join(expected_names)
        lines.append(f"{fortran_file}: {module_name}::{routine} -> expected api export(s): {expected}")
    return "\n".join(lines)


def _expected_api_names(fortran_routine: str) -> tuple[str, ...]:
    """Map a Fortran public routine name to the expected Python API export(s)."""
    lower_name = fortran_routine.lower()
    if lower_name.endswith("_api"):
        return (lower_name[:-4],)
    return (lower_name,)


def test_full_clubb_core_cross_module_api_coverage():
    """Report any cross-module CLUBB_core public routine usage that lacks a Python API export."""
    gaps = sorted(_collect_all_clubb_core_cross_module_api_gaps())
    assert not gaps, _format_gap_report(gaps)


if __name__ == "__main__":
    gaps = sorted(_collect_all_clubb_core_cross_module_api_gaps())
    if gaps:
        print(_format_gap_report(gaps))
        raise SystemExit(1)
    print("No CLUBB_core cross-module API coverage gaps found.")
    raise SystemExit(0)
