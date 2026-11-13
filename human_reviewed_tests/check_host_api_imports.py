#!/usr/bin/env python3
"""Check that host models access CLUBB and SILHS through their API modules.

Run from the repository root, using an existing host checkout:
    python3 human_reviewed_tests/check_host_api_imports.py /path/to/host
    python3 human_reviewed_tests/check_host_api_imports.py /path/to/sam \
        -exclude_dir SRC/CLUBB SRC/SILHS

Usage:
    check_host_api_imports(clubb_source, host_source, exclude_dirs=()): return
        (file count, private-module count, violations). Each violation is a
        (path, starting line number, module name) tuple.
    main(argv=None): print violations and a summary; return 0 for a complete
        passing scan, 1 for forbidden imports, or 2 for incomplete input.

The module inventory comes from actual MODULE declarations in the selected
CLUBB checkout's src/CLUBB_core and src/SILHS directories. Only clubb_api_module
and silhs_api_module may be imported directly by a host. Other host-owned
modules and intrinsic modules are outside this check.

The scanner reads .f/.for/.ftn fixed-form and .f90/.f95/.f03/.f08 free-form
files, ignoring .git directories. Exclusions name exact directories relative
to the host root, typically its vendored CLUBB
and SILHS sources. Files are read without modification; legacy non-UTF8
comment bytes are preserved. Comments, strings, continuations and semicolon
statements are normalized by utilities/CLUBBStandardsCheck.py. All preprocessor
branches are scanned; include files are not expanded. A pass establishes the
source-level import boundary, not argument compatibility or numerical agreement.

The default CLUBB checkout contains this script; -clubb_source selects another.
Git/checkout failures belong to the caller. Missing or empty source trees and
read failures are errors, never passing scans. Jenkins supplies SAM, CAM and WRF checkouts;
this script has no host-specific paths and needs only the Python standard library.
"""

import argparse
import os
from pathlib import Path
import re
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
if __name__ == "__main__":
    sys.path.insert(0, str(REPO_ROOT))

from utilities.CLUBBStandardsCheck import iter_logical_statements

FIXED_SUFFIXES = {".f", ".for", ".ftn"}
SOURCE_SUFFIXES = FIXED_SUFFIXES | {".f90", ".f95", ".f03", ".f08"}
API_MODULES = {"clubb_api_module", "silhs_api_module"}
MODULE = re.compile(r"module\s+([a-z]\w*)", re.IGNORECASE)
USE = re.compile(
    r"use(?:\s*,\s*(intrinsic|non_intrinsic)\s*::|\s*::|\s+)\s*([a-z]\w*)\b",
    re.IGNORECASE,
)
FIXED_USE = re.compile(r"use(?:(?:,(intrinsic|non_intrinsic))?::)?([a-z]\w*)(?:,.*)?", re.IGNORECASE)


def _source_files(root, excluded=()):
    """Yield sorted Fortran source paths, pruning exact excluded directories."""
    if root.is_file():
        if root.suffix.lower() not in SOURCE_SUFFIXES:
            raise ValueError(f"Not a Fortran source file: {root}")
        yield root
        return
    if not root.is_dir():
        raise ValueError(f"Source path does not exist: {root}")

    def walk_error(error):
        """Propagate directory read errors instead of silently omitting source."""
        raise error

    for directory, dirs, files in os.walk(root, onerror=walk_error):
        directory = Path(directory)
        dirs[:] = sorted(name for name in dirs if name != ".git"
                         and (directory / name).relative_to(root) not in excluded)
        for name in sorted(files):
            path = directory / name
            if path.suffix.lower() in SOURCE_SUFFIXES:
                yield path


def _statements(path):
    """Read source statements with original line numbers and source-form handling."""
    with path.open(encoding="utf-8", errors="surrogateescape") as source:
        yield from iter_logical_statements(source, fixed_form=path.suffix.lower() in FIXED_SUFFIXES)


def _private_modules(clubb_source):
    """Discover core/SILHS modules and require both public APIs and private source."""
    modules = set()
    for name in ("CLUBB_core", "SILHS"):
        found = set()
        for path in _source_files(clubb_source / "src" / name):
            for _, _, statement in _statements(path):
                match = MODULE.fullmatch(statement)
                if match:
                    found.add(match[1].lower())
        if not found:
            raise ValueError(f"No module declarations found in {clubb_source / 'src' / name}")
        modules.update(found)
    if not API_MODULES <= modules or not modules - API_MODULES:
        raise ValueError(f"Incomplete CLUBB/SILHS module inventory in {clubb_source}")
    return modules - API_MODULES


def check_host_api_imports(clubb_source, host_source, exclude_dirs=()):
    """Return scan counts and forbidden imports; raise on incomplete source input."""
    clubb_source, host_source = Path(clubb_source).resolve(), Path(host_source).resolve()
    excluded = {Path(name) for name in exclude_dirs}
    for path in excluded:
        if path.is_absolute() or ".." in path.parts or not (host_source / path).is_dir():
            raise ValueError(f"Exclusion must name an existing directory relative to {host_source}: {path}")
        if path == Path("."):
            raise ValueError("Cannot exclude the entire host source tree")
    private = _private_modules(clubb_source)
    sources = list(_source_files(host_source, excluded))
    if not sources:
        raise ValueError(f"No Fortran source files to check in {host_source}")
    violations = []
    for path in sources:
        for line, _, statement in _statements(path):
            if path.suffix.lower() in FIXED_SUFFIXES:
                # Fixed-form blanks are insignificant, even within an identifier.
                match = FIXED_USE.fullmatch(re.sub(r"\s+", "", statement))
            else:
                match = USE.match(statement)
            if match and (match[1] or "").lower() != "intrinsic" and match[2].lower() in private:
                violations.append((path, line, match[2].lower()))
    return len(sources), len(private), violations


def main(argv=None):
    """Check CLI-selected source trees and return pass, violation or input-error status."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                     add_help=False, allow_abbrev=False)
    parser.add_argument("-h", "-help", action="help")
    parser.add_argument("host_source", type=Path, help="Host source directory or one Fortran file.")
    parser.add_argument("-clubb_source", type=Path, default=REPO_ROOT,
                        help="CLUBB checkout providing module declarations (default: this checkout).")
    parser.add_argument("-exclude_dir", nargs="+", default=[], metavar="PATH",
                        help="Exact directory paths relative to host_source to omit.")
    args = parser.parse_args(argv)
    try:
        files, modules, violations = check_host_api_imports(
            args.clubb_source, args.host_source, args.exclude_dir)
    except (OSError, ValueError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 2
    for path, line, module in violations:
        print(f"{path}:{line}: forbidden import of {module}; use the CLUBB/SILHS API")
    print(f"Checked {files} Fortran files against {modules} private modules: "
          f"{len(violations)} forbidden imports.")
    return 1 if violations else 0


if __name__ == "__main__":
    sys.exit(main())
