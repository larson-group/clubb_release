#!/usr/bin/env python3
"""Add async/wait clauses to selected OpenACC directives."""

from __future__ import annotations

import argparse
import re
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


# Directories converted by default. Paths are relative to the CLUBB repo root.
# Each directory is scanned non-recursively for *.F90 files, matching the old
# convert_to_async.bash behavior.
DIRECTORIES_TO_CONVERT = [
    "src",
    "src/CLUBB_core",
    "src/Radiation",
    "src/Benchmark_cases",
    "src/Microphys",
    "src/SILHS",
]

FORTRAN_FILE_PATTERN = "*.F90"

DEFAULT_WORKERS = 12
ASYNC_CLAUSE = "async(1)"
WAIT_CLAUSE = "wait"

REPO_ROOT = Path(__file__).resolve().parents[1]
ACC_DIRECTIVE_RE = re.compile(r"^\s*!\$acc\b", re.IGNORECASE)
ASYNC_1_RE = re.compile(r"(?i)(?<![\w])async\s*\(\s*1\s*\)(?![\w])")
WAIT_RE = re.compile(r"(?i)(?<![\w])wait(?:\s*\([^)]*\))?(?![\w])")


def read_file(file_path: Path) -> list[str]:
    """Read a source file as lines so directive spans can be rewritten in place."""
    with file_path.open("r") as file:
        return file.readlines()


def write_file(file_path: Path, lines: list[str]) -> None:
    """Write converted lines back without changing the surrounding file workflow."""
    with file_path.open("w") as file:
        file.writelines(lines)


def find_directives(lines: list[str]) -> list[tuple[int, int]]:
    """Group continued OpenACC directives into single conversion units."""
    directive_spans = []
    start_idx = None

    for idx, line in enumerate(lines):
        if ACC_DIRECTIVE_RE.match(line):
            if start_idx is None:
                start_idx = idx

        if start_idx is not None and not line.rstrip().endswith("&"):
            directive_spans.append((start_idx, idx + 1))
            start_idx = None

    if start_idx is not None:
        directive_spans.append((start_idx, len(lines)))

    return directive_spans


def contains_keyword(directive_lines: list[str], keyword: str) -> bool:
    """Check directive words using the script's coarse OpenACC token model."""
    return any(keyword.lower() in token.lower() for token in " ".join(directive_lines).split())


def contains_clause(directive_lines: list[str], clause: str) -> bool:
    """Detect an exact clause before adding a duplicate synchronization marker."""
    return clause.lower() in (token.lower() for token in " ".join(directive_lines).split())


def clean_directive_line(line: str) -> str:
    """Normalize whitespace left behind by removing clauses from a directive."""
    newline = "\n" if line.endswith("\n") else ""
    body = line[:-1] if newline else line
    body = re.sub(r"[ \t]+(&\s*)$", r" \1", body)
    body = re.sub(r"[ \t]+$", "", body)
    return body + newline


def remove_clause(directive_lines: list[str], clause_re: re.Pattern[str]) -> list[str]:
    """Strip clauses that conflict with the directive's assigned synchronization role."""
    return [clean_directive_line(clause_re.sub("", line)) for line in directive_lines]


def add_clause(directive_lines: list[str], clause: str) -> list[str]:
    """Append a clause to the directive tail, where OpenACC continuation is complete."""
    if contains_clause(directive_lines, clause):
        return directive_lines

    updated_lines = directive_lines.copy()
    last_line = updated_lines[-1]
    newline = "\n" if last_line.endswith("\n") else ""
    body = last_line[:-1] if newline else last_line
    updated_lines[-1] = f"{body.rstrip()} {clause}{newline}"
    return updated_lines


def wait_directive_after(directive_lines: list[str]) -> list[str]:
    """Insert a standalone wait that preserves the directive block's indentation."""
    indent = re.match(r"^(\s*)", directive_lines[0]).group(1)
    return directive_lines + [f"{indent}!$acc wait\n"]


def process_directive(
    directive_lines: list[str],
    pending_wait_end_keyword: str | None = None,
) -> tuple[list[str], str | None]:
    """Modify one OpenACC directive when async/wait rules apply."""
    # Sequential loops are intentionally synchronous regions. Leaving async on
    # them can change ordering assumptions in scalar or dependency-sensitive code.
    if contains_keyword(directive_lines, "seq"):
        directive_lines = remove_clause(directive_lines, ASYNC_1_RE)
        return directive_lines, pending_wait_end_keyword

    # Some constructs need a wait after the whole OpenACC region, not on the
    # opening directive. This closes the pending region once its matching end
    # directive appears.
    if (
        pending_wait_end_keyword is not None
        and contains_keyword(directive_lines, "end")
        and contains_keyword(directive_lines, pending_wait_end_keyword)
    ):
        return wait_directive_after(directive_lines), None

    # Reductions may run asynchronously, but the reduced scalar is often read by
    # nearby host control flow. Waiting after the region keeps that contract.
    if (    contains_keyword(directive_lines, "reduction")
        and contains_keyword(directive_lines, "parallel")):
        directive_lines = remove_clause(directive_lines, WAIT_RE)
        directive_lines = add_clause(directive_lines, ASYNC_CLAUSE)
        return directive_lines, "parallel"

    # host_data often wraps CUDA library calls such as cuRAND. The library call
    # is not itself an OpenACC async kernel, so we wait after the wrapper before
    # later OpenACC work consumes the device data it produced.
    if (
        contains_keyword(directive_lines, "host_data")
        and not contains_keyword(directive_lines, "end")
    ):
        return directive_lines, "host_data"

    # Data motion directives define visibility boundaries between host and
    # device. Keeping them synchronous avoids stale host reads and incomplete
    # device copies at CPU/GPU handoff points.
    if (   contains_keyword(directive_lines, "copy")
        or contains_keyword(directive_lines, "copyout")
        or contains_keyword(directive_lines, "exit")
        or contains_keyword(directive_lines, "update")
    ):
        directive_lines = remove_clause(directive_lines, ASYNC_1_RE)
        directive_lines = add_clause(directive_lines, WAIT_CLAUSE)
        return directive_lines, pending_wait_end_keyword

    # Structural directives, routine metadata, and already-async regions are
    # treated as explicit author intent. The converter only changes directives
    # where the synchronization policy is clear.
    if (   contains_keyword(directive_lines, "end")
        or contains_keyword(directive_lines, "declare")
        or contains_keyword(directive_lines, "routine")
        or contains_keyword(directive_lines, "host_data")
        or contains_keyword(directive_lines, "async")
    ):
        # leave anything with these unchanged
        return directive_lines, pending_wait_end_keyword

    # Compute/data lifetime regions are the main target: putting these on one
    # async queue lets adjacent GPU work overlap while preserving queue order.
    if (   contains_keyword(directive_lines, "enter")
        or contains_keyword(directive_lines, "data")
        or (
                contains_keyword(directive_lines, "parallel")
            and contains_keyword(directive_lines, "loop")
        )
    ):
        directive_lines = remove_clause(directive_lines, WAIT_RE)
        directive_lines = add_clause(directive_lines, ASYNC_CLAUSE)

    return directive_lines, pending_wait_end_keyword


def convert_lines(original_lines: list[str]) -> list[str]:
    """Convert one file's lines while preserving all non-directive text."""
    modified_lines = []
    pending_wait_end_keyword = None
    previous_end = 0

    for start, end in find_directives(original_lines):
        original_directive = original_lines[start:end]
        modified_directive, pending_wait_end_keyword = process_directive(
            original_directive,
            pending_wait_end_keyword,
        )
        modified_lines.extend(original_lines[previous_end:start])
        modified_lines.extend(modified_directive)
        previous_end = end

    modified_lines.extend(original_lines[previous_end:])

    return modified_lines


def process_file(file_path: Path) -> bool:
    """Convert a single file and report whether its contents changed."""
    original_lines = read_file(file_path)
    modified_lines = convert_lines(original_lines)
    changed = modified_lines != original_lines

    if changed:
        write_file(file_path, modified_lines)

    return changed


def repo_path(path: str) -> Path:
    """Resolve command-line paths relative to the CLUBB repository root."""
    candidate = Path(path)
    if candidate.is_absolute():
        return candidate
    return REPO_ROOT / candidate


def files_in_directory(directory: Path) -> list[Path]:
    """List the Fortran sources directly owned by one configured directory."""
    return sorted(path for path in directory.glob(FORTRAN_FILE_PATTERN) if path.is_file())


def unique_paths(paths: list[Path]) -> list[Path]:
    """Deduplicate files so overlapping inputs do not race in the worker pool."""
    seen = set()
    unique = []

    for path in paths:
        resolved = path.resolve()
        if resolved not in seen:
            seen.add(resolved)
            unique.append(path)

    return unique


def paths_to_files(paths: list[str]) -> tuple[list[Path], list[Path]]:
    """Expand user-provided files/directories into the conversion work list."""
    files_to_convert = []
    missing_paths = []

    for path_string in paths:
        path = repo_path(path_string)
        if path.is_dir():
            files_to_convert.extend(files_in_directory(path))
        elif path.is_file():
            files_to_convert.append(path)
        else:
            missing_paths.append(path)

    return unique_paths(files_to_convert), missing_paths


def display_path(path: Path) -> Path:
    """Show repo-relative paths in user-facing output when possible."""
    try:
        return path.relative_to(REPO_ROOT)
    except ValueError:
        return path


def parse_args() -> argparse.Namespace:
    """Define the small CLI used by Jenkins and local conversion checks."""
    parser = argparse.ArgumentParser(
        description="Add async/wait clauses to configured OpenACC Fortran files.",
    )
    parser.add_argument(
        "paths",
        nargs="*",
        help="Optional files or directories to convert instead of DIRECTORIES_TO_CONVERT.",
    )
    parser.add_argument(
        "-nproc",
        type=int,
        default=DEFAULT_WORKERS,
        help=f"Number of files to process in parallel. Default: {DEFAULT_WORKERS}.",
    )
    return parser.parse_args()


def main() -> int:
    """Coordinate path expansion, parallel conversion, and concise reporting."""
    args = parse_args()
    configured_paths = args.paths if args.paths else DIRECTORIES_TO_CONVERT
    files_to_convert, missing_paths = paths_to_files(configured_paths)

    if missing_paths:
        for path in missing_paths:
            print(f"Error: path does not exist: {path}")
        return 1

    if not files_to_convert:
        print("No Fortran files found.")
        return 0

    workers = max(1, args.nproc)
    with ThreadPoolExecutor(max_workers=workers) as executor:
        changed = list(
            executor.map(
                process_file,
                files_to_convert,
            ),
        )

    changed_files = [path for path, did_change in zip(files_to_convert, changed) if did_change]

    for path in changed_files:
        print(f"Updated: {display_path(path)}")

    if not changed_files:
        print("No files changed.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
