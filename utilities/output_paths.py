"""Small shared conventions for repository-owned runtime output paths."""

from __future__ import annotations

from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_ROOT = REPO_ROOT / "output"


def resolve_output_dir(value: str | Path | None) -> Path:
    """Resolve a user-facing output name without creating or deleting it.

    A bare relative name is rooted below the repository's ``output/``
    directory.  Existing explicit ``output/...`` spellings remain valid, as
    do absolute paths for advanced external workflows.  Relative traversal is
    rejected so a convenient output label cannot escape the output tree.
    """
    if value is None or not str(value).strip():
        return OUTPUT_ROOT

    requested = Path(str(value).strip()).expanduser()
    if requested.is_absolute():
        return requested
    if ".." in requested.parts:
        raise ValueError("relative output directories may not contain '..'")
    if requested.parts and requested.parts[0] == "output":
        return REPO_ROOT / requested
    return OUTPUT_ROOT / requested
