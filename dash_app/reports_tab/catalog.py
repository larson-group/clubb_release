"""Find PDF, Markdown and plain-text documents for the dashboard's Documents tab.

Discovery is recursive under doc/, excluding hidden paths.
Paths relative to doc/ identify documents; no publication index is required.
Use discover_documents for the list, resolve_document before reading a selection,
and catalog_token to avoid refreshing an unchanged list. resolve_linked_file
confines linked-file reads to visible doc/ files or Git-tracked repository files.
"""

from dataclasses import dataclass
from hashlib import sha256
import os
import subprocess
from pathlib import Path
from urllib.parse import quote

DOC_ROOT = Path(__file__).resolve().parents[2] / "doc"
DOCUMENT_TYPES = {".pdf": "PDFs", ".md": "Markdown", ".txt": "Text"}


@dataclass(frozen=True)
class Document:
    """A displayable file and its version, identified relative to doc/."""

    path: str
    version: str

    @property
    def label(self):
        """Show the full relative path without a redundant format prefix."""
        return self.path

    @property
    def url(self):
        """Address this file version through the document-only Flask route."""
        return f"/documents/{quote(self.path, safe='/')}?v={self.version}"


def _resolve_file(relative_path, root, suffixes):
    """Confine a regular file to visible paths inside the requested root."""
    path = Path(str(relative_path or ""))
    if path.is_absolute() or any(part.startswith('.') for part in path.parts):
        return None
    if suffixes is not None and path.suffix.lower() not in suffixes:
        return None
    base = Path(root) if root is not None else DOC_ROOT
    candidate = base / path
    try:
        resolved = candidate.resolve().relative_to(base.resolve())
        if any(part.startswith('.') for part in resolved.parts):
            return None
        return candidate if candidate.is_file() else None
    except (OSError, RuntimeError, ValueError):
        return None


def resolve_document(relative_path, root=None):
    """Return a PDF, Markdown or plain-text file inside doc/, including READMEs."""
    return _resolve_file(relative_path, root, set(DOCUMENT_TYPES))


def resolve_image(relative_path, root=None):
    """Resolve a local Markdown image without exposing scripts or other files."""
    return _resolve_file(relative_path, root, {'.png', '.jpg', '.jpeg', '.gif', '.webp', '.svg'})


def resolve_linked_file(relative_path):
    """Resolve a visible doc/ file or Git-tracked repository file for read-only links."""
    root = DOC_ROOT.parent
    candidate = _resolve_file(relative_path, root, None)
    if candidate is None:
        return None
    resolved = candidate.resolve()
    if resolved.is_relative_to(DOC_ROOT.resolve()):
        return candidate
    try:
        tracked = subprocess.run(
            ['git', '--literal-pathspecs', 'ls-files', '--error-unmatch', '--', resolved.relative_to(root.resolve()).as_posix()],
            cwd=root, capture_output=True, timeout=3, check=False,
        )
    except (OSError, subprocess.TimeoutExpired):
        return None
    return candidate if tracked.returncode == 0 else None


def discover_documents(root=None, query=""):
    """List matching PDF/Markdown/plain-text paths, without opening or interpreting other files."""
    base = Path(root) if root is not None else DOC_ROOT
    words = str(query or "").casefold().split()
    documents = []
    for directory, folders, files in os.walk(base):
        folders[:] = sorted(name for name in folders if not name.startswith('.'))
        for name in sorted(files):
            relative = (Path(directory) / name).relative_to(base).as_posix()
            if not all(word in relative.casefold() for word in words):
                continue
            path = resolve_document(relative, base)
            if path is None:
                continue
            try:
                stat = path.stat()
            except OSError:
                continue
            documents.append(Document(relative, f"{stat.st_mtime_ns}-{stat.st_ctime_ns}-{stat.st_size}-{stat.st_ino}"))
    type_order = {suffix: index for index, suffix in enumerate(DOCUMENT_TYPES)}
    return tuple(sorted(documents, key=lambda document: (
        type_order[Path(document.path).suffix.lower()], document.path.casefold()
    )))


def catalog_token(documents):
    """Identify a list and its file versions so unchanged polls do no UI work."""
    return sha256(repr(documents).encode()).hexdigest()[:16]
