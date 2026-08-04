"""Validated discovery for static report bundles under ``doc/reports``."""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path
import re
from typing import Any
from urllib.parse import quote


REPO_ROOT = Path(__file__).resolve().parents[2]
REPORTS_ROOT = REPO_ROOT / "doc" / "reports"
INDEX_PATH = REPORTS_ROOT / "index.json"
INDEX_SCHEMA_VERSION = 1
_REPORT_ID_RE = re.compile(r"^[a-z0-9][a-z0-9-]{0,95}$")


@dataclass(frozen=True)
class StaticReport:
    """One published, browser-displayable report bundle."""

    report_id: str
    title: str
    summary: str
    created_at: str
    tags: tuple[str, ...]
    entrypoint: str = "report.html"
    source_revision: str = ""

    @property
    def tab_value(self) -> str:
        return report_tab_value(self.report_id)

    @property
    def url(self) -> str:
        return report_url(self.report_id, self.entrypoint)


def report_tab_value(report_id: str) -> str:
    """Return the stable Dash selection value for a published report."""
    return f"static-report-{report_id}"


def report_url(report_id: str, entrypoint: str = "report.html") -> str:
    """Return the app route for an already-validated bundle asset."""
    return f"/static-reports/{quote(report_id, safe='')}/{quote(entrypoint, safe='/')}"


def _valid_report_id(value: Any) -> str | None:
    report_id = str(value or "").strip()
    return report_id if _REPORT_ID_RE.fullmatch(report_id) else None


def _valid_entrypoint(value: Any) -> str | None:
    entrypoint = str(value or "report.html").strip()
    candidate = Path(entrypoint)
    if not entrypoint or candidate.is_absolute() or ".." in candidate.parts:
        return None
    return entrypoint


def _read_index(root: Path) -> dict[str, Any]:
    try:
        payload = json.loads((root / "index.json").read_text(encoding="utf-8"))
    except (OSError, TypeError, ValueError):
        return {"schema_version": INDEX_SCHEMA_VERSION, "reports": []}
    if not isinstance(payload, dict) or not isinstance(payload.get("reports"), list):
        return {"schema_version": INDEX_SCHEMA_VERSION, "reports": []}
    return payload


def discover_reports(root: Path | str | None = None) -> tuple[StaticReport, ...]:
    """Read only complete, indexed report bundles; never import report code."""
    reports_root = Path(root) if root is not None else REPORTS_ROOT
    reports: list[StaticReport] = []
    seen: set[str] = set()
    for raw in _read_index(reports_root).get("reports", []):
        if not isinstance(raw, dict):
            continue
        report_id = _valid_report_id(raw.get("id"))
        entrypoint = _valid_entrypoint(raw.get("entrypoint"))
        title = str(raw.get("title") or "").strip()
        if not report_id or not entrypoint or not title or report_id in seen:
            continue
        asset = resolve_report_asset(report_id, entrypoint, root=reports_root)
        if asset is None:
            continue
        seen.add(report_id)
        tags = tuple(
            str(tag).strip()
            for tag in raw.get("tags", [])
            if isinstance(tag, str) and str(tag).strip()
        )
        reports.append(
            StaticReport(
                report_id=report_id,
                title=title,
                summary=str(raw.get("summary") or "").strip(),
                created_at=str(raw.get("created_at") or "").strip(),
                tags=tags,
                entrypoint=entrypoint,
                source_revision=str(raw.get("source_revision") or "").strip(),
            )
        )
    return tuple(reports)


def catalog_token(reports: tuple[StaticReport, ...] | list[StaticReport]) -> str:
    """Return a stable change token without using report file mtimes."""
    values = [
        {
            "id": report.report_id,
            "title": report.title,
            "summary": report.summary,
            "created_at": report.created_at,
            "tags": list(report.tags),
            "entrypoint": report.entrypoint,
            "source_revision": report.source_revision,
        }
        for report in reports
    ]
    encoded = json.dumps(values, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return sha256(encoded).hexdigest()[:16]


def resolve_report_asset(
    report_id: str,
    asset_path: str,
    *,
    root: Path | str | None = None,
) -> Path | None:
    """Return a regular asset inside a bundle, rejecting path traversal."""
    safe_id = _valid_report_id(report_id)
    safe_asset = _valid_entrypoint(asset_path)
    if not safe_id or not safe_asset:
        return None
    reports_root = Path(root) if root is not None else REPORTS_ROOT
    bundle = (reports_root / safe_id).resolve()
    candidate = (bundle / safe_asset).resolve()
    try:
        candidate.relative_to(bundle)
    except ValueError:
        return None
    return candidate if candidate.is_file() else None


def report_by_id(report_id: str, root: Path | str | None = None) -> StaticReport | None:
    """Find one indexed report by report id."""
    return next((report for report in discover_reports(root) if report.report_id == report_id), None)
