"""Small stdlib-first builder for atomic, self-contained static report bundles."""

from __future__ import annotations

from datetime import datetime, timezone
from html import escape
import json
import os
from pathlib import Path
import re
import shutil
from typing import Iterable
from uuid import uuid4

from .catalog import INDEX_SCHEMA_VERSION, REPORTS_ROOT, StaticReport, _valid_report_id


_SAFE_FILE_RE = re.compile(r"[^A-Za-z0-9_.-]+")
MAX_RASTER_DIMENSION = 1200
_RASTER_SUFFIXES = {".jpg", ".jpeg", ".png"}


def _atomic_write_json(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{uuid4().hex}.tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def _safe_filename(value: str, fallback: str) -> str:
    name = _SAFE_FILE_RE.sub("-", Path(value).name).strip(".-")
    return name or fallback


def cap_raster_image(
    source: Path | str,
    destination: Path | str,
    *,
    max_dimension: int = MAX_RASTER_DIMENSION,
) -> bool:
    """Copy a report raster while bounding its largest pixel dimension.

    Returns whether a resize occurred. SVG remains vector evidence and is
    intentionally outside this policy. Keeping this operation here makes the
    report bundle itself—not a particular plotting script—the enforcement
    boundary for all future static reports.
    """
    source_path = Path(source)
    destination_path = Path(destination)
    if source_path.suffix.lower() not in _RASTER_SUFFIXES:
        shutil.copy2(source_path, destination_path)
        return False
    if max_dimension < 1:
        raise ValueError("max_dimension must be positive")

    try:
        from PIL import Image
    except ModuleNotFoundError as exc:
        raise RuntimeError("Pillow is required to publish raster report figures") from exc

    with Image.open(source_path) as image:
        width, height = image.size
        if max(width, height) <= max_dimension:
            shutil.copy2(source_path, destination_path)
            return False
        resized = image.copy()
        resized.thumbnail((max_dimension, max_dimension), Image.Resampling.LANCZOS)
        image_format = image.format or ("JPEG" if source_path.suffix.lower() in {".jpg", ".jpeg"} else "PNG")

    if image_format.upper() == "JPEG" and resized.mode not in {"RGB", "L"}:
        resized = resized.convert("RGB")
    temporary = destination_path.with_name(f".{destination_path.name}.{uuid4().hex}.tmp")
    save_options = {"optimize": True}
    if image_format.upper() == "JPEG":
        save_options["quality"] = 92
    try:
        resized.save(temporary, format=image_format, **save_options)
        shutil.copystat(source_path, temporary, follow_symlinks=False)
        os.replace(temporary, destination_path)
    finally:
        if temporary.exists():
            temporary.unlink()
    return True


class ReportBuilder:
    """Build a polished report in staging, then publish it with one atomic move.

    The builder intentionally writes only beneath ``doc/reports``.  It creates
    plain HTML, figures, snippets, and JSON; publishing one report therefore
    never edits Dash source files or triggers its Python reloader.
    """

    def __init__(
        self,
        report_id: str,
        title: str,
        *,
        summary: str,
        tags: Iterable[str] = (),
        replace: bool = False,
        source_revision: str = "",
        root: Path | str | None = None,
    ):
        safe_id = _valid_report_id(report_id)
        if not safe_id:
            raise ValueError("report_id must use lowercase letters, digits, and hyphens")
        if not str(title).strip() or not str(summary).strip():
            raise ValueError("title and summary are required")
        self.root = Path(root) if root is not None else REPORTS_ROOT
        self.report_id = safe_id
        self.title = str(title).strip()
        self.summary = str(summary).strip()
        self.tags = tuple(dict.fromkeys(str(tag).strip() for tag in tags if str(tag).strip()))
        self.replace = bool(replace)
        self.source_revision = str(source_revision).strip()
        self.created_at = datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")
        self.stage = self.root / ".staging" / f"{self.report_id}-{uuid4().hex}"
        (self.stage / "figures").mkdir(parents=True, exist_ok=False)
        (self.stage / "snippets").mkdir(parents=True, exist_ok=True)
        self._sections: list[str] = []
        self._asset_counts: dict[str, int] = {}

    def heading(self, title: str, *, level: int = 2) -> "ReportBuilder":
        level = min(4, max(2, int(level)))
        self._sections.append(f"<h{level}>{escape(str(title))}</h{level}>")
        return self

    def paragraph(self, text: str) -> "ReportBuilder":
        self._sections.append(f"<p>{escape(str(text))}</p>")
        return self

    def callout(self, title: str, text: str, *, tone: str = "note") -> "ReportBuilder":
        safe_tone = tone if tone in {"note", "success", "warning"} else "note"
        self._sections.append(
            f'<section class="callout callout-{safe_tone}"><h3>{escape(str(title))}</h3><p>{escape(str(text))}</p></section>'
        )
        return self

    def metrics(self, values: Iterable[tuple[str, str, str]]) -> "ReportBuilder":
        cards = "".join(
            f'<div class="metric"><span>{escape(str(label))}</span><strong>{escape(str(value))}</strong><small>{escape(str(detail))}</small></div>'
            for label, value, detail in values
        )
        self._sections.append(f'<div class="metrics">{cards}</div>')
        return self

    def figure(self, source: Path | str, *, caption: str, alt: str = "") -> "ReportBuilder":
        source_path = Path(source)
        if not source_path.is_file():
            raise FileNotFoundError(source_path)
        filename = self._next_asset_name(_safe_filename(source_path.name, "figure"))
        destination = self.stage / "figures" / filename
        cap_raster_image(source_path, destination)
        self._sections.append(
            '<figure><img src="figures/{src}" alt="{alt}"><figcaption>{caption}</figcaption></figure>'.format(
                src=escape(filename, quote=True), alt=escape(alt or caption, quote=True), caption=escape(caption)
            )
        )
        return self

    def code(self, text: str, *, language: str = "text", caption: str = "") -> "ReportBuilder":
        heading = f"<figcaption>{escape(caption)}</figcaption>" if caption else ""
        self._sections.append(
            f'<figure class="code-figure">{heading}<pre><code class="language-{escape(language, quote=True)}">{escape(str(text))}</code></pre></figure>'
        )
        return self

    def code_file(self, source: Path | str, *, language: str = "text", caption: str = "") -> "ReportBuilder":
        source_path = Path(source)
        text = source_path.read_text(encoding="utf-8")
        filename = self._next_asset_name(_safe_filename(source_path.name, "snippet.txt"))
        shutil.copy2(source_path, self.stage / "snippets" / filename)
        return self.code(text, language=language, caption=caption or source_path.name)

    def equation(self, latex: str, *, caption: str = "") -> "ReportBuilder":
        """Render a MathText equation to portable SVG without browser JavaScript."""
        os.environ.setdefault("MPLCONFIGDIR", "/tmp/clubb_dash_report_math")
        from matplotlib.backends.backend_svg import FigureCanvasSVG
        from matplotlib.figure import Figure

        expression = str(latex).strip()
        if not expression:
            raise ValueError("equation cannot be empty")
        if not (expression.startswith("$") and expression.endswith("$")):
            expression = f"${expression}$"
        filename = self._next_asset_name("equation.svg")
        path = self.stage / "figures" / filename
        figure = Figure(figsize=(0.1, 0.1), facecolor="none")
        FigureCanvasSVG(figure)
        figure.text(0, 0, expression, color="#e8edf5", fontsize=17)
        figure.savefig(path, format="svg", transparent=True, bbox_inches="tight", pad_inches=0.06)
        self._sections.append(
            '<figure class="equation"><img src="figures/{src}" alt="{alt}">{caption}</figure>'.format(
                src=escape(filename, quote=True),
                alt=escape(latex, quote=True),
                caption=f"<figcaption>{escape(caption)}</figcaption>" if caption else "",
            )
        )
        return self

    def publish(self) -> StaticReport:
        """Make the completed report discoverable in a single atomic publication."""
        destination = self.root / self.report_id
        if destination.exists() and not self.replace:
            raise FileExistsError(
                f"report id already exists: {self.report_id}; pass replace=True to update it"
            )
        self.root.mkdir(parents=True, exist_ok=True)
        self._write_bundle_files()
        previous = None
        if destination.exists():
            previous = self.root / ".staging" / f"{self.report_id}-previous-{uuid4().hex}"
            os.replace(destination, previous)
        try:
            os.replace(self.stage, destination)
        except OSError:
            if previous is not None and previous.exists():
                os.replace(previous, destination)
            raise
        if previous is not None:
            shutil.rmtree(previous)
        index_path = self.root / "index.json"
        try:
            current = json.loads(index_path.read_text(encoding="utf-8"))
        except (OSError, TypeError, ValueError):
            current = {"schema_version": INDEX_SCHEMA_VERSION, "reports": []}
        rows = [
            row
            for row in current.get("reports", [])
            if isinstance(row, dict) and str(row.get("id") or "").strip() != self.report_id
        ]
        rows.insert(0, self._manifest())
        _atomic_write_json(index_path, {"schema_version": INDEX_SCHEMA_VERSION, "reports": rows})
        return StaticReport(
            self.report_id,
            self.title,
            self.summary,
            self.created_at,
            self.tags,
            source_revision=self.source_revision,
        )

    def _next_asset_name(self, filename: str) -> str:
        stem, suffix = Path(filename).stem, Path(filename).suffix
        count = self._asset_counts.get(filename, 0)
        self._asset_counts[filename] = count + 1
        return filename if count == 0 else f"{stem}-{count + 1}{suffix}"

    def _manifest(self) -> dict:
        return {
            "schema_version": INDEX_SCHEMA_VERSION,
            "id": self.report_id,
            "title": self.title,
            "summary": self.summary,
            "created_at": self.created_at,
            "tags": list(self.tags),
            "entrypoint": "report.html",
            "source_revision": self.source_revision,
        }

    def _write_bundle_files(self) -> None:
        _atomic_write_json(self.stage / "manifest.json", self._manifest())
        (self.stage / "report.html").write_text(
            _render_document(self.title, self.summary, self.created_at, self.tags, self.source_revision, self._sections),
            encoding="utf-8",
        )


def _render_document(
    title: str,
    summary: str,
    created_at: str,
    tags: tuple[str, ...],
    source_revision: str,
    sections: list[str],
) -> str:
    tag_html = "".join(f"<span class=\"tag\">{escape(tag)}</span>" for tag in tags)
    revision = f"<span>Revision: {escape(source_revision)}</span>" if source_revision else ""
    body = "\n".join(sections) or "<p>No analysis sections were added.</p>"
    return f"""<!doctype html>
<html lang=\"en\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">
<title>{escape(title)}</title><style>
:root {{ color-scheme: dark; --bg:#0d1525; --panel:#162237; --line:#30425f; --ink:#e8edf5; --muted:#a9b8cd; --accent:#38bdf8; --gold:#facc15; }}
* {{ box-sizing:border-box; }} body {{ margin:0; background:var(--bg); color:var(--ink); font:16px/1.65 Inter,ui-sans-serif,system-ui,sans-serif; }}
main {{ width:100%; max-width:none; margin:0; padding:48px clamp(22px,5vw,72px) 72px; }}
header {{ border-bottom:1px solid var(--line); padding-bottom:28px; margin-bottom:32px; }} h1 {{ margin:0 0 12px; font:700 clamp(2rem,5vw,3.4rem)/1.08 Georgia,serif; }} h2 {{ margin:40px 0 12px; font:700 1.7rem/1.2 Georgia,serif; }} h3 {{ margin:0 0 6px; }} p {{ color:#d9e1ed; }} .summary {{ max-width:800px; color:var(--muted); font-size:1.12rem; }} .meta {{ display:flex; flex-wrap:wrap; gap:9px 14px; color:var(--muted); font-size:.85rem; }} .tag {{ border:1px solid #356182; color:#bfe7ff; padding:2px 9px; border-radius:999px; font-weight:700; }}
.callout {{ border-left:4px solid var(--accent); padding:16px 19px; margin:24px 0; background:var(--panel); border-radius:0 10px 10px 0; }} .callout-success {{ border-color:#34d399; }} .callout-warning {{ border-color:var(--gold); }} .callout p {{ margin:0; color:var(--muted); }}
.metrics {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(170px,1fr)); gap:14px; margin:22px 0; }} .metric {{ border:1px solid var(--line); border-radius:11px; background:var(--panel); padding:16px; }} .metric span,.metric small {{ display:block; color:var(--muted); font-size:.82rem; }} .metric strong {{ display:block; margin:5px 0; color:var(--gold); font:700 1.6rem Georgia,serif; }}
figure {{ margin:26px 0; border:1px solid var(--line); border-radius:12px; background:#101a2c; overflow:hidden; }} figure img {{ display:block; width:100%; max-width:100%; background:#fff; }} figcaption {{ padding:11px 14px; color:var(--muted); font-size:.9rem; }} .equation {{ display:inline-block; min-width:250px; padding:18px; }} .equation img {{ width:auto; max-width:100%; background:transparent; }} .code-figure {{ background:#08111f; }} pre {{ margin:0; padding:18px; overflow:auto; color:#d9e8f7; font:13px/1.55 ui-monospace,SFMono-Regular,Consolas,monospace; }}
footer {{ margin-top:50px; padding-top:18px; border-top:1px solid var(--line); color:var(--muted); font-size:.85rem; }} @media print {{ body {{ background:#fff; color:#111; }} main {{ max-width:none; }} .metric,figure,.callout {{ break-inside:avoid; }} }}
</style></head><body><main><header><h1>{escape(title)}</h1><p class=\"summary\">{escape(summary)}</p><div class=\"meta\"><span>Published {escape(created_at)}</span>{tag_html}{revision}</div></header>{body}<footer>Static CLUBB investigation report</footer></main></body></html>"""
