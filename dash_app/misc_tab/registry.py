"""Registration for the self-contained pages within the Misc tab."""

from __future__ import annotations

from dataclasses import dataclass, field
from importlib import import_module
import pkgutil
from typing import Callable, Iterable

from dash.development.base_component import Component


LayoutBuilder = Callable[[], Component]
CallbackRegistrar = Callable[[object], None]


@dataclass(frozen=True)
class SubtabSpec:
    """Metadata and entrypoints for one self-contained Misc subtab."""

    slug: str
    title: str
    summary: str
    build_layout: LayoutBuilder
    category: str = "Investigation"
    updated: str = ""
    tags: tuple[str, ...] = field(default_factory=tuple)
    order: int = 100
    register_callbacks: CallbackRegistrar | None = None
    page_value: str | None = None


_SUBTABS: dict[str, SubtabSpec] = {}
_DISCOVERED = False


def register_subtab(subtab: SubtabSpec) -> SubtabSpec:
    """Register one Misc subtab, rejecting ambiguous duplicate slugs."""

    slug = subtab.slug.strip()
    if not slug or slug != subtab.slug:
        raise ValueError("A report slug must be a non-empty, trimmed string.")
    previous = _SUBTABS.get(slug)
    if previous is not None:
        # ``dash_app/app.py`` supports script-style imports while tests and
        # services use package imports.  The same subtab may therefore be
        # loaded through both module paths in one interpreter.  A slug is the
        # stable identity; keep the first registration rather than making the
        # import style observable to users.
        return previous
    _SUBTABS[slug] = subtab
    return subtab


def discover_subtabs() -> tuple[SubtabSpec, ...]:
    """Import each direct Misc subtab package and return its registered specs."""

    global _DISCOVERED
    if not _DISCOVERED:
        package = import_module(__package__)
        for module in pkgutil.iter_modules(package.__path__):
            if module.ispkg and not module.name.startswith("_"):
                import_module(f"{package.__name__}.{module.name}")
        _DISCOVERED = True
    return tuple(
        sorted(_SUBTABS.values(), key=lambda subtab: (subtab.order, subtab.title.lower()))
    )


def register_subtab_callbacks(app, subtabs: Iterable[SubtabSpec] | None = None) -> None:
    """Register optional callbacks owned by discovered Misc subtabs."""

    for subtab in subtabs if subtabs is not None else discover_subtabs():
        if subtab.register_callbacks is not None:
            subtab.register_callbacks(app)
