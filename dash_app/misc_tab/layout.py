"""Left-rail catalog that combines research reports and small diagnostics."""

from __future__ import annotations

from collections.abc import Sequence

from dash import dcc, html

from dash_app.misc_tab.registry import SubtabSpec, discover_subtabs


def subtab_page_value(subtab: SubtabSpec) -> str:
    """Return the stable persisted page value for one Misc subtab."""
    return subtab.page_value or f"misc-{subtab.slug}"


def build_layout(subtabs: Sequence[SubtabSpec] | None = None):
    """Build the persistent left-side page rail for miscellaneous tools."""
    specs = tuple(subtabs) if subtabs is not None else discover_subtabs()
    first_value = subtab_page_value(specs[0]) if specs else "misc-empty"
    pages = [
        dcc.Tab(
            label=subtab.title,
            value=subtab_page_value(subtab),
            className="misc-directory-card",
            selected_className="misc-directory-card-selected",
            children=subtab.build_layout(),
        )
        for subtab in specs
    ]

    return html.Div(
        [
            html.Div(
                [
                    html.Div("Tool shelf", className="misc-directory-eyebrow"),
                    html.H2("Misc"),
                    html.P(
                        "Interactive investigations and focused diagnostics.",
                        className="misc-directory-intro",
                    ),
                ],
                className="misc-directory-header",
            ),
            dcc.Tabs(
                id="misc-pages",
                value=first_value,
                vertical=True,
                persistence=True,
                persistence_type="local",
                parent_className="misc-tabs-shell",
                className="misc-directory",
                content_className="misc-report-content",
                children=pages,
            ),
        ],
        # Existing laboratory-specific CSS stays scoped beneath this root;
        # Misc owns the enclosing navigation and all interactive pages.
        className="misc-labs-root misc-tab-root",
    )
