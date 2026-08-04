"""Dash layout for the reload-safe static report library."""

from __future__ import annotations

from dash import dcc, html

from .catalog import StaticReport, catalog_token, discover_reports


EMPTY_VALUE = "reports-library"
POLL_INTERVAL_MS = 3000


def report_tabs(reports: tuple[StaticReport, ...] | list[StaticReport]):
    """Build navigation only; the iframe is kept separate to avoid reloads on polling."""
    if reports:
        return [
            dcc.Tab(
                label=report.title,
                value=report.tab_value,
                className="reports-directory-card",
                selected_className="reports-directory-card-selected",
            )
            for report in reports
        ]
    return [
        dcc.Tab(
            label="Getting started",
            value=EMPTY_VALUE,
            className="reports-directory-card",
            selected_className="reports-directory-card-selected",
        )
    ]


def selected_value(reports: tuple[StaticReport, ...] | list[StaticReport], current: str | None = None) -> str:
    values = {report.tab_value for report in reports}
    if current in values:
        return str(current)
    return reports[0].tab_value if reports else EMPTY_VALUE


def build_layout(reports: tuple[StaticReport, ...] | None = None):
    """Build a left directory and stable iframe viewer for static report bundles."""
    report_specs = reports if reports is not None else discover_reports()
    current = selected_value(report_specs)
    first_url = report_specs[0].url if report_specs else "about:blank"
    return html.Div(
        [
            dcc.Interval(id="reports-catalog-poll", interval=POLL_INTERVAL_MS, n_intervals=0),
            dcc.Store(id="reports-catalog-token", data=catalog_token(report_specs)),
            html.Aside(
                [
                    html.Div(
                        [
                            html.Div("Evidence library", className="reports-directory-eyebrow"),
                            html.H2("Reports"),
                            html.P(
                                "Published investigations remain available while new bundles appear here.",
                                className="reports-directory-intro",
                            ),
                        ],
                        className="reports-directory-header",
                    ),
                    dcc.Tabs(
                        id="reports-pages",
                        value=current,
                        vertical=True,
                        persistence=True,
                        persistence_type="local",
                        className="reports-directory",
                        children=report_tabs(report_specs),
                    ),
                ],
                className="reports-directory-shell",
            ),
            html.Section(
                [
                    html.Div(
                        [
                            html.H2("No static reports yet"),
                            html.P(
                                "Publish a completed bundle under doc/reports; this page will add it automatically without restarting Dash."
                            ),
                        ],
                        id="reports-empty-state",
                        className="reports-empty-state" if not report_specs else "reports-empty-state reports-empty-state-hidden",
                    ),
                    html.Iframe(
                        id="reports-frame",
                        src=first_url,
                        title="Static CLUBB report",
                        className="reports-frame" if report_specs else "reports-frame reports-frame-hidden",
                        sandbox="allow-same-origin",
                    ),
                ],
                className="reports-view",
            ),
        ],
        className="reports-tab-root",
    )
