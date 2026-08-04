"""Entrypoint and callbacks for the static Reports dashboard tab."""

from __future__ import annotations

from dash import Input, Output, State, dcc, no_update
from urllib.parse import quote

from .catalog import catalog_token, discover_reports, report_by_id, report_tab_value
from .layout import build_layout, report_tabs, selected_value


def _report_id_from_value(value: str | None) -> str:
    prefix = "static-report-"
    selected = str(value or "")
    return selected.removeprefix(prefix) if selected.startswith(prefix) else ""


def _opened_report_url(report, request_id) -> str:
    """Give each explicit agent handoff a distinct iframe URL.

    A report can already be selected when an agent republishes a corrected
    artifact.  Updating ``dcc.Tabs.value`` to that same value does not trigger
    the normal selection callback, so attach the durable activity id as a
    harmless query parameter and explicitly reload the static iframe.
    """
    token = quote(str(request_id or "open"), safe="")
    return f"{report.url}?open={token}"


def build_tab(app):
    """Build the permanent tab once; report publication only changes data files."""

    @app.callback(
        Output("reports-pages", "children"),
        Output("reports-catalog-token", "data"),
        Output("reports-pages", "value", allow_duplicate=True),
        Input("reports-catalog-poll", "n_intervals"),
        State("reports-catalog-token", "data"),
        State("reports-pages", "value"),
        prevent_initial_call=True,
    )
    def refresh_report_catalog(_poll_count, previous_token, current_value):
        reports = discover_reports()
        token = catalog_token(reports)
        if token == previous_token:
            return no_update, no_update, no_update
        return report_tabs(reports), token, selected_value(reports, current_value)

    @app.callback(
        Output("reports-frame", "src"),
        Output("reports-frame", "className"),
        Output("reports-empty-state", "className"),
        Input("reports-pages", "value"),
    )
    def display_selected_report(value):
        report = report_by_id(_report_id_from_value(value))
        if report is None:
            return "about:blank", "reports-frame reports-frame-hidden", "reports-empty-state"
        return report.url, "reports-frame", "reports-empty-state reports-empty-state-hidden"

    @app.callback(
        Output("reports-pages", "value", allow_duplicate=True),
        Output("reports-frame", "src", allow_duplicate=True),
        Input("dashboard-request", "data"),
        prevent_initial_call=True,
    )
    def select_agent_report(request):
        if (request or {}).get("tab") != "reports" or (request or {}).get("operation") != "open_report":
            return no_update, no_update
        report_id = str((request or {}).get("report_id") or "")
        report = report_by_id(report_id)
        if report is None:
            return no_update, no_update
        return report_tab_value(report.report_id), _opened_report_url(report, (request or {}).get("id"))

    return dcc.Tab(label="Reports", value="reports", children=build_layout())
