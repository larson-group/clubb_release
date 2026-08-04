"""Small, reusable building blocks for dashboard reports and notes."""

from __future__ import annotations

from collections.abc import Iterable, Sequence

from dash import html


def report_header(
    title: str,
    summary: str,
    *,
    eyebrow: str | None = None,
    badges: Iterable[str] = (),
):
    """Build a compact report heading with optional provenance/status badges."""

    badge_nodes = [
        html.Span(str(badge), className="shared-report-badge")
        for badge in badges
        if str(badge).strip()
    ]
    return html.Header(
        [
            html.Div(eyebrow, className="shared-report-eyebrow") if eyebrow else None,
            html.H1(title),
            html.P(summary, className="shared-report-lead"),
            html.Div(badge_nodes, className="shared-report-badges")
            if badge_nodes
            else None,
        ],
        className="shared-report-header",
    )


def report_section(
    title: str,
    children,
    *,
    intro: str | None = None,
    class_name: str = "",
    collapsible: bool = False,
    open_by_default: bool = False,
):
    """Wrap one report topic in a consistently styled, reusable section."""

    if not isinstance(children, (list, tuple)):
        children = [children]
    classes = "shared-report-section"
    if class_name:
        classes = f"{classes} {class_name}"
    if collapsible:
        return html.Details(
            [
                html.Summary(
                    [
                        html.Span(title),
                        html.Span(
                            "Expand / collapse",
                            className="shared-report-section-toggle-label",
                        ),
                    ],
                    className="shared-report-section-summary",
                ),
                html.Div(
                    [
                        html.P(intro, className="shared-report-section-intro")
                        if intro
                        else None,
                        *children,
                    ],
                    className="shared-report-section-body",
                ),
            ],
            open=open_by_default,
            className=f"{classes} shared-report-section-collapsible",
        )
    return html.Section(
        [
            html.H2(title),
            html.P(intro, className="shared-report-section-intro") if intro else None,
            *children,
        ],
        className=classes,
    )


def fact_card(
    label: str,
    value: str,
    detail: str = "",
    *,
    tone: str = "neutral",
):
    """Build a concise labeled fact for a report summary grid."""

    tone_name = tone if tone in {"neutral", "good", "warning", "info"} else "neutral"
    return html.Div(
        [
            html.Div(label, className="shared-report-fact-label"),
            html.Div(value, className="shared-report-fact-value"),
            html.Div(detail, className="shared-report-fact-detail") if detail else None,
        ],
        className=f"shared-report-fact shared-report-fact-{tone_name}",
    )


def fact_grid(items: Sequence[dict[str, str]]):
    """Render dictionaries accepted by :func:`fact_card` as a responsive grid."""

    return html.Div(
        [
            fact_card(
                item["label"],
                item["value"],
                item.get("detail", ""),
                tone=item.get("tone", "neutral"),
            )
            for item in items
        ],
        className="shared-report-fact-grid",
    )


def empty_state(title: str, message: str):
    """Show an intentional placeholder where a report will later attach data."""

    return html.Div(
        [
            html.Div("○", className="shared-report-empty-icon", **{"aria-hidden": "true"}),
            html.H3(title),
            html.P(message),
        ],
        className="shared-report-empty",
    )
