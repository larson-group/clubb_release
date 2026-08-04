"""Reusable overlay notecards for explanatory text, diagnostics, and logs."""

from __future__ import annotations

from dash import html


NOTECARD_SIZES = frozenset({"small", "medium", "large", "full"})
OVERLAY_STYLE = {
    "position": "fixed",
    "inset": 0,
    "zIndex": 2000,
    "display": "flex",
    "alignItems": "center",
    "justifyContent": "center",
    "padding": "4vh 4vw",
    "boxSizing": "border-box",
    "background": "rgba(2, 6, 23, 0.58)",
}
PANEL_STYLE = {
    "display": "flex",
    "flexDirection": "column",
    "maxWidth": "calc(100vw - 32px)",
    "maxHeight": "calc(100vh - 32px)",
    "border": "1px solid var(--shared-notecard-border, #475569)",
    "borderRadius": "10px",
    "overflow": "hidden",
    "boxShadow": "0 24px 80px rgba(0, 0, 0, 0.5)",
    "background": "var(--shared-notecard-background, #0f172a)",
    "color": "var(--shared-notecard-text, #e2e8f0)",
}
SIZE_STYLES = {
    "small": {"width": "min(520px, 92vw)", "maxHeight": "72vh"},
    "medium": {"width": "min(760px, 92vw)", "maxHeight": "82vh"},
    "large": {"width": "min(1080px, 94vw)", "maxHeight": "88vh"},
    "full": {"width": "90vw", "height": "90vh"},
}
HEADER_STYLE = {
    "display": "flex",
    "alignItems": "center",
    "justifyContent": "space-between",
    "gap": "16px",
    "flex": "0 0 auto",
    "borderBottom": "1px solid var(--shared-notecard-border, #475569)",
    "padding": "14px 16px",
}
BODY_STYLE = {
    "flex": "1 1 auto",
    "minHeight": 0,
    "overflow": "auto",
    "padding": "18px 20px 22px",
    "fontSize": "14px",
    "lineHeight": 1.5,
}


def notecard(
    title,
    body,
    close_button_id,
    *,
    subtitle=None,
    size="medium",
    body_class_name="",
    panel_class_name="",
):
    """Return an accessible themed overlay notecard with configurable content."""
    normalized_size = str(size or "medium").lower()
    if normalized_size not in NOTECARD_SIZES:
        raise ValueError(f"Unsupported notecard size: {size}")
    body_classes = " ".join(
        value for value in ("shared-notecard-body", body_class_name) if value
    )
    panel_classes = " ".join(
        value
        for value in (
            "shared-notecard-panel",
            f"shared-notecard-size-{normalized_size}",
            panel_class_name,
        )
        if value
    )
    return html.Div(
        html.Div(
            [
                html.Div(
                    [
                        html.Div(
                            [
                                html.Div(
                                    title,
                                    className="shared-notecard-title",
                                    style={
                                        "fontSize": "19px",
                                        "fontWeight": 750,
                                        "lineHeight": 1.25,
                                    },
                                ),
                                html.Div(
                                    subtitle,
                                    className="shared-notecard-subtitle",
                                    style={
                                        "marginTop": "4px",
                                        "overflow": "hidden",
                                        "textOverflow": "ellipsis",
                                        "whiteSpace": "nowrap",
                                        "fontFamily": "ui-monospace, SFMono-Regular, Menlo, Consolas, monospace",
                                        "fontSize": "12px",
                                    },
                                )
                                if subtitle
                                else None,
                            ],
                            className="shared-notecard-header-text",
                            style={"minWidth": 0},
                        ),
                        html.Button(
                            "Close",
                            id=close_button_id,
                            type="button",
                            n_clicks=0,
                            className="shared-notecard-close",
                            title="Close notecard",
                            style={
                                "flex": "0 0 auto",
                                "border": "1px solid var(--shared-notecard-border, #64748b)",
                                "borderRadius": "6px",
                                "cursor": "pointer",
                                "fontWeight": 700,
                                "padding": "8px 12px",
                                "background": "var(--shared-notecard-header, #1e293b)",
                                "color": "var(--shared-notecard-text, #e2e8f0)",
                            },
                        ),
                    ],
                    className="shared-notecard-header",
                    style=HEADER_STYLE,
                ),
                html.Div(body, className=body_classes, style=BODY_STYLE),
            ],
            className=panel_classes,
            role="dialog",
            style={**PANEL_STYLE, **SIZE_STYLES[normalized_size]},
            **{"aria-modal": "true", "aria-label": str(title)},
        ),
        className="shared-notecard-overlay",
        style=OVERLAY_STYLE,
    )


def information_body(overview, sections=()):
    """Build consistently structured overview-and-nuance notecard content."""
    children = [
        html.Section(
            [
                html.H3("Overview", className="shared-notecard-section-title"),
                html.P(overview, className="shared-notecard-overview"),
            ],
            className="shared-notecard-section shared-notecard-overview-section",
        )
    ]
    for section in sections or ():
        heading = section.get("heading")
        paragraphs = section.get("paragraphs") or ()
        bullets = section.get("bullets") or ()
        section_children = []
        if heading:
            section_children.append(
                html.H3(heading, className="shared-notecard-section-title")
            )
        section_children.extend(html.P(paragraph) for paragraph in paragraphs)
        if bullets:
            section_children.append(html.Ul([html.Li(item) for item in bullets]))
        children.append(
            html.Section(section_children, className="shared-notecard-section")
        )
    return html.Div(children, className="shared-information-content")
