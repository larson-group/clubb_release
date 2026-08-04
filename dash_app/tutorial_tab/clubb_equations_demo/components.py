"""Dash component builders for the CLUBB equation guide."""

from __future__ import annotations

import re

from dash import dcc, html

from .content import (
    ABBREVIATIONS,
    EQUATION_GROUPS,
    NOTATION,
    OWNERSHIP_META,
    PROCESS_META,
    TERM_OCCURRENCES,
)


DEFAULT_TERM = "wp2--3--bp"


def _contains_fraction(latex):
    return any(marker in latex for marker in (r"\frac", r"\dfrac", r"\tfrac"))


def _term_min_width(latex):
    """Estimate a readable card width for one unbroken MathJax expression.

    MathJax content does not wrap like prose.  Raw LaTeX length is therefore a
    useful conservative proxy: commands make fractions and nested derivatives
    visually wider, which is exactly where a larger card is needed.
    """
    visual = re.sub(r"\\[A-Za-z]+", "x", latex)
    visual = re.sub(r"[{}_^]", "", visual)
    visual = re.sub(r"\\[,;! ]|\s+", "", visual)
    fraction_count = sum(latex.count(marker) for marker in (r"\frac", r"\dfrac", r"\tfrac"))
    width = 7.0 + 0.42 * len(visual) + 0.7 * fraction_count
    return float(min(max(width, 10.0), 30.0))


def math_fragment(latex, class_name=""):
    """Render one compact MathJax fragment without paragraph margins."""
    has_fraction = _contains_fraction(latex)
    rendered_latex = rf"\displaystyle {latex}" if has_fraction else latex
    return dcc.Markdown(
        f"${rendered_latex}$",
        mathjax=True,
        className=" ".join(
            value
            for value in (
                "eq-math",
                "eq-math-fraction" if has_fraction else "",
                class_name,
            )
            if value
        ),
    )


def ownership_legend():
    return html.Div(
        [
            html.Div(
                [
                    html.Span(className=f"eq-legend-dot eq-owner-{key}"),
                    html.Span(label),
                ],
                className="eq-legend-item",
                title=description,
            )
            for key, (label, description) in OWNERSHIP_META.items()
        ],
        className="eq-ownership-legend",
        **{"aria-label": "Official CLUBB equation color legend"},
    )


def flow_map():
    stages = (
        ("1", "Host state", "Mean winds, rₜ, θₗ and forcing", "host"),
        ("2", "Prognosed moments", "Variances, covariances, fluxes and w′³", "clubb_internal"),
        ("3", "PDF diagnosis", "Component weights, means, widths and correlations", "pdf"),
        ("4", "Closed terms", "Higher moments, cloud properties and buoyancy", "pdf"),
        ("5", "Advance", "PDF and classical closures update CLUBB moments", "classical"),
    )
    children = []
    for index, (number, title, description, owner) in enumerate(stages):
        if index:
            children.append(html.Div("→", className="eq-flow-arrow", **{"aria-hidden": "true"}))
        children.append(
            html.Div(
                [
                    html.Span(number, className=f"eq-flow-number eq-owner-{owner}"),
                    html.Div([html.Strong(title), html.Small(description)]),
                ],
                className="eq-flow-stage",
            )
        )
    return html.Div(children, className="eq-flow-map")


def inventory_strip():
    columns = (
        (
            "Host/grid means",
            ("r̄ₜ", "θ̄ₗ", "ū", "v̄"),
            "The resolved environment CLUBB mixes.",
            "host",
        ),
        (
            "Core CLUBB moments",
            ("w′²", "w′³", "rₜ′²", "θₗ′²", "rₜ′θₗ′", "w′rₜ′", "w′θₗ′"),
            "Advanced in time; these constrain the PDF.",
            "clubb_internal",
        ),
        (
            "PDF diagnostics",
            ("w′⁴", "w′²rₜ′", "w′r_c′", "cloud fraction", "r̄_c", "w′θᵥ′"),
            "Needed by budgets but not independently prognosed.",
            "pdf",
        ),
    )
    return html.Div(
        [
            html.Div(
                [
                    html.Span(title, className=f"eq-inventory-label eq-owner-{owner}"),
                    html.Div(
                        [html.Strong(symbol, className="eq-inventory-symbol") for symbol in symbols],
                        className="eq-inventory-symbols",
                    ),
                    html.Small(description),
                ],
                className="eq-inventory-card",
            )
            for title, symbols, description, owner in columns
        ],
        className="eq-inventory-strip",
    )


def dependency_chain():
    """Show one concrete closure-to-tendency path relevant to this checkout."""
    nodes = (
        ("PDF geometry", "places cloud relative to w", "pdf"),
        ("w′r_c′", "cloud-water transport", "pdf"),
        ("w′θᵥ′", "buoyancy flux", "pdf"),
        ("wp2_bp", "buoyancy production", "pdf"),
        ("∂w′²/∂t", "vertical-turbulence tendency", "clubb_internal"),
    )
    children = [html.Strong("Follow one dependency:", className="eq-chain-title")]
    for index, (label, description, owner) in enumerate(nodes):
        if index:
            children.append(html.Span("→", className="eq-chain-arrow"))
        children.append(
            html.Div(
                [html.Span(label, className=f"eq-owner-{owner}"), html.Small(description)],
                className="eq-chain-node",
            )
        )
    return html.Div(children, className="eq-dependency-chain")


def term_button(equation_id, position, item, *, selected=False):
    occurrence = f"{equation_id}--{position}--{item['key']}"
    classes = [
        "eq-term",
        f"eq-owner-{item['ownership']}",
        f"eq-process-{item['process']}",
    ]
    if selected:
        classes.append("eq-term-selected")
    if _contains_fraction(item["latex"]):
        classes.append("eq-term-fraction")
    return html.Button(
        [
            math_fragment(item["latex"]),
            html.Span(item["code"], className="eq-budget-code"),
        ],
        id={"type": "equation-guide-term", "index": occurrence},
        n_clicks=0,
        type="button",
        className=" ".join(classes),
        style={"--eq-term-min": f"{_term_min_width(item['latex']):.1f}rem"},
        title=f"{item['title']} — click for explanation",
        **{"aria-label": f"Explain {item['title']} in {TERM_OCCURRENCES[occurrence]['equation_title']}"},
    )


def equation_card(equation):
    lhs_owner = equation["lhs_owner"]
    return html.Article(
        [
            html.Header(
                [
                    html.Div(
                        [
                            html.H3(equation["title"]),
                            html.P(equation["summary"]),
                        ]
                    ),
                    html.Code(equation["source"], className="eq-source-chip"),
                ],
                className="eq-card-header",
            ),
            html.Div(
                [
                    html.Div(
                        math_fragment(equation["lhs"]),
                        className=f"eq-lhs eq-owner-{lhs_owner}",
                        title=OWNERSHIP_META[lhs_owner][1],
                    ),
                    html.Span("=", className="eq-equals", **{"aria-hidden": "true"}),
                    html.Div(
                        [
                            term_button(
                                equation["id"],
                                position,
                                item,
                                selected=(
                                    f"{equation['id']}--{position}--{item['key']}"
                                    == DEFAULT_TERM
                                ),
                            )
                            for position, item in enumerate(equation["terms"])
                        ],
                        className="eq-rhs",
                    ),
                ],
                className="eq-expression",
            ),
            html.P(equation["note"], className="eq-card-note"),
        ],
        className="eq-card",
        id=f"equation-{equation['id']}",
    )


def equation_group(group):
    return html.Div(
        [
            html.Div(
                [
                    html.Span(group["eyebrow"], className="eq-section-eyebrow"),
                    html.P(group["intro"]),
                ],
                className="eq-group-intro",
            ),
            html.Div(
                [equation_card(equation) for equation in group["equations"]],
                className="eq-card-stack",
            ),
        ],
        className="eq-group",
    )


def equation_tabs():
    return dcc.Tabs(
        id="equation-guide-groups",
        value=EQUATION_GROUPS[0]["value"],
        persistence=True,
        persistence_type="local",
        className="eq-family-tabs",
        parent_className="eq-family-tabs-parent",
        content_className="eq-family-content",
        children=[
            dcc.Tab(
                label=group["label"],
                value=group["value"],
                className="eq-family-tab",
                selected_className="eq-family-tab-selected",
                children=equation_group(group),
            )
            for group in EQUATION_GROUPS
        ],
    )


def term_inspector(metadata=None):
    metadata = metadata or TERM_OCCURRENCES[DEFAULT_TERM]
    process_label, process_description = PROCESS_META[metadata["process"]]
    owner_label, owner_description = OWNERSHIP_META[metadata["ownership"]]
    return html.Div(
        [
            html.Div(
                [
                    html.Span(process_label, className=f"eq-inspector-badge eq-process-{metadata['process']}"),
                    html.Span(owner_label, className=f"eq-inspector-badge eq-owner-{metadata['ownership']}"),
                ],
                className="eq-inspector-badges",
            ),
            html.H2(metadata["title"]),
            html.Div(math_fragment(metadata["latex"]), className="eq-inspector-formula"),
            html.Div(
                [html.Span("Budget label"), html.Code(metadata["code"])],
                className="eq-inspector-code",
            ),
            html.Section([html.H3("What it means"), html.P(metadata["meaning"])]),
            html.Section([html.H3("How CLUBB gets it"), html.P(metadata["obtained"])]),
            html.Section([html.H3("Why you care"), html.P(metadata["why"])]),
            html.Details(
                [
                    html.Summary("Implementation details"),
                    html.Dl(
                        [
                            html.Dt("Appears in"),
                            html.Dd(metadata["equation_title"]),
                            html.Dt("Units"),
                            html.Dd(metadata["units"]),
                            html.Dt("Equation source"),
                            html.Dd(metadata["source"]),
                            html.Dt("Equation card"),
                            html.Dd(metadata["equation_source"]),
                            html.Dt("Color meaning"),
                            html.Dd(owner_description),
                            html.Dt("Physical family"),
                            html.Dd(process_description),
                        ]
                    ),
                ],
                className="eq-inspector-details",
            ),
        ],
        className="eq-inspector-inner",
    )


def reference_drawer():
    return html.Div(
        [
            html.Details(
                [
                    html.Summary("Notation"),
                    html.Div(
                        [
                            html.Div([math_fragment(symbol), html.Span(description)])
                            for symbol, description in NOTATION
                        ],
                        className="eq-reference-list",
                    ),
                ],
                open=True,
            ),
            html.Details(
                [
                    html.Summary("Budget abbreviations"),
                    html.Div(
                        [html.Div([html.Code(code), html.Span(description)]) for code, description in ABBREVIATIONS],
                        className="eq-reference-list eq-abbreviation-list",
                    ),
                ]
            ),
            html.Details(
                [
                    html.Summary("Continuous equations versus the code"),
                    html.P(
                        "These cards show the continuous closed equations. The actual advance mixes explicit and implicit pieces and may add surface forcing, clipping, positive-definite hole filling, flux limiting, and host-model tendencies according to the active flags. A saved before/after difference therefore need not equal only the terms inside one advance routine."
                    ),
                ]
            ),
        ],
        className="eq-reference-drawer",
    )
