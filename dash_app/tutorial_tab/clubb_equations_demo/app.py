"""Interactive equation-guide layout and callbacks."""

from __future__ import annotations

import argparse
from pathlib import Path

from dash import ALL, Dash, Input, Output, ctx, html

from .components import (
    DEFAULT_TERM,
    _contains_fraction,
    dependency_chain,
    equation_tabs,
    flow_map,
    inventory_strip,
    ownership_legend,
    reference_drawer,
    term_inspector,
)
from .content import TERM_OCCURRENCES


def build_layout():
    return html.Div(
        [
            html.Header(
                [
                    html.Div(
                        [
                            html.Span("Interactive quick reference", className="eq-eyebrow"),
                            html.H1("How CLUBB's equations fit together"),
                            html.P(
                                "Follow a moment from its prognostic budget, through the closure it needs, to the cloud and buoyancy diagnostics that feed back on turbulence. Click any colored term for a plain-language explanation."
                            ),
                        ]
                    ),
                    html.A(
                        "Open the full equation document ↗",
                        href="https://arxiv.org/pdf/1711.03675",
                        target="_blank",
                        rel="noreferrer",
                        className="eq-source-link",
                    ),
                ],
                className="eq-page-header",
            ),
            html.Section(
                [
                    html.Div(
                        [
                            html.H2("The CLUBB loop in one glance"),
                            html.P("Moments constrain a PDF; the PDF supplies missing moments and cloud diagnostics; all closures return to the budgets."),
                        ],
                        className="eq-section-heading",
                    ),
                    flow_map(),
                    inventory_strip(),
                    dependency_chain(),
                ],
                className="eq-overview-panel",
            ),
            html.Section(
                [
                    html.Div(
                        [
                            html.H2("Equation explorer"),
                            html.P("Colors follow the official CLUBB equation-document ownership convention. Physical-process labels appear in the inspector."),
                        ],
                        className="eq-section-heading",
                    ),
                    ownership_legend(),
                    html.Div(
                        [
                            html.Main(equation_tabs(), className="eq-explorer-main"),
                            html.Aside(
                                term_inspector(),
                                id="equation-term-inspector",
                                className="eq-term-inspector",
                                **{"aria-live": "polite"},
                            ),
                        ],
                        className="eq-explorer-grid",
                    ),
                ],
                className="eq-explorer-panel",
            ),
            html.Section(
                [
                    html.Div(
                        [
                            html.H2("Keep beside you while reading the code"),
                            html.P("Notation, budget suffixes, and the caveat that explains most apparent budget mismatches."),
                        ],
                        className="eq-section-heading",
                    ),
                    reference_drawer(),
                ],
                className="eq-reference-panel",
            ),
            html.Footer(
                [
                    html.Span("Primary sources: "),
                    html.Code("doc/CLUBBeqns.tex"),
                    html.Span(" and current routines in "),
                    html.Code("src/CLUBB_core/"),
                    html.Span(". The arXiv document describes a particular revision; active flags and newer source paths can change details."),
                ],
                className="eq-page-footer",
            ),
        ],
        className="eq-guide-page",
    )


def _button_class(occurrence, selected):
    metadata = TERM_OCCURRENCES[occurrence]
    classes = [
        "eq-term",
        f"eq-owner-{metadata['ownership']}",
        f"eq-process-{metadata['process']}",
    ]
    if occurrence == selected:
        classes.append("eq-term-selected")
    if _contains_fraction(metadata["latex"]):
        classes.append("eq-term-fraction")
    return " ".join(classes)


def register_callbacks(app):
    @app.callback(
        Output("equation-term-inspector", "children"),
        Output({"type": "equation-guide-term", "index": ALL}, "className"),
        Input({"type": "equation-guide-term", "index": ALL}, "n_clicks"),
    )
    def _explain_term(_clicks):
        input_group = ctx.inputs_list[0] if ctx.inputs_list else []
        visible_occurrences = [
            entry["id"]["index"]
            for entry in input_group
            if isinstance(entry.get("id"), dict)
            and entry["id"].get("index") in TERM_OCCURRENCES
        ]
        triggered = ctx.triggered_id
        selected = (
            DEFAULT_TERM
            if DEFAULT_TERM in visible_occurrences
            else visible_occurrences[0] if visible_occurrences else DEFAULT_TERM
        )
        if isinstance(triggered, dict):
            candidate = triggered.get("index")
            if candidate in TERM_OCCURRENCES:
                selected = candidate
        classes = [
            _button_class(occurrence, selected) for occurrence in visible_occurrences
        ]
        return term_inspector(TERM_OCCURRENCES[selected]), classes


def create_app():
    app = Dash(
        __name__,
        suppress_callback_exceptions=True,
        title="CLUBB Equation Guide",
        assets_folder=str(Path(__file__).resolve().parents[2] / "assets"),
    )
    app.layout = build_layout()
    register_callbacks(app)
    return app


def main():
    parser = argparse.ArgumentParser(description="Launch the standalone CLUBB equation guide.")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=23406)
    parser.add_argument("-debug", action="store_true")
    args = parser.parse_args()
    create_app().run(host=args.host, port=args.port, debug=args.debug)


if __name__ == "__main__":
    main()
