"""Welcome page for the dashboard tutorial section."""

from dash import html


def build_page():
    """Return the short introduction shown when the Tutorial tab opens."""
    return html.Div(
        [
            html.Div("CLUBB tutorials", className="tutorial-eyebrow"),
            html.H1("Welcome"),
            html.P(
                "These tutorials are small, interactive explanations of CLUBB. "
                "They are meant for people who recognize the model variables but "
                "want to understand how its closures turn those variables into behavior.",
                className="tutorial-lead",
            ),
            html.Div(
                [
                    html.Div(
                        [
                            html.H2("A good way to begin"),
                            html.Ol(
                                [
                                    html.Li("Open CLUBB Equations and learn what the model advances versus diagnoses."),
                                    html.Li("Click the turbulent-advection and buoyancy terms in the w′² budget."),
                                    html.Li("Follow the PDF → w′r_c′ → w′θᵥ′ → w′² dependency chain."),
                                    html.Li("Then use the ADG1 two-Gaussian explorer to see the PDF assumptions geometrically."),
                                ],
                                className="tutorial-steps",
                            ),
                            html.Button(
                                "Start the equations guide",
                                id="tutorial-start-equations",
                                n_clicks=0,
                                className="tutorial-start-button",
                            ),
                            html.Button(
                                "Open the ADG1 explorer",
                                id="tutorial-start-adg1-explorer",
                                n_clicks=0,
                                className="tutorial-start-button tutorial-start-button-secondary",
                            ),
                        ],
                        className="tutorial-welcome-card tutorial-welcome-primary",
                    ),
                    html.Div(
                        [
                            html.H2("What to look for"),
                            html.P(
                                "The important comparison is not which contour looks nicest. "
                                "Watch which component properties are allowed to change, which "
                                "moments remain fixed, and when a physical constraint forces an adjustment."
                            ),
                            html.P(
                                "Each lesson begins with the simplest useful controls. More detailed "
                                "information is kept in expandable sections."
                            ),
                        ],
                        className="tutorial-welcome-card",
                    ),
                ],
                className="tutorial-welcome-grid",
            ),
        ],
        className="tutorial-welcome-page",
    )
