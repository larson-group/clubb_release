"""Dash layout for the fill-holes viewer."""

from __future__ import annotations

from pathlib import Path

from dash import dcc, html

from dash_app.services.profiles import discover_output_directories
from dash_app.shared.components import styled_dropdown

from .analysis import MOMENTS, THRESHOLD_NAMES, THRESHOLDS


THRESHOLD_UNITS = {
    "wp2": "m^2/s^2",
    "rtp2": "(kg/kg)^2",
    "thlp2": "K^2",
    "up2": "m^2/s^2",
    "vp2": "m^2/s^2",
    "rtm": "kg/kg",
    "thlm": "K",
}


REPO_ROOT = Path(__file__).resolve().parents[3]


def _catalog():
    return discover_output_directories(REPO_ROOT / "output")


def _folder_options(records):
    return [
        {"label": record["label"], "value": record["path"]}
        for record in records
        if record.get("available")
    ]


def _case_options(record):
    return [
        {"label": case_name, "value": path}
        for case_name, path in (record or {}).get("cases", {}).items()
    ]


def _threshold_panel():
    rows = [
        html.Div(
            [
                html.Span(MOMENTS[moment]["label"], className="hf-threshold-label"),
                html.Code(
                    f"= {THRESHOLD_NAMES[moment]} = "
                    f"{THRESHOLDS[moment]:.3g} {THRESHOLD_UNITS[moment]}",
                    className="hf-threshold-value",
                ),
            ],
            className="hf-threshold-row",
        )
        for moment in MOMENTS
    ]
    return html.Div(
        [
            html.Div("Fill thresholds", className="hf-threshold-title"),
            html.P(
                [
                    "Hardcoded in the viewer; keep these values synchronized with "
                    "the corresponding tolerances in ",
                    html.Code("src/CLUBB_core/constants_clubb.F90"),
                    ".",
                ],
                className="hf-threshold-note",
            ),
            html.Div(rows, className="hf-threshold-list"),
        ],
        id="hf-threshold-panel",
        className="hf-threshold-panel",
    )


def build_layout():
    records = _catalog()
    default_record = next((record for record in records if record.get("case_count")), None)
    default_folder = default_record["path"] if default_record else None
    cases = _case_options(default_record)
    default_case = cases[0]["value"] if cases else None

    return html.Div(
        [
            dcc.Store(id="hf-output-catalog", data=records),
            dcc.Interval(id="hf-playback-timer", interval=500, disabled=True),
            html.Div(
                [
                    html.Main(
                        [
                            html.Div(
                                [
                                    html.Div("Numerical diagnostic", className="hf-kicker"),
                                    html.H1("Hole-filling snapshots"),
                                    html.P(
                                        "Inspect the CLUBB moments immediately before and after each "
                                        "vertical hole-fill call. Red markers identify levels changed "
                                        "by the algorithm; the timeline exposes unresolved holes."
                                    ),
                                ],
                                className="hf-hero",
                            ),
                            html.Div(
                                [
                                    dcc.Graph(
                                        id="hf-profile-figure",
                                        config={"displaylogo": False},
                                    ),
                                    dcc.Graph(
                                        id="hf-summary-figure",
                                        config={"displaylogo": False},
                                    ),
                                ],
                                className="hf-primary-grid",
                            ),
                            dcc.Graph(id="hf-map-figure", config={"displaylogo": False}),
                            html.Div(
                                "Before and after fields are direct call-site snapshots. "
                                "Thresholds use CLUBB's fixed hole-filling tolerances. "
                                "When the stats interval spans multiple model steps, each NetCDF record "
                                "contains the corresponding stats-window average.",
                                className="hf-note",
                            ),
                        ],
                        className="hf-main-pane",
                    ),
                    html.Div(className="hf-pane-divider"),
                    html.Aside(
                        [
                            html.Div(
                                [
                                    html.Div("Display controls", className="hf-kicker"),
                                    html.H2("Snapshot view"),
                                    html.P("Choose the output state and how its field profile is displayed."),
                                ],
                                className="hf-sidebar-heading",
                            ),
                            html.Div(
                                [
                                    html.Label("Output folder"),
                                    styled_dropdown(
                                        id="hf-output-folder",
                                        options=_folder_options(records),
                                        value=default_folder,
                                        clearable=False,
                                    ),
                                ],
                                className="hf-control",
                            ),
                            html.Div(
                                [
                                    html.Label("Case"),
                                    styled_dropdown(
                                        id="hf-case-file",
                                        options=cases,
                                        value=default_case,
                                        clearable=False,
                                    ),
                                ],
                                className="hf-control",
                            ),
                            html.Button(
                                "Refresh outputs",
                                id="hf-refresh-outputs",
                                n_clicks=0,
                                className="hf-sidebar-button",
                            ),
                            html.Div(id="hf-file-status", className="hf-status"),
                            html.Div(className="hf-sidebar-rule"),
                            html.Div(
                                [
                                    html.Label("CLUBB moment"),
                                    styled_dropdown(
                                        id="hf-moment",
                                        options=[],
                                        value=None,
                                        clearable=False,
                                    ),
                                ],
                                className="hf-control",
                            ),
                            html.Div(
                                [
                                    html.Label("Grid column"),
                                    styled_dropdown(
                                        id="hf-column",
                                        options=[],
                                        value=None,
                                        clearable=False,
                                    ),
                                ],
                                className="hf-control",
                            ),
                            html.Div(
                                [
                                    html.Label("Field magnitude scale"),
                                    dcc.RadioItems(
                                        id="hf-field-scale",
                                        options=[
                                            {"label": "Log", "value": "log"},
                                            {"label": "Linear", "value": "linear"},
                                        ],
                                        value="log",
                                        inline=True,
                                        className="hf-radio-row",
                                    ),
                                ],
                                className="hf-control",
                            ),
                            _threshold_panel(),
                            html.Div(
                                [
                                    html.Div(
                                        [
                                            html.Label("Playback speed"),
                                            dcc.RadioItems(
                                                id="hf-speed",
                                                options=[
                                                    {"label": "0.5×", "value": 1000},
                                                    {"label": "1×", "value": 500},
                                                    {"label": "2×", "value": 250},
                                                ],
                                                value=500,
                                                inline=True,
                                                className="hf-radio-row",
                                            ),
                                        ],
                                        className="hf-control hf-speed-control",
                                    ),
                                    html.Button(
                                        "Play",
                                        id="hf-play",
                                        n_clicks=0,
                                        className="hf-sidebar-button hf-play-button",
                                    ),
                                ],
                                className="hf-playback-row",
                            ),
                            html.Div(className="hf-sidebar-rule"),
                            html.Div(
                                [
                                    html.Label("NetCDF record / timestep"),
                                    dcc.Slider(
                                        id="hf-record",
                                        min=0,
                                        max=1,
                                        step=1,
                                        value=0,
                                        marks={0: "0", 1: "1"},
                                        updatemode="drag",
                                        tooltip={
                                            "placement": "bottom",
                                            "always_visible": False,
                                        },
                                    ),
                                ],
                                className="hf-control hf-record-control",
                            ),
                            html.Div(id="hf-time-label", className="hf-time-label"),
                        ],
                        className="hf-right-pane",
                    ),
                ],
                className="hf-workspace",
            ),
        ],
        className="hf-viewer",
    )
