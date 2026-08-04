"""Dash layout for the mixing-length parcel trajectory explorer."""

from __future__ import annotations

from pathlib import Path

from dash import dcc, html

from dash_app.shared.components import styled_dropdown

from .analysis import discover_netcdf_file_records

REPO_ROOT = Path(__file__).resolve().parents[3]


def _file_options(records):
    return [
        {
            "label": str(record["label"]),
            "value": str(record["path"]),
        }
        for record in records
    ]


def build_layout():
    records = discover_netcdf_file_records(REPO_ROOT)
    paths = [Path(record["path"]) for record in records]
    default_path = (
        next(
            (
                path
                for path in paths
                if path == REPO_ROOT / "output" / "dash_default" / "arm_stats.nc"
            ),
            paths[0] if paths else None,
        )
    )
    return html.Div(
        [
            html.Div(
                [
                    html.Main(
                        [
                            html.Div(
                                [
                                    html.Div(
                                        "Interactive diagnostic",
                                        className="mlt-kicker",
                                    ),
                                    html.H1(
                                        "Mixing-length parcel trajectories"
                                    ),
                                    html.P(
                                        "At one CLUBB output record, launch a "
                                        "parcel from every thermodynamic level "
                                        "and follow its remaining energy K upward "
                                        "and downward. The profile panel shows "
                                        "how those paths become L↑, L↓, and the "
                                        "final geometric-mean Lscale."
                                    ),
                                ],
                                className="mlt-hero",
                            ),
                            html.Div(
                                id="mlt-metrics",
                                className="mlt-metric-grid",
                            ),
                            html.Div(
                                [
                                    dcc.Graph(
                                        id="mlt-upward-figure",
                                        className="mlt-graph",
                                        config={"displaylogo": False},
                                    ),
                                    dcc.Graph(
                                        id="mlt-downward-figure",
                                        className="mlt-graph",
                                        config={"displaylogo": False},
                                    ),
                                ],
                                className="mlt-spaghetti-grid",
                            ),
                            html.Div(
                                [
                                    html.Div(
                                        "Parcel buoyancy",
                                        className="mlt-kicker",
                                    ),
                                    html.H2(
                                        "Buoyancy along each diagnosed path"
                                    ),
                                    html.P(
                                        "The same upward and downward parcels "
                                        "are colored by launch height. Positive "
                                        "buoyancy adds kinetic energy; negative "
                                        "buoyancy removes it."
                                    ),
                                ],
                                className="mlt-section-heading",
                            ),
                            html.Div(
                                [
                                    dcc.Graph(
                                        id="mlt-upward-buoyancy-figure",
                                        className="mlt-graph",
                                        config={"displaylogo": False},
                                    ),
                                    dcc.Graph(
                                        id="mlt-downward-buoyancy-figure",
                                        className="mlt-graph",
                                        config={"displaylogo": False},
                                    ),
                                ],
                                className="mlt-spaghetti-grid",
                            ),
                            html.Div(
                                [
                                    html.Div(
                                        "Entraining parcel state",
                                        className="mlt-kicker",
                                    ),
                                    html.H2(
                                        "Upward parcel properties versus altitude"
                                    ),
                                    html.P(
                                        "Parcel rₜ and θₗ follow the same "
                                        "layerwise entrainment used by the "
                                        "buoyancy calculation. Parcel θᵥ is "
                                        "shown against grid-mean thvm and is "
                                        "exactly consistent with the diagnosed "
                                        "parcel buoyancy."
                                    ),
                                ],
                                className="mlt-section-heading",
                            ),
                            html.Div(
                                [
                                    dcc.Graph(
                                        id="mlt-parcel-thv-figure",
                                        className="mlt-graph",
                                        config={"displaylogo": False},
                                    ),
                                    dcc.Graph(
                                        id="mlt-parcel-rt-figure",
                                        className="mlt-graph",
                                        config={"displaylogo": False},
                                    ),
                                    dcc.Graph(
                                        id="mlt-parcel-thl-figure",
                                        className="mlt-graph",
                                        config={"displaylogo": False},
                                    ),
                                ],
                                className="mlt-parcel-state-grid",
                            ),
                            dcc.Graph(
                                id="mlt-profile-figure",
                                className="mlt-graph",
                                config={"displaylogo": False},
                            ),
                            html.Div(
                                [
                                    html.Strong("Replica scope. "),
                                    "This follows the ascending-grid standalone "
                                    "SCM path with Flatau saturation "
                                    "(saturation_formula=3), including the "
                                    "0.1 m internal floor, nonlocal endpoint "
                                    "envelopes, the surface lmin floor, and the "
                                    "final top-boundary copy. The NetCDF must "
                                    "contain the pre-advance thlm_old/rtm_old "
                                    "fields.",
                                ],
                                className="mlt-note",
                            ),
                        ],
                        className="mlt-main-pane",
                    ),
                    html.Div(className="mlt-pane-divider"),
                    html.Aside(
                        [
                            html.Div(
                                [
                                    html.Div(
                                        "Display controls",
                                        className="mlt-kicker",
                                    ),
                                    html.H2("Trajectory view"),
                                    html.P(
                                        "Choose the source state used for all "
                                        "three plots."
                                    ),
                                ],
                                className="mlt-sidebar-heading",
                            ),
                            html.Div(
                                [
                                    html.Label("NetCDF file"),
                                    styled_dropdown(
                                        id="mlt-file",
                                        options=_file_options(records),
                                        value=(
                                            str(default_path)
                                            if default_path
                                            else None
                                        ),
                                        clearable=False,
                                    ),
                                ],
                                className="mlt-control",
                            ),
                            html.Div(
                                [
                                    html.Label("Add file by path"),
                                    dcc.Input(
                                        id="mlt-custom-path",
                                        type="text",
                                        placeholder="/path/to/case_stats.nc",
                                    ),
                                    html.Div(
                                        [
                                            html.Button(
                                                "Load",
                                                id="mlt-load-path",
                                                n_clicks=0,
                                            ),
                                            html.Button(
                                                "Refresh",
                                                id="mlt-refresh-files",
                                                n_clicks=0,
                                            ),
                                        ],
                                        className="mlt-path-buttons",
                                    ),
                                ],
                                className="mlt-control",
                            ),
                            html.Div(
                                [
                                    html.Label("Grid column"),
                                    styled_dropdown(
                                        id="mlt-column",
                                        options=[
                                            {"label": "0", "value": 0}
                                        ],
                                        value=0,
                                        clearable=False,
                                    ),
                                ],
                                className="mlt-control",
                            ),
                            html.Div(
                                [
                                    html.Label(
                                        [
                                            "Parcel entrainment μ",
                                            html.Span(
                                                id="mlt-mu-label",
                                                className="mlt-slider-value",
                                            ),
                                        ]
                                    ),
                                    dcc.Slider(
                                        id="mlt-mu",
                                        min=1.0e-5,
                                        max=3.0e-3,
                                        step=1.0e-5,
                                        value=1.0e-3,
                                        marks={
                                            1.0e-5: "0.00001",
                                            1.0e-3: "0.001",
                                            2.0e-3: "0.002",
                                            3.0e-3: "0.003",
                                        },
                                        updatemode="mouseup",
                                        tooltip={
                                            "placement": "bottom",
                                            "always_visible": False,
                                        },
                                    ),
                                    html.Div(
                                        "Initialized from the selected NetCDF "
                                        "column. The orange Fortran reference "
                                        "retains the file's original μ.",
                                        className="mlt-control-help",
                                    ),
                                ],
                                className="mlt-control mlt-mu-control",
                            ),
                            html.Div(className="mlt-sidebar-rule"),
                            html.Div(
                                [
                                    html.Label(
                                        "NetCDF record / timestep"
                                    ),
                                    dcc.Slider(
                                        id="mlt-record",
                                        min=0,
                                        max=1,
                                        step=1,
                                        value=0,
                                        marks={0: "0", 1: "1"},
                                        updatemode="mouseup",
                                        tooltip={
                                            "placement": "bottom",
                                            "always_visible": False,
                                        },
                                    ),
                                ],
                                className="mlt-control mlt-time-slider",
                            ),
                            html.Div(
                                id="mlt-time-label",
                                className="mlt-time-label",
                            ),
                            html.Div(
                                id="mlt-file-status",
                                className="mlt-status",
                            ),
                        ],
                        className="mlt-right-pane",
                    ),
                ],
                className="mlt-workspace",
            ),
        ],
        className="mlt-page",
    )
