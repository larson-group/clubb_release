"""Build controls and plot containers for coherent mixing-length diagnostics.

Usage:
  ``build_layout()`` returns the Affine Coherence Mixing Length page for Dash's
  Misc registry. ``callbacks.register_callbacks(app)`` supplies its behavior.

Understanding:
  File choices come from discoverable NetCDF output, preferring
  output/dash_default/arm_stats.nc when present and otherwise the newest file.
  The page selects one zero-based record/column and sets entrainment [1/m],
  altitude ceiling [m], energy multiplier, PDF launch fraction, directional
  weight, and aggregation. It displays five figures without launching a model
  run. IDs use the amlt- prefix; plots are filled after selection callbacks.
"""

from __future__ import annotations

from pathlib import Path

from dash import dcc, html

from dash_app.shared.components import styled_dropdown

from .analysis import discover_netcdf_file_records

REPO_ROOT = Path(__file__).resolve().parents[3]


def _directory_label(directory):
    """Return a repository-relative directory label, or its full path when outside the checkout."""
    path = Path(directory)
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def _directory_options(records):
    """Return distinct output-directory choices sorted by newest file, then directory label."""
    modified = {}
    for record in records:
        directory = str(record["directory"])
        modified[directory] = max(
            modified.get(directory, 0.0), float(record.get("modified", 0.0))
        )
    directories = sorted(
        modified,
        key=lambda directory: (-modified[directory], _directory_label(directory)),
    )
    return [
        {"label": _directory_label(directory), "value": directory}
        for directory in directories
    ]


def _case_options(records, directory):
    """Return case/path choices in one directory, sorted by case name and path."""
    return [
        {"label": str(record["case"]), "value": str(record["path"])}
        for record in sorted(
            (
                record
                for record in records
                if str(record["directory"]) == str(directory)
            ),
            key=lambda record: (str(record["case"]), str(record["path"])),
        )
    ]


def _graph(identifier):
    """Return a consistently styled empty graph with the supplied component ID."""
    return dcc.Graph(
        id=identifier,
        className="mlt-graph",
        config={"displaylogo": False},
    )


def _section(kicker, title, description):
    """Return a section heading and explanatory paragraph for the supplied title/text."""
    return html.Div(
        [
            html.Div(kicker, className="mlt-kicker"),
            html.H2(title),
            html.P(description),
        ],
        className="mlt-section-heading",
    )


def build_layout():
    """Return the saved-column controls and five diagnostic graph containers.

    Discover output files, select the preferred ARM file or newest available
    file, and initialize the catalog. With no files, leave selection empty;
    callbacks load metadata and populate plots when a file is chosen.
    """
    records = discover_netcdf_file_records(REPO_ROOT)
    preferred = REPO_ROOT / "output" / "dash_default" / "arm_stats.nc"
    default_record = next(
        (
            record
            for record in records
            if Path(record["path"]) == preferred
        ),
        max(records, key=lambda record: float(record.get("modified", 0.0)))
        if records
        else None,
    )
    default_path = str(default_record["path"]) if default_record else None
    default_directory = (
        str(default_record["directory"]) if default_record else None
    )

    return html.Div(
        [
            dcc.Store(
                id="amlt-file-catalog",
                data={"records": records, "preferred_path": default_path},
            ),
            html.Div(
                [
                    html.Main(
                        [
                            html.Div(
                                [
                                    html.Div(
                                        "Coherence formulation laboratory",
                                        className="mlt-kicker",
                                    ),
                                    html.H1("Affine coherence mixing length"),
                                    html.P(
                                        "Follow air parcels through a saved model "
                                        "column and see how energy barriers and "
                                        "launch settings change the mixing length."
                                    ),
                                ],
                                className="mlt-hero",
                            ),
                            html.Div(id="amlt-metrics", className="mlt-metric-grid"),
                            _section(
                                "Affine parcel calculation",
                                "Dense parcel energy and buoyancy",
                                "Every cell represents a parcel launched at the "
                                "horizontal altitude and evaluated at the vertical "
                                "altitude. Signed-log colors preserve both weak "
                                "structure and the large negative-energy tail.",
                            ),
                            _graph("amlt-parcel-matrix-figure"),
                            _section(
                                "Launch-target physics",
                                "Energy barrier and survival weight",
                                "Launch altitude runs horizontally and target "
                                "altitude vertically. The barrier is the greatest "
                                "adverse endpoint work encountered along each path; "
                                "coherence decreases continuously as that barrier grows.",
                            ),
                            _graph("amlt-matrix-figure"),
                            _section(
                                "Directional reduction",
                                "Selected overwrite comparison",
                                "Each panel compares the inverse-square "
                                "method using the selected "
                                "overwrite, alongside saved CLUBB output.",
                            ),
                            _graph("amlt-profile-figure"),
                            _section(
                                "Temporal behavior",
                                "Transition smoothness",
                                "The selected record is compared with nearby "
                                "timesteps at the most rapidly changing saved-Fortran "
                                "level and across the complete profile.",
                            ),
                            _graph("amlt-temporal-figure"),
                            _section(
                                "Coherent launch population",
                                "Coherence sums by target level",
                                "Upward, downward, and total sums show how many "
                                "coherent launch contributions reach each target "
                                "level.",
                            ),
                            _graph("amlt-coherence-sum-figure"),
                        ],
                        className="mlt-main-pane",
                    ),
                    html.Div(className="mlt-pane-divider"),
                    html.Aside(
                        [
                            html.Div(
                                [
                                    html.Div("Display controls", className="mlt-kicker"),
                                    html.H2("Coherence matrix view"),
                                    html.P(
                                        "Choose a saved model state, truncate the "
                                        "calculation if desired, and vary the current "
                                        "mixing-length parameters."
                                    ),
                                ],
                                className="mlt-sidebar-heading",
                            ),
                            html.Div(
                                [
                                    html.Label("Output directory"),
                                    styled_dropdown(
                                        id="amlt-output-directory",
                                        options=_directory_options(records),
                                        value=default_directory,
                                        clearable=False,
                                    ),
                                ],
                                className="mlt-control",
                            ),
                            html.Div(
                                [
                                    html.Label("Case"),
                                    styled_dropdown(
                                        id="amlt-case",
                                        options=_case_options(
                                            records, default_directory
                                        ),
                                        value=default_path,
                                        clearable=False,
                                    ),
                                ],
                                className="mlt-control",
                            ),
                            html.Div(
                                [
                                    html.Label("Add file by path"),
                                    dcc.Input(
                                        id="amlt-custom-path",
                                        type="text",
                                        placeholder="/path/to/case_stats.nc",
                                    ),
                                    html.Div(
                                        [
                                            html.Button(
                                                "Load",
                                                id="amlt-load-path",
                                                n_clicks=0,
                                            ),
                                            html.Button(
                                                "Refresh",
                                                id="amlt-refresh-files",
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
                                        id="amlt-column",
                                        options=[{"label": "0", "value": 0}],
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
                                                id="amlt-mu-label",
                                                className="mlt-slider-value",
                                            ),
                                        ]
                                    ),
                                    dcc.Slider(
                                        id="amlt-mu",
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
                                    ),
                                    html.Div(
                                        "Initialized from Lscale_mu_entrain in the "
                                        "NetCDF parameter table.",
                                        className="mlt-control-help",
                                    ),
                                ],
                                className="mlt-control mlt-mu-control",
                            ),
                            html.Div(
                                [
                                    html.Label("Overwrite method"),
                                    dcc.RadioItems(
                                        id="amlt-overwrite-mode",
                                        options=[
                                            {
                                                "label": "No overwrite",
                                                "value": "no_overwrite",
                                            },
                                            {
                                                "label": "Legacy max overwrite",
                                                "value": "legacy_max",
                                            },
                                            {
                                                "label": "Coherence averaged",
                                                "value": "coherence_averaged",
                                            },
                                        ],
                                        value="coherence_averaged",
                                        labelStyle={"display": "block"},
                                    ),
                                    html.Div(
                                        "Changes how calculated trajectories are "
                                        "combined; saved Fortran output is unchanged.",
                                        className="mlt-control-help",
                                    ),
                                ],
                                className="mlt-control",
                            ),
                            html.Div(
                                [
                                    html.Label(
                                        [
                                            "Maximum calculated altitude",
                                            html.Span(
                                                id="amlt-altitude-max-label",
                                                className="mlt-slider-value",
                                            ),
                                        ]
                                    ),
                                    dcc.Slider(
                                        id="amlt-altitude-max",
                                        min=500.0,
                                        max=20000.0,
                                        step=100.0,
                                        value=3500.0,
                                        marks={
                                            500.0: "0.5 km",
                                            2000.0: "2 km",
                                            3500.0: "3.5 km",
                                            10000.0: "10 km",
                                            20000.0: "20 km",
                                        },
                                        updatemode="mouseup",
                                    ),
                                    html.Div(
                                        "Sets a new upper boundary for the "
                                        "calculation, rather than zooming the plot.",
                                        className="mlt-control-help",
                                    ),
                                ],
                                className="mlt-control mlt-mu-control",
                            ),
                            html.Div(
                                [
                                    html.Label(
                                        [
                                            "Launch-energy TKE multiplier",
                                            html.Span(
                                                id="amlt-tke-multiplier-label",
                                                className="mlt-slider-value",
                                            ),
                                        ]
                                    ),
                                    dcc.Input(
                                        id="amlt-tke-multiplier",
                                        type="number",
                                        min=0.0,
                                        step="any",
                                        value=1.0,
                                        debounce=True,
                                    ),
                                    html.Div(
                                        "Mean launch energy divided by local TKE. "
                                        "Nonnegative, with no upper limit; 1 uses TKE.",
                                        className="mlt-control-help",
                                    ),
                                ],
                                className="mlt-control mlt-mu-control",
                            ),
                            html.Div(
                                [
                                    html.Label(
                                        [
                                            "PDF launch fraction",
                                            html.Span(
                                                id="amlt-pdf-launch-frac-label",
                                                className="mlt-slider-value",
                                            ),
                                        ]
                                    ),
                                    dcc.Slider(
                                        id="amlt-pdf-launch-frac",
                                        min=0.0,
                                        max=1.0,
                                        step=0.01,
                                        value=0.0,
                                        marks={
                                            0.0: "0 mean",
                                            0.25: "0.25",
                                            0.5: "0.50",
                                            0.75: "0.75",
                                            1.0: "1 PDF",
                                        },
                                        updatemode="mouseup",
                                    ),
                                    html.Div(
                                        "Interpolates directly from the mean state to "
                                        "PDF component 1 upward and component 2 "
                                        "downward; mixture fraction is not applied.",
                                        className="mlt-control-help",
                                    ),
                                ],
                                className="mlt-control mlt-mu-control",
                            ),
                            html.Div(
                                [
                                    html.Label(
                                        [
                                            "Upward mixing-length weight",
                                            html.Span(
                                                id="amlt-up-weight-label",
                                                className="mlt-slider-value",
                                            ),
                                        ]
                                    ),
                                    dcc.Slider(
                                        id="amlt-up-weight",
                                        min=0.0,
                                        max=1.0,
                                        step=0.01,
                                        value=0.5,
                                        marks={
                                            0.0: "down",
                                            0.25: "0.25",
                                            0.5: "equal",
                                            0.75: "0.75",
                                            1.0: "up",
                                        },
                                        updatemode="mouseup",
                                    ),
                                    html.Div(
                                        "Exponent on Lscale_up in the weighted "
                                        "geometric directional combination.",
                                        className="mlt-control-help",
                                    ),
                                ],
                                className="mlt-control mlt-mu-control",
                            ),
                            html.Div(className="mlt-sidebar-rule"),
                            html.Div(
                                [
                                    html.Label("NetCDF record / timestep"),
                                    dcc.Slider(
                                        id="amlt-record",
                                        min=0,
                                        max=1,
                                        step=1,
                                        value=0,
                                        marks={0: "0", 1: "1"},
                                        updatemode="mouseup",
                                    ),
                                ],
                                className="mlt-control mlt-time-slider",
                            ),
                            html.Div(id="amlt-time-label", className="mlt-time-label"),
                            html.Div(id="amlt-file-status", className="mlt-status"),
                        ],
                        className="mlt-right-pane",
                    ),
                ],
                className="mlt-workspace",
            )
        ],
        className="mlt-page",
    )
