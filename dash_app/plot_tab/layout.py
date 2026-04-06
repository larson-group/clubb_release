from dash import dcc, html

from .plot_types import shared
from .plot_types.registry import PLOT_TYPES
from .state import DEFAULT_OUTPUT_DIR, DEFAULT_OUTPUT_ENTRY, DEFAULT_PLAYBACK_INTERVAL_S


CASE_BUTTON_STYLE = {
    "padding": "6px 10px",
    "margin": "4px",
    "borderRadius": "4px",
    "fontSize": "16px",
}
DIR_ENTRY_STYLE = {"width": "420px", "fontSize": "14pt"}
MODE_RADIO_LABEL_STYLE = {"display": "inline-block", "marginRight": "16px", "fontSize": "14pt"}
SECTION_HEADING_STYLE = {"fontSize": "16pt"}
RIGHT_PANE_STYLE = {
    "paddingLeft": "16px",
    "height": "calc(100vh - 96px)",
    "minHeight": 0,
    "overflowY": "auto",
    "overflowX": "auto",
}
GRID_STYLE = {
    "display": "grid",
    "gridTemplateColumns": "repeat(auto-fill, minmax(420px, 1fr))",
    "gridAutoFlow": "dense",
    "gap": "16px",
    "marginTop": "16px",
}


def case_button(name, available, selected=False):
    """Render one case-selection button with the current availability styling."""
    return html.Button(
        name,
        id={"type": "plots-case-button", "name": name},
        disabled=not available,
        style={
            **CASE_BUTTON_STYLE,
            "backgroundColor": "#2563eb" if available else "#c9c9c9",
            "color": "#ffffff" if available else "#5f5f5f",
            "border": "2px solid #f59e0b" if selected else "1px solid transparent",
            "cursor": "pointer" if available else "not-allowed",
        },
    )


def benchmark_button(source, label, available, selected=False):
    """Render one benchmark toggle button with case-button-like styling."""
    return html.Button(
        label,
        id={"type": "plots-benchmark-button", "source": source},
        n_clicks=0,
        n_clicks_timestamp=-1,
        disabled=not available,
        style={
            **CASE_BUTTON_STYLE,
            "backgroundColor": "#2563eb" if available else "#c9c9c9",
            "color": "#ffffff" if available else "#5f5f5f",
            "border": "2px solid #f59e0b" if selected else "1px solid transparent",
            "cursor": "pointer" if available else "not-allowed",
        },
    )


def directory_entry(index, value):
    """Render one editable directory row in the header controls."""
    return html.Div(
        [
            html.Label(f"Directory {index + 1}", style={"marginRight": "8px", "whiteSpace": "nowrap", "fontSize": "14pt"}),
            dcc.Input(
                id={"type": "plots-dir-entry", "index": index},
                type="text",
                value=value,
                style=DIR_ENTRY_STYLE,
            ),
            html.Button(
                "Remove",
                id={"type": "plots-remove-dir", "index": index},
                n_clicks=0,
                style={"marginLeft": "8px", "fontSize": "14pt"},
            ),
        ],
        style={"marginBottom": "10px", "display": "flex", "alignItems": "center", "gap": "0px"},
    )


def render_plot_card(plot_id, plot_state, case_data):
    """Render a single plot card and splice in any plot-specific auxiliary nodes."""
    plot_type = plot_state.get("plot_type")
    module = PLOT_TYPES[plot_type]
    card = module.render_card(plot_id, plot_state, {"case_data": case_data})
    aux = list(module.auxiliary_components(plot_id))
    if aux:
        card.children[2:2] = aux
    return card


def add_plot_controls_card():
    """Render the trailing grid card that exposes the add-plot actions."""
    button_style = {"width": "70%", "fontSize": "16pt"}
    return html.Div(
        [
            html.Div("", style={"fontWeight": "600", "fontSize": "16px", "marginBottom": "10px"}),
            html.Div(
                [
                    html.Button("Add plot", id="plots-add-profile", n_clicks=0, style=button_style),
                    html.Button("Add budget plot", id="plots-add-budget", n_clicks=0, style={**button_style, "marginTop": "10px"}),
                    html.Button("Add time-series plot", id="plots-add-timeseries", n_clicks=0, style={**button_style, "marginTop": "10px"}),
                    html.Button("Add time-height plot", id="plots-add-timeheight", n_clicks=0, style={**button_style, "marginTop": "10px"}),
                    html.Button("Add subcol plot", id="plots-add-subcolumn", n_clicks=0, style={**button_style, "marginTop": "10px"}),
                ],
                style={"display": "flex", "flexDirection": "column", "alignItems": "center", "justifyContent": "center", "minHeight": "100%"},
            ),
        ],
        id="plots-add-card",
        key="plots-add-card",
        style={
            "padding": "10px",
            "alignSelf": "stretch",
            "display": "flex",
            "flexDirection": "column",
            "justifyContent": "center",
            "background": "transparent",
            "border": "none",
            "boxShadow": "none",
        },
    )


def render_plot_grid(plot_order, plot_state, case_data):
    """Render the ordered plot cards followed by the add-controls card."""
    children = []
    for plot_id in plot_order or []:
        state = (plot_state or {}).get(str(plot_id))
        if state is None:
            continue
        children.append(render_plot_card(int(plot_id), state, case_data))
    children.append(add_plot_controls_card())
    return children


def child_id(child):
    """Extract a Dash child id from either a component instance or serialized dict."""
    if isinstance(child, dict):
        return child.get("props", {}).get("id")
    return getattr(child, "id", None)


def _directory_case_selector():
    """Build the combined directory/case selection header block."""
    return html.Div(
        [
            html.Div(
                [
                    html.Div(
                        [
                            html.Div(id="plots-dir-list"),
                            html.Div(
                                [
                                    html.Button("Add directory", id="plots-add-dir-row", n_clicks=0, style={"fontSize": "14pt"}),
                                    html.Button("Refresh", id="plots-refresh-cases", n_clicks=0, style={"fontSize": "14pt", "marginLeft": "8px"}),
                                ],
                                style={"marginTop": "8px", "display": "flex", "alignItems": "center"},
                            ),
                        ],
                        style={"display": "flex", "flexDirection": "column", "alignItems": "flex-start", "flex": "1 1 auto"},
                    ),
                    html.Div(style={"width": "1px", "backgroundColor": "#d0d0d0", "alignSelf": "stretch", "margin": "0 16px"}),
                    html.Div(
                        [
                            html.Div("Benchmarks", style={"fontWeight": "600", "marginBottom": "8px"}),
                            html.Div(id="plots-benchmark-button-container"),
                        ],
                        style={"display": "flex", "flexDirection": "column", "justifyContent": "center", "minWidth": "220px"},
                    ),
                ],
                style={"padding": "12px", "display": "flex", "alignItems": "stretch"},
            ),
            html.Div(style={"width": "1px", "backgroundColor": "#d0d0d0", "alignSelf": "stretch"}),
            html.Div(
                [
                    html.Div("Cases", style={"fontWeight": "600", "marginBottom": "8px"}),
                    html.Div(id="plots-case-button-container"),
                ],
                style={"padding": "12px", "minHeight": "100%"},
            ),
        ],
        style={
            "display": "grid",
            "gridTemplateColumns": "1fr 1px 1fr",
            "gap": "16px",
            "alignItems": "stretch",
            "padding": "12px",
            "border": "1px solid #d0d0d0",
            "borderRadius": "6px",
        },
    )


def _plots_stores(initial_state):
    """Build the stores and interval used to coordinate the plots tab."""
    return [
        dcc.Store(id="plots-output-dir-entries", data=[DEFAULT_OUTPUT_ENTRY]),
        dcc.Store(id="plots-output-dirs", data=[DEFAULT_OUTPUT_DIR]),
        dcc.Store(id="plots-case-data", data=initial_state["case_data"]),
        dcc.Store(id="plots-enabled-benchmarks", data=initial_state["enabled_benchmarks"]),
        dcc.Store(id="plots-plot-order", data=initial_state["plot_order"]),
        dcc.Store(id="plots-plot-state", data=initial_state["plot_state"]),
        dcc.Store(id="plots-next-id", data=initial_state["next_id"]),
        dcc.Store(id="plots-last-add-ts", data=0),
        dcc.Store(id="plots-param-data", data=None),
        dcc.Store(id="plots-param-names", data=None),
        dcc.Store(id="plots-selected-column", data=initial_state["selected_column"]),
        dcc.Store(id="plots-playback", data={"playing": False, "interval_s": DEFAULT_PLAYBACK_INTERVAL_S, "inflight": False, "target_point": None}),
        dcc.Interval(id="plots-playback-interval", interval=int(DEFAULT_PLAYBACK_INTERVAL_S * 1000), disabled=True, n_intervals=0),
    ]


def _left_pane(initial_state):
    """Build the grid pane that holds plot cards and the add-controls card."""
    return html.Div(
        [
            html.Div(
                id="plots-plot-container",
                children=render_plot_grid(initial_state["plot_order"], initial_state["plot_state"], initial_state["case_data"]),
                style=GRID_STYLE,
            ),
        ],
        id="plots-left-pane",
        className="plots-left-pane",
    )


def _time_section(initial_state):
    """Build the time controls section in the right-hand UI pane."""
    return [
        html.Div(id="plots-time-heading", className="run-settings-heading", children="Time", style=SECTION_HEADING_STYLE),
        dcc.RadioItems(
            id="plots-time-mode",
            options=[{"label": "Average Range", "value": "range"}, {"label": "Single Time", "value": "point"}],
            value=initial_state["time_mode"],
            labelStyle=MODE_RADIO_LABEL_STYLE,
            style={"marginBottom": "10px", "textAlign": "center"},
        ),
        html.Div(
            dcc.RangeSlider(
                id="plots-global-time-range",
                min=1,
                max=initial_state["time_slider_max"],
                value=initial_state["time_range"],
                step=1,
                allowCross=False,
                marks=initial_state["time_marks"],
                tooltip={"always_visible": True, "placement": "bottom"},
            ),
            id="plots-global-time-range-wrapper",
            className="plots-slider-block",
        ),
        html.Div(
            [
                html.Div(
                    [
                        html.Button("<<", id="plots-playback-slower", n_clicks=0, style={"minWidth": "42px", "fontSize": "16pt"}),
                        html.Button("Play (1s)", id="plots-playback-toggle", n_clicks=0, style={"minWidth": "110px", "fontSize": "16pt"}),
                        html.Button(">>", id="plots-playback-faster", n_clicks=0, style={"minWidth": "42px", "fontSize": "16pt"}),
                    ],
                    id="plots-playback-controls",
                    style={"display": "flex", "justifyContent": "center", "alignItems": "center", "gap": "10px", "marginBottom": "16px"},
                ),
                dcc.Slider(
                    id="plots-global-time-point",
                    min=1,
                    max=initial_state["time_point_max"],
                    value=initial_state["time_point"],
                    step=1,
                    marks=initial_state["time_point_marks"],
                    tooltip={"always_visible": True, "placement": "bottom"},
                    included=False,
                ),
            ],
            id="plots-global-time-point-wrapper",
            className="plots-slider-block",
            style={"display": "none"},
        ),
    ]


def _height_section(initial_state):
    """Build the global height controls section in the right-hand UI pane."""
    return [
        html.Div(id="plots-height-heading", className="run-settings-heading", children="Height", style=SECTION_HEADING_STYLE),
        html.Div(
            dcc.RangeSlider(
                id="plots-global-height-range",
                min=initial_state["height_min"],
                max=initial_state["height_max"],
                value=initial_state["height_range"],
                step=initial_state["height_step"],
                allowCross=False,
                marks=initial_state["height_marks"],
                tooltip={"always_visible": True, "placement": "bottom"},
            ),
            className="plots-slider-block",
        ),
        html.Div(style={"marginTop": "16px", "marginBottom": "24px"}),
    ]


def _column_section(initial_state):
    """Build the column-mode and parameter controls section."""
    return [
        html.Div("Columns", id="plots-column-heading", className="run-settings-heading", style=SECTION_HEADING_STYLE),
        dcc.RadioItems(
            id="plots-column-mode",
            options=[{"label": "Single", "value": "single"}, {"label": "Overplot", "value": "all"}],
            value=initial_state["column_mode"],
            labelStyle=MODE_RADIO_LABEL_STYLE,
            style={"marginBottom": "10px", "textAlign": "center"},
        ),
        html.Div(id="plots-param-panel", style={"overflowY": "auto", "maxHeight": "none", "paddingRight": "8px", "marginTop": "10px"}),
    ]


def _right_pane(initial_state):
    """Build the full controls pane shown to the right of the plots grid."""
    return html.Div(
        _time_section(initial_state) + _height_section(initial_state) + _column_section(initial_state),
        id="plots-right-pane",
        className="plots-right-pane",
        style=RIGHT_PANE_STYLE,
    )


def build_layout(initial_state):
    """Assemble the full static plots-tab layout from the provided initial state."""
    return html.Div(
        [
            html.Div([_directory_case_selector()], style={"marginBottom": "12px"}),
            *_plots_stores(initial_state),
            html.Div(
                [_left_pane(initial_state), html.Div(id="plots-pane-divider", className="plots-pane-divider"), _right_pane(initial_state)],
                id="plots-tab-layout",
                className="plots-tab-layout",
                style={"display": "grid", "gridTemplateColumns": "minmax(0,1fr) 8px 540px", "gap": "16px", "marginTop": "16px"},
            ),
        ],
        style={"padding": "10px"},
    )
