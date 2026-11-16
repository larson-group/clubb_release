from dash import dcc, html

from .plot_types import shared
from .case_cache import compact_case_data
from .plot_types.registry import PLOT_TYPES
from .state import DEFAULT_PLAYBACK_INTERVAL_S
from .output_chooser import build_output_chooser


CASE_BUTTON_STYLE = {
    "padding": "6px 10px",
    "margin": "4px",
    "borderRadius": "4px",
    "fontSize": "16px",
}
MODE_RADIO_LABEL_STYLE = {"display": "inline-block", "marginRight": "16px", "fontSize": "14pt"}
SECTION_HEADING_STYLE = {"fontSize": "16pt"}
RIGHT_PANE_STYLE = {
    "paddingLeft": "16px",
    "alignSelf": "start",
    "minHeight": 0,
    "boxSizing": "border-box",
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


def benchmark_button(source, label, available, selected=False, acquisition=None):
    """Show source selection, an animated preparation indicator, or its acquisition error."""
    state = (acquisition or {}).get("state")
    loading = state == "loading"
    error = (acquisition or {}).get("error") if state == "error" else None
    button = html.Button(
        [html.Span(className="plots-les-spinner", **{"aria-hidden": "true"}), label,
         html.Span(" · Loading…", className="plots-les-status")] if loading else label,
        id={"type": "plots-benchmark-button", "source": source},
        n_clicks=0,
        n_clicks_timestamp=-1,
        disabled=not available,
        title=error or ("Finding or downloading LES data; plots remain usable." if loading else label),
        **{"aria-busy": "true" if loading else "false"},
        style={
            **CASE_BUTTON_STYLE,
            "backgroundColor": "#2563eb" if available else "#c9c9c9",
            "color": "#ffffff" if available else "#5f5f5f",
            "border": "2px solid #f59e0b" if selected else "1px solid transparent",
            "cursor": "pointer" if available else "not-allowed",
        },
    )
    if error:
        return html.Div([button, html.Div(
            ["LES unavailable: ", error, " Toggle off and on to retry."],
            className="plots-les-error", role="alert",
        )])
    return button


def render_plot_card(plot_id, plot_state, case_data):
    """Render a single plot card and splice in any plot-specific auxiliary nodes."""
    plot_type = plot_state.get("plot_type")
    module = PLOT_TYPES[plot_type]
    card = module.render_card(plot_id, plot_state, {"case_data": case_data})
    aux = list(module.auxiliary_components(plot_id))
    if aux:
        card.children.extend(aux)
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
                    html.Button("Add Custom", id="plots-add-custom", n_clicks=0, style={**button_style, "marginTop": "10px"}),
                    html.Button("Add PDF contours", id="plots-add-pdf-contour", n_clicks=0, style={**button_style, "marginTop": "10px"}),
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
    seen_plot_ids = set()
    for plot_id in plot_order or []:
        normalized_id = int(plot_id)
        # A stale browser callback can transiently repeat an id in persisted
        # order state.  Never mount duplicate pattern-matching component ids.
        if normalized_id in seen_plot_ids:
            continue
        seen_plot_ids.add(normalized_id)
        state = (plot_state or {}).get(str(plot_id))
        if state is None:
            continue
        children.append(render_plot_card(normalized_id, state, case_data))
    children.append(add_plot_controls_card())
    return children


def child_id(child):
    """Extract a Dash child id from either a component instance or serialized dict."""
    if isinstance(child, dict):
        return child.get("props", {}).get("id")
    return getattr(child, "id", None)


def _initial_case_buttons(initial_state):
    """Render initial case buttons before callback hydration."""
    selected_name = ((initial_state or {}).get("case_data") or {}).get("name")
    available_names = ((initial_state or {}).get("case_data") or {}).get("available_cases") or []
    if not available_names:
        return [html.Div("No cases found in the active outputs.")]
    return [case_button(name, True, selected=(name == selected_name)) for name in available_names]


def _pyplotgen_controls():
    """Keep static gallery export and its help-driven options in one header category."""
    return html.Div(
        [
            html.Div([
                html.Span('PyPlotGen'),
                html.Button('?', id='plots-pyplotgen-help-open', n_clicks=0,
                            title='About PyPlotGen and its command-line options',
                            **{'aria-label': 'PyPlotGen help'}),
            ], className='plots-pyplotgen-heading'),
            dcc.Dropdown(id='plots-pyplotgen-options', options=[], value=['-les'],
                         multi=True, searchable=True, closeOnSelect=False,
                         placeholder='Search options…', className='clubb-dropdown'),
            html.Small('Type an option and press Enter; repeat to add more.'),
            html.Div(id='plots-pyplotgen-option-values', className='plots-pyplotgen-option-values'),
            html.Div(id='plots-pyplotgen-help-status', role='status'),
            html.Button(
                html.Span('Generate PyPlotGen', id='plots-pyplotgen-run-label'),
                id='plots-pyplotgen-run', n_clicks=0, className='plots-pyplotgen-button',
                title='Generate a static PyPlotGen gallery from the active outputs',
            ),
            html.Div(id='plots-pyplotgen-status', className='plots-pyplotgen-status'),
        ],
        className='plots-pyplotgen-control',
    )


def _directory_case_selector(initial_state):
    """Build the combined directory/case selection header block."""
    return html.Div(
        [
            html.Div(
                [
                    build_output_chooser(),
                    html.Div(className="plots-output-benchmark-divider"),
                    html.Div(
                        [
                            html.Div("Benchmarks", style={"fontWeight": "600", "marginBottom": "8px"}),
                            html.Div(id="plots-benchmark-button-container"),
                        ],
                        className="plots-benchmark-column",
                    ),
                    html.Div(className="plots-output-benchmark-divider"),
                    _pyplotgen_controls(),
                ],
                className="plots-output-benchmark-row",
            ),
            html.Div(style={"width": "1px", "backgroundColor": "#d0d0d0", "alignSelf": "stretch"}),
            html.Div(
                [
                    html.Div("Cases", style={"fontWeight": "600", "marginBottom": "8px"}),
                    html.Div(id="plots-case-button-container", children=_initial_case_buttons(initial_state)),
                    html.Div(id="plots-load-status", role="status"),
                ],
                style={"padding": "12px", "minHeight": "100%"},
            ),
        ],
        className="plots-directory-case-selector",
        style={
            "padding": "12px",
            "border": "1px solid #d0d0d0",
            "borderRadius": "6px",
        },
    )


def _plots_stores(initial_state):
    """Build plot/control state; the output chooser supplies its own stores."""
    return [
        dcc.Store(id="plots-case-data", data=compact_case_data(initial_state["case_data"])),
        dcc.Store(id="plots-case-selection", data=None),
        dcc.Store(id="plots-enabled-benchmarks", data=initial_state["enabled_benchmarks"]),
        dcc.Store(id="plots-ready-benchmarks", data=None),
        dcc.Store(id="plots-benchmark-sam", data=None),
        dcc.Store(id="plots-benchmark-coamps", data=None),
        dcc.Store(id="plots-plot-order", data=initial_state["plot_order"]),
        dcc.Store(id="plots-plot-state", data=initial_state["plot_state"]),
        dcc.Store(id="plots-next-id", data=initial_state["next_id"]),
        dcc.Store(id="plots-instance-snapshot", data=None),
        dcc.Store(id="plots-last-add-ts", data=0),
        dcc.Store(id="plots-param-data", data=None),
        dcc.Store(id="plots-param-names", data=None),
        dcc.Store(id="plots-column-filters", data={"indices": None, "filters": {}}),
        dcc.Store(id="plots-selected-column", data=initial_state["selected_column"]),
        dcc.Store(id="plots-time-override", data=None),
        dcc.Store(id="plots-pyplotgen-action", data=None),
        dcc.Store(id="plots-pyplotgen-help", data=None),
        dcc.Store(id="plots-pyplotgen-popup", data=None),
        dcc.Store(id="plots-pyplotgen-opened-run", data=None),
        dcc.Store(id="plots-playback", data={"playing": False, "interval_s": DEFAULT_PLAYBACK_INTERVAL_S, "inflight": False, "target_point": None}),
        dcc.Interval(id="plots-pyplotgen-progress-interval", interval=500, disabled=True, n_intervals=0),
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
        html.Div(
            [
                html.Button("Loss window", id="plots-use-loss-window", n_clicks=0, style={**CASE_BUTTON_STYLE, "backgroundColor": "#2563eb", "color": "#ffffff"}),
                html.Button("Pyplotgen window", id="plots-use-pyplotgen-window", n_clicks=0, style={**CASE_BUTTON_STYLE, "backgroundColor": "#334155", "color": "#ffffff"}),
            ],
            style={"display": "flex", "justifyContent": "center", "flexWrap": "wrap", "gap": "8px", "marginBottom": "10px"},
        ),
        html.Div(
            [
                html.Div(id="plots-time-start-label", children="Start time", style={"fontSize": "13pt", "fontWeight": "600", "marginBottom": "4px"}),
                dcc.Slider(
                    id="plots-global-time-point",
                    min=initial_state["time_point_min"],
                    max=initial_state["time_point_max"],
                    value=initial_state["time_point"],
                    step=initial_state["time_point_step"],
                    marks=initial_state["time_point_marks"],
                    tooltip={"always_visible": True, "placement": "bottom", "transform": "formatPlotSeconds"},
                    included=False,
                    dots=True,
                ),
            ],
            id="plots-global-time-point-wrapper",
            className="plots-slider-block plots-time-slider-block",
        ),
        html.Div(
            [
                html.Div(id="plots-time-average-label", children="Average Length", style={"fontSize": "13pt", "fontWeight": "600", "marginBottom": "4px"}),
                dcc.Slider(
                    id="plots-global-time-range",
                    min=initial_state["time_slider_min"],
                    max=initial_state["time_slider_max"],
                    value=initial_state["time_range"],
                    step=initial_state["time_slider_step"],
                    marks=initial_state["time_marks"],
                    tooltip={"always_visible": True, "placement": "bottom", "transform": "formatPlotMinutes"},
                    included=False,
                ),
            ],
            id="plots-global-time-range-wrapper",
            className="plots-slider-block plots-time-slider-block",
        ),
        html.Div(
            [
                html.Button("<<", id="plots-playback-slower", n_clicks=0, style={"minWidth": "42px", "fontSize": "16pt"}),
                html.Button("Play (1s)", id="plots-playback-toggle", n_clicks=0, style={"minWidth": "110px", "fontSize": "16pt"}),
                html.Button(">>", id="plots-playback-faster", n_clicks=0, style={"minWidth": "42px", "fontSize": "16pt"}),
            ],
            id="plots-playback-controls",
            style={"display": "flex", "justifyContent": "center", "alignItems": "center", "gap": "10px", "marginBottom": "16px"},
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
        html.Div(id="plots-param-status", role="status"),
        html.Div(id="plots-param-panel", style={"paddingRight": "8px", "marginTop": "10px"}),
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
            html.Div([_directory_case_selector(initial_state)], style={"marginBottom": "12px"}),
            *_plots_stores(initial_state),
            html.Div(id="plots-help-modal"),
            html.Div(id="plots-pyplotgen-help-modal"),
            html.Div(
                [_left_pane(initial_state), html.Div(id="plots-pane-divider", className="plots-pane-divider"), _right_pane(initial_state)],
                id="plots-tab-layout",
                className="plots-tab-layout",
                style={
                    "display": "grid",
                    "gridTemplateColumns": "minmax(0,1fr) 8px minmax(360px, min(540px, 34vw))",
                    "gap": "16px",
                    "marginTop": "16px",
                    "alignItems": "start",
                },
            ),
        ],
        style={"padding": "10px"},
    )
