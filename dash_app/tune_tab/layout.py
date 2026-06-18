"""Static layout builders for the tuning tab."""

from __future__ import annotations

from dash import dcc, html

from tuner.taylor_metrics import DEFAULT_AGGREGATION_MODE, DEFAULT_LOSS_MODE


def action_button_style(color, disabled=False):
    """Return the shared style for tuning-tab action buttons."""
    style = {
        "backgroundColor": color,
        "color": "#ffffff",
        "border": "none",
        "padding": "10px 16px",
        "margin": "4px",
        "borderRadius": "6px",
        "cursor": "pointer",
        "fontSize": "14px",
        "fontWeight": "600",
    }
    if disabled:
        style.update({"backgroundColor": "#9ca3af", "color": "#f3f4f6", "cursor": "not-allowed"})
    return style


def mode_button_style(selected=False, disabled=False):
    """Return the style for tuning-mode selector buttons."""
    style = action_button_style("#2563eb" if selected else "#374151", disabled=disabled)
    style.update({
        "border": "2px solid #f59e0b" if selected else "1px solid transparent",
        "minWidth": "96px",
    })
    return style


def mode_options_block_style(visible: bool) -> dict:
    """Return show/hide style for a mode-options block."""
    return {"display": "block" if visible else "none"}


def build_top_controls(initial_data):
    """Render mode selection, mode options, and run controls."""
    selected_mode = initial_data.get("strategy_mode")
    selected_loss_mode = initial_data.get("loss_mode", DEFAULT_LOSS_MODE)
    selected_aggregation_mode = initial_data.get("aggregation_mode", DEFAULT_AGGREGATION_MODE)
    return html.Div(
        [
            html.Div(
                [
                    html.Div("Mode", className="tune-section-title"),
                    html.Div(
                        [
                            html.Button(
                                "Random",
                                id="tune-mode-random",
                                n_clicks=0,
                                style=mode_button_style(selected=selected_mode == "random"),
                            ),
                            html.Button(
                                "Resolve",
                                id="tune-mode-resolve",
                                n_clicks=0,
                                style=mode_button_style(selected=selected_mode == "resolve"),
                            ),
                        ],
                        style={"display": "flex", "gap": "8px", "flexWrap": "wrap"},
                    ),
                    html.Div(
                        "Select a mode.",
                        id="tune-no-mode-options",
                        style={**mode_options_block_style(selected_mode is None), "marginTop": "10px"},
                    ),
                    html.Div(
                        [
                            html.Label("Max Samples", htmlFor="tune-random-max-samples"),
                            dcc.Input(
                                id="tune-random-max-samples",
                                type="number",
                                min=1,
                                step=1,
                                value=initial_data["random_max_samples"],
                                debounce=True,
                                style={"width": "140px", "marginLeft": "8px"},
                            ),
                        ],
                        id="tune-random-options",
                        style={**mode_options_block_style(selected_mode == "random"), "marginTop": "10px"},
                    ),
                    html.Div(
                        [
                            html.Label("Spacing", htmlFor="tune-resolve-spacing"),
                            dcc.Input(
                                id="tune-resolve-spacing",
                                type="number",
                                min=0,
                                step="any",
                                value=initial_data["resolve_spacing"],
                                debounce=True,
                                style={"width": "140px", "marginLeft": "8px"},
                            ),
                            html.Div(
                                id="tune-resolve-total-samples",
                                style={"marginTop": "8px", "opacity": "0.85"},
                            ),
                        ],
                        id="tune-resolve-options",
                        style={**mode_options_block_style(selected_mode == "resolve"), "marginTop": "10px"},
                    ),
                ],
                style={"padding": "12px"},
            ),
            html.Div(style={"width": "1px", "backgroundColor": "#d0d0d0", "alignSelf": "stretch"}),
            html.Div(
                [
                    html.Div("Loss Mode", className="tune-section-title"),
                    html.Div(
                        [
                            html.Button(
                                "Scaled RMSE",
                                id="tune-loss-mode-scaled-rmse",
                                n_clicks=0,
                                title="Use the scaled profile mismatch returned by the loss driver.",
                                style=mode_button_style(selected=selected_loss_mode == "scaled_rmse"),
                            ),
                            html.Button(
                                "Centered RMS + Bias",
                                id="tune-loss-mode-centered-rmse-bias",
                                n_clicks=0,
                                title="Legacy Python smart loss: centered_rmse_norm + abs(bias_norm).",
                                style=mode_button_style(selected=selected_loss_mode == "centered_rmse_bias"),
                            ),
                            html.Button(
                                "Taylor Components",
                                id="tune-loss-mode-taylor-components",
                                n_clicks=0,
                                title="Centered RMSE plus mild bias, correlation, and std-ratio penalties.",
                                style=mode_button_style(selected=selected_loss_mode == "taylor_components"),
                            ),
                            html.Button(
                                "Squared",
                                id="tune-loss-mode-taylor-components-squared",
                                n_clicks=0,
                                title="Squared Taylor components; rejects extreme component failures more strongly.",
                                style=mode_button_style(selected=selected_loss_mode == "taylor_components_squared"),
                            ),
                            html.Button(
                                "Shape First",
                                id="tune-loss-mode-shape-first",
                                n_clicks=0,
                                title="More direct weight on vertical-shape correlation and variability.",
                                style=mode_button_style(selected=selected_loss_mode == "shape_first"),
                            ),
                            html.Button(
                                "Bias Light",
                                id="tune-loss-mode-bias-light-taylor",
                                n_clicks=0,
                                title="Prioritizes centered vertical structure while retaining light bias pressure.",
                                style=mode_button_style(selected=selected_loss_mode == "bias_light_taylor"),
                            ),
                            html.Button(
                                "Decomposed Taylor",
                                id="tune-loss-mode-decomposed-taylor",
                                n_clicks=0,
                                title="abs(bias_norm) + (1 - R) + abs(log(std_ratio)); omits centered RMSE.",
                                style=mode_button_style(selected=selected_loss_mode == "decomposed_taylor"),
                            ),
                        ],
                        style={"display": "flex", "gap": "8px", "flexWrap": "wrap"},
                    ),
                ],
                style={"padding": "12px"},
            ),
            html.Div(style={"width": "1px", "backgroundColor": "#d0d0d0", "alignSelf": "stretch"}),
            html.Div(
                [
                    html.Div("Aggregation", className="tune-section-title"),
                    html.Div(
                        [
                            html.Button(
                                "Mean/Max",
                                id="tune-aggregation-mean-max",
                                n_clicks=0,
                                title="Legacy aggregation: 0.6 weighted mean + 0.4 max active loss.",
                                style=mode_button_style(selected=selected_aggregation_mode == "mean_max"),
                            ),
                            html.Button(
                                "Worst Quantile",
                                id="tune-aggregation-mean-worst-quantile",
                                n_clicks=0,
                                title="0.7 weighted mean + 0.3 weighted mean of the worst 25% active losses.",
                                style=mode_button_style(selected=selected_aggregation_mode == "mean_worst_quantile"),
                            ),
                        ],
                        style={"display": "flex", "gap": "8px", "flexWrap": "wrap"},
                    ),
                ],
                style={"padding": "12px"},
            ),
            html.Div(style={"width": "1px", "backgroundColor": "#d0d0d0", "alignSelf": "stretch"}),
            html.Div(
                [
                    html.Div("Run Tuner", className="tune-section-title"),
                    html.Button(
                        "Start",
                        id="tune-start-button",
                        n_clicks=0,
                        style=action_button_style("#16a34a", disabled=True),
                        disabled=True,
                    ),
                    html.Button(
                        "Stop",
                        id="tune-stop-button",
                        n_clicks=0,
                        style=action_button_style("#b91c1c", disabled=True),
                        disabled=True,
                    ),
                ],
                style={"padding": "12px", "display": "flex", "flexDirection": "column", "alignItems": "flex-start"},
            ),
            html.Div(style={"width": "1px", "backgroundColor": "#d0d0d0", "alignSelf": "stretch"}),
            html.Div(
                [
                    html.Div("Run best results", className="tune-section-title"),
                    html.Button(
                        "Run loss window",
                        id={"type": "tune-loss-run-button", "action": "window"},
                        n_clicks=0,
                        title="Run every listed parameter set using the selected case time window.",
                        style=action_button_style("#2563eb", disabled=True),
                        disabled=True,
                    ),
                    html.Button(
                        "Run complete",
                        id={"type": "tune-loss-run-button", "action": "complete"},
                        n_clicks=0,
                        title="Run every listed parameter set with default full case settings.",
                        style=action_button_style("#2563eb", disabled=True),
                        disabled=True,
                    ),
                ],
                style={"padding": "12px", "display": "flex", "flexDirection": "column", "alignItems": "flex-start"},
            ),
        ],
        className="tune-header-controls",
        style={
            "display": "grid",
            "gridTemplateColumns": "1fr 1px 1fr 1px 1fr 1px 1fr 1px 1fr",
            "gap": "16px",
            "alignItems": "stretch",
            "padding": "12px",
            "border": "1px solid #d0d0d0",
            "borderRadius": "6px",
        },
    )


def build_param_range_row(row, tunable_names):
    """Render one parameter-range row."""
    row_id = row.get("id")
    options = [{"label": name, "value": name} for name in tunable_names]
    return html.Div(
        [
            dcc.Dropdown(
                id={"type": "tune-range-param", "index": row_id},
                options=options,
                value=row.get("param", "") or None,
                placeholder="parameter",
                clearable=True,
                searchable=True,
                style={"minWidth": "170px", "flex": "2 1 170px"},
            ),
            dcc.Input(
                id={"type": "tune-range-min", "index": row_id},
                type="text",
                value=row.get("min", ""),
                placeholder="min",
                style={"width": "90px", "flex": "0 0 90px"},
            ),
            dcc.Input(
                id={"type": "tune-range-max", "index": row_id},
                type="text",
                value=row.get("max", ""),
                placeholder="max",
                style={"width": "90px", "flex": "0 0 90px"},
            ),
            html.Button(
                "Remove",
                id={"type": "tune-range-remove", "index": row_id},
                n_clicks=0,
                style=action_button_style("#6b7280"),
            ),
        ],
        className="tune-range-row",
        style={"display": "flex", "flexWrap": "wrap", "gap": "8px", "alignItems": "center", "marginBottom": "8px"},
    )


def build_case_config_row(row, case_names):
    """Render one editable per-case tuning setup section."""
    row_id = row.get("id")
    options = case_names if case_names and isinstance(case_names[0], dict) else [{"label": name, "value": name} for name in case_names]
    return html.Div(
        [
            html.Div(
                [
                    dcc.Dropdown(
                        id={"type": "tune-case-name", "index": row_id},
                        options=options,
                        value=row.get("case_name") or None,
                        placeholder="case",
                        clearable=True,
                        searchable=True,
                        style={"minWidth": "170px", "flex": "1 1 170px"},
                    ),
                    html.Button(
                        "Remove",
                        id={"type": "tune-case-remove", "index": row_id},
                        n_clicks=0,
                        style=action_button_style("#6b7280"),
                    ),
                ],
                style={"display": "flex", "gap": "8px", "alignItems": "center", "marginBottom": "8px"},
            ),
            html.Div(
                [
                    html.Div(
                        [
                            html.Label("Start", className="tune-case-input-label"),
                            dcc.Input(
                                id={"type": "tune-case-time-start", "index": row_id},
                                type="number",
                                value=row.get("time_start", ""),
                                placeholder="seconds",
                                debounce=True,
                                style={"width": "112px"},
                            ),
                        ],
                        className="tune-case-input-group",
                    ),
                    html.Div(
                        [
                            html.Label("End", className="tune-case-input-label"),
                            dcc.Input(
                                id={"type": "tune-case-time-end", "index": row_id},
                                type="number",
                                value=row.get("time_end", ""),
                                placeholder="seconds",
                                debounce=True,
                                style={"width": "112px"},
                            ),
                        ],
                        className="tune-case-input-group",
                    ),
                    html.Div(
                        [
                            html.Label(
                                "Average length",
                                className="tune-case-input-label",
                            ),
                            dcc.Input(
                                id={"type": "tune-case-average-time", "index": row_id},
                                type="number",
                                min=1,
                                step=1,
                                value=row.get("average_time_seconds", ""),
                                placeholder="seconds",
                                debounce=True,
                                className="tune-average-time-input",
                                style={"width": "112px"},
                            ),
                        ],
                        className="tune-case-input-group",
                    ),
                ],
                style={"display": "flex", "flexWrap": "wrap", "gap": "8px"},
            ),
        ],
        className="tune-case-config-row",
        style={
            "border": "1px solid #d1d5db",
            "borderRadius": "6px",
            "padding": "10px",
            "marginBottom": "10px",
        },
    )


def build_results_placeholder():
    """Render the empty-state message shown before any tuning results exist."""
    return html.Div("No tuning results yet.", className="tune-empty-message", style={"fontStyle": "italic", "padding": "8px 0"})


def build_layout(initial_data):
    """Assemble the full static tuning-tab layout."""
    initial_param_rows = list(initial_data["initial_param_rows"])
    initial_case_rows = list(initial_data["initial_case_rows"])
    return html.Div(
        [
            build_top_controls(initial_data),
            dcc.Store(id="tune-case-data", data=initial_data["case_data"]),
            dcc.Store(id="tune-tunable-names", data=initial_data["tunable_names"]),
            dcc.Store(id="tune-tunable-default-ranges", data=initial_data["tunable_default_ranges"]),
            dcc.Store(id="tune-strategy-mode", data=initial_data["strategy_mode"]),
            dcc.Store(id="tune-loss-mode", data=initial_data["loss_mode"]),
            dcc.Store(id="tune-aggregation-mode", data=initial_data["aggregation_mode"]),
            dcc.Store(id="tune-case-next-id", data=len(initial_case_rows)),
            dcc.Store(id="tune-case-row-order", data=[row["id"] for row in initial_case_rows]),
            dcc.Store(id="tune-range-next-id", data=len(initial_param_rows)),
            dcc.Store(id="tune-range-row-order", data=[row["id"] for row in initial_param_rows]),
            dcc.Store(id="tune-active-job", data={}),
            dcc.Store(id="tune-status", data=initial_data["status"]),
            dcc.Store(id="tune-top-results", data=[]),
            dcc.Store(id="tune-best-results", data=[]),
            dcc.Store(id="tune-best-results-by-case", data={}),
            dcc.Store(id="tune-diagnostics-signature", data=""),
            dcc.Store(id="tune-loss-runs", data={}),
            dcc.Interval(id="tune-interval", interval=500, disabled=True),
            dcc.Interval(id="tune-loss-run-interval", interval=1000, disabled=True),
            html.Div(
                [
                    html.Div(
                        [
                            html.Div(id="tune-results-summary", className="tune-status-card", children=initial_data["status_text"]),
                            html.Div(
                                [
                                    html.Div(
                                        dcc.Graph(
                                            id="tune-taylor-diagram",
                                            className="tune-taylor-graph",
                                            config={"responsive": True, "displaylogo": False},
                                            style={"width": "100%", "minWidth": 0, "height": "430px"},
                                        ),
                                        className="tune-diagnostics-taylor",
                                    ),
                                    html.Div(
                                        [
                                            dcc.Checklist(
                                                id="tune-parameter-box-groups",
                                                options=[{"label": "Aggregate", "value": "aggregate"}],
                                                value=["aggregate"],
                                                inline=True,
                                                inputStyle={"marginRight": "4px"},
                                                labelStyle={"marginRight": "14px", "fontSize": "13px", "fontWeight": "600"},
                                                style={"padding": "0 8px 4px 8px"},
                                            ),
                                            dcc.Graph(
                                                id="tune-parameter-box-plot",
                                                className="tune-parameter-box-graph",
                                                config={"responsive": True, "displaylogo": False},
                                                style={"width": "100%", "minWidth": 0, "height": "398px"},
                                            ),
                                        ],
                                        className="tune-diagnostics-params",
                                    ),
                                ],
                                id="tune-taylor-container",
                                className="tune-taylor-container tune-diagnostics-row",
                            ),
                            html.Div(id="tune-results-container", className="tune-results-container", children=build_results_placeholder()),
                        ],
                        id="tune-left-pane",
                        className="tune-left-pane",
                        style={"height": "calc(100vh - 158px)", "minHeight": 0, "overflowY": "auto", "overflowX": "auto"},
                    ),
                    html.Div(id="tune-pane-divider", className="tune-pane-divider"),
                    html.Div(
                        [
                            html.H4("Tuning", className="tune-settings-heading"),
                            html.Div("Cases", className="tune-section-title"),
                            html.Div(
                                [
                                    build_case_config_row(row, initial_data["cases"])
                                    for row in initial_case_rows
                                ],
                                id="tune-case-rows",
                            ),
                            html.Div(
                                [
                                    html.Button(
                                        "Add case",
                                        id="tune-case-add",
                                        n_clicks=0,
                                        style=action_button_style("#111827"),
                                    )
                                ],
                                style={"display": "flex", "gap": "8px", "marginBottom": "8px"},
                            ),
                            html.Div("Fields", className="tune-section-title", style={"marginTop": "14px"}),
                            dcc.Dropdown(
                                id="tune-field-selector",
                                options=[{"label": field_name, "value": field_name} for field_name in initial_data["field_options"]],
                                value=initial_data["selected_fields"],
                                multi=True,
                                clearable=False,
                                searchable=True,
                            ),
                            html.Div("Batch Size", className="tune-section-title", style={"marginTop": "14px"}),
                            dcc.Input(
                                id="tune-batch-size",
                                type="number",
                                min=1,
                                step=1,
                                value=initial_data["batch_size"],
                                debounce=True,
                                style={"width": "120px"},
                            ),
                            html.Div("Max Workers", className="tune-section-title", style={"marginTop": "14px"}),
                            dcc.Input(
                                id="tune-max-workers",
                                type="number",
                                min=1,
                                step=1,
                                value=initial_data["max_workers"],
                                debounce=True,
                                style={"width": "120px"},
                            ),
                            html.H4("Parameter Ranges", className="tune-settings-heading"),
                            html.Div(
                                "Search modes vary only the parameters listed below across the supplied ranges.",
                                style={"marginBottom": "8px", "opacity": "0.85"},
                            ),
                            html.Div(
                                [build_param_range_row(row, initial_data["tunable_names"]) for row in initial_param_rows],
                                id="tune-range-rows",
                            ),
                            html.Div(
                                [
                                    html.Button(
                                        "Add parameter",
                                        id="tune-range-add",
                                        n_clicks=0,
                                        style=action_button_style("#111827"),
                                    )
                                ],
                                style={"display": "flex", "gap": "8px", "marginBottom": "8px"},
                            ),
                            html.Div(id="tune-loss-run-message", className="tune-info-message", style={"marginTop": "10px"}),
                            html.Div(id="tune-validation-message", className="tune-validation-message", style={"marginTop": "10px"}),
                        ],
                        id="tune-right-pane",
                        className="tune-right-pane",
                        style={"paddingLeft": "16px", "height": "calc(100vh - 158px)", "minHeight": 0, "overflowY": "auto", "overflowX": "auto"},
                    ),
                ],
                id="tune-tab-layout",
                className="tune-tab-layout",
                style={
                    "display": "grid",
                    "gridTemplateColumns": f"minmax(0,1fr) 8px {initial_data['right_pane_width_px']}px",
                    "gap": "16px",
                    "padding": "10px",
                    "marginTop": "16px",
                    "overflowX": "auto",
                },
            ),
        ]
    )
