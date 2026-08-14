"""Static layout builders for the tuning tab."""

from __future__ import annotations

from dash import dcc, html

from dash_app.shared import styled_dropdown
from dash_app.shared.selected_build import selected_build_badge
from tuner.taylor_metrics import (
    DEFAULT_AGGREGATION_WEIGHTS,
    DEFAULT_LOSS_MODE,
    DEFAULT_TIME_WINDOW_AGGREGATION_SCOPE,
)
from tuner.presets import list_presets

DEFAULT_AVERAGE_TIME_SECONDS = 3600


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


def config_button_style(selected=False, disabled=False):
    """Return the style for tunable-config selector buttons."""
    style = mode_button_style(selected=selected, disabled=disabled)
    style.update({"minWidth": "116px"})
    return style


def mode_options_block_style(visible: bool) -> dict:
    """Return show/hide style for a mode-options block."""
    return {"display": "block" if visible else "none"}


def build_top_controls(initial_data):
    """Render mode selection, mode options, and run controls."""
    selected_mode = initial_data.get("strategy_mode")
    selected_loss_mode = initial_data.get("loss_mode", DEFAULT_LOSS_MODE)
    selected_aggregation_scope = initial_data.get(
        "time_window_aggregation_scope", DEFAULT_TIME_WINDOW_AGGREGATION_SCOPE
    )
    aggregation_weights = initial_data.get("aggregation_weights", DEFAULT_AGGREGATION_WEIGHTS)
    return html.Div(
        [
            html.Div(
                [
                    dcc.Input(
                        id="tune-workspace-name",
                        type="text",
                        value="",
                        debounce=True,
                        placeholder="Workspace name",
                        className="clubb-input tune-workspace-name-input",
                        style={
                            "width": "100%",
                            "minWidth": "260px",
                            "fontSize": "20px",
                            "fontWeight": "650",
                            "padding": "10px 12px",
                            "boxSizing": "border-box",
                        },
                    ),
                    html.Div(
                        [
                            html.Button(
                                "Workspaces",
                                id="tune-workspace-change",
                                n_clicks=0,
                                className="tune-workspace-change",
                                style=action_button_style("#374151"),
                            ),
                            dcc.Dropdown(
                                id="tune-workspace-revision-selector",
                                options=[],
                                value=None,
                                clearable=False,
                                searchable=False,
                                placeholder="Revision",
                                style={"minWidth": "155px", "flex": "1 1 155px"},
                            ),
                        ],
                        style={"display": "flex", "gap": "8px", "alignItems": "center", "flexWrap": "wrap", "marginTop": "10px"},
                    ),
                    html.Div(
                        html.Button(
                            "New revision",
                            id="tune-workspace-new-revision",
                            n_clicks=0,
                            style=action_button_style("#7c3aed", disabled=True),
                            disabled=True,
                        ),
                        style={"display": "flex", "gap": "8px", "flexWrap": "wrap", "marginTop": "6px"},
                    ),
                    html.Div(
                        id="tune-workspace-summary",
                        className="tune-info-message",
                        style={"marginTop": "8px", "minHeight": "18px", "fontSize": "12px"},
                    ),
                ],
                className="tune-workspace-controls",
                style={"padding": "12px", "display": "flex", "flexDirection": "column", "alignItems": "stretch"},
            ),
            html.Div(style={"width": "1px", "backgroundColor": "#d0d0d0", "alignSelf": "stretch"}),
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
                            html.Button(
                                "SimAnn",
                                id="tune-mode-simann",
                                n_clicks=0,
                                style=mode_button_style(selected=selected_mode == "simann"),
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
                                type="text",
                                value=initial_data["random_max_samples"],
                                debounce=True,
                                className="clubb-input",
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
                                type="text",
                                value=initial_data["resolve_spacing"],
                                debounce=True,
                                className="clubb-input",
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
                    html.Div(
                        [
                            html.Label("Max Iterations", htmlFor="tune-simann-max-iters"),
                            dcc.Input(
                                id="tune-simann-max-iters",
                                type="number",
                                min=1,
                                step=1,
                                value=initial_data["simann_max_iters"],
                                debounce=True,
                                style={"width": "120px", "marginLeft": "8px", "marginRight": "12px"},
                            ),
                            html.Label("Initial Temp", htmlFor="tune-simann-initial-temp"),
                            dcc.Input(
                                id="tune-simann-initial-temp",
                                type="number",
                                min=0,
                                step="any",
                                value=initial_data["simann_initial_temp"],
                                debounce=True,
                                style={"width": "100px", "marginLeft": "8px", "marginRight": "12px"},
                            ),
                            html.Label("Final Temp", htmlFor="tune-simann-final-temp"),
                            dcc.Input(
                                id="tune-simann-final-temp",
                                type="number",
                                min=0,
                                step="any",
                                value=initial_data["simann_final_temp"],
                                debounce=True,
                                style={"width": "110px", "marginLeft": "8px"},
                            ),
                        ],
                        id="tune-simann-options",
                        style={**mode_options_block_style(selected_mode == "simann"), "marginTop": "10px"},
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
                                title="L = scaled_rmse.",
                                style=mode_button_style(selected=selected_loss_mode == "scaled_rmse"),
                            ),
                            html.Button(
                                "Centered RMS + Bias",
                                id="tune-loss-mode-centered-rmse-bias",
                                n_clicks=0,
                                title="L = centered_rmse_norm + |bias_norm|.",
                                style=mode_button_style(selected=selected_loss_mode == "centered_rmse_bias"),
                            ),
                            html.Button(
                                "Taylor Components",
                                id="tune-loss-mode-taylor-components",
                                n_clicks=0,
                                title="L = centered_rmse_norm + 0.5|bias_norm| + 0.25(1-R) + 0.25|log(std_ratio)|.",
                                style=mode_button_style(selected=selected_loss_mode == "taylor_components"),
                            ),
                            html.Button(
                                "Squared",
                                id="tune-loss-mode-taylor-components-squared",
                                n_clicks=0,
                                title="L = centered_rmse_norm² + 0.5|bias_norm|² + 0.25(1-R)² + 0.25|log(std_ratio)|².",
                                style=mode_button_style(selected=selected_loss_mode == "taylor_components_squared"),
                            ),
                            html.Button(
                                "Shape First",
                                id="tune-loss-mode-shape-first",
                                n_clicks=0,
                                title="L = 0.5 centered_rmse_norm + 0.5|bias_norm| + (1-R) + 0.5|log(std_ratio)|.",
                                style=mode_button_style(selected=selected_loss_mode == "shape_first"),
                            ),
                            html.Button(
                                "Bias Light",
                                id="tune-loss-mode-bias-light-taylor",
                                n_clicks=0,
                                title="L = centered_rmse_norm + 0.25|bias_norm| + 0.5(1-R) + 0.5|log(std_ratio)|.",
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
                    html.Div("Time-window aggregation", className="tune-section-title"),
                    html.Div(
                        [
                            html.Button(
                                "Overall",
                                id="tune-aggregation-overall",
                                n_clicks=0,
                                title="Pool all active case/field/time-window losses, then rank them best to worst.",
                                style=mode_button_style(selected=selected_aggregation_scope == "overall"),
                            ),
                            html.Button(
                                "By case",
                                id="tune-aggregation-by-case",
                                n_clicks=0,
                                title="Rank each case's active time-window losses separately, then take the case-weighted mean.",
                                style=mode_button_style(selected=selected_aggregation_scope == "by_case"),
                            ),
                        ],
                        style={"display": "flex", "gap": "8px", "flexWrap": "wrap"},
                    ),
                    html.Div(
                        [
                            html.Span("Best → worst weights", style={"fontSize": "12px", "fontWeight": "600"}),
                            *[
                                dcc.Input(
                                    id=f"tune-aggregation-weight-{index + 1}",
                                    type="number",
                                    min=0,
                                    step="any",
                                    value=aggregation_weights[index] if index < len(aggregation_weights) else DEFAULT_AGGREGATION_WEIGHTS[index],
                                    debounce=True,
                                    style={"width": "68px"},
                                )
                                for index, label in enumerate(("Best 25%", "Lower-middle 25%", "Upper-middle 25%", "Worst 25%"))
                            ],
                        ],
                        style={"display": "flex", "alignItems": "center", "gap": "6px", "marginTop": "10px", "flexWrap": "wrap"},
                    ),
                ],
                style={"padding": "12px"},
            ),
            html.Div(style={"width": "1px", "backgroundColor": "#d0d0d0", "alignSelf": "stretch"}),
            html.Div(
                [
                    html.Div("Run Tuner", className="tune-section-title"),
                    html.Div(
                        [
                            html.Div(
                                [
                                    html.Label("Batch Size", htmlFor="tune-batch-size"),
                                    dcc.Input(
                                        id="tune-batch-size",
                                        type="number",
                                        min=1,
                                        step=1,
                                        value=initial_data["batch_size"],
                                        debounce=True,
                                        style={"width": "96px"},
                                    ),
                                ],
                                style={"display": "flex", "flexDirection": "column", "gap": "4px"},
                            ),
                            html.Div(
                                [
                                    html.Label("Max Workers", htmlFor="tune-max-workers"),
                                    dcc.Input(
                                        id="tune-max-workers",
                                        type="number",
                                        min=1,
                                        step=1,
                                        value=initial_data["max_workers"],
                                        debounce=True,
                                        style={"width": "96px"},
                                    ),
                                ],
                                style={"display": "flex", "flexDirection": "column", "gap": "4px"},
                            ),
                        ],
                        style={"display": "flex", "gap": "10px", "flexWrap": "wrap", "marginBottom": "8px"},
                    ),
                    html.Div(
                        [
                            selected_build_badge("tune-selected-build-badge"),
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
                            html.Button(
                                "Reset",
                                id="tune-reset-button",
                                n_clicks=0,
                                style=action_button_style("#b45309", disabled=True),
                                disabled=True,
                            ),
                        ],
                        style={"display": "flex", "gap": "8px", "flexWrap": "wrap"},
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
            "gridTemplateColumns": (
                "minmax(300px, 1.35fr) 1px minmax(250px, 1fr) 1px "
                "minmax(310px, 1.25fr) 1px minmax(205px, 0.8fr) 1px "
                "minmax(190px, 0.75fr) 1px minmax(190px, 0.75fr)"
            ),
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
    targets = [str(target).strip() for target in row.get("targets", [row.get("param", "")]) if str(target).strip()]
    if not targets:
        targets = [str(row.get("param", "")).strip()]
    linked_label = " = ".join(targets) if len(targets) > 1 else ""
    options = [{"label": name, "value": name} for name in tunable_names]
    return html.Div(
        [
            dcc.Store(id={"type": "tune-range-targets", "index": row_id}, data=targets),
            dcc.Dropdown(
                id={"type": "tune-range-param", "index": row_id},
                options=options,
                value=row.get("param", "") or None,
                placeholder="parameter",
                clearable=True,
                searchable=True,
                className="clubb-dropdown",
                style={"minWidth": "170px", "flex": "2 1 170px"},
            ),
            html.Span(
                linked_label,
                id={"type": "tune-range-link-label", "index": row_id},
                className="tune-range-link-label",
                title="This one sampled value is applied to each linked physical CLUBB parameter.",
                style={"fontSize": "12px", "opacity": 0.8, "whiteSpace": "nowrap", "display": "inline-block" if linked_label else "none"},
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


def build_preset_buttons():
    """Return compact preset actions from the shared tuner definition file."""
    buttons = []
    for preset in list_presets():
        buttons.append(
            html.Button(
                preset["label"],
                id={"type": "tune-preset", "name": preset["name"]},
                n_clicks=0,
                title=preset["description"],
                style=action_button_style("#0f766e"),
            )
        )
    return html.Div(
        buttons,
        id="tune-preset-buttons",
        style={"display": "flex", "flexWrap": "wrap", "gap": "4px", "marginBottom": "6px"},
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
                        className="clubb-dropdown",
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
                                type="text",
                                value=row.get("time_start", ""),
                                placeholder="seconds",
                                debounce=True,
                                className="clubb-input",
                                style={"width": "68px"},
                            ),
                        ],
                        className="tune-case-input-group",
                    ),
                    html.Div(
                        [
                            html.Label("End", className="tune-case-input-label"),
                            dcc.Input(
                                id={"type": "tune-case-time-end", "index": row_id},
                                type="text",
                                value=row.get("time_end", ""),
                                placeholder="seconds",
                                debounce=True,
                                className="clubb-input",
                                style={"width": "68px"},
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
                                type="text",
                                value=row.get("average_time_seconds", DEFAULT_AVERAGE_TIME_SECONDS),
                                placeholder="seconds",
                                debounce=True,
                                className="clubb-input tune-average-time-input",
                                style={"width": "68px"},
                            ),
                        ],
                        className="tune-case-input-group",
                    ),
                    html.Div(
                        [
                            html.Label("Altitude min", className="tune-case-input-label"),
                            dcc.Input(
                                id={"type": "tune-case-altitude-min", "index": row_id},
                                type="text",
                                value=row.get("altitude_min", ""),
                                placeholder="m",
                                debounce=True,
                                className="clubb-input",
                                style={"width": "68px"},
                            ),
                        ],
                        className="tune-case-input-group",
                    ),
                    html.Div(
                        [
                            html.Label("Altitude max", className="tune-case-input-label"),
                            dcc.Input(
                                id={"type": "tune-case-altitude-max", "index": row_id},
                                type="text",
                                value=row.get("altitude_max", ""),
                                placeholder="m",
                                debounce=True,
                                className="clubb-input",
                                style={"width": "68px"},
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


def build_config_buttons(configs, selected_config):
    """Render selectable tunable config buttons."""
    buttons = []
    for config in configs or []:
        value = str(config.get("value", "")).strip()
        if not value:
            continue
        label = str(config.get("label") or value)
        buttons.append(
            html.Button(
                label,
                id={"type": "tune-config-button", "name": value},
                n_clicks=0,
                title=f"Use input/parameter_and_flag_configs/{value}",
                style=config_button_style(selected=value == selected_config),
            )
        )
    if buttons:
        return html.Div(buttons, style={"display": "flex", "gap": "8px", "flexWrap": "wrap"})
    return html.Div("No complete tunable configs found.", className="tune-empty-message")


def build_results_placeholder():
    """Render the empty-state message shown before any tuning results exist."""
    return html.Div("No tuning results yet.", className="tune-empty-message", style={"fontStyle": "italic", "padding": "8px 0"})


def landscape_control(label, control):
    """Wrap one Landscape selector in the compact labelled-control pattern."""
    return html.Div(
        [html.Span(label, className="tune-landscape-control-label"), control],
        className="tune-landscape-control-group",
    )


def build_layout(initial_data):
    """Assemble the full static tuning-tab layout."""
    initial_param_rows = list(initial_data["initial_param_rows"])
    initial_case_rows = list(initial_data["initial_case_rows"])
    return html.Div(
        [
            html.Div(
                build_top_controls(initial_data),
                id="tune-top-controls",
                style={"opacity": 0.45, "pointerEvents": "none"},
            ),
            dcc.Store(id="tune-case-data", data=initial_data["case_data"]),
            dcc.Store(id="tune-tunable-configs", data=initial_data["tunable_configs"]),
            dcc.Store(id="tune-selected-config", data=initial_data["selected_config"]),
            dcc.Store(id="tune-tunable-names", data=initial_data["tunable_names"]),
            dcc.Store(id="tune-settings-resolution", data={}),
            dcc.Store(id="tune-tunable-default-ranges", data=initial_data["tunable_default_ranges"]),
            dcc.Store(id="tune-strategy-mode", data=initial_data["strategy_mode"]),
            dcc.Store(id="tune-loss-mode", data=initial_data["loss_mode"]),
            dcc.Store(id="tune-aggregation-mode", data=initial_data["aggregation_mode"]),
            dcc.Store(id="tune-time-window-aggregation-scope", data=initial_data.get("time_window_aggregation_scope", DEFAULT_TIME_WINDOW_AGGREGATION_SCOPE)),
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
            dcc.Store(id="tune-landscape-signature", data=""),
            dcc.Store(id="tune-runtime-loss-signature", data=""),
            dcc.Store(id="tune-loss-runs", data={}),
            dcc.Store(id="tune-workspace-list", data=[]),
            dcc.Store(id="tune-workspace-activity", data=[]),
            dcc.Store(id="tune-workspace-selection", data={}, storage_type="local"),
            dcc.Store(id="tune-workspace-refresh", data=0),
            dcc.Store(id="tune-workspace-delete-target", data=None),
            dcc.Store(id="tune-reset-prompt", data=False),
            dcc.ConfirmDialog(id="tune-workspace-delete-confirm", message="Delete this saved Tune workspace and all of its revisions?"),
            html.Div(
                html.Div(
                    [
                        html.Div("Reset this started revision?", style={"fontWeight": "700", "marginBottom": "6px"}),
                        html.Div("Choose whether to erase this revision's saved data or preserve it and branch a fresh editable revision.", className="tune-info-message"),
                        html.Div(
                            [
                                html.Button("Cancel", id="tune-reset-cancel", n_clicks=0, style=action_button_style("#475569")),
                                html.Button("Delete data", id="tune-reset-delete-data", n_clicks=0, style=action_button_style("#b91c1c")),
                                html.Button("New revision", id="tune-reset-new-revision", n_clicks=0, style=action_button_style("#7c3aed")),
                            ],
                            style={"display": "flex", "gap": "8px", "justifyContent": "flex-end"},
                        ),
                    ],
                    className="tune-status-card",
                    style={"maxWidth": "520px", "margin": "10vh auto", "padding": "18px"},
                ),
                id="tune-reset-dialog",
                style={"display": "none", "position": "fixed", "zIndex": 2500, "inset": 0, "backgroundColor": "rgba(0, 0, 0, 0.55)", "padding": "24px"},
            ),
            # A one-second cadence keeps the live status responsive while
            # avoiding two full multi-graph Dash render passes every second.
            dcc.Interval(id="tune-interval", interval=1000, disabled=True),
            # Workspace status is deliberately lightweight and independent of
            # the selected revision.  A Tune worker may continue in another
            # workspace while this tab is inspecting an inactive one.
            dcc.Interval(id="tune-workspace-status-interval", interval=5000, n_intervals=0),
            dcc.Interval(id="tune-loss-run-interval", interval=1000, disabled=True),
            html.Div(
                [
                    html.Div(
                        [
                            html.Div(
                                [
                                    html.Div(
                                        [
                                            html.H4("Saved Tune workspaces", className="tune-settings-heading", style={"margin": "0 0 10px"}),
                                        ],
                                        style={"display": "flex", "alignItems": "center"},
                                    ),
                                    html.Div(
                                        html.Div(
                                            html.Button("New workspace", id="tune-workspace-new", n_clicks=0, style=action_button_style("#2563eb")),
                                            className="tune-status-card",
                                            style={"display": "flex", "alignItems": "center", "gap": "10px", "marginBottom": "8px"},
                                        ),
                                        id="tune-workspace-picker-list",
                                    ),
                                ],
                                id="tune-workspace-picker",
                                className="tune-workspace-picker",
                            ),
                            html.Div(
                                [
                            html.Div(
                                [
                                    html.Div(
                                        [
                                            html.Div(
                                                [
                                                    html.Span("Taylor subwindow display", className="tune-diagnostics-label"),
                                                    dcc.RadioItems(
                                                        id="tune-taylor-window-display",
                                                        options=[
                                                            {"label": "Worst", "value": "worst"},
                                                            {"label": "Average", "value": "average"},
                                                            {"label": "Best", "value": "best"},
                                                        ],
                                                        value="average",
                                                        inline=True,
                                                        className="tune-segmented-control",
                                                    ),
                                                ],
                                                className="tune-diagnostics-control-row",
                                            ),
                                            dcc.Graph(
                                                id="tune-taylor-diagram",
                                                className="tune-taylor-graph",
                                                config={"responsive": True, "displaylogo": False},
                                                style={"width": "100%", "minWidth": 0, "height": "430px"},
                                            ),
                                        ],
                                        className="tune-diagnostics-taylor",
                                    ),
                                    html.Div(
                                        [
                                            dcc.RadioItems(
                                                id="tune-parameter-box-groups",
                                                className="tune-segmented-control tune-parameter-box-groups",
                                                options=[
                                                    {"label": "Aggregate", "value": "aggregate"},
                                                    {"label": "All", "value": "all"},
                                                ],
                                                value="aggregate",
                                                inline=True,
                                            ),
                                            dcc.Graph(
                                                id="tune-parameter-box-plot",
                                                className="tune-parameter-box-graph",
                                                config={"responsive": True, "displaylogo": False},
                                                style={"width": "100%", "minWidth": 0, "height": "430px"},
                                            ),
                                        ],
                                        className="tune-diagnostics-params",
                                    ),
                                ],
                                id="tune-taylor-container",
                                className="tune-taylor-container tune-diagnostics-row",
                            ),
                            html.Div(
                                [
                                    html.Div(
                                        [
                                            html.Div(
                                                [
                                                    html.Div("Explore parameter response", className="tune-landscape-heading"),
                                                    html.Div(
                                                        "Each sample is one tested parameter set. Color shows the selected loss or diagnostic; choose a field to expose its response and pairwise effects.",
                                                        className="tune-landscape-caption",
                                                    ),
                                                ],
                                                className="tune-landscape-intro",
                                            ),
                                            landscape_control("X parameter", styled_dropdown(
                                                id="tune-landscape-x-param",
                                                options=[],
                                                value=None,
                                                clearable=False,
                                                searchable=True,
                                                placeholder="x parameter",
                                                className="tune-landscape-control",
                                            )),
                                            landscape_control("Y parameter", styled_dropdown(
                                                id="tune-landscape-y-param",
                                                options=[],
                                                value=None,
                                                clearable=False,
                                                searchable=True,
                                                placeholder="y parameter",
                                                className="tune-landscape-control",
                                            )),
                                            landscape_control("Metric", styled_dropdown(
                                                id="tune-landscape-metric",
                                                options=[],
                                                value=None,
                                                clearable=False,
                                                searchable=True,
                                                placeholder="metric",
                                                className="tune-landscape-control",
                                            )),
                                            landscape_control("Summary", styled_dropdown(
                                                id="tune-landscape-aggregation",
                                                options=[
                                                    {"label": "Mean", "value": "mean"},
                                                    {"label": "Median", "value": "median"},
                                                    {"label": "Min", "value": "min"},
                                                    {"label": "P90", "value": "p90"},
                                                    {"label": "Max", "value": "max"},
                                                ],
                                                value="mean",
                                                clearable=False,
                                                searchable=False,
                                                className="tune-landscape-control",
                                            )),
                                            landscape_control("Display", styled_dropdown(
                                                id="tune-landscape-mode",
                                                options=[
                                                    {"label": "Samples", "value": "samples"},
                                                    {"label": "Binned heatmap", "value": "binned"},
                                                ],
                                                value="samples",
                                                clearable=False,
                                                searchable=False,
                                                className="tune-landscape-control",
                                            )),
                                            landscape_control("Case", styled_dropdown(
                                                id="tune-landscape-case",
                                                options=[{"label": "All cases", "value": "__all__"}],
                                                value=None,
                                                clearable=False,
                                                searchable=True,
                                                className="tune-landscape-control",
                                            )),
                                            landscape_control("Field", styled_dropdown(
                                                id="tune-landscape-field",
                                                options=[{"label": "All fields", "value": "__all__"}],
                                                value=None,
                                                clearable=False,
                                                searchable=True,
                                                className="tune-landscape-control",
                                            )),
                                            landscape_control("Time window", styled_dropdown(
                                                id="tune-landscape-window",
                                                options=[{"label": "All windows", "value": "__all__"}],
                                                value="__all__",
                                                clearable=False,
                                                searchable=False,
                                                className="tune-landscape-control",
                                            )),
                                        ],
                                        className="tune-landscape-toolbar",
                                    ),
                                    html.Div(
                                        [
                                            dcc.Graph(
                                                id="tune-landscape-plot",
                                                className="tune-landscape-graph",
                                                config={"responsive": True, "displaylogo": False},
                                                style={"width": "100%", "minWidth": 0, "height": "430px"},
                                            ),
                                            dcc.Graph(
                                                id="tune-parameter-correlation-plot",
                                                className="tune-landscape-graph",
                                                config={"responsive": True, "displaylogo": False},
                                                style={"width": "100%", "minWidth": 0, "height": "430px"},
                                            ),
                                        ],
                                        className="tune-landscape-row",
                                    ),
                                    html.Div(
                                        [
                                            dcc.Graph(
                                                id="tune-field-sensitivity-plot",
                                                className="tune-landscape-graph",
                                                config={"responsive": True, "displaylogo": False},
                                                style={"width": "100%", "minWidth": 0, "height": "430px"},
                                            ),
                                            dcc.Graph(
                                                id="tune-field-interaction-plot",
                                                className="tune-landscape-graph",
                                                config={"responsive": True, "displaylogo": False},
                                                style={"width": "100%", "minWidth": 0, "height": "430px"},
                                            ),
                                        ],
                                        className="tune-field-sensitivity-row",
                                        style={
                                            "display": "grid",
                                            "gridTemplateColumns": "minmax(0, 1fr) minmax(0, 1fr)",
                                            "gap": "12px",
                                            "alignItems": "stretch",
                                            "marginTop": "12px",
                                        },
                                    ),
                                ],
                                id="tune-landscape-container",
                                className="tune-landscape-container",
                            ),
                            html.Div(id="tune-results-container", className="tune-results-container", children=build_results_placeholder()),
                                ],
                                id="tune-results-view",
                                style={"display": "none"},
                            ),
                        ],
                        id="tune-left-pane",
                        className="tune-left-pane",
                        style={"height": "calc(100vh - 158px)", "minHeight": 0, "overflowY": "auto", "overflowX": "auto"},
                    ),
                    html.Div(id="tune-pane-divider", className="tune-pane-divider"),
                    html.Div(
                        [
                            dcc.Tabs(
                                id="tune-right-tabs",
                                value="configure",
                                className="tune-right-tabs",
                                children=[
                                    dcc.Tab(
                                        label="Configure",
                                        value="configure",
                                        className="tune-right-tab",
                                        selected_className="tune-right-tab--selected",
                                        children=html.Div(
                                            [
                                                html.H4("Tuning", className="tune-settings-heading"),
                                                html.Div("Presets", className="tune-section-title"),
                                                build_preset_buttons(),
                                                html.Div(
                                                    "A preset replaces this editable draft with its cases, fields, parameter coordinates, and any required override. Ranges come from Fortran hard bounds plus config defaults.",
                                                    className="tune-info-message",
                                                    style={"marginBottom": "12px"},
                                                ),
                                                html.Div("Cases", className="tune-section-title"),
                                                html.Div(
                                                    [build_case_config_row(row, initial_data["cases"]) for row in initial_case_rows],
                                                    id="tune-case-rows",
                                                ),
                                                html.Div(
                                                    html.Button("Add case", id="tune-case-add", n_clicks=0, style=action_button_style("#111827")),
                                                    style={"display": "flex", "gap": "8px", "marginBottom": "8px"},
                                                ),
                                                html.Div("Fields", className="tune-section-title", style={"marginTop": "14px"}),
                                                dcc.Dropdown(
                                                    id="tune-field-selector",
                                                    options=[{"label": field_name, "value": field_name} for field_name in initial_data["field_options"]],
                                                    value=initial_data["selected_fields"], multi=True, clearable=False, searchable=True,
                                                    className="clubb-dropdown",
                                                ),
                                                html.H4("Parameter Ranges", className="tune-settings-heading"),
                                                html.Div("Search modes vary only the parameters listed below across the supplied ranges.", style={"marginBottom": "8px", "opacity": "0.85"}),
                                                html.Div(
                                                    [build_param_range_row(row, initial_data["tunable_names"]) for row in initial_param_rows],
                                                    id="tune-range-rows",
                                                ),
                                                html.Div(
                                                    html.Button("Add parameter", id="tune-range-add", n_clicks=0, style=action_button_style("#111827")),
                                                    style={"display": "flex", "gap": "8px", "marginBottom": "8px"},
                                                ),
                                                html.H4("Config", className="tune-settings-heading"),
                                                html.Div(build_config_buttons(initial_data["tunable_configs"], initial_data["selected_config"]), id="tune-config-buttons", style={"marginBottom": "8px"}),
                                                html.Label("SCM override", htmlFor="tune-scm-override"),
                                                dcc.Input(id="tune-scm-override", type="text", value=initial_data["scm_override"], placeholder="iiPDF_type=1", debounce=True, style={"width": "100%", "marginTop": "4px"}),
                                                html.Div(id="tune-settings-resolution-note", className="tune-info-message", style={"display": "none", "marginTop": "5px", "marginBottom": "5px"}),
                                                html.Div(
                                                    "Comma-separated key=value pairs, passed as -override to tuning and result reruns. A pasted '-override …' prefix is also accepted.",
                                                    className="tune-info-message", style={"marginTop": "4px", "marginBottom": "8px"},
                                                ),
                                            ],
                                            id="tune-configure-content",
                                        ),
                                    ),
                                    dcc.Tab(
                                        label="Runtime",
                                        value="runtime",
                                        className="tune-right-tab",
                                        selected_className="tune-right-tab--selected",
                                        children=html.Div(
                                            [
                                                html.Div(
                                                    "Start a Tune run to see its state and case queues.",
                                                    id="tune-runtime-status",
                                                    className="tune-status-card tune-runtime-status-card",
                                                ),
                                                dcc.Graph(
                                                    id="tune-runtime-best-loss-graph",
                                                    className="tune-runtime-loss-graph",
                                                    config={"displayModeBar": False, "responsive": True},
                                                    style={"width": "100%", "height": "310px"},
                                                ),
                                                html.Div(id="tune-runtime-error", className="tune-runtime-error"),
                                                html.Div(id="tune-loss-run-message", className="tune-runtime-message tune-info-message"),
                                                html.Div(id="tune-validation-message", className="tune-runtime-message tune-validation-message"),
                                            ],
                                            className="tune-runtime-pane",
                                        ),
                                    ),
                                ],
                            ),
                        ],
                        id="tune-right-pane",
                        className="tune-right-pane",
                        # Keep controls in document flow so their full height
                        # uses the browser scrollbar rather than a nested pane
                        # scrollbar beside the Tune controls.
                        style={"paddingLeft": "16px", "height": "auto", "minHeight": 0, "overflow": "visible", "opacity": 0.45, "pointerEvents": "none"},
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
