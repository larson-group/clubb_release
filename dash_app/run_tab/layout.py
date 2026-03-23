"""Static layout builders and UI styling helpers for the run tab."""

from dash import dcc, html

from .state import MAX_RUN_PROCS


def field_style(changed, disabled=False):
    """Return the field row style, highlighting modified or disabled values when needed."""
    base = {"padding": "2px 6px", "display": "flex", "alignItems": "center", "gap": "10px", "borderRadius": "4px"}
    if changed:
        base.update({"outline": "2px solid #f59e0b"})
    if disabled:
        base.update({"opacity": "0.5", "filter": "grayscale(0.25)"})
    return base


def case_button_style(color, selected=False):
    """Return the style for one case button, including selected outline state."""
    border_color = "#f59e0b" if selected else "transparent"
    border_width = "3px" if selected else "1px"
    style = {
        "color": "#ffffff",
        "border": f"{border_width} solid {border_color}",
        "padding": "6px 10px",
        "margin": "4px",
        "borderRadius": "4px",
        "cursor": "pointer",
        "boxSizing": "border-box",
        "opacity": "1",
        "fontSize": "16px",
    }
    style["background"] = color
    style["backgroundColor"] = "#dc2626" if "gradient" in color else color
    return style


def run_action_button_style(color, disabled=False):
    """Return the common style for run-tab action buttons."""
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


def stats_button_style(selected=False):
    """Return the style for one stats-file selection button."""
    return {
        "backgroundColor": "#0ea5e9" if selected else "#e2e8f0",
        "color": "#0f172a" if selected else "#334155",
        "border": "1px solid #94a3b8",
        "padding": "6px 10px",
        "margin": "4px",
        "borderRadius": "4px",
        "cursor": "pointer",
        "fontWeight": "600" if selected else "500",
    }


def case_dom_id(case_name):
    """Return a DOM-safe identifier suffix for a case name."""
    import re

    return re.sub(r"[^A-Za-z0-9_-]", "_", case_name)


def compute_width_hints(all_config_names):
    """Compute label and pane widths from the longest config parameter name."""
    max_name_len = max((len(name) for name in all_config_names), default=16)
    label_width_px = max(170, min(560, int(5 + max_name_len * 7.0)))
    value_width_px = 65
    right_pane_width_px = max(340, min(960, label_width_px + value_width_px + 80))
    return label_width_px, value_width_px, right_pane_width_px


def build_select_actions(case_groups):
    """Render the case selection helper buttons."""
    return html.Div(
        [
            html.Button("Deselect", id="run-deselect", n_clicks=0, style=run_action_button_style("#6b7280")),
            html.Button("Select all", id={"type": "run-group-button", "name": "all"}, n_clicks=0, disabled=not case_groups["all"], className="run-select-group-button", style=run_action_button_style("#111827", disabled=not case_groups["all"])),
            html.Button("Select standard", id={"type": "run-group-button", "name": "standard"}, n_clicks=0, disabled=not case_groups["standard"], className="run-select-group-button", style=run_action_button_style("#111827", disabled=not case_groups["standard"])),
            html.Button("Select priority", id={"type": "run-group-button", "name": "priority"}, n_clicks=0, disabled=not case_groups["priority"], className="run-select-group-button", style=run_action_button_style("#111827", disabled=not case_groups["priority"])),
            html.Button("Select minimum", id={"type": "run-group-button", "name": "minimum"}, n_clicks=0, disabled=not case_groups["minimum"], className="run-select-group-button", style=run_action_button_style("#111827", disabled=not case_groups["minimum"])),
            html.Button("Select short", id={"type": "run-group-button", "name": "short"}, n_clicks=0, disabled=not case_groups["short"], className="run-select-group-button", style=run_action_button_style("#111827", disabled=not case_groups["short"])),
        ],
        className="run-select-actions",
        style={"display": "flex", "flexWrap": "wrap", "gap": "8px", "marginBottom": "6px"},
    )


def build_case_buttons(cases):
    """Render all case buttons with the default unselected style."""
    return [
        html.Button(case_name, id={"type": "run-case-button", "name": case_name}, n_clicks=0, style=case_button_style("#2563eb", False))
        for case_name in cases
    ]


def build_stats_buttons(stats_files, default_stats_name, no_stats_name):
    """Render stats-file selection buttons, including the synthetic none option."""
    buttons = [
        html.Button(stats_name, id={"type": "run-stats-button", "name": stats_name}, n_clicks=0, style=stats_button_style(stats_name == default_stats_name))
        for stats_name in stats_files
    ]
    if no_stats_name not in stats_files:
        buttons.append(
            html.Button("none", id={"type": "run-stats-button", "name": no_stats_name}, n_clicks=0, style=stats_button_style(no_stats_name == default_stats_name))
        )
    return buttons


def build_optional_args_section():
    """Render the optional run_scm.py argument inputs."""
    return html.Div(
        [
            html.Div("Optional run args:", className="run-section-title"),
            html.Div(
                [
                    dcc.Input(id="run-opt-max-iters", type="text", value="", placeholder="max_iters", style={"width": "130px"}),
                    dcc.Input(id="run-opt-debug", type="text", value="", placeholder="debug", style={"width": "130px"}),
                    dcc.Input(id="run-opt-dt-main", type="text", value="", placeholder="dt_main", style={"width": "130px"}),
                    dcc.Input(id="run-opt-dt-rad", type="text", value="", placeholder="dt_rad", style={"width": "130px"}),
                    dcc.Input(id="run-opt-tout", type="text", value="", placeholder="tout", style={"width": "130px"}),
                ],
                style={"display": "flex", "flexWrap": "wrap", "gap": "8px", "marginTop": "4px"},
            ),
            html.Div(
                [
                    html.Label("Output dir:", style={"whiteSpace": "nowrap", "alignSelf": "center"}),
                    dcc.Input(id="run-opt-out-dir", type="text", value="output", placeholder="output", style={"width": "200px"}),
                ],
                style={"display": "flex", "gap": "8px", "marginTop": "6px"},
            ),
        ],
        style={"marginTop": "6px"},
    )


def build_multicol_row(row, tunable_names):
    """Render one multicol hypergrid specification row."""
    row_id = row.get("id")
    options = [{"label": name, "value": name} for name in tunable_names]
    return html.Div(
        [
            dcc.Dropdown(
                id={"type": "run-hr-param", "index": row_id},
                options=options,
                value=row.get("param", "") or None,
                placeholder="parameter",
                clearable=True,
                searchable=True,
                style={"minWidth": "170px", "flex": "2 1 170px"},
            ),
            dcc.Input(
                id={"type": "run-hr-min", "index": row_id},
                type="text",
                value=row.get("min", ""),
                placeholder="min",
                style={"width": "90px", "flex": "0 0 90px"},
            ),
            dcc.Input(
                id={"type": "run-hr-max", "index": row_id},
                type="text",
                value=row.get("max", ""),
                placeholder="max",
                style={"width": "90px", "flex": "0 0 90px"},
            ),
            dcc.Input(
                id={"type": "run-hr-npoints", "index": row_id},
                type="text",
                value=row.get("npoints", ""),
                placeholder="npoints",
                style={"width": "110px", "flex": "0 0 110px"},
            ),
            html.Button(
                "Remove",
                id={"type": "run-hr-remove", "index": row_id},
                n_clicks=0,
                style=run_action_button_style("#6b7280"),
            ),
        ],
        className="run-multicol-row",
        style={"display": "flex", "flexWrap": "wrap", "gap": "8px", "alignItems": "center", "marginBottom": "8px"},
    )


def build_multicol_section(tunable_names):
    """Render the multicol hypergrid controls shown above the parameter editors."""
    return [
        html.H4("Multicol", className="run-settings-heading"),
        html.Div(
            "Add parameter ranges to generate a custom hypergrid passed to run_scm.py with -hr.",
            style={"marginBottom": "8px", "opacity": "0.85"},
        ),
        html.Div(
            [build_multicol_row({"id": 0, "param": "", "min": "", "max": "", "npoints": ""}, tunable_names)],
            id="run-multicol-rows",
        ),
        html.Div(
            [
                html.Button(
                    "Add parameter",
                    id="run-multicol-add",
                    n_clicks=0,
                    style=run_action_button_style("#111827"),
                )
            ],
            style={"display": "flex", "gap": "8px", "marginBottom": "8px"},
        ),
        html.Div(
            "Format sent to run_scm.py: PARAM/MIN:MAX/NPOINTS,...",
            style={"marginBottom": "14px", "opacity": "0.75", "fontSize": "13px"},
        ),
    ]


def build_run_action_section():
    """Render the primary run/cancel/clear action buttons."""
    return html.Div(
        [
            html.Button("Run selected", id="run-button", n_clicks=0, className="run-button-run-selected", style=run_action_button_style("#111827")),
            html.Button("Cancel runs", id="run-cancel", n_clicks=0, style=run_action_button_style("#b91c1c")),
            html.Button("Clear", id="run-clear", n_clicks=0, style=run_action_button_style("#374151")),
            dcc.Input(
                id="run-max-tasks",
                type="number",
                min=1,
                step=1,
                debounce=True,
                placeholder=str(MAX_RUN_PROCS),
                style={
                    "width": "130px",
                    "padding": "10px 12px",
                    "margin": "4px",
                    "borderRadius": "6px",
                    "border": "1px solid #9ca3af",
                    "fontSize": "14px",
                },
            ),
        ],
        className="run-action-buttons",
        style={"display": "flex", "flexWrap": "wrap", "gap": "8px", "marginTop": "6px"},
    )


def build_console_shell():
    """Render the empty console container shell."""
    return html.Div(id="run-console-container", className="run-console-container", style={"display": "flex", "flexDirection": "column", "gap": "10px"})


def build_left_header(case_groups, case_buttons, stats_buttons):
    """Render the left header block with selections and action controls."""
    return html.Div(
        [
            build_select_actions(case_groups),
            html.Div("Cases:", className="run-section-title"),
            html.Div(case_buttons, className="run-case-buttons"),
            html.Div([html.Div("Stats file:", className="run-section-title"), html.Div(stats_buttons, className="run-stats-buttons")], className="run-stats-section", style={"marginTop": "6px"}),
            build_optional_args_section(),
            build_run_action_section(),
        ],
        className="run-left-header",
        style={"marginBottom": "10px"},
    )


def build_param_input(entry, label_width_px, display_value):
    """Render one editable numeric/text parameter row."""
    return html.Div(
        [
            html.Label(entry["name"], style={"whiteSpace": "nowrap", "minWidth": f"{label_width_px}px"}),
            dcc.Input(id={"type": "run-param", "file": entry["file"], "name": entry["name"]}, type="text", value=display_value, debounce=True, style={"flex": "1 1 auto", "minWidth": "0"}),
        ],
        id={"type": "run-param-container", "file": entry["file"], "name": entry["name"]},
        style=field_style(False),
        className="run-param-container run-param-row--default",
    )


def build_flag_controls(flag_bools, is_true_func):
    """Render the boolean flag checklist controls."""
    controls = []
    for entry in flag_bools:
        controls.append(
            html.Div(
                [
                    dcc.Checklist(
                        id={"type": "run-flag", "name": entry["name"]},
                        className="run-flag-checklist",
                        options=[{"label": entry["name"], "value": "on"}],
                        value=["on"] if is_true_func(entry["value"]) else [],
                        labelStyle={"display": "inline-flex", "alignItems": "center", "gap": "6px"},
                    )
                ],
                id={"type": "run-flag-container", "name": entry["name"]},
                className="run-param-container run-flag-container",
                style=field_style(False),
            )
        )
    return controls


def build_flag_value_section(flag_params, label_width_px, normalize_numeric_display):
    """Render editable non-boolean flag values."""
    if not flag_params:
        return []
    return [
        html.H4("Flag vals", className="run-settings-heading"),
        html.Div([build_param_input({"file": "flags", **entry}, label_width_px, normalize_numeric_display(entry["value"])) for entry in flag_params], className="run-param-list"),
    ]


def build_flags_section(flag_controls):
    """Render boolean flag controls."""
    return [html.H4("Flags", className="run-settings-heading"), html.Div(flag_controls, className="run-param-list")]


def build_tunable_section(tunable_entries, label_width_px, normalize_numeric_display):
    """Render tunable parameter inputs."""
    if not tunable_entries:
        return []
    return [
        html.H4("Tunables", className="run-settings-heading"),
        html.Div([build_param_input({"file": "tunable", **entry}, label_width_px, normalize_numeric_display(entry["value"])) for entry in tunable_entries], className="run-param-list"),
    ]


def build_silhs_section(silhs_entries, label_width_px, normalize_numeric_display):
    """Render SILHS parameter inputs."""
    if not silhs_entries:
        return []
    return [
        html.H4("SILHS", className="run-settings-heading"),
        html.Div([build_param_input({"file": "silhs", **entry}, label_width_px, normalize_numeric_display(entry["value"])) for entry in silhs_entries], className="run-param-list"),
    ]


def build_param_sections(flag_params, flag_controls, tunable_entries, silhs_entries, label_width_px, normalize_numeric_display):
    """Build the full right-pane parameter section list."""
    sections = []
    sections.extend(build_flag_value_section(flag_params, label_width_px, normalize_numeric_display))
    sections.extend(build_flags_section(flag_controls))
    sections.extend(build_tunable_section(tunable_entries, label_width_px, normalize_numeric_display))
    sections.extend(build_silhs_section(silhs_entries, label_width_px, normalize_numeric_display))
    return sections


def build_layout(initial_data):
    """Assemble the full static run-tab layout from precomputed initial metadata."""
    return html.Div(
        [
            dcc.Store(id="run-defaults", data=initial_data["defaults"]),
            dcc.Store(id="run-defaults-by-key", data=initial_data["defaults_by_key"]),
            dcc.Store(id="run-flag-names", data=initial_data["flag_names"]),
            dcc.Store(id="run-param-meta", data=initial_data["param_meta"]),
            dcc.Store(id="run-tunable-names", data=initial_data["tunable_names"]),
            dcc.Store(id="run-multicol-next-id", data=1),
            dcc.Store(id="run-multicol-row-order", data=[0]),
            dcc.Store(id="run-selected-cases", data=[]),
            dcc.Store(id="run-selected-stats-file", data=initial_data["default_stats_name"]),
            dcc.Store(id="run-completed-cases", data=[]),
            dcc.Store(id="run-failed-cases", data=[]),
            dcc.Store(id="run-running-cases", data={}),
            dcc.Store(id="run-queued-cases", data=[]),
            dcc.Store(id="run-max-tasks-active", data=MAX_RUN_PROCS),
            dcc.Store(id="run-case-logs", data={}),
            dcc.Store(id="run-case-commands", data={}),
            dcc.Store(id="run-case-runtimes", data={}),
            dcc.Store(id="run-log-offsets", data={}),
            dcc.Store(id="run-case-order", data=[]),
            dcc.Interval(id="run-interval", interval=500, disabled=True),
            html.Div([build_left_header(initial_data["case_groups"], initial_data["case_buttons"], initial_data["stats_buttons"]), build_console_shell()], className="run-left-pane"),
            html.Div(id="run-pane-divider", className="run-pane-divider"),
            html.Div(build_multicol_section(initial_data["tunable_names"]) + initial_data["param_sections"], id="run-right-pane", className="run-right-pane", style={"paddingLeft": "16px", "height": "calc(100vh - 96px)", "minHeight": 0, "overflowY": "auto", "overflowX": "auto"}),
        ],
        id="run-tab-layout",
        className="run-tab-layout",
        style={"display": "grid", "gridTemplateColumns": f"minmax(0,1fr) 8px {initial_data['right_pane_width_px']}px", "gap": "16px", "padding": "10px", "overflowX": "auto"},
    )
