"""Static layout builders and UI styling helpers for the run tab."""

from dash import dcc, html

from dash_app.persistence import WORKSPACE_TOKEN
from dash_app.compile_tab.build_selector import build_selector_trigger

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


def run_config_button_style(selected=False, disabled=False):
    """Return the style for one tunable-config selection button."""
    style = run_action_button_style("#2563eb" if selected else "#374151", disabled=disabled)
    style.update(
        {
            "border": "3px solid #f59e0b" if selected else "1px solid transparent",
            "minWidth": "174px",
            "padding": "15px 24px",
            "fontSize": "16px",
            "boxSizing": "border-box",
        }
    )
    return style


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
            html.Button("Deselect", id="run-deselect", n_clicks=0, className="run-select-action", style=run_action_button_style("#6b7280")),
            html.Button("Select all", id={"type": "run-group-button", "name": "all"}, n_clicks=0, disabled=not case_groups["all"], className="run-select-group-button", style=run_action_button_style("#111827", disabled=not case_groups["all"])),
            html.Button("Select standard", id={"type": "run-group-button", "name": "standard"}, n_clicks=0, disabled=not case_groups["standard"], className="run-select-group-button", style=run_action_button_style("#111827", disabled=not case_groups["standard"])),
            html.Button("Select priority", id={"type": "run-group-button", "name": "priority"}, n_clicks=0, disabled=not case_groups["priority"], className="run-select-group-button", style=run_action_button_style("#111827", disabled=not case_groups["priority"])),
            html.Button("Select minimum", id={"type": "run-group-button", "name": "minimum"}, n_clicks=0, disabled=not case_groups["minimum"], className="run-select-group-button", style=run_action_button_style("#111827", disabled=not case_groups["minimum"])),
            html.Button("Select short", id={"type": "run-group-button", "name": "short"}, n_clicks=0, disabled=not case_groups["short"], className="run-select-group-button", style=run_action_button_style("#111827", disabled=not case_groups["short"])),
        ],
        className="run-select-actions",
    )


def build_case_buttons(cases):
    """Render all case buttons with the default unselected style."""
    return [
        html.Button(
            case_name,
            id={"type": "run-case-button", "name": case_name},
            n_clicks=0,
            className="run-case-button",
            style=case_button_style("#2563eb", False),
            **{"data-case-name": case_name},
        )
        for case_name in cases
    ]


def build_stats_buttons(stats_files, default_stats_name, no_stats_name):
    """Render stats-file selection buttons, including the synthetic none option."""
    buttons = [
        html.Button(stats_name, id={"type": "run-stats-button", "name": stats_name}, n_clicks=0, className="run-stats-button", style=stats_button_style(stats_name == default_stats_name))
        for stats_name in stats_files
    ]
    if no_stats_name not in stats_files:
        buttons.append(
            html.Button("none", id={"type": "run-stats-button", "name": no_stats_name}, n_clicks=0, className="run-stats-button", style=stats_button_style(no_stats_name == default_stats_name))
        )
    return buttons


def build_run_config_buttons(configs, selected_config):
    """Render selectable run-tab tunable config buttons."""
    buttons = []
    for config in configs or []:
        value = str(config.get("value", "")).strip()
        if not value:
            continue
        label = str(config.get("label") or value)
        buttons.append(
            html.Button(
                label,
                id={"type": "run-config-button", "name": value},
                n_clicks=0,
                title=f"Use input/parameter_and_flag_configs/{value}",
                style=run_config_button_style(selected=value == selected_config),
            )
        )
    save_style = run_config_button_style(disabled=not buttons)
    save_style.update(
        {
            "backgroundColor": "#0f766e" if buttons else "#9ca3af",
            "border": "2px dashed #5eead4" if buttons else "2px dashed #cbd5e1",
        }
    )
    buttons.append(
        html.Button(
            "Save config",
            id="run-config-save",
            n_clicks=0,
            disabled=not buttons,
            title="Save the current settings as a named config",
            style=save_style,
        )
    )
    return html.Div(buttons, style={"display": "flex", "gap": "12px", "flexWrap": "wrap"})


def build_config_save_dialog():
    """Render the name/note form used to save the current Run configuration."""
    return html.Div(
        html.Div(
            [
                html.Div(
                    [
                        html.Div(
                            "Save configuration",
                            id="run-config-save-dialog-title",
                            className="shared-notecard-title",
                        ),
                        html.Div(
                            "Clone the selected config with the current settings",
                            className="run-config-save-subtitle",
                        ),
                    ],
                    className="shared-notecard-header run-config-save-header",
                ),
                html.Div(
                    [
                        html.Div(
                            [
                                html.Label("Config name", htmlFor="run-config-save-name"),
                                dcc.Input(
                                    id="run-config-save-name",
                                    type="text",
                                    value="",
                                    placeholder="example: adg2_experiment",
                                    className="run-config-save-input",
                                ),
                            ],
                            className="run-config-save-field",
                        ),
                        html.Div(
                            [
                                html.Label("Note", htmlFor="run-config-save-note"),
                                dcc.Textarea(
                                    id="run-config-save-note",
                                    value="",
                                    placeholder="Optional context, purpose, or follow-up notes…",
                                    className="run-config-save-note",
                                ),
                                html.Div(
                                    "Saved in the new config's README.md.",
                                    className="run-config-save-help",
                                ),
                            ],
                            className="run-config-save-field",
                        ),
                        html.Div(id="run-config-save-feedback"),
                    ],
                    className="shared-notecard-body run-config-save-form",
                ),
                html.Div(
                    [
                        html.Button(
                            "Cancel",
                            id="run-config-save-cancel",
                            type="button",
                            n_clicks=0,
                            className="run-config-save-cancel",
                        ),
                        html.Button(
                            "Save config",
                            id="run-config-save-submit",
                            type="button",
                            n_clicks=0,
                            className="run-config-save-submit",
                        ),
                    ],
                    className="run-config-save-footer",
                ),
            ],
            className="shared-notecard-panel shared-notecard-size-small run-config-save-panel",
            role="dialog",
            **{
                "aria-modal": "true",
                "aria-labelledby": "run-config-save-dialog-title",
            },
        ),
        id="run-config-save-modal",
        className="shared-notecard-overlay run-config-save-modal--hidden",
    )


def build_optional_args_section():
    """Render the optional run_scm.py argument inputs."""
    return html.Div(
        [
            html.Div(
                [
                    html.Div("Execution options", className="run-setup-section-title"),
                    html.Div("Optional limits and CLI settings", className="run-setup-section-note"),
                ],
                className="run-setup-section-heading",
            ),
            html.Div(
                [
                    dcc.Input(id="run-opt-max-iters", type="text", value="", placeholder="max_iters", className="run-quick-input", style={"width": "130px"}),
                    dcc.Input(id="run-opt-debug", type="text", value="", placeholder="debug", className="run-quick-input", style={"width": "130px"}),
                    dcc.Input(id="run-opt-dt-main", type="text", value="", placeholder="dt_main", className="run-quick-input", style={"width": "130px"}),
                    dcc.Input(id="run-opt-dt-rad", type="text", value="", placeholder="dt_rad", className="run-quick-input", style={"width": "130px"}),
                    dcc.Input(id="run-opt-tout", type="text", value="", placeholder="tout", className="run-quick-input", style={"width": "130px"}),
                ],
                className="run-option-input-row",
            ),
            html.Div(
                [
                    html.Label("Output", style={"whiteSpace": "nowrap", "alignSelf": "center"}),
                    dcc.Input(
                        id="run-opt-out-dir",
                        type="text",
                        value="dash_default",
                        placeholder="test1 → output/test1",
                        className="run-quick-input run-output-dir-input",
                        style={"width": "200px"},
                    ),
                    html.Span(
                        id="run-output-dir-warning",
                        className="run-output-dir-warning run-output-dir-warning--hidden",
                    ),
                ],
                className="run-option-detail-row",
            ),
            html.Div(
                [
                    html.Label("Extra args", htmlFor="run-opt-extra-args", style={"whiteSpace": "nowrap", "alignSelf": "center"}),
                    dcc.Input(
                        id="run-opt-extra-args",
                        type="text",
                        value="",
                        debounce=True,
                        placeholder="-max_iters 10 -dt_main 60",
                        className="run-quick-input",
                        style={"minWidth": "280px", "flex": "1 1 360px"},
                    ),
                ],
                className="run-option-detail-row",
            ),
        ],
        className="run-setup-section run-optional-args",
    )


def build_multicol_row(row, tunable_names):
    """Render one multicol hypergrid specification row."""
    row_id = row.get("id")
    selected_param = str(row.get("param", "") or "").strip()
    options = [{"label": name, "value": name} for name in tunable_names]
    if selected_param and selected_param not in set(tunable_names or []):
        options.insert(
            0,
            {
                "label": f"{selected_param} (not in config)",
                "value": selected_param,
            },
        )
    return html.Div(
        [
            dcc.Dropdown(
                id={"type": "run-hr-param", "index": row_id},
                options=options,
                value=selected_param or None,
                placeholder="parameter",
                clearable=True,
                searchable=True,
                persistence=WORKSPACE_TOKEN,
                persistence_type="local",
                className="clubb-dropdown",
                style={"minWidth": "170px", "flex": "2 1 170px"},
            ),
            dcc.Input(
                id={"type": "run-hr-min", "index": row_id},
                type="text",
                value=row.get("min", ""),
                placeholder="min",
                persistence=WORKSPACE_TOKEN,
                persistence_type="local",
                style={"width": "90px", "flex": "0 0 90px"},
            ),
            dcc.Input(
                id={"type": "run-hr-max", "index": row_id},
                type="text",
                value=row.get("max", ""),
                placeholder="max",
                persistence=WORKSPACE_TOKEN,
                persistence_type="local",
                style={"width": "90px", "flex": "0 0 90px"},
            ),
            dcc.Input(
                id={"type": "run-hr-npoints", "index": row_id},
                type="text",
                value=row.get("npoints", ""),
                placeholder="npoints",
                persistence=WORKSPACE_TOKEN,
                persistence_type="local",
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


def build_multicol_section(tunable_names, rows=None):
    """Render the multicol hypergrid controls shown above the parameter editors."""
    row_data = list(rows or [])
    return [
        html.H4("Multicol", className="run-settings-heading"),
        html.Div(
            "Add parameter ranges to generate a custom hypergrid passed to run_scm.py with -multicol.",
            style={"marginBottom": "8px", "opacity": "0.85"},
        ),
        html.Div(
            [build_multicol_row(row, tunable_names) for row in row_data],
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


def build_run_limit_control(label, control_id, value, disabled=False, title=None):
    """Render one compact numeric run-limit control."""
    return html.Div(
        [
            html.Label(label, htmlFor=control_id, className="run-limit-label"),
            dcc.Input(
                id=control_id,
                type="text",
                value=str(value),
                debounce=True,
                persistence=WORKSPACE_TOKEN,
                persistence_type="local",
                disabled=disabled,
                className="clubb-input run-limit-input",
            ),
        ],
        className="run-limit-control",
        title=title,
    )


def build_run_action_section():
    """Render the primary run/cancel/clear action buttons."""
    return html.Div(
        [
            html.Div(
                [
                    html.Div("Launch", className="run-launch-title"),
                    html.Div("Selected cases run with the settings above.", className="run-launch-note"),
                ],
                className="run-launch-intro",
            ),
            html.Div(
                [
                    build_selector_trigger("run-selected-build-badge"),
                    html.Button("Run selected", id="run-button", n_clicks=0, className="run-button-run-selected", style=run_action_button_style("#111827")),
                    html.Button("Cancel runs", id="run-cancel", n_clicks=0, disabled=False, style=run_action_button_style("#b91c1c")),
                    html.Button("Clear", id="run-clear", n_clicks=0, style=run_action_button_style("#374151")),
                ],
                className="run-primary-actions",
            ),
            html.Div(
                [
                    build_run_limit_control("Workers", "run-max-tasks", MAX_RUN_PROCS),
                    build_run_limit_control("Batch size", "run-batch-size", 8, disabled=True, title="Used only for multicol runs."),
                ],
                className="run-limit-panel",
            ),
        ],
        className="run-action-bar",
    )


def build_console_shell(_cases):
    """Render the browser-owned container for broker-discovered runs."""
    return html.Div(
        [
            html.Div("No runs yet.", id="run-console-empty", className="run-empty-message"),
        ],
        id="run-console-container",
        className="run-console-container",
        style={
            "display": "flex",
            "flexDirection": "column",
            "gap": "10px",
            "minHeight": "220px",
            "overflowAnchor": "none",
        },
    )


def build_left_header(case_groups, case_buttons, stats_buttons):
    """Render the left header block with selections and action controls."""
    return html.Div(
        [
            html.Div(
                [
                    html.Div(
                        [
                            html.Div("Cases", className="run-setup-section-title"),
                            html.Div("Quick-select a group or choose individual cases.", className="run-setup-section-note"),
                        ],
                        className="run-setup-section-heading",
                    ),
                    build_select_actions(case_groups),
                    html.Div(case_buttons, className="run-case-buttons"),
                ],
                className="run-setup-section run-cases-section",
            ),
            html.Div(
                [
                    html.Div(
                        [
                            html.Div("Statistics", className="run-setup-section-title"),
                            html.Div("Select the requested CLUBB statistics definition.", className="run-setup-section-note"),
                        ],
                        className="run-setup-section-heading",
                    ),
                    html.Div(stats_buttons, className="run-stats-buttons"),
                ],
                className="run-setup-section run-stats-section",
            ),
            build_optional_args_section(),
            build_run_action_section(),
        ],
        className="run-left-header",
    )


def build_param_input(entry, label_width_px, value_width_px, display_value):
    """Render one editable numeric/text parameter row."""
    return html.Div(
        [
            html.Label(
                entry["name"],
                title=entry["name"],
                style={
                    "whiteSpace": "nowrap",
                    "flex": f"1 1 {label_width_px}px",
                    "width": f"{label_width_px}px",
                    "minWidth": "0",
                    "overflow": "hidden",
                    "textOverflow": "ellipsis",
                },
            ),
            dcc.Input(
                id={"type": "run-param", "file": entry["file"], "name": entry["name"]},
                type="text",
                value=display_value,
                disabled=bool(entry.get("disabled", False)),
                debounce=True,
                persistence=WORKSPACE_TOKEN,
                persistence_type="local",
                className="run-param-value-input",
                style={
                    "flex": f"0 0 {value_width_px}px",
                    "width": f"{value_width_px}px",
                    "minWidth": f"{value_width_px}px",
                    "boxSizing": "border-box",
                },
            ),
        ],
        id={"type": "run-param-container", "file": entry["file"], "name": entry["name"]},
        style=field_style(False),
        className="run-param-container run-param-row--default",
    )


def linked_parameter_key(members):
    """Return the stable UI key for one equal-value parameter group."""
    return "=".join(str(member) for member in members)


def build_linked_tunable_input(entries, members, value_width_px, normalize_numeric_display):
    """Render a single visible value with hidden physical member inputs.

    The hidden inputs retain the existing callback IDs and namelist plumbing;
    a callback mirrors the visible logical value to every member.  This keeps
    the Run request explicit while making an equality constraint obvious.
    """
    by_name = {str(entry["name"]): entry for entry in entries}
    group_entries = [by_name[name] for name in members if name in by_name]
    if len(group_entries) != len(members):
        return []
    group_key = linked_parameter_key(members)
    display_value = normalize_numeric_display(group_entries[0]["value"])
    # Keep a hidden physical row for each member.  The existing ALL-pattern
    # style callback intentionally sees every ``run-param`` input and every
    # matching ``run-param-container``.  Preserving that one-to-one shape
    # avoids a linked group changing the callback's wildcard output length.
    shadow_inputs = [
        html.Div(
            dcc.Input(
                id={"type": "run-param", "file": "tunable", "name": entry["name"]},
                type="text",
                value=normalize_numeric_display(entry["value"]),
                disabled=bool(entry.get("disabled", False)),
                debounce=True,
                persistence=WORKSPACE_TOKEN,
                persistence_type="local",
                className="run-linked-param-shadow",
            ),
            id={"type": "run-param-container", "file": "tunable", "name": entry["name"]},
            className="run-param-container run-param-row--default",
            style={"display": "none"},
        )
        for entry in group_entries
    ]
    return html.Div(
        [
            html.Div(
                [html.Span(name, className="run-linked-param-name") for name in members],
                className="run-linked-param-names",
                title="These parameters are constrained to the same value by CLUBB.",
            ),
            html.Div("linked", className="run-linked-param-badge"),
            dcc.Input(
                id={"type": "run-linked-param", "group": group_key},
                type="text",
                value=display_value,
                debounce=True,
                persistence=WORKSPACE_TOKEN,
                persistence_type="local",
                className="run-param-value-input run-linked-param-value",
                style={"width": f"{value_width_px}px", "minWidth": f"{value_width_px}px", "boxSizing": "border-box"},
            ),
            *shadow_inputs,
        ],
        id={"type": "run-linked-param-container", "group": group_key},
        className="run-linked-param-container run-param-row--default",
        # A linked control has one visible row per member plus its shared
        # value.  Reserving compact grid tracks lets ordinary controls fill
        # the other columns instead of stretching to this card's height.
        style={"gridRow": f"span {max(2, len(members))}"},
    )


def build_flag_control(entry, is_true_func):
    """Render one physical boolean flag control, retaining its stable ID."""
    return html.Div(
        [
            dcc.Checklist(
                id={"type": "run-flag", "name": entry["name"]},
                className="run-flag-checklist",
                options=[{"label": entry["name"], "value": "on"}],
                value=["on"] if is_true_func(entry["value"]) else [],
                persistence=WORKSPACE_TOKEN,
                persistence_type="local",
                labelStyle={"display": "inline-flex", "alignItems": "center", "gap": "6px"},
            )
        ],
        id={"type": "run-flag-container", "name": entry["name"]},
        className="run-param-container run-flag-container run-param-row--default",
    )


def build_flag_controls(flag_bools, is_true_func, relationships=None):
    """Render boolean flags, grouping declared related options together."""
    by_name = {str(entry["name"]): entry for entry in flag_bools}
    relation_by_member = {
        member: relation
        for relation in (relationships or [])
        for member in relation.get("members") or []
    }
    rendered_relationships: set[tuple[str, ...]] = set()
    controls = []
    for entry in flag_bools:
        name = str(entry["name"])
        relation = relation_by_member.get(name)
        if not relation:
            controls.append(build_flag_control(entry, is_true_func))
            continue
        members = tuple(str(member) for member in relation.get("members") or [])
        if members in rendered_relationships:
            continue
        rendered_relationships.add(members)
        related_entries = [by_name[member] for member in members if member in by_name]
        controls.append(
            html.Div(
                [
                    html.Div(
                        [
                            html.Span(relation.get("label", "linked"), className="run-linked-flag-badge"),
                            html.Span(relation.get("description", "Related CLUBB settings."), className="run-linked-flag-description"),
                        ],
                        className="run-linked-flag-heading",
                    ),
                    html.Div(
                        [build_flag_control(related, is_true_func) for related in related_entries],
                        className="run-linked-flag-members",
                    ),
                ],
                className="run-linked-flag-group",
                # One track for the relationship note and one for every
                # physical checkbox.  This makes the group an atomic card
                # while preserving dense packing of the surrounding flags.
                style={"gridRow": f"span {max(2, len(related_entries) + 1)}"},
            )
        )
    return controls


def build_flag_value_section(flag_params, label_width_px, value_width_px, normalize_numeric_display):
    """Render editable non-boolean flag values."""
    if not flag_params:
        return []
    return [
        html.H4("Flag vals", className="run-settings-heading"),
        html.Div(
            [
                build_param_input({"file": "flags", **entry}, label_width_px, value_width_px, normalize_numeric_display(entry["value"]))
                for entry in flag_params
            ],
            className="run-param-list",
        ),
    ]


def build_flags_section(flag_controls):
    """Render boolean flag controls."""
    return [html.H4("Flags", className="run-settings-heading"), html.Div(flag_controls, className="run-param-list")]


def build_tunable_controls(tunable_entries, linked_groups, label_width_px, value_width_px, normalize_numeric_display):
    """Build logical tunable controls, collapsing active equal-value groups."""
    active_groups = [tuple(group) for group in (linked_groups or [])]
    group_by_member = {member: group for group in active_groups for member in group}
    rendered_groups: set[tuple[str, ...]] = set()
    controls = []
    for entry in tunable_entries:
        name = str(entry["name"])
        group = group_by_member.get(name)
        if group:
            if group in rendered_groups:
                continue
            rendered_groups.add(group)
            controls.append(build_linked_tunable_input(tunable_entries, group, value_width_px, normalize_numeric_display))
            continue
        controls.append(build_param_input({"file": "tunable", **entry}, label_width_px, value_width_px, normalize_numeric_display(entry["value"])))
    return controls


def build_tunable_section(tunable_entries, linked_groups, label_width_px, value_width_px, normalize_numeric_display):
    """Render tunable parameter inputs."""
    if not tunable_entries:
        return []
    return [
        html.H4("Tunables", className="run-settings-heading"),
        html.Div(
            build_tunable_controls(tunable_entries, linked_groups, label_width_px, value_width_px, normalize_numeric_display),
            id="run-tunable-controls",
            className="run-param-list",
        ),
    ]


def build_silhs_section(silhs_entries, label_width_px, value_width_px, normalize_numeric_display):
    """Render SILHS parameter inputs."""
    if not silhs_entries:
        return []
    return [
        html.H4("SILHS", className="run-settings-heading"),
        html.Div(
            [
                build_param_input({"file": "silhs", **entry}, label_width_px, value_width_px, normalize_numeric_display(entry["value"]))
                for entry in silhs_entries
            ],
            className="run-param-list",
        ),
    ]


def build_param_sections(flag_params, flag_controls, tunable_entries, linked_groups, silhs_entries, label_width_px, value_width_px, normalize_numeric_display):
    """Build the full right-pane parameter section list."""
    sections = []
    sections.extend(build_flag_value_section(flag_params, label_width_px, value_width_px, normalize_numeric_display))
    sections.extend(build_flags_section(flag_controls))
    sections.extend(build_tunable_section(tunable_entries, linked_groups, label_width_px, value_width_px, normalize_numeric_display))
    sections.extend(build_silhs_section(silhs_entries, label_width_px, value_width_px, normalize_numeric_display))
    return sections


def build_right_pane(initial_data):
    """Render the run-tab right pane for the selected tunable config."""
    return [
        html.H4("Config", className="run-settings-heading"),
        html.Div(
            build_run_config_buttons(
                initial_data["tunable_configs"],
                initial_data["selected_config"],
            ),
            id="run-config-buttons",
            style={"marginBottom": "8px"},
        ),
        html.Div(id="run-settings-resolution-note", className="run-settings-resolution-note"),
    ] + build_multicol_section(
        initial_data["tunable_names"],
        initial_data.get("multicol_rows"),
    ) + initial_data["param_sections"]


def build_layout(initial_data):
    """Assemble the full static run-tab layout from precomputed initial metadata."""
    multicol_rows = initial_data.get("multicol_rows")
    if multicol_rows is None:
        multicol_rows = []
    initial_data = {**initial_data, "multicol_rows": multicol_rows}
    row_ids = [row["id"] for row in multicol_rows]
    return html.Div(
        [
            dcc.Store(id="run-settings-schema", data=initial_data.get("settings_schema") or {}),
            dcc.Store(id="run-param-meta", data=initial_data["param_meta"]),
            dcc.Store(id="run-tunable-names", data=initial_data["tunable_names"]),
            dcc.Store(id="run-tunable-default-ranges", data=initial_data["tunable_default_ranges"]),
            dcc.Store(id="run-linked-parameter-groups", data=initial_data.get("linked_parameter_groups") or []),
            dcc.Store(id="run-tunable-configs", data=initial_data["tunable_configs"]),
            dcc.Store(id="run-selected-config", data=initial_data["selected_config"]),
            # Unlike the selected config, this describes the controls that are
            # actually on screen and must never be restored from the browser.
            dcc.Store(id="run-rendered-config", data=initial_data["selected_config"]),
            dcc.Store(id="run-settings-resolution", data=initial_data.get("settings_resolution") or {}),
            # A config button deliberately discards the saved per-field Run
            # values before the server rebuilds the matching controls.
            dcc.Store(id="run-config-reset-signal"),
            dcc.Store(id="run-config-save-request"),
            dcc.Store(id="run-config-save-overwrite"),
            build_config_save_dialog(),
            html.Div(
                id="run-config-save-status",
                hidden=True,
                style={
                    "position": "fixed",
                    "right": "24px",
                    "bottom": "24px",
                    "zIndex": "2000",
                    "padding": "10px 14px",
                    "borderRadius": "6px",
                    "backgroundColor": "#0f766e",
                    "color": "white",
                },
            ),
            # Ephemeral acknowledgement for the browser-side removal of
            # transient Run local-storage entries after Clear.
            dcc.Store(id="run-clear-persistence-signal"),
            dcc.Store(
                id="run-multicol-rows-state",
                data=multicol_rows,
            ),
            dcc.Store(id="run-multicol-next-id", data=max(row_ids, default=-1) + 1),
            dcc.Store(id="run-multicol-row-order", data=row_ids),
            dcc.Store(id="run-selected-cases", data=[]),
            dcc.Store(id="run-selected-stats-file", data=initial_data["default_stats_name"]),
            dcc.Store(id="run-resolved-output-dir"),
            dcc.Store(id="run-action-result"),
            dcc.Store(id="run-ui-render-signal"),
            html.Div([build_left_header(initial_data["case_groups"], initial_data["case_buttons"], initial_data["stats_buttons"]), build_console_shell(initial_data["cases"])], className="run-left-pane"),
            html.Div(id="run-pane-divider", className="run-pane-divider"),
            # The settings pane intentionally grows with its controls.  It is
            # part of the document flow, so the browser's main scrollbar—not
            # a narrow nested pane—owns the reading position.
            html.Div(
                build_right_pane(initial_data),
                id="run-right-pane",
                className="run-right-pane",
                style={"paddingLeft": "16px", "paddingRight": "16px"},
            ),
        ],
        id="run-tab-layout",
        className="run-tab-layout",
        style={"display": "grid", "gridTemplateColumns": f"minmax(0,1fr) 8px {initial_data['right_pane_width_px']}px", "gap": "16px", "padding": "10px"},
    )
