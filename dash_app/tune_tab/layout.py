"""Static layout builders for the tuning tab."""

from __future__ import annotations

from pathlib import Path
import re
from urllib.parse import urljoin

from dash import dcc, html

from dash_app.shared import styled_dropdown, information_body, notecard
from dash_app.shared.provenance import source_provenance
from dash_app.compile_tab.build_selector import build_selector_trigger
from tuner.taylor_metrics import (
    DEFAULT_AGGREGATION_WEIGHTS,
    DEFAULT_LOSS_WEIGHTS,
)
from tuner.presets import list_presets
from tuner.algorithm_registry import ALGORITHMS, get_algorithm
from .algorithm_controls import build_algorithm_controls
from .state import REPO_ROOT

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


def mode_button_style(selected=False, disabled=False, incompatible=False):
    """Style tuning selectors and highlight methods incompatible with the runtime."""
    style = action_button_style("#2563eb" if selected else "#374151", disabled=disabled)
    style.update({
        "border": "2px solid #f59e0b" if selected else "1px solid transparent",
        "minWidth": "96px",
    })
    if incompatible:
        style.update({"backgroundColor": "#991b1b" if selected else "#7f1d1d",
                      "border": "2px solid #f87171"})
    return style


def config_button_style(selected=False, disabled=False):
    """Return the style for tunable-config selector buttons."""
    style = mode_button_style(selected=selected, disabled=disabled)
    style.update({"minWidth": "116px"})
    return style



def settings_heading(title, section):
    """Pair a settings title with an accessible explanatory notecard button."""
    return html.Div([
        html.Div(title, className="tune-section-title"),
        html.Button("?", id={"type": "tune-settings-help-open", "section": section},
                    type="button", n_clicks=0, className="tune-settings-help-button",
                    title=f"About {title}", **{"aria-label": f"About {title}"}),
    ], className="tune-settings-heading-row")


def setting_row(label, control, *, wide=False):
    """Align one labelled control with the other entries in its settings panel."""
    return html.Div([html.Label(label, htmlFor=control.id), control],
                    className="tune-setting-row" + (" tune-setting-row--wide" if wide else ""))


def settings_help_notecard(section):
    """Explain the selected setting using the shared themed notecard component."""
    if section == "algorithms":
        guide = Path(REPO_ROOT) / "tuner" / "algorithms" / "README.md"
        markdown = guide.read_text(encoding="utf-8")
        revision = source_provenance(REPO_ROOT)["commit"]
        base_url = f"https://github.com/larson-group/clubb/blob/{revision}/tuner/algorithms/README.md"
        # Resolve repository links against the committed guide instead of the Dash URL.
        markdown = re.sub(r"(\[[^\]]+\]\()([^)]+)(\))",
                          lambda match: match[1] + urljoin(base_url, match[2]) + match[3], markdown)
        return notecard(
            "Tuning methods",
            dcc.Markdown(markdown, link_target="_blank", className="tune-methods-guide"),
            {"type": "tune-settings-help-close", "section": section}, size="large",
        )
    specs = {
        "gradient": (
            "Gradients",
            "Gradient-based algorithms use JAX derivatives of the same combined scalar loss used to rank each parameter set. Other algorithms use loss values only.",
            [
                {"heading": "Differentiation mode", "bullets": [
                    "Auto selects the differentiation mode for the requested parameters.",
                    "Forward propagates parameter directions through the model; it is useful for a small number of parameters.",
                    "Backward (reverse) propagates the scalar loss back to the requested parameters, using timestep checkpoints.",
                ]},
                {"heading": "What is differentiated", "paragraphs": [
                    "Field, time-window and across-case aggregation are included. Linked physical parameter derivatives are combined into their shared tuning coordinate.",
                    "The selected algorithm enables gradients automatically. Adam (grad) requires the JAX backend; value-only algorithms do not compute optional gradients from this form.",
                ]},
            ],
        ),
        "loss": (
            "Taylor RMSE weights",
            "These weights control three contributions to each model-versus-reference profile error before fields and time windows are combined into a case loss.",
            [
                {"heading": "The three contributions", "bullets": [
                    "Mean bias: difference between the model and reference profile means, normalized by the reference standard deviation.",
                    "Standard deviation: departure of the model-to-reference standard deviation ratio from 1.",
                    "Correlation: departure of the centered profiles' correlation from 1, scaled by twice the standard deviation ratio.",
                ]},
                {"heading": "Weighted profile RMSE", "paragraphs": [
                    html.Span([
                        "RMSE = √[w", html.Sub("mean"), " b", html.Sup("2"),
                        " + w", html.Sub("std"), " (s − 1)", html.Sup("2"),
                        " + w", html.Sub("corr"), " 2s(1 − r)]",
                    ], style={"display": "block", "fontFamily": "serif", "fontSize": "18px"}),
                    "Here b = (μ_model − μ_ref) / σ_ref is bias_norm, s = σ_model / σ_ref is std_ratio, and r is correlation. μ and σ are the mean and population standard deviation over the compared heights, after time averaging.",
                    "The weights w_mean, w_std and w_corr are the Mean bias, Standard deviation and Correlation entries, respectively. Each multiplies its squared-error contribution inside the square root; the weights are not normalized by their sum.",
                    "This form uses the reported Taylor metrics when both profiles have nonzero variance. For flat profiles, correlation or std_ratio can be undefined; the loss uses the equivalent moment-based calculation so the error remains defined.",
                ]},
                {"heading": "Choosing weights", "paragraphs": [
                    "Weights multiply squared-error contributions before the square root. All weights must be strictly positive and remain fixed during the run.",
                    "Weights of 1, 1 and 1 recover RMSE normalized by the reference profile's standard deviation. Flat reference profiles use the existing unit normalization scale.",
                    "Increasing a weight makes that type of error more important. Increasing all three by the same factor rescales the loss rather than changing their relative importance.",
                ]},
            ],
        ),
        "aggregation": (
            "Across-case aggregation",
            "For each parameter set, the tuner sorts the active case losses from best (smallest) to worst and groups them into four approximately equal parts.",
            [
                {"heading": "Best-to-worst weights", "paragraphs": [
                    "The four entries weight the mean loss of each group. They are relative weights and are normalized automatically; they do not change how many cases belong to a group.",
                    "A larger worst-group weight puts more pressure on difficult cases. A larger best-group weight favors the already well-performing cases. Group membership can change as parameters change.",
                    "Weights must be nonnegative, with a positive weight on at least one nonempty group. Empty groups are omitted and the remaining weights are renormalized. With one case, give the best group a positive weight; the score remains that case loss.",
                ]},
                {"heading": "Case weights and saved revisions", "paragraphs": [
                    "Individual case weights are applied within each group; zero-weight cases are excluded. These four controls apply to quantile-weighted aggregation. Older saved revisions keep their stored aggregation policy.",
                ]},
            ],
        ),
    }
    title, overview, sections = specs[section]
    return notecard(title, information_body(overview, sections),
                    {"type": "tune-settings-help-close", "section": section}, size="medium")


def build_top_controls(initial_data):
    """Render mode selection, mode options, and run controls."""
    selected_mode = initial_data.get("strategy_mode")
    selected_loss_weights = initial_data.get("loss_weights", DEFAULT_LOSS_WEIGHTS)
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
                            "minWidth": 0,
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
                style={"display": "flex", "flexDirection": "column", "alignItems": "stretch"},
            ),
            html.Div(
                [
                    settings_heading("Algorithm", "algorithms"),
                    html.Div([
                        html.Button(get_algorithm(name).label,
                                    id={"type": "tune-algorithm-button", "name": name},
                                    n_clicks=0, style=mode_button_style(selected=selected_mode == name))
                        for name in ALGORITHMS
                    ], style={"display": "flex", "gap": "8px", "flexWrap": "wrap"}),
                    *build_algorithm_controls(selected_mode, initial_data.get("strategy_options")),
                    html.Div(id="tune-algorithm-summary", style={"marginTop": "8px", "opacity": "0.85"}),
                ],
                className="tune-algorithm-controls",
            ),
            html.Div(
                [
                    html.Section([
                        settings_heading("Gradients", "gradient"),
                        setting_row("Mode", styled_dropdown(
                            id="tune-gradient-mode", options=[
                                {"label": "Auto", "value": "auto"},
                                {"label": "Forward", "value": "forward"},
                                {"label": "Backward (reverse)", "value": "reverse"},
                            ], value=initial_data.get("gradient_mode", "auto"), clearable=False,
                            disabled=not get_algorithm(selected_mode or "random").requires_gradients,
                        ), wide=True),
                        html.Small("Enabled automatically by gradient-based algorithms.",
                                   className="tune-setting-caption"),
                    ], className="tune-settings-block"),
                    html.Section([
                        settings_heading("Taylor RMSE weights", "loss"),
                        html.Div([
                            setting_row(label, dcc.Input(
                                id=f"tune-loss-weight-{name}", type="number", min=0,
                                step="any", value=selected_loss_weights[name], debounce=True,
                                className="clubb-input tune-setting-input",
                            )) for name, label in (("mean", "Mean bias"), ("std", "Standard deviation"), ("corr", "Correlation"))
                        ], className="tune-settings-fields"),
                        html.Small("Positive weights; 1 / 1 / 1 gives normalized RMSE.",
                                   className="tune-setting-caption"),
                    ], className="tune-settings-block"),
                ], className="tune-gradient-controls",
            ),
            html.Div(
                [
                    settings_heading("Across-case aggregation", "aggregation"),
                    html.Small("Relative importance, best to worst.", className="tune-setting-caption"),
                    html.Div([
                        setting_row(label, dcc.Input(
                            id=f"tune-aggregation-weight-{index + 1}", type="number", min=0,
                            step="any", debounce=True,
                            value=aggregation_weights[index] if index < len(aggregation_weights) else DEFAULT_AGGREGATION_WEIGHTS[index],
                            className="clubb-input tune-setting-input",
                        )) for index, label in enumerate(("Best 25%", "Lower-middle 25%", "Upper-middle 25%", "Worst 25%"))
                    ], className="tune-settings-fields"),
                ], className="tune-aggregation-controls",
            ),
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
                            build_selector_trigger("tune-selected-build-badge"),
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
                className="tune-run-controls",
                style={"display": "flex", "flexDirection": "column", "alignItems": "flex-start"},
            ),
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
                    html.Button(
                        "Save tuned config",
                        id="tune-config-save-open",
                        n_clicks=0,
                        title="Save one of the top tuning results as a named CLUBB configuration.",
                        style=action_button_style("#0f766e", disabled=True),
                        disabled=True,
                    ),
                    html.Div(
                        id="tune-config-save-status",
                        className="tune-info-message tune-config-save-status",
                    ),
                ],
                className="tune-result-controls",
                style={"display": "flex", "flexDirection": "column", "alignItems": "flex-start"},
            ),
        ],
        className="tune-header-controls",
    )


def build_param_range_row(row, tunable_names):
    """Render one parameter-range row."""
    row_id = row.get("id")
    raw_targets = row.get("targets")
    if raw_targets is None:
        raw_targets = [row.get("param", "")]
    targets = [str(target or "").strip() for target in raw_targets]
    if row.get("locked") and len(targets) < 2:
        targets.extend([""] * (2 - len(targets)))
    if not targets:
        targets = [""]
    options = [{"label": name, "value": name} for name in tunable_names]

    members = []
    for member_index, target in enumerate(targets):
        member = [
            dcc.Dropdown(
                id={"type": "tune-range-member", "row": row_id, "member": member_index},
                options=options,
                value=target or None,
                placeholder="parameter",
                clearable=True,
                searchable=True,
                className="clubb-dropdown",
                style={"minWidth": "170px", "flex": "1 1 170px"},
            )
        ]
        if len(targets) > 1:
            member.append(
                html.Button(
                    "Remove",
                    id={"type": "tune-range-member-remove", "row": row_id, "member": member_index},
                    n_clicks=0,
                    style=action_button_style("#6b7280"),
                )
            )
        members.append(
            html.Div(
                member,
                style={"display": "flex", "gap": "8px", "alignItems": "center"},
            )
        )

    member_block = html.Div(members, style={"display": "flex", "flexDirection": "column", "gap": "7px"})
    if len(targets) > 1:
        member_block = html.Div(
            [
                html.Div(
                    style={
                        "width": "12px",
                        "alignSelf": "stretch",
                        "borderLeft": "3px solid #38bdf8",
                        "borderTop": "3px solid #38bdf8",
                        "borderBottom": "3px solid #38bdf8",
                        "borderRadius": "5px 0 0 5px",
                    },
                    title="These parameters share one sampled value.",
                ),
                html.Div(
                    [
                        member_block,
                        html.Button(
                            "+ Add to group",
                            id={"type": "tune-range-member-add", "index": row_id},
                            n_clicks=0,
                            style={**action_button_style("#0369a1"), "alignSelf": "flex-start"},
                        ),
                    ],
                    style={"display": "flex", "flexDirection": "column", "gap": "7px"},
                ),
            ],
            style={"display": "flex", "gap": "8px", "alignItems": "stretch"},
        )

    return html.Div(
        [
            member_block,
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
        style={
            "display": "flex",
            "flexWrap": "wrap",
            "gap": "8px",
            "alignItems": "center",
            "marginBottom": "10px",
            "padding": "8px",
            "border": "1px solid #334155" if len(targets) > 1 else "none",
            "borderRadius": "6px",
        },
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
                title=f"Use the {value} parameter and flag presets",
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


def build_tuned_config_save_dialog():
    """Render the Tune result selector and named-config save form."""
    return html.Div(
        html.Div(
            [
                html.Div(
                    [
                        html.Div("Save tuned configuration", className="shared-notecard-title"),
                        html.Div(
                            "Clone the run's base config, apply its SCM override, then apply one retained result.",
                            className="tune-config-save-subtitle",
                        ),
                    ],
                    className="shared-notecard-header tune-config-save-header",
                ),
                html.Div(
                    [
                        html.Div(
                            [
                                html.Label("Tuning result", htmlFor="tune-config-save-result"),
                                dcc.Dropdown(
                                    id="tune-config-save-result",
                                    options=[],
                                    value=None,
                                    clearable=False,
                                    searchable=True,
                                    className="clubb-dropdown",
                                ),
                            ],
                            className="tune-config-save-field",
                        ),
                        html.Div(
                            [
                                html.Label("Config name", htmlFor="tune-config-save-name"),
                                dcc.Input(
                                    id="tune-config-save-name",
                                    type="text",
                                    value="",
                                    placeholder="example: tuned_arm_2026",
                                    className="tune-config-save-input",
                                ),
                            ],
                            className="tune-config-save-field",
                        ),
                        html.Div(id="tune-config-save-feedback"),
                    ],
                    className="shared-notecard-body tune-config-save-form",
                ),
                html.Div(
                    [
                        html.Button(
                            "Cancel",
                            id="tune-config-save-cancel",
                            n_clicks=0,
                            className="tune-config-save-cancel",
                        ),
                        html.Button(
                            "Save config",
                            id="tune-config-save-submit",
                            n_clicks=0,
                            className="tune-config-save-submit",
                        ),
                    ],
                    className="tune-config-save-footer",
                ),
            ],
            className="shared-notecard-panel shared-notecard-size-small tune-config-save-panel",
            role="dialog",
            **{"aria-modal": "true"},
        ),
        id="tune-config-save-modal",
        className="shared-notecard-overlay tune-config-save-modal--hidden",
    )


def build_layout(initial_data):
    """Assemble the full static tuning-tab layout."""
    initial_param_rows = list(initial_data["initial_param_rows"])
    initial_case_rows = list(initial_data["initial_case_rows"])
    return html.Div(
        [
            dcc.Store(id="tune-config-save-overwrite", data=None),
            build_tuned_config_save_dialog(),
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
            dcc.Store(id="tune-loss-weights", data=initial_data["loss_weights"]),
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
            dcc.Store(id="tune-landscape-signature", data=""),
            html.Div(id="tune-settings-help-modal"),
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
                                                    dcc.RadioItems(
                                                        id="tune-parameter-box-scale",
                                                        className="tune-segmented-control tune-parameter-box-scale",
                                                        options=[
                                                            {"label": "Range normalized", "value": "normalized"},
                                                            {"label": "Unscaled", "value": "unscaled"},
                                                        ],
                                                        value="normalized",
                                                        inline=True,
                                                    ),
                                                ],
                                                className="tune-parameter-box-controls",
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
                                                    [
                                                        html.Button("Add parameter", id="tune-range-add", n_clicks=0, style=action_button_style("#111827")),
                                                        html.Button(
                                                            "Add locked parameters",
                                                            id="tune-range-add-locked",
                                                            n_clicks=0,
                                                            style=action_button_style("#0369a1"),
                                                        ),
                                                    ],
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
