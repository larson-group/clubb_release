"""Static layout for the Profile benchmark and comparison workspace."""

from __future__ import annotations

from dash import dcc, html

from dash_app.shared.components import styled_dropdown
from dash_app.compile_tab.build_selector import build_selector_trigger

from .figures import METRICS, X_AXES, empty_profile_figure
from .library import profile_option


def _control_label(label: str, tooltip: str):
    return html.Label(
        [
            html.Span(label),
            html.Span("?", className="profile-help-icon", title=tooltip, **{"aria-label": tooltip}),
        ],
        className="profile-field-label",
        title=tooltip,
    )


def _field(
    label: str,
    control,
    *,
    tooltip: str,
    help_text: str = "",
    compact: bool = False,
    class_name: str = "",
):
    children = [_control_label(label, tooltip), control]
    if help_text:
        children.append(html.Div(help_text, className="profile-field-help"))
    classes = ["profile-field"]
    if compact:
        classes.append("profile-field-compact")
    if class_name:
        classes.append(class_name)
    return html.Div(children, className=" ".join(classes), title=tooltip)


def _number_input(component_id: str, value, *, minimum=None, maximum=None, step=1):
    return dcc.Input(
        id=component_id,
        type="number",
        value=value,
        min=minimum,
        max=maximum,
        step=step,
        debounce=True,
        className="profile-input",
    )


def _section_header(kicker: str, title: str, note: str):
    return html.Div(
        [
            html.Div(kicker, className="profile-section-kicker"),
            html.Div(
                [
                    html.H3(title, className="profile-section-title"),
                    html.P(note, className="profile-section-note"),
                ],
                className="profile-section-heading-copy",
            ),
        ],
        className="profile-section-heading",
    )


def _benchmark_group(number: str, title: str, note: str, children, *, class_name: str = ""):
    classes = " ".join(part for part in ("profile-benchmark-group", class_name) if part)
    return html.Section(
        [
            html.Div(
                [
                    html.Span(number, className="profile-benchmark-group-number", **{"aria-hidden": "true"}),
                    html.Div(
                        [
                            html.H3(title, className="profile-benchmark-group-title"),
                            html.P(note, className="profile-benchmark-group-note"),
                        ]
                    ),
                ],
                className="profile-benchmark-group-heading",
            ),
            html.Div(children, className="profile-benchmark-group-fields"),
        ],
        className=classes,
    )


def _profile_record_title(record: dict[str, object]) -> str:
    parts = [profile_option(record)["label"]]
    started = str(record.get("started_utc") or "").replace("T", " ")[:16]
    revision = str(record.get("git_revision") or "")[:8]
    status = str(record.get("status") or "")
    if started:
        parts.append(started)
    if revision:
        parts.append(revision)
    if status:
        parts.append(status)
    return " · ".join(parts)


def _profile_delete_button(record, confirmation=None):
    """Render the same guarded delete action used by the Plot output chooser."""
    run_id = str(record.get("run_id") or "")
    confirming = bool(run_id) and (confirmation or {}).get("run_id") == run_id
    return html.Button(
        "Confirm" if confirming else "×",
        id={"type": "profile-delete-run", "run_id": run_id},
        n_clicks=0,
        n_clicks_timestamp=-1,
        disabled=not run_id,
        title=(
            "Click again within 3 seconds to permanently delete this profile"
            if confirming
            else "Permanently delete this stored profile"
        ),
        className=(
            "plots-output-dir-delete plots-output-dir-delete--confirm"
            if confirming
            else "plots-output-dir-delete"
        ),
    )


def available_profile_buttons(
    records,
    selected_ids=None,
    expanded: bool = False,
    delete_confirmation=None,
):
    """Render the three newest or all currently unselected profiles."""
    selected = {str(run_id) for run_id in (selected_ids or [])}
    available = [
        record
        for record in (records or [])
        if str(record.get("run_id") or "") and str(record.get("run_id")) not in selected
    ]
    if not available:
        return [html.Div("All profiles selected", className="profile-available-empty")]
    visible = available if expanded else available[:3]
    buttons = []
    for record in visible:
        add_button = html.Button(
            profile_option(record)["label"],
            id={"type": "profile-add-run", "run_id": str(record["run_id"])},
            n_clicks=0,
            n_clicks_timestamp=-1,
            title=f"Add {_profile_record_title(record)}",
            className="profile-available-button",
        )
        buttons.append(
            html.Div(
                [add_button, _profile_delete_button(record, delete_confirmation)],
                className="profile-available-row",
            )
            if expanded
            else add_button
        )
    return buttons


def active_profile_items(records, selected_ids, expanded: bool = False, delete_confirmation=None):
    """Render selected profiles as compact, removable comparison pills."""
    by_id = {str(record.get("run_id") or ""): record for record in (records or [])}
    items = []
    for raw_run_id in selected_ids or []:
        run_id = str(raw_run_id or "")
        record = by_id.get(run_id)
        if not run_id or record is None:
            continue
        label = profile_option(record)["label"]
        items.append(
            html.Div(
                [
                    html.Span(label, className="profile-active-label"),
                    *(
                        [_profile_delete_button(record, delete_confirmation)]
                        if expanded
                        else []
                    ),
                    html.Button(
                        "×",
                        id={"type": "profile-remove-run", "run_id": run_id},
                        n_clicks=0,
                        n_clicks_timestamp=-1,
                        title=f"Remove {label} from comparison",
                        **{"aria-label": f"Remove {label} from comparison"},
                        className="profile-active-remove",
                    ),
                ],
                title=_profile_record_title(record),
                className="profile-active-item",
            )
        )
    return items or [html.Div("No profiles selected", className="profile-active-empty")]


def _chart_card(component_id: str, title: str, message: str, controls=None):
    body = [
        dcc.Graph(
            id=component_id,
            figure=empty_profile_figure(message),
            config={"displaylogo": False, "responsive": True},
            className="profile-chart-graph",
        )
    ]
    if controls:
        body.append(html.Aside(controls, className="profile-chart-controls"))
    return html.Section(
        [
            html.Div(title, className="profile-chart-heading profile-chart-title"),
            html.Div(body, className="profile-chart-body"),
        ],
        id=f"{component_id}-card",
        className="profile-chart-card",
    )


def _benchmark_controls(initial_state: dict[str, object]):
    cases = list(initial_state.get("cases") or [])
    configs = [
        {"label": item.get("label", item.get("value")), "value": item.get("value")}
        for item in (initial_state.get("configs") or [])
    ]
    runner_options = [
        _field(
            "Tunable configuration",
            html.Div(
                [
                    dcc.Input(
                        id="profile-config",
                        type="text",
                        value=initial_state.get("default_config") or "default",
                        list="profile-config-options",
                        debounce=True,
                        className="profile-input",
                    ),
                    html.Datalist(
                        id="profile-config-options",
                        children=[
                            html.Option(value=option.get("value"), label=option.get("label"))
                            for option in configs
                        ],
                    ),
                ]
            ),
            tooltip="Value forwarded with -config. It chooses a named parameter_and_flag_configs directory for the benchmark.",
            compact=True,
        ),
        _field(
            "Executable",
            html.Div(
                [
                    dcc.Input(
                        id="profile-executable",
                        type="text",
                        value="",
                        list="profile-executable-options",
                        placeholder="Auto, or enter an executable path",
                        debounce=True,
                        className="profile-input",
                    ),
                    html.Datalist(
                        id="profile-executable-options",
                        children=[
                            html.Option(value=option.get("value"), label=option.get("label"))
                            for option in (initial_state.get("executables") or [])
                            if option.get("value")
                        ],
                    ),
                ]
            ),
            tooltip="Executable forwarded with -exe. Leave blank to use run_scm.py's selected or latest CLUBB installation.",
            help_text="Auto uses run_scm.py's selected/latest install.",
            compact=True,
        ),
        _field(
            "Namelist overrides",
            dcc.Input(
                id="profile-override",
                type="text",
                value="",
                placeholder="C2=2.0,l_predict_upwp_vpwp=false",
                debounce=True,
                className="profile-input",
            ),
            tooltip="Comma-separated namelist assignments forwarded with -override and applied by run_scm.py.",
            compact=True,
        ),
        _field(
            "Additional run_scm.py arguments",
            dcc.Input(
                id="profile-extra-args",
                type="text",
                value="",
                placeholder="-max_iters 100 -dt_main 60",
                debounce=True,
                className="profile-input",
            ),
            tooltip="Additional current or future run_scm.py options. Shell-style quoting is parsed, but the command is launched without a shell. Workload and output options are managed by the benchmark.",
            help_text="Shell-style quoting is supported; arguments are not run through a shell.",
            compact=True,
        ),
    ]
    return html.Section(
        [
            html.Div(
                [
                    _benchmark_group(
                        "01",
                        "Profile identity",
                        "Choose what to run and how this result will be named.",
                        [
                            _field(
                                "Case",
                                styled_dropdown(
                                    id="profile-case",
                                    options=[{"label": case, "value": case} for case in cases],
                                    value=initial_state.get("default_case"),
                                    clearable=False,
                                    searchable=True,
                                    className="profile-dropdown",
                                ),
                                tooltip="Choose the CLUBB SCM case to benchmark. Case settings determine the simulated time, timestep, and vertical grid unless overridden through run_scm.py arguments.",
                                compact=True,
                            ),
                            _field(
                                "Benchmark label",
                                html.Div(
                                    [
                                        dcc.Input(
                                            id="profile-name",
                                            type="text",
                                            value="",
                                            placeholder="Defaults to the case name",
                                            debounce=False,
                                            className="profile-input",
                                        ),
                                        html.Div(
                                            id="profile-name-status",
                                            className="profile-name-status",
                                            **{"aria-live": "polite"},
                                        ),
                                    ],
                                    id="profile-name-control",
                                    className="profile-name-control",
                                ),
                                tooltip="Choose the saved profile name. Reusing an existing name replaces that profile only after you confirm the overwrite.",
                                compact=True,
                            ),
                        ],
                        class_name="profile-benchmark-identity",
                    ),
                    _benchmark_group(
                        "02",
                        "Workload matrix",
                        "Sweep process counts against per-process CLUBB batch sizes.",
                        [
                            _field(
                                "Processes",
                                dcc.Input(
                                    id="profile-processes",
                                    type="text",
                                    value="1",
                                    debounce=True,
                                    className="profile-input",
                                ),
                                tooltip="Comma-separated counts of independent Python-launched CLUBB processes. Every process in a group runs concurrently.",
                                help_text="Example: 1,2,4,8",
                                compact=True,
                            ),
                            _field(
                                "Batch size / process",
                                dcc.Input(
                                    id="profile-columns",
                                    type="text",
                                    value="1,2,4,8",
                                    debounce=True,
                                    className="profile-input",
                                ),
                                tooltip="Comma-separated CLUBB runtime batch sizes assigned to each process. Overall batch size equals processes multiplied by this value.",
                                help_text="Overall batch = processes × this value.",
                                compact=True,
                            ),
                        ],
                        class_name="profile-benchmark-workload",
                    ),
                    _benchmark_group(
                        "03",
                        "Sampling & resilience",
                        "Control repeatability and what happens when a workload fails.",
                        [
                            html.Div(
                                [
                                    _field(
                                        "Warmups",
                                        _number_input("profile-warmups", 1, minimum=0),
                                        tooltip="Number of warm-up executions at each workload. Warmup timer rows are stored but ignored by the default plots.",
                                        compact=True,
                                    ),
                                    _field(
                                        "Repetitions",
                                        _number_input("profile-repeats", 3, minimum=1),
                                        tooltip="Number of measured executions at each process/batch-size point. Repetitions provide variability and error-bar information.",
                                        compact=True,
                                    ),
                                    _field(
                                        "Timeout",
                                        _number_input("profile-timeout", None, minimum=0.001, step="any"),
                                        tooltip="Maximum wall time in seconds allowed for an entire concurrent process group. Leave blank to wait without a timeout.",
                                        help_text="Seconds; blank waits indefinitely.",
                                        compact=True,
                                    ),
                                ],
                                className="profile-benchmark-sampling-grid",
                            ),
                            _field(
                                "On error",
                                dcc.Checklist(
                                    id="profile-continue-on-error",
                                    options=[{"label": "Keep the sweep running after a failed workload", "value": "continue"}],
                                    value=[],
                                    className="profile-checklist",
                                ),
                                tooltip="When enabled, preserve failed points and continue with the remaining workload matrix instead of stopping at the first error.",
                                compact=True,
                                class_name="profile-benchmark-errors",
                            ),
                        ],
                        class_name="profile-benchmark-sampling",
                    ),
                    _benchmark_group(
                        "04",
                        "Choose compiler & run",
                        "Select the CLUBB build, then launch or stop the timing sweep.",
                        [
                            build_selector_trigger("profile-selected-build-badge"),
                            html.Div(
                                [
                                    html.Button(
                                        "Run benchmark",
                                        id="profile-start",
                                        n_clicks=0,
                                        className="profile-button profile-button-start",
                                        title="Launch the configured timing sweep through the durable dashboard broker.",
                                    ),
                                    html.Button(
                                        "Cancel",
                                        id="profile-cancel",
                                        n_clicks=0,
                                        disabled=True,
                                        className="profile-button profile-button-cancel",
                                        title="Gracefully interrupt the active timing sweep and preserve all completed CSV rows and artifacts.",
                                    ),
                                ],
                                className="profile-benchmark-actions",
                            ),
                            html.Div(
                                [
                                    html.Div(
                                        [
                                            html.Span("Results", className="profile-results-label"),
                                            html.Div(
                                                "Idle",
                                                id="profile-status",
                                                className="profile-status profile-status-idle",
                                            ),
                                        ],
                                        className="profile-results-state",
                                    ),
                                    html.Div(id="profile-result-summary", className="profile-result-summary"),
                                ],
                                className="profile-status-row",
                            ),
                        ],
                        class_name="profile-benchmark-launch",
                    ),
                ],
                className="profile-benchmark-grid",
            ),
            html.Details(
                [
                    html.Summary(
                        "Advanced runner options",
                        title="Show configuration, executable, namelist overrides, and additional run_scm.py arguments.",
                    ),
                    html.P(
                        "Configuration, executable, and overrides forwarded to run_scripts/run_scm.py.",
                        className="profile-benchmark-runner-copy",
                    ),
                    html.Div(runner_options, className="profile-benchmark-runner-grid"),
                    html.Div(
                        [
                            html.Div(
                                [
                                    html.Span("Command"),
                                    html.Span(
                                        "?",
                                        className="profile-help-icon",
                                        title="Exact shell-escaped command that the dashboard will launch. It is shown for reproducibility and inspection.",
                                    ),
                                ],
                                className="profile-command-label",
                            ),
                            html.Pre(id="profile-command-preview", className="profile-command-preview"),
                        ],
                        className="profile-benchmark-command",
                    ),
                ],
                className="profile-benchmark-runner",
            ),
        ],
        className="profile-benchmark-panel",
    )


def _comparison_controls():
    return [
        _field(
            "Timer",
            styled_dropdown(
                id="profile-result-timer",
                options=[{"label": "advance_clubb_to_end", "value": "advance_clubb_to_end"}],
                value="advance_clubb_to_end",
                clearable=False,
                searchable=True,
                className="profile-dropdown",
            ),
            tooltip="Select the named CLUBB timer used by performance, baseline, and variability plots. It also identifies the critical process for cost decomposition.",
        ),
        _field(
            "Metric",
            styled_dropdown(
                id="profile-result-metric",
                options=[{"label": label, "value": value} for value, label in METRICS.items()],
                value="throughput_columns_per_second",
                clearable=False,
                searchable=False,
                className="profile-dropdown",
            ),
            tooltip="Choose the value shown by the performance and baseline comparison plots.",
        ),
        html.Div(
            [
                _field(
                    "Horizontal axis",
                    styled_dropdown(
                        id="profile-x-axis",
                        options=[{"label": label, "value": value} for value, label in X_AXES.items()],
                        value="total_columns",
                        clearable=False,
                        searchable=False,
                        className="profile-dropdown",
                    ),
                    tooltip="Use overall batch size to compare the full concurrent workload, or batch size per process to compare each CLUBB instance.",
                ),
                _field(
                    "Axis scale",
                    styled_dropdown(
                        id="profile-x-scale",
                        options=[{"label": "Linear", "value": "linear"}, {"label": "Log", "value": "log"}],
                        value="linear",
                        clearable=False,
                        searchable=False,
                        className="profile-dropdown",
                    ),
                    tooltip="Display the workload axis using a linear or logarithmic scale on the performance and baseline plots.",
                ),
            ],
            className="profile-two-column",
        ),
        _field(
            "Detail profile",
            styled_dropdown(
                id="profile-detail-run",
                options=[],
                value=None,
                clearable=False,
                searchable=True,
                className="profile-dropdown",
            ),
            tooltip="Choose the single stored profile used by the cost-decomposition and process-variability plots.",
        ),
    ]


def build_controls():
    output_tooltip = (
        "Choose the directory containing compact profile folders. New benchmarks create "
        "one named profile folder here."
    )
    profiles_tooltip = (
        "Select one or more stored timing profiles from this library. Selected profiles are "
        "overlaid using stable colors and made available to every comparison plot."
    )
    return html.Div(
        [
            html.Section(
                [
                    html.Div(
                        [
                            _field(
                                "Profile library directory",
                                dcc.Input(
                                    id="profile-output",
                                    type="text",
                                    value="output/timing",
                                    debounce=True,
                                    className="profile-input",
                                ),
                                tooltip=output_tooltip,
                                compact=True,
                            ),
                            html.Button(
                                "Load library",
                                id="profile-library-refresh",
                                n_clicks=0,
                                className="profile-library-refresh-button",
                                title="Load or refresh compact profile directories from the specified library.",
                            ),
                        ],
                        className="profile-library-picker",
                        title=output_tooltip,
                    ),
                    html.Div(
                        [
                            html.Button(
                                [
                                    html.Span("Available profiles"),
                                    html.Span("▾", className="profile-picker-chevron"),
                                ],
                                id="profile-picker-toggle",
                                n_clicks=0,
                                n_clicks_timestamp=-1,
                                className="profile-picker-header",
                                title="Show or hide all unselected profiles in this library.",
                            ),
                            html.Div(
                                html.Div(id="profile-available-list", className="profile-available-list"),
                                className="profile-picker-body",
                            ),
                        ],
                        id="profile-picker-menu",
                        className="profile-picker-menu",
                        title=profiles_tooltip,
                    ),
                    html.Div(
                        [
                            html.Div(
                                [
                                    html.Span("Profiles to compare"),
                                    html.Span(
                                        "?",
                                        className="profile-help-icon",
                                        title=profiles_tooltip,
                                        **{"aria-label": profiles_tooltip},
                                    ),
                                ],
                                className="profile-active-heading",
                            ),
                            html.Div(id="profile-active-list", className="profile-active-list"),
                        ],
                        className="profile-active-tray",
                    ),
                    html.Div(
                        styled_dropdown(
                            id="profile-selected-runs",
                            options=[],
                            value=[],
                            multi=True,
                            clearable=True,
                            searchable=True,
                            className="profile-dropdown",
                        ),
                        className="profile-selection-state",
                    ),
                    html.Div(
                        [
                            html.Div(
                                dcc.Upload(
                                    id="profile-library-import",
                                    children=html.Button(
                                        "Import profile",
                                        className="profile-tool-button",
                                        title="Import a complete compact profile ZIP into the active profile library.",
                                    ),
                                    accept=".zip,application/zip",
                                    multiple=False,
                                ),
                                title="Import a complete compact profile ZIP into the active profile library.",
                            ),
                            html.Button(
                                "Export selected",
                                id="profile-library-export",
                                n_clicks=0,
                                className="profile-tool-button",
                                title="Download the selected profiles as a portable ZIP containing manifests, CSVs, native timings, and representative logs.",
                            ),
                            dcc.Download(id="profile-library-download"),
                        ],
                        className="profile-library-actions profile-library-actions-two",
                    ),
                    html.Div(id="profile-library-status", className="profile-library-status"),
                    html.Div(id="profile-export-status", className="profile-library-status"),
                ],
                className="profile-sidebar-section profile-output-section",
            ),
            html.Section(_comparison_controls(), className="profile-sidebar-section profile-view-section"),
        ],
        className="profile-controls-panel",
    )


def _results_panel():
    relative_controls = [
        _section_header(
            "PLOT OPTIONS",
            "Baseline",
            "Relative changes require an exact workload match.",
        ),
        _field(
            "Baseline profile",
            styled_dropdown(
                id="profile-baseline-run",
                options=[],
                value=None,
                clearable=True,
                searchable=True,
                placeholder="Choose baseline",
                className="profile-dropdown",
            ),
            tooltip="Choose the reference profile. Each comparison point is shown as its percent change from the same process count and per-process batch-size point in this profile.",
        ),
    ]
    decomposition_controls = [
        _section_header(
            "PLOT OPTIONS",
            "Cost detail",
            "Small components are grouped into Other.",
        ),
        _field(
            "Top cost timers",
            _number_input("profile-top-timers", 8, minimum=1, maximum=20),
            tooltip="Number of largest exclusive timer components to show individually. Remaining exclusive cost is grouped into Other.",
        ),
        _field(
            "Normalization",
            dcc.Checklist(
                id="profile-cost-percentage",
                options=[{"label": "Show percentage of total", "value": "percentage"}],
                value=[],
                className="profile-checklist",
            ),
            tooltip="Normalize every workload stack to 100 percent instead of showing absolute exclusive timer seconds.",
        ),
    ]
    return html.Div(
        [
            html.Div(
                [
                    _chart_card(
                        "profile-graph",
                        "Performance",
                        "Select or run a profile to see timing results.",
                    ),
                    _chart_card(
                        "profile-relative-graph",
                        "Baseline comparison",
                        "Select a baseline and comparison profile.",
                        relative_controls,
                    ),
                    _chart_card(
                        "profile-decomposition-graph",
                        "Cost decomposition",
                        "Select a detail profile.",
                        decomposition_controls,
                    ),
                    _chart_card(
                        "profile-variability-graph",
                        "Process variability",
                        "Select a detail profile.",
                    ),
                ],
                className="profile-plot-grid",
            ),
            html.Details(
                [
                    html.Summary(
                        "Active timing log",
                        title="Show or hide stdout and stderr from the currently running or most recent timing command.",
                    ),
                    html.Pre("No profile runs yet.", id="profile-log", className="profile-log"),
                ],
                className="profile-log-panel",
            ),
        ],
        className="profile-results-panel",
    )


def _overwrite_modal():
    return html.Div(
        html.Div(
            [
                html.Div(
                    [
                        html.Div("!", className="profile-overwrite-icon", **{"aria-hidden": "true"}),
                        html.Div(
                            [
                                html.Div("Benchmark label already exists", className="profile-overwrite-title"),
                                html.Div(
                                    id="profile-overwrite-message",
                                    className="profile-overwrite-message",
                                ),
                            ]
                        ),
                    ],
                    className="profile-overwrite-heading",
                ),
                html.Label("Benchmark label", htmlFor="profile-overwrite-name", className="profile-overwrite-label"),
                dcc.Input(
                    id="profile-overwrite-name",
                    type="text",
                    value="",
                    debounce=False,
                    className="profile-input profile-overwrite-input",
                ),
                html.Div(
                    [
                        html.Button(
                            "Overwrite",
                            id="profile-overwrite-button",
                            type="button",
                            n_clicks=0,
                            className="profile-overwrite-button profile-overwrite-button-danger",
                        ),
                        html.Button(
                            "Rename",
                            id="profile-rename-button",
                            type="button",
                            n_clicks=0,
                            disabled=True,
                            className="profile-overwrite-button profile-overwrite-button-primary",
                            title="Enter a different benchmark label to rename and run.",
                        ),
                        html.Button(
                            "Cancel",
                            id="profile-overwrite-cancel-button",
                            type="button",
                            n_clicks=0,
                            className="profile-overwrite-button profile-overwrite-button-cancel",
                        ),
                    ],
                    className="profile-overwrite-actions",
                ),
            ],
            className="profile-overwrite-panel",
            role="dialog",
            **{"aria-modal": "true", "aria-labelledby": "profile-overwrite-message"},
        ),
        id="profile-overwrite-modal",
        className="profile-overwrite-modal profile-overwrite-modal-hidden",
    )


def build_layout(initial_state: dict[str, object]):
    return html.Div(
        [
            dcc.Store(id="profile-action", data={}),
            dcc.Store(id="profile-job", data={}),
            dcc.Store(id="profile-active-results", data={}),
            dcc.Store(id="profile-catalog", data=[]),
            dcc.Store(id="profile-results", data=[]),
            dcc.Store(id="profile-library-action", data={}),
            dcc.Store(id="profile-selection-library", data=""),
            dcc.Store(id="profile-picker-expanded", data=False),
            dcc.Store(id="profile-delete-confirm", data=None),
            dcc.Store(id="profile-pending-run", data={}),
            _overwrite_modal(),
            dcc.Interval(id="profile-interval", interval=1000, n_intervals=0, disabled=True),
            dcc.Interval(
                id="profile-delete-expiry",
                interval=250,
                n_intervals=0,
                disabled=True,
            ),
            _benchmark_controls(initial_state),
            html.Div(
                [_results_panel(), build_controls()],
                className="profile-workspace",
            ),
        ],
        id="profile-tab-layout",
        className="profile-tab-layout",
    )
