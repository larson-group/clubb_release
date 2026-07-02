"""Static layout builders for the compile tab."""

from __future__ import annotations

from dash import dcc, html


def action_button_style(color, disabled=False):
    """Return the common compile-tab action button style."""
    style = {
        "backgroundColor": color,
        "color": "#ffffff",
        "border": "none",
        "padding": "10px 16px",
        "borderRadius": "6px",
        "cursor": "pointer",
        "fontSize": "14px",
        "fontWeight": "600",
    }
    if disabled:
        style.update({"backgroundColor": "#9ca3af", "color": "#f3f4f6", "cursor": "not-allowed"})
    return style


def build_option_checklist():
    """Render build feature toggles."""
    return dcc.Checklist(
        id="compile-feature-flags",
        options=[
            {"label": "Debug", "value": "debug"},
            {"label": "Run CTest", "value": "run_tests"},
            {"label": "Python/F2PY", "value": "python"},
            {"label": "Fresh build", "value": "fresh"},
            {"label": "OpenMP", "value": "openmp"},
            {"label": "Tuning", "value": "tuning"},
            {"label": "GPTL", "value": "gptl"},
            {"label": "Disable NetCDF", "value": "disable_netcdf"},
            {"label": "Disable SILHS", "value": "disable_silhs"},
        ],
        value=[],
        className="compile-checklist",
        labelStyle={"display": "inline-flex", "alignItems": "center", "gap": "6px", "marginRight": "14px", "marginBottom": "8px"},
    )


def build_controls():
    """Render compile controls."""
    return html.Div(
        [
            html.Div(
                [
                    html.Label("Environment", className="compile-field-label"),
                    dcc.Dropdown(id="compile-env-select", clearable=False, className="clubb-dropdown compile-dropdown"),
                ],
                id="compile-native-env-field",
                className="compile-field",
            ),
            html.Div(
                [
                    html.Div("Module Stack", className="compile-section-title"),
                    html.Div(
                        [
                            html.Label("Compiler module", className="compile-field-label"),
                            dcc.Dropdown(
                                id="compile-lmod-compiler",
                                clearable=False,
                                searchable=True,
                                className="clubb-dropdown compile-dropdown",
                            ),
                        ],
                        className="compile-field",
                    ),
                    html.Div(
                        [
                            html.Label("NetCDF module", className="compile-field-label"),
                            dcc.Dropdown(
                                id="compile-lmod-netcdf",
                                clearable=True,
                                searchable=True,
                                placeholder="optional NetCDF module",
                                className="clubb-dropdown compile-dropdown",
                            ),
                        ],
                        className="compile-field compile-field-inline",
                    ),
                    html.Div(
                        [
                            html.Label("Extra modules", className="compile-field-label"),
                            dcc.Dropdown(
                                id="compile-lmod-extra",
                                multi=True,
                                clearable=True,
                                searchable=True,
                                placeholder="optional extra modules",
                                className="clubb-dropdown compile-dropdown",
                            ),
                        ],
                        className="compile-field compile-field-inline",
                    ),
                    html.Div(id="compile-lmod-stack-summary", className="compile-lmod-stack-summary"),
                ],
                id="compile-lmod-panel",
                className="compile-lmod-panel",
            ),
            html.Div(
                [
                    html.Label("Toolchain", className="compile-field-label"),
                    dcc.Dropdown(id="compile-toolchain-select", clearable=False, className="clubb-dropdown compile-dropdown"),
                ],
                className="compile-field",
            ),
            html.Div(
                [
                    html.Label("Precision", className="compile-field-label"),
                    dcc.Dropdown(
                        id="compile-precision-select",
                        clearable=False,
                        value="double",
                        options=[
                            {"label": "double", "value": "double"},
                            {"label": "single", "value": "single"},
                            {"label": "quad", "value": "quad"},
                        ],
                        className="clubb-dropdown compile-dropdown",
                    ),
                ],
                className="compile-field compile-field-inline",
            ),
            html.Div(
                [
                    html.Label("GPU", className="compile-field-label"),
                    dcc.Dropdown(
                        id="compile-gpu-select",
                        clearable=False,
                        value="none",
                        options=[
                            {"label": "none", "value": "none"},
                            {"label": "OpenACC", "value": "openacc"},
                            {"label": "OpenMP target", "value": "openmp"},
                        ],
                        className="clubb-dropdown compile-dropdown",
                    ),
                ],
                className="compile-field compile-field-inline",
            ),
            html.Div([html.Div("Options", className="compile-section-title"), build_option_checklist()], className="compile-field"),
            html.Div(
                [
                    html.Label("Extra CMake args", className="compile-field-label"),
                    dcc.Input(id="compile-extra-args", type="text", value="", placeholder="-DSOME_OPTION=VALUE", className="clubb-input compile-text-input"),
                ],
                className="compile-field compile-extra-args-field",
            ),
            html.Div(id="compile-warnings", className="compile-warnings"),
            html.Div(
                [
                    html.Button("Compile", id="compile-start", n_clicks=0, className="compile-button-start", style=action_button_style("#16a34a")),
                    html.Button("Cancel", id="compile-cancel", n_clicks=0, style=action_button_style("#b91c1c")),
                    html.Button("Clear", id="compile-clear", n_clicks=0, style=action_button_style("#374151")),
                ],
                className="compile-action-row",
            ),
            html.Div(
                [
                    html.Div(id="compile-command-preview", className="compile-command-text"),
                    dcc.Clipboard(id="compile-command-copy", title="Copy compile command", className="compile-copy-command-button"),
                ],
                className="compile-command-row",
            ),
        ],
        className="compile-controls",
    )


def build_status_panel():
    """Render the right-hand discovery/build status panel."""
    return html.Div(
        [
            html.Div("Detected", className="compile-section-heading"),
            html.Div(id="compile-detection-summary", className="compile-detection-summary"),
            html.Div("Existing Builds", className="compile-section-heading"),
            html.Div(
                [
                    html.Button("Recompile All", id="compile-rebuild-all", n_clicks=0, className="compile-build-toolbar-button"),
                    html.Button("Refresh", id="compile-refresh", n_clicks=0, className="compile-build-toolbar-button"),
                    html.Span(id="compile-build-status-summary", className="compile-muted"),
                ],
                className="compile-build-toolbar",
            ),
            html.Div(id="compile-build-select-message", className="compile-build-select-message"),
            html.Div(id="compile-build-list", className="compile-build-list"),
        ],
        className="compile-status-panel",
    )


def build_console():
    """Render compile output console."""
    return html.Div(
        [
            html.Div(id="compile-job-summary", className="compile-job-summary"),
            html.Pre(id="compile-console", className="compile-console", children="No compile jobs yet."),
        ],
        className="compile-console-panel",
    )


def build_layout(initial_data):
    """Assemble the compile tab layout."""
    return html.Div(
        [
            dcc.Store(id="compile-discovery", data=initial_data),
            dcc.Store(id="compile-build-statuses", data={}),
            dcc.Store(id="compile-build-failures", data={}),
            dcc.Store(id="compile-build-delete-target", data=None),
            dcc.Store(id="compile-job", data={}),
            dcc.Store(id="compile-log", data=""),
            dcc.Store(id="compile-log-offset", data=0),
            dcc.Interval(id="compile-interval", interval=500, disabled=True),
            dcc.Interval(id="compile-build-status-interval", interval=30000),
            html.Div(
                [build_controls(), build_console()],
                className="compile-left-pane",
            ),
            build_status_panel(),
        ],
        id="compile-tab-layout",
        className="compile-tab-layout",
    )
