"""Callbacks for rendering run logs and console panels."""

import time

from dash import Input, Output, State, dcc, html

from .layout import case_dom_id
from .runtime import extract_progress, format_eta, format_runtime


def register_console_callbacks(app, cases):
    """Register callbacks that render the live per-case console panels."""

    @app.callback(
        Output("run-console-container", "children"),
        Output("run-case-order", "data"),
        Input("run-case-logs", "data"),
        Input("run-selected-cases", "data"),
        Input("run-completed-cases", "data"),
        Input("run-failed-cases", "data"),
        Input("run-running-cases", "data"),
        Input("run-case-runtimes", "data"),
        Input("run-case-commands", "data"),
        State("run-case-order", "data"),
    )
    def render_consoles(case_logs, selected_cases, completed_cases, failed_cases, running_cases, case_runtimes, case_commands, case_order):
        """Render log consoles in stable order, keeping active runs pinned to the top."""
        logs = case_logs or {}
        runtimes = case_runtimes or {}
        commands = case_commands or {}
        running_data = running_cases or {}
        now = time.time()
        selected = set(selected_cases or [])
        completed = set(completed_cases or [])
        failed = set(failed_cases or [])
        running = set(running_data.keys())
        candidates = [
            case
            for case in cases
            if case in logs or case in selected or case in completed or case in failed or case in running
        ]
        if not candidates:
            return html.Div("No runs yet.", className="run-empty-message", style={"fontStyle": "italic", "padding": "8px 0"}), []

        candidate_set = set(candidates)
        prev_order = [name for name in (case_order or []) if name in candidate_set]
        for name in candidates:
            if name not in prev_order:
                prev_order.append(name)
        running_now = [name for name in prev_order if name in running]
        not_running_now = [name for name in prev_order if name not in running]
        display_cases = running_now + not_running_now

        panels = []
        for case_name in display_cases:
            if case_name in running:
                color = "linear-gradient(90deg, #dc2626, #16a34a)"
                status = "running"
            elif case_name in failed:
                color = "#dc2626"
                status = "failed"
            elif case_name in completed:
                color = "#16a34a"
                status = "completed"
            elif case_name in selected:
                color = "#f59e0b"
                status = "selected"
            else:
                color = "#2563eb"
                status = "pending"
            if case_name in running:
                status_style = {
                    "backgroundImage": color,
                    "WebkitBackgroundClip": "text",
                    "WebkitTextFillColor": "transparent",
                    "fontSize": "14px",
                }
            else:
                status_style = {"color": color, "fontSize": "14px"}

            if case_name in running:
                start_time = running_data.get(case_name, {}).get("start_time")
                runtime_txt = format_runtime(now - float(start_time)) if start_time is not None else ""
            else:
                runtime_txt = format_runtime(runtimes.get(case_name))

            progress = extract_progress(logs.get(case_name, ""))
            if progress is None:
                if case_name in completed or case_name in failed:
                    progress = 1.0
                else:
                    progress = 0.0

            eta_txt = ""
            if case_name in running and progress and progress > 0:
                start_time = running_data.get(case_name, {}).get("start_time")
                if start_time is not None:
                    elapsed = now - float(start_time)
                    remaining = max(0.0, elapsed * (1.0 / progress - 1.0))
                    eta_txt = format_eta(remaining)

            summary_children = [
                html.Span(case_name, style={"fontWeight": "600", "marginRight": "8px"}),
                html.Span(status, style=status_style),
            ]
            if runtime_txt:
                summary_children.append(
                    html.Span(
                        f"{'elapsed' if case_name in running else 'runtime'}: {runtime_txt}",
                        className="run-muted-text",
                        style={"marginLeft": "10px", "fontSize": "14px"},
                    )
                )

            progress_bar = html.Div(
                [
                    html.Span(f"eta {eta_txt}" if eta_txt else "eta --", className="run-muted-text", style={"fontSize": "14px"}),
                    html.Div(
                        [
                            html.Div(
                                className="run-progress-fill",
                                style={
                                    "width": f"{int(progress * 100)}%",
                                    "height": "16px",
                                    "background": "#16a34a",
                                    "borderRadius": "999px",
                                    "transition": "width 0.3s ease",
                                },
                            )
                        ],
                        className="run-progress-track",
                        style={"width": "240px", "height": "16px", "borderRadius": "999px", "overflow": "hidden"},
                    ),
                ],
                style={"display": "flex", "alignItems": "center", "gap": "8px", "marginLeft": "auto", "flexShrink": 0},
            )
            summary = html.Summary(
                [
                    html.Div(summary_children, className="run-summary-left", style={"display": "flex", "alignItems": "center", "flexWrap": "wrap", "gap": "6px"}),
                    progress_bar,
                ],
                className="run-console-summary",
                style={"cursor": "pointer", "display": "flex", "alignItems": "center", "gap": "12px", "fontSize": "14px"},
            )
            panel_children = [summary]
            command_text = commands.get(case_name)
            if command_text:
                panel_children.append(
                    html.Div(
                        [
                            html.Div(command_text, className="run-command-text"),
                            html.Span("Copy run command", className="run-copy-command-label"),
                            dcc.Clipboard(content=command_text, title="Copy run command", className="run-copy-command-button"),
                        ],
                        className="run-command-row",
                    )
                )
            panel_children.append(
                html.Pre(
                    logs.get(case_name, "Not run yet."),
                    id=f"run-console-{case_dom_id(case_name)}",
                    className="run-console run-console-active" if case_name in running else "run-console",
                    style={"padding": "12px", "borderRadius": "6px", "maxHeight": "300px", "overflowY": "auto", "marginTop": "6px"},
                )
            )
            panels.append(
                html.Details(
                    panel_children,
                    open=case_name in running,
                    className="run-console-panel",
                    style={"borderRadius": "6px", "padding": "8px 10px"},
                )
            )
        return panels, display_cases
