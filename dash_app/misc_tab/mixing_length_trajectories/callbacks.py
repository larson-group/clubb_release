"""Callbacks for the mixing-length parcel trajectory explorer."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

from dash import Input, Output, State, ctx, html, no_update

from .analysis import (
    compute_record,
    discover_netcdf_files,
    inspect_dataset,
    load_dataset_record,
    profile_metrics,
)
from .figures import (
    buoyancy_trajectory_figure,
    length_profile_figure,
    parcel_state_figure,
    trajectory_figure,
)

REPO_ROOT = Path(__file__).resolve().parents[3]


def _options(paths):
    options = []
    for path in paths:
        try:
            label = str(path.relative_to(REPO_ROOT))
        except ValueError:
            label = str(path)
        options.append({"label": label, "value": str(path)})
    return options


def _slider_marks(times):
    count = len(times)
    if count <= 1:
        return {0: "0"}
    indices = sorted(
        {0, count - 1, *(round(i * (count - 1) / 5) for i in range(1, 5))}
    )
    return {
        index: f"{index}\n({times[index] / 60.0:.0f}m)"
        for index in indices
    }


def _metric_card(label, metrics, color):
    return html.Div(
        [
            html.Div(label, className="mlt-metric-label"),
            html.Div(
                f"{metrics['rmse']:.3g} m",
                className="mlt-metric-value",
                style={"color": color},
            ),
            html.Div(
                f"max |Δ| {metrics['max_abs']:.3g} m · "
                f"r {metrics['correlation']:.6f}",
                className="mlt-metric-detail",
            ),
        ],
        className="mlt-metric-card",
    )


def register_callbacks(app):
    @app.callback(
        Output("mlt-file", "options"),
        Output("mlt-file", "value"),
        Input("mlt-load-path", "n_clicks"),
        Input("mlt-refresh-files", "n_clicks"),
        State("mlt-custom-path", "value"),
        State("mlt-file", "value"),
        prevent_initial_call=True,
    )
    def update_file_choices(_load_clicks, _refresh_clicks, custom_path, current):
        paths = discover_netcdf_files(REPO_ROOT)
        selected = current
        if ctx.triggered_id == "mlt-load-path" and custom_path:
            candidate = Path(custom_path).expanduser().resolve()
            inspect_dataset(candidate)
            if candidate not in paths:
                paths.append(candidate)
            selected = str(candidate)
        elif current not in {str(path) for path in paths}:
            selected = str(paths[0]) if paths else None
        return _options(sorted(paths)), selected

    @app.callback(
        Output("mlt-record", "max"),
        Output("mlt-record", "marks"),
        Output("mlt-record", "value"),
        Output("mlt-column", "options"),
        Output("mlt-column", "value"),
        Output("mlt-file-status", "children"),
        Input("mlt-file", "value"),
    )
    def update_file_metadata(path):
        if not path:
            return 1, {0: "0", 1: "1"}, 0, [], None, "Choose a file."
        try:
            metadata = inspect_dataset(path)
        except Exception as error:
            return (
                1,
                {0: "0", 1: "1"},
                0,
                [],
                None,
                html.Span(str(error), className="mlt-error"),
            )
        default_record = min(612, metadata["record_count"] - 1)
        columns = [
            {"label": str(index), "value": index}
            for index in range(metadata["column_count"])
        ]
        status = (
            f"{metadata['path']} · {metadata['record_count']} records · "
            f"{metadata['column_count']} column(s)"
        )
        return (
            metadata["record_count"] - 1,
            _slider_marks(metadata["times"]),
            default_record,
            columns,
            0,
            status,
        )

    @app.callback(
        Output("mlt-mu", "max"),
        Output("mlt-mu", "marks"),
        Output("mlt-mu", "value"),
        Input("mlt-file", "value"),
        Input("mlt-column", "value"),
    )
    def update_mu_control(path, column):
        if path is None or column is None:
            return 3.0e-3, {}, 1.0e-3
        try:
            record = load_dataset_record(path, 0, column)
        except (IndexError, OSError, RuntimeError, TypeError, ValueError):
            # File metadata owns the visible validation message. This callback
            # can race that validation when a newly-created NetCDF is empty.
            return 3.0e-3, {}, 1.0e-3
        slider_max = max(3.0e-3, 2.0 * record.mu)
        marks = {
            1.0e-5: "0.00001",
            record.mu: f"{record.mu:.4g} file",
            slider_max: f"{slider_max:.4g}",
        }
        return slider_max, marks, record.mu

    @app.callback(
        Output("mlt-mu-label", "children"),
        Input("mlt-mu", "value"),
    )
    def update_mu_label(mu):
        return "" if mu is None else f"{float(mu):.6g} m⁻¹"

    @app.callback(
        Output("mlt-upward-figure", "figure"),
        Output("mlt-downward-figure", "figure"),
        Output("mlt-upward-buoyancy-figure", "figure"),
        Output("mlt-downward-buoyancy-figure", "figure"),
        Output("mlt-parcel-thv-figure", "figure"),
        Output("mlt-parcel-rt-figure", "figure"),
        Output("mlt-parcel-thl-figure", "figure"),
        Output("mlt-profile-figure", "figure"),
        Output("mlt-metrics", "children"),
        Output("mlt-time-label", "children"),
        Input("mlt-file", "value"),
        Input("mlt-record", "value"),
        Input("mlt-column", "value"),
        Input("mlt-mu", "value"),
    )
    def update_diagnostics(path, record_index, column, mu):
        if (
            path is None
            or record_index is None
            or column is None
            or mu is None
        ):
            return (
                no_update,
                no_update,
                no_update,
                no_update,
                no_update,
                no_update,
                no_update,
                no_update,
                [],
                "",
            )
        try:
            record = load_dataset_record(path, record_index, column)
        except (IndexError, OSError, RuntimeError, TypeError, ValueError) as error:
            message = f"Could not load this record: {error}"
            return (
                no_update,
                no_update,
                no_update,
                no_update,
                no_update,
                no_update,
                no_update,
                no_update,
                html.Span(message, className="mlt-error"),
                message,
            )
        file_mu = record.mu
        record = replace(record, mu=float(mu))
        result = compute_record(record)
        metrics_up = profile_metrics(result.lscale_up, record.fortran_up)
        metrics_down = profile_metrics(
            result.lscale_down, record.fortran_down
        )
        metrics_final = profile_metrics(
            result.lscale, record.fortran_lscale
        )
        metrics = [
            _metric_card("Upward RMSE", metrics_up, "#38bdf8"),
            _metric_card("Downward RMSE", metrics_down, "#a78bfa"),
            _metric_card("Final Lscale RMSE", metrics_final, "#34d399"),
            html.Div(
                [
                    html.Div("Calculation inputs", className="mlt-metric-label"),
                    html.Div(
                        f"μ={record.mu:.4g} m⁻¹ · lmin={record.lmin:.3g} m",
                        className="mlt-input-value",
                    ),
                    html.Div(
                        (
                            "Using the NetCDF clubb_params value."
                            if abs(record.mu - file_mu)
                            <= 1.0e-12 * max(abs(file_mu), 1.0)
                            else f"Slider override; file μ={file_mu:.4g} m⁻¹."
                        ),
                        className="mlt-metric-detail",
                    ),
                ],
                className="mlt-metric-card",
            ),
        ]
        time_label = (
            f"record {record.record} · model time "
            f"{record.time_seconds / 60.0:.2f} min · column {record.column}"
        )
        return (
            trajectory_figure(result, "up"),
            trajectory_figure(result, "down"),
            buoyancy_trajectory_figure(result, "up"),
            buoyancy_trajectory_figure(result, "down"),
            parcel_state_figure(record, result, "thv"),
            parcel_state_figure(record, result, "rt"),
            parcel_state_figure(record, result, "thl"),
            length_profile_figure(record, result),
            metrics,
            time_label,
        )
