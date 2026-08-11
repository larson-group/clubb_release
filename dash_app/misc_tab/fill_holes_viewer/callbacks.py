"""Callbacks for the fill-holes viewer."""

from __future__ import annotations

from pathlib import Path

from dash import Input, Output, State, html, no_update

from dash_app.services.profiles import discover_output_directories

from .analysis import MOMENTS, inspect_dataset, load_moment
from .figures import empty_figure, holes_over_time_figure, profile_figure, summary_figure


REPO_ROOT = Path(__file__).resolve().parents[3]


def _folder_options(records):
    return [
        {"label": record["label"], "value": record["path"]}
        for record in records
        if record.get("available")
    ]


def _selected_record(records, folder):
    return next((record for record in records or [] if record.get("path") == folder), None)


def _case_options(record):
    return [
        {"label": name, "value": path}
        for name, path in (record or {}).get("cases", {}).items()
    ]


def _slider_marks(times):
    count = len(times)
    if count <= 1:
        return {0: "0"}
    indices = sorted({0, count - 1, *(round(i * (count - 1) / 4) for i in range(1, 4))})
    return {index: f"{times[index] / 60.0:.0f}m" for index in indices}


def register_callbacks(app):
    @app.callback(
        Output("hf-output-catalog", "data"),
        Input("hf-refresh-outputs", "n_clicks"),
        prevent_initial_call=True,
    )
    def refresh_outputs(_clicks):
        return discover_output_directories(REPO_ROOT / "output")

    @app.callback(
        Output("hf-output-folder", "options"),
        Output("hf-output-folder", "value"),
        Output("hf-case-file", "options"),
        Output("hf-case-file", "value"),
        Input("hf-output-catalog", "data"),
        Input("hf-output-folder", "value"),
        State("hf-case-file", "value"),
    )
    def update_sources(records, folder, current_case):
        records = records or []
        valid_folders = {record["path"] for record in records if record.get("available")}
        selected_folder = folder if folder in valid_folders else next(
            (record["path"] for record in records if record.get("available")),
            None,
        )
        record = _selected_record(records, selected_folder)
        case_options = _case_options(record)
        valid_cases = {option["value"] for option in case_options}
        selected_case = current_case if current_case in valid_cases else (
            case_options[0]["value"] if case_options else None
        )
        return _folder_options(records), selected_folder, case_options, selected_case

    @app.callback(
        Output("hf-moment", "options"),
        Output("hf-moment", "value"),
        Output("hf-column", "options"),
        Output("hf-column", "value"),
        Output("hf-record", "max"),
        Output("hf-record", "marks"),
        Output("hf-record", "value"),
        Output("hf-file-status", "children"),
        Input("hf-case-file", "value"),
        State("hf-moment", "value"),
    )
    def update_file_metadata(path, current_moment):
        if not path:
            return [], None, [], None, 1, {0: "0", 1: "1"}, 0, "Choose an output folder and case."
        try:
            metadata = inspect_dataset(path)
        except Exception as error:
            return (
                [], None, [], None, 1, {0: "0", 1: "1"}, 0,
                html.Span(str(error), className="hf-error"),
            )
        moments = metadata["moments"]
        moment_options = [
            {"label": f"{MOMENTS[name]['label']} ({name})", "value": name}
            for name in moments
        ]
        selected_moment = current_moment if current_moment in moments else moments[0]
        columns = [{"label": str(index), "value": index} for index in range(metadata["column_count"])]
        status = (
            f"{metadata['path']} · {metadata['record_count']} records · "
            f"{len(moments)} hole-filled CLUBB moments"
        )
        return (
            moment_options,
            selected_moment,
            columns,
            0,
            metadata["record_count"] - 1,
            _slider_marks(metadata["times"]),
            0,
            status,
        )

    @app.callback(
        Output("hf-playback-timer", "disabled"),
        Output("hf-play", "children"),
        Input("hf-play", "n_clicks"),
    )
    def toggle_playback(clicks):
        playing = bool((clicks or 0) % 2)
        return not playing, "Pause" if playing else "Play"

    @app.callback(
        Output("hf-playback-timer", "interval"),
        Input("hf-speed", "value"),
    )
    def update_speed(milliseconds):
        return int(milliseconds or 500)

    @app.callback(
        Output("hf-record", "value", allow_duplicate=True),
        Input("hf-playback-timer", "n_intervals"),
        State("hf-record", "value"),
        State("hf-record", "max"),
        prevent_initial_call=True,
    )
    def advance_record(_ticks, record, maximum):
        if record is None or maximum is None:
            return no_update
        return (int(record) + 1) % (int(maximum) + 1)

    @app.callback(
        Output("hf-profile-figure", "figure"),
        Output("hf-map-figure", "figure"),
        Output("hf-summary-figure", "figure"),
        Output("hf-time-label", "children"),
        Input("hf-case-file", "value"),
        Input("hf-moment", "value"),
        Input("hf-record", "value"),
        Input("hf-column", "value"),
        Input("hf-field-scale", "value"),
    )
    def update_figures(path, moment, record, column, field_scale):
        if path is None or moment is None or record is None or column is None:
            blank = empty_figure("Choose a stats file with dedicated hole-filling diagnostics.")
            return blank, blank, blank, ""
        try:
            metadata = inspect_dataset(path)
            data = load_moment(path, moment)
            record_index = min(max(int(record), 0), len(data["times"]) - 1)
            column_index = min(max(int(column), 0), data["before"].shape[2] - 1)
            return (
                profile_figure(
                    data,
                    moment,
                    record_index,
                    column_index,
                    field_scale,
                ),
                holes_over_time_figure(data, moment, column_index),
                summary_figure(path, metadata["moments"], column_index),
                f"record {record_index} · model time {data['times'][record_index] / 60.0:.2f} min",
            )
        except Exception as error:
            blank = empty_figure(str(error))
            return blank, blank, blank, ""
