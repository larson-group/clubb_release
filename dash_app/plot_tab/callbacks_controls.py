from dash import ALL, Input, Output, State, callback_context, no_update

from .plot_types.budget_plot import PLOT as budget_plot
from .plot_types.profile_plot import PLOT as profile_plot
from .plot_types.shared import average_length_label, physical_time_window_label, start_time_label, time_start_max_for_duration
from .plot_types.subcolumn_plot import PLOT as subcolumn_plot
from .state import (
    DEFAULT_PLAYBACK_INTERVAL_S,
    MAX_PLAYBACK_INTERVAL_S,
    MIN_PLAYBACK_INTERVAL_S,
    PLAYBACK_INTERVAL_STEP_S,
    normalize_playback_interval,
)


def register_control_callbacks(app):
    """Register time, height, and playback callbacks for the plots tab."""
    @app.callback(
        Output("plots-time-heading", "children"),
        Output("plots-time-start-label", "children"),
        Output("plots-time-average-label", "children"),
        Input("plots-case-data", "data"),
        Input("plots-time-mode", "value"),
        Input("plots-global-time-range", "value"),
        Input("plots-global-time-point", "value"),
    )
    def update_time_label(case_data, time_mode, time_range, time_point):
        """Show the active time window or point directly in the section header."""
        if not case_data:
            return "Time", "Start time", "Average Length"
        start_value = time_point if time_point is not None else case_data.get("default_time_start_seconds")
        average_value = time_range if time_range is not None else case_data.get("default_time_duration_minutes")
        return (
            physical_time_window_label(case_data, start_value, average_value),
            start_time_label(start_value),
            average_length_label(average_value),
        )

    @app.callback(
        Output("plots-height-heading", "children"),
        Input("plots-case-data", "data"),
        Input("plots-global-height-range", "value"),
    )
    def update_height_label(case_data, height_range):
        """Show the active height window directly in the section header."""
        if not case_data:
            return "Height"
        active_range = height_range or case_data.get("default_height_range") or [0.0, 1.0]
        return f"Height: {float(active_range[0]):g} - {float(active_range[1]):g}"

    @app.callback(
        Output("plots-global-time-range-wrapper", "style"),
        Output("plots-global-time-point-wrapper", "style"),
        Input("plots-time-mode", "value"),
    )
    def toggle_time_mode(time_mode):
        """Keep both time sliders visible; the hidden mode store preserves callback compatibility."""
        return {"display": "block"}, {"display": "block"}

    @app.callback(
        Output("plots-global-time-point", "step", allow_duplicate=True),
        Output("plots-global-time-point", "max", allow_duplicate=True),
        Output("plots-global-time-point", "value", allow_duplicate=True),
        Input("plots-global-time-range", "value"),
        State("plots-case-data", "data"),
        State("plots-global-time-point", "value"),
        prevent_initial_call=True,
    )
    def sync_start_control_to_average_length(average_minutes, case_data, current_start):
        """Move start time in chunks that match the selected averaging window."""
        step = max(1.0e-6, float(average_minutes or 1.0)) * 60.0
        if not case_data:
            return step, no_update, no_update
        start_min = float(case_data.get("time_slider_start_min_seconds") or 0.0)
        start_max = time_start_max_for_duration(case_data, average_minutes)
        active_start = float(current_start if current_start is not None else case_data.get("default_time_start_seconds", start_min))
        return step, start_max, max(start_min, min(active_start, start_max))

    @app.callback(
        Output("plots-global-time-point", "value", allow_duplicate=True),
        Output("plots-global-time-range", "value", allow_duplicate=True),
        Input("plots-use-loss-window", "n_clicks"),
        Input("plots-use-pyplotgen-window", "n_clicks"),
        State("plots-case-data", "data"),
        prevent_initial_call=True,
    )
    def apply_time_window_preset(_loss_clicks, _pyplotgen_clicks, case_data):
        """Jump the time sliders to the loss or pyplotgen averaging window."""
        if not case_data:
            return no_update, no_update
        trigger = callback_context.triggered_id
        if trigger == "plots-use-loss-window":
            start = case_data.get("loss_time_start_seconds")
            duration = case_data.get("loss_time_duration_minutes")
        elif trigger == "plots-use-pyplotgen-window":
            start = case_data.get("pyplotgen_time_start_seconds")
            duration = case_data.get("pyplotgen_time_duration_minutes")
        else:
            return no_update, no_update
        if start is None or duration is None:
            return no_update, no_update
        return float(start), float(duration)

    @app.callback(
        Output("plots-playback", "data"),
        Input("plots-playback-toggle", "n_clicks"),
        Input("plots-playback-slower", "n_clicks"),
        Input("plots-playback-faster", "n_clicks"),
        Input("plots-time-mode", "value"),
        Input("plots-case-data", "data"),
        State("plots-playback", "data"),
        prevent_initial_call=True,
    )
    def update_playback(_toggle_clicks, _slower_clicks, _faster_clicks, time_mode, case_data, playback):
        """Update playback state in response to transport controls or mode changes."""
        current = dict(playback or {"playing": False, "interval_s": DEFAULT_PLAYBACK_INTERVAL_S, "inflight": False, "target_point": None})
        current["interval_s"] = normalize_playback_interval(current.get("interval_s", DEFAULT_PLAYBACK_INTERVAL_S))
        trigger = callback_context.triggered_id
        if trigger in {"plots-time-mode", "plots-case-data"}:
            if not case_data:
                current["playing"] = False
                current["inflight"] = False
                current["target_point"] = None
            return current
        if not case_data:
            current["playing"] = False
            current["inflight"] = False
            current["target_point"] = None
            return current
        if trigger == "plots-playback-toggle":
            current["playing"] = not bool(current.get("playing"))
            current["inflight"] = False
            current["target_point"] = None
        elif trigger == "plots-playback-slower":
            current["interval_s"] = normalize_playback_interval(current["interval_s"] + PLAYBACK_INTERVAL_STEP_S)
            current["inflight"] = False
            current["target_point"] = None
        elif trigger == "plots-playback-faster":
            current["interval_s"] = normalize_playback_interval(current["interval_s"] - PLAYBACK_INTERVAL_STEP_S)
            current["inflight"] = False
            current["target_point"] = None
        return current

    @app.callback(
        Output("plots-playback-interval", "disabled"),
        Output("plots-playback-interval", "interval"),
        Output("plots-playback-toggle", "children"),
        Output("plots-playback-toggle", "style"),
        Output("plots-playback-toggle", "disabled"),
        Output("plots-playback-slower", "disabled"),
        Output("plots-playback-faster", "disabled"),
        Input("plots-playback", "data"),
        Input("plots-time-mode", "value"),
        Input("plots-case-data", "data"),
    )
    def sync_playback_ui(playback, time_mode, case_data):
        """Translate playback state into interval timing and button presentation."""
        state = dict(playback or {"playing": False, "interval_s": DEFAULT_PLAYBACK_INTERVAL_S, "inflight": False})
        interval_s = normalize_playback_interval(state.get("interval_s", DEFAULT_PLAYBACK_INTERVAL_S))
        enabled = bool(case_data)
        playing = enabled and bool(state.get("playing"))
        inflight = bool(state.get("inflight"))
        button_style = {
            "minWidth": "110px",
            "backgroundColor": "#dc2626" if playing else "#16a34a",
            "color": "#ffffff",
            "border": "none",
            "fontSize": "16pt",
        }
        return (
            (not playing) or inflight,
            int(interval_s * 1000),
            f"{'Pause' if playing else 'Play'} ({interval_s:g}s)",
            button_style,
            not enabled,
            (not enabled) or interval_s >= MAX_PLAYBACK_INTERVAL_S,
            (not enabled) or interval_s <= MIN_PLAYBACK_INTERVAL_S,
        )

    @app.callback(
        Output("plots-playback", "data", allow_duplicate=True),
        Output("plots-global-time-point", "value", allow_duplicate=True),
        Input("plots-playback-interval", "n_intervals"),
        State("plots-playback", "data"),
        State("plots-time-mode", "value"),
        State("plots-global-time-range", "value"),
        State("plots-global-time-point", "value"),
        State("plots-global-time-point", "min"),
        State("plots-global-time-point", "max"),
        State("plots-case-data", "data"),
        prevent_initial_call=True,
    )
    def advance_time_point(_n_intervals, playback, time_mode, average_minutes, time_point, time_point_min, time_point_max, case_data):
        """Advance one playback frame when the current frame is not still rendering."""
        if not playback or not playback.get("playing") or playback.get("inflight"):
            return no_update, no_update
        min_point = float(time_point_min if time_point_min is not None else 0.0)
        max_point = max(min_point, float(time_point_max if time_point_max is not None else min_point))
        step = max(1.0e-6, float(average_minutes or 1.0)) * 60.0
        current_point = max(min_point, min(float(time_point if time_point is not None else min_point), max_point))
        next_point = min_point if current_point >= max_point else min(current_point + step, max_point)
        updated_playback = dict(playback)
        updated_playback["inflight"] = True
        updated_playback["target_point"] = next_point
        return updated_playback, next_point

    @app.callback(
        Output("plots-playback", "data", allow_duplicate=True),
        Input({"type": "budget-render-signal", "index": ALL}, "children"),
        Input({"type": "profile-render-signal", "index": ALL}, "children"),
        Input({"type": "subcolumn-render-signal", "index": ALL}, "children"),
        Input("plots-time-mode", "value"),
        State({"type": "budget-render-signal", "index": ALL}, "id"),
        State({"type": "profile-render-signal", "index": ALL}, "id"),
        State({"type": "subcolumn-render-signal", "index": ALL}, "id"),
        State("plots-plot-order", "data"),
        State("plots-plot-state", "data"),
        State("plots-playback", "data"),
        prevent_initial_call=True,
    )
    def unlock_playback(budget_signals, profile_signals, subcolumn_signals, time_mode, budget_ids, profile_ids, subcolumn_ids, plot_order, plot_state, playback):
        """Release the playback lock once all time-dependent plots finish a frame."""
        current = dict(playback or {})
        if not current.get("inflight"):
            return no_update
        target_point = current.get("target_point")
        if target_point is None:
            current["inflight"] = False
            return current
        time_dependent_types = {budget_plot.plot_type_id, profile_plot.plot_type_id, subcolumn_plot.plot_type_id}
        relevant_ids = [
            int(plot_id)
            for plot_id in (plot_order or [])
            if (plot_state or {}).get(str(plot_id), {}).get("plot_type") in time_dependent_types
        ]
        if not relevant_ids:
            current["inflight"] = False
            current["target_point"] = None
            return current
        signal_by_id = {}
        for ids, signals in (
            (budget_ids, budget_signals),
            (profile_ids, profile_signals),
            (subcolumn_ids, subcolumn_signals),
        ):
            for meta, value in zip(ids or [], signals or []):
                if not isinstance(meta, dict):
                    continue
                signal_by_id[int(meta.get("index"))] = value
        if all(signal_by_id.get(plot_id) == target_point for plot_id in relevant_ids):
            current["inflight"] = False
            current["target_point"] = None
            return current
        return no_update
