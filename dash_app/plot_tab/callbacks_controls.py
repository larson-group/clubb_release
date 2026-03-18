from dash import ALL, Input, Output, State, callback_context, no_update

from .plot_types.budget_plot import PLOT as budget_plot
from .plot_types.profile_plot import PLOT as profile_plot
from .plot_types.shared import time_label_for_point, time_label_for_range
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
        Input("plots-case-data", "data"),
        Input("plots-time-mode", "value"),
        Input("plots-global-time-range", "value"),
        Input("plots-global-time-point", "value"),
    )
    def update_time_label(case_data, time_mode, time_range, time_point):
        """Show the active time window or point directly in the section header."""
        if not case_data:
            return "Time"
        time_seconds = case_data.get("time_seconds")
        if time_mode == "point":
            return time_label_for_point(time_seconds, time_point or 1)
        return time_label_for_range(time_seconds, time_range or case_data.get("default_time_range"))

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
        """Swap the visible time control between range and single-point modes."""
        if time_mode == "point":
            return {"display": "none"}, {"display": "block"}
        return {"display": "block"}, {"display": "none"}

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
            if time_mode != "point" or not case_data:
                current["playing"] = False
                current["inflight"] = False
                current["target_point"] = None
            return current
        if time_mode != "point" or not case_data:
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
        enabled = bool(case_data) and time_mode == "point"
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
        State("plots-global-time-point", "value"),
        State("plots-global-time-point", "max"),
        prevent_initial_call=True,
    )
    def advance_time_point(_n_intervals, playback, time_mode, time_point, time_point_max):
        """Advance one playback frame when the current frame is not still rendering."""
        if time_mode != "point" or not playback or not playback.get("playing") or playback.get("inflight"):
            return no_update, no_update
        max_point = max(1, int(time_point_max or 1))
        current_point = max(1, min(int(time_point or 1), max_point))
        next_point = 1 if current_point >= max_point else current_point + 1
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
        if time_mode != "point":
            current["inflight"] = False
            current["target_point"] = None
            return current
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
