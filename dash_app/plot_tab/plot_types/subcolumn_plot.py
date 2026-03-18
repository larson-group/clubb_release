import threading

import numpy as np
import plotly.graph_objects as go
from dash import Input, Output, MATCH, Patch, State, callback_context, html

from . import shared
from .base_plot import BasePlotType


class SubcolumnPlotType(BasePlotType):
    def __init__(self):
        super().__init__(
            plot_type_id="subcolumn",
            default_vars=shared.SUBCOLUMN_DEFAULT_VARS,
            var_input_type="subcolumn-var",
            graph_type="subcolumn-graph",
            case_data_var_key="subcolumn_vars",
            subtitle="Uses the global time controls. In all-column mode, it plots every available line from every column.",
        )
        self._render_state_lock = threading.Lock()
        self._rendered_plot_ids = set()

    def get_variable_options(self, collection, case_data):
        return case_data.get("subcolumn_vars") or shared.list_subcolumn_vars(collection)

    def auxiliary_components(self, plot_id):
        return [
            html.Div(
                id={"type": "subcolumn-note", "index": plot_id},
                style={"marginBottom": "8px", "fontSize": "12px", "color": "#555"},
            )
        ]

    def _mark_full_render(self, plot_id):
        with self._render_state_lock:
            self._rendered_plot_ids.add(int(plot_id))

    def _has_full_render(self, plot_id):
        with self._render_state_lock:
            return int(plot_id) in self._rendered_plot_ids

    def clear_render_state(self):
        with self._render_state_lock:
            self._rendered_plot_ids.clear()

    def _build_note(self, has_true_subcolumns, column_mode, note_suffix):
        if has_true_subcolumns and column_mode == "all":
            return f"All columns and all subcolumns, {note_suffix}"
        if has_true_subcolumns:
            return f"True subcolumn lines, {note_suffix}"
        if column_mode == "all":
            return f"Fallback: no explicit subcolumn dimension in file, plotting all columns, {note_suffix}"
        return f"Fallback: no explicit subcolumn dimension in file, {note_suffix}"

    def _trace_bundle(self, files, base_name, case_data, global_context):
        path, resolved_name, meta = shared.subcolumn_source_for_paths(files, base_name)
        if path is None:
            return None
        slider_range = global_context.get("time_range") or case_data.get("default_time_range") or [1, max(int(case_data.get("time_len") or 1), 1)]
        time_mode = global_context.get("time_mode") or "range"
        time_point = global_context.get("time_point") or slider_range[0]
        time_len = max(int(meta["time_len"] or 1), 1)
        if time_mode == "point":
            point_idx = shared.slider_value_to_index(time_point, time_len)
            time_indices = [point_idx, point_idx]
            note_suffix = "single time"
        else:
            time_indices = shared.slider_range_to_indices(slider_range, time_len)
            note_suffix = "time-averaged range"
        column_mode = global_context.get("column_mode") or "single"
        col_index = int(global_context.get("selected_column") or 0)
        shared.ensure_subcolumn_plot_data(path, base_name)
        z_vals, profiles, labels, has_true_subcolumns = shared.extract_subcolumn_profiles_for_path(
            path,
            base_name,
            time_indices,
            col_index=col_index,
            column_mode=column_mode,
        )
        if z_vals is None:
            return None
        info = meta["var_info"].get(resolved_name, {})
        units = info.get("units", "") if resolved_name else ""
        var_long_name = info.get("long_name", "") if resolved_name else ""
        if time_mode == "point":
            x_range = shared.subcolumn_x_range_for_path(path, base_name, col_index=col_index, column_mode=column_mode)
        else:
            x_range = shared.padded_data_range(profiles)
        return {
            "z_units": meta["dim_units"].get(info.get("z_dim"), ""),
            "resolved_name": resolved_name,
            "z_vals": z_vals,
            "profiles": profiles,
            "labels": labels,
            "has_true_subcolumns": has_true_subcolumns,
            "units": units,
            "var_long_name": var_long_name,
            "column_mode": column_mode,
            "note": self._build_note(has_true_subcolumns, column_mode, note_suffix),
            "x_range": x_range,
        }

    def build_figure(self, state, global_context):
        case_data = global_context.get("case_data") or {}
        theme_name = global_context.get("theme_name")
        files = case_data.get("files") or []
        base_name = state.get("var")
        if not files or not base_name:
            return shared.make_empty_figure("Select a case and variable.", theme_name), ""
        if case_data.get("compare_mode"):
            return shared.make_empty_figure("Subcolumn plots are disabled in compare mode.", theme_name), ""
        trace_bundle = self._trace_bundle(files, base_name, case_data, global_context)
        if trace_bundle is None:
            return shared.make_empty_figure(f"No subcolumn-compatible variable found for {base_name}.", theme_name), ""
        plot_height = shared.figure_height_for_size(global_context.get("size"))
        fig = go.Figure()
        for idx in range(trace_bundle["profiles"].shape[1]):
            fig.add_trace(
                go.Scatter(
                    x=trace_bundle["profiles"][:, idx],
                    y=trace_bundle["z_vals"],
                    mode="lines",
                    name=trace_bundle["labels"][idx],
                    line={"width": 1.5 if trace_bundle["profiles"].shape[1] > 1 else 2.0},
                    opacity=0.9,
                )
            )
        shared.apply_figure_chrome(
            fig,
            title=shared.format_plot_title(base_name, trace_bundle["var_long_name"]),
            showlegend=False,
            height=plot_height,
            uirevision=shared.figure_uirevision(self.plot_type_id, case_data, base_name, trace_bundle["column_mode"]),
        )
        fig.update_layout(
            xaxis_title=shared.format_axis_title(base_name, trace_bundle["units"]),
            yaxis_title=shared.format_height_axis_title(trace_bundle["z_units"]),
        )
        if trace_bundle["x_range"]:
            fig.update_xaxes(range=trace_bundle["x_range"])
        height_range = global_context.get("height_range") or case_data.get("default_height_range")
        if height_range:
            shared.apply_axis_bounds(fig, "y", shared.padded_range(height_range[0], height_range[1]))
        shared.apply_plot_theme(fig, theme_name)
        return fig, trace_bundle["note"]

    def build_patch(self, state, global_context):
        case_data = global_context.get("case_data") or {}
        files = case_data.get("files") or []
        base_name = state.get("var")
        if (
            not files
            or not base_name
            or case_data.get("compare_mode")
            or (global_context.get("time_mode") or "range") != "point"
        ):
            return None, None
        trace_bundle = self._trace_bundle(files, base_name, case_data, global_context)
        if trace_bundle is None:
            return None, None
        patch = Patch()
        for idx in range(trace_bundle["profiles"].shape[1]):
            patch["data"][idx]["x"] = np.asarray(trace_bundle["profiles"][:, idx], dtype=float).tolist()
        if trace_bundle["x_range"]:
            patch["layout"]["xaxis"]["range"] = list(trace_bundle["x_range"])
        return patch, trace_bundle["note"]

    def register_callbacks(self, app):
        @app.callback(
            Output(self.graph_id(MATCH), "figure"),
            Output({"type": "subcolumn-note", "index": MATCH}, "children"),
            Output(self.render_signal_id(MATCH), "children"),
            Input(self.var_input_id(MATCH), "value"),
            Input("plots-case-data", "data"),
            Input("plots-time-mode", "value"),
            Input("plots-global-time-range", "value"),
            Input("plots-global-time-point", "value"),
            Input("plots-global-height-range", "value"),
            Input("plots-selected-column", "data"),
            Input("plots-column-mode", "value"),
            Input("theme-store", "data"),
            Input(self.size_store_id(MATCH), "data"),
            State(self.graph_id(MATCH), "id"),
        )
        def _update_subcolumn_graph(var_name, case_data, time_mode, time_range, time_point, height_range, selected_column, column_mode, theme_name, size_store_value, graph_id):
            plot_id = int((graph_id or {}).get("index", -1))
            size_value = shared.normalize_plot_size(size_store_value)
            signal = int(time_point) if time_mode == "point" and time_point is not None else ""
            if callback_context.triggered_id == "plots-global-time-point" and plot_id >= 0 and self._has_full_render(plot_id):
                patch, note = self.build_patch(
                    {"var": var_name, "size": size_value},
                    {
                        "case_data": case_data,
                        "time_mode": time_mode,
                        "time_range": time_range,
                        "time_point": time_point,
                        "height_range": height_range,
                        "selected_column": selected_column,
                        "column_mode": column_mode,
                        "size": size_value,
                        "theme_name": theme_name,
                    },
                )
                if patch is not None:
                    return patch, note, signal
            fig, note = self.build_figure(
                {"var": var_name, "size": size_value},
                {
                    "case_data": case_data,
                    "time_mode": time_mode,
                    "time_range": time_range,
                    "time_point": time_point,
                    "height_range": height_range,
                    "selected_column": selected_column,
                    "column_mode": column_mode,
                    "size": size_value,
                    "theme_name": theme_name,
                },
            )
            if plot_id >= 0:
                self._mark_full_render(plot_id)
            return fig, note, signal


PLOT = SubcolumnPlotType()
