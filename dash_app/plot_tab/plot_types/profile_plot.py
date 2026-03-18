import threading

import numpy as np
import plotly.graph_objects as go
from dash import Input, MATCH, Output, Patch, State, callback_context

from .. import benchmark_overlay
from . import shared
from .base_plot import BasePlotType


class ProfilePlotType(BasePlotType):
    def __init__(self):
        super().__init__(
            plot_type_id="profile",
            default_vars=shared.DEFAULT_PROFILE_VARS,
            var_input_type="profile-var",
            graph_type="profile-graph",
            case_data_var_key="profile_vars",
            subtitle="Uses the global time controls and the global column controls.",
        )
        self._render_state_lock = threading.Lock()
        self._rendered_plot_ids = set()

    def get_variable_options(self, collection, case_data):
        require_all = bool(case_data.get("compare_mode"))
        return case_data.get("profile_vars") or shared.list_profile_vars(collection, require_all=require_all)

    def supports_case(self, case_data, collection):
        return bool(self.get_variable_options(collection, case_data))

    def _time_indices(self, global_context):
        case_data = global_context.get("case_data") or {}
        time_len = max(int(case_data.get("time_len") or 1), 1)
        slider_range = global_context.get("time_range") or case_data.get("default_time_range") or [1, time_len]
        time_mode = global_context.get("time_mode") or "range"
        time_point = global_context.get("time_point") or slider_range[0]
        if time_mode == "point":
            point_idx = shared.slider_value_to_index(time_point, time_len)
            return [point_idx, point_idx], time_mode, time_point
        return shared.slider_range_to_indices(slider_range, time_len), time_mode, time_point

    def _mark_full_render(self, plot_id):
        with self._render_state_lock:
            self._rendered_plot_ids.add(int(plot_id))

    def _has_full_render(self, plot_id):
        with self._render_state_lock:
            return int(plot_id) in self._rendered_plot_ids

    def clear_render_state(self):
        with self._render_state_lock:
            self._rendered_plot_ids.clear()

    def _point_mode_x_range(self, bounds_list):
        finite_bounds = [bounds for bounds in bounds_list if bounds is not None]
        if not finite_bounds:
            return None
        low = min(bounds[0] for bounds in finite_bounds)
        high = max(bounds[1] for bounds in finite_bounds)
        return shared.padded_range(low, high)

    def _single_trace_specs(self, files, var_name, case_data, global_context):
        path, _meta = shared.dataset_info_for_var(files, var_name)
        if path is None:
            return None
        time_indices, time_mode, _time_point = self._time_indices(global_context)
        col_index = int(global_context.get("selected_column") or 0)
        column_mode = global_context.get("column_mode") or "single"
        trace_specs = []
        x_values = []
        point_bounds = []
        enabled_sources = benchmark_overlay.sanitize_enabled_sources(case_data, global_context.get("enabled_benchmarks"))
        for source_name in enabled_sources:
            benchmark_profile = benchmark_overlay.extract_benchmark_profile(
                case_data,
                source_name,
                var_name,
                global_context.get("time_mode") or "range",
                global_context.get("time_range") or case_data.get("default_time_range"),
                global_context.get("time_point"),
            )
            if benchmark_profile is None:
                continue
            trace_specs.append(
                {
                    "x": benchmark_profile["profile"],
                    "y": benchmark_profile["z_values"],
                    "name": benchmark_profile["label"],
                    "line": benchmark_profile["line"],
                    "opacity": None,
                    "showlegend": True,
                }
            )
            x_values.extend([value for value in benchmark_profile["profile"] if value is not None])
            point_bounds.append(shared._finite_bounds(benchmark_profile["profile"]))
        extracted = shared.extract_time_avg_profile_for_path(
            path,
            var_name,
            time_indices,
            col_index=col_index,
            column_mode=column_mode,
        )
        if extracted is None:
            return None
        point_bounds.append(extracted["bounds"])
        for idx, profile in enumerate(extracted["profiles"]):
            line = {"width": 1.5 if column_mode == "all" else 2.0}
            if column_mode == "all":
                line["color"] = shared.column_line_color(idx)
            else:
                line["color"] = shared.source_line_color(0)
            trace_specs.append(
                {
                    "x": profile,
                    "y": extracted["z_values"],
                    "name": extracted["labels"][idx],
                    "line": line,
                    "opacity": 0.85 if column_mode == "all" else None,
                    "showlegend": False,
                }
            )
            x_values.extend([value for value in profile if value is not None])
        x_range = self._point_mode_x_range(point_bounds) if time_mode == "point" else shared.padded_data_range(x_values)
        return {
            "trace_specs": trace_specs,
            "x_range": x_range,
            "z_units": extracted["z_units"],
            "var_units": extracted["units"],
            "var_long_name": extracted["long_name"],
            "column_mode": column_mode,
            "showlegend": any(spec["showlegend"] for spec in trace_specs),
            "enabled_benchmarks": tuple(enabled_sources),
        }

    def _build_single_figure(self, files, var_name, case_data, global_context, theme_name):
        trace_bundle = self._single_trace_specs(files, var_name, case_data, global_context)
        if trace_bundle is None:
            return shared.make_empty_figure(f"{var_name} is not available for this case.", theme_name)
        fig = go.Figure()
        for spec in trace_bundle["trace_specs"]:
            fig.add_trace(
                go.Scatter(
                    x=spec["x"],
                    y=spec["y"],
                    mode="lines",
                    name=spec["name"],
                    line=spec["line"],
                    opacity=spec["opacity"],
                    showlegend=spec["showlegend"],
                )
            )
        if not fig.data:
            return shared.make_empty_figure(f"{var_name} is not available for this case.", theme_name)
        plot_height = shared.figure_height_for_size(global_context.get("size"))
        shared.apply_figure_chrome(
            fig,
            title=shared.format_plot_title(var_name, trace_bundle["var_long_name"]),
            showlegend=trace_bundle["showlegend"],
            height=plot_height,
            uirevision=shared.figure_uirevision(
                self.plot_type_id,
                case_data,
                var_name,
                trace_bundle["column_mode"],
                ",".join(trace_bundle["enabled_benchmarks"]),
            ),
        )
        fig.update_layout(
            xaxis_title=shared.format_axis_title(var_name, trace_bundle["var_units"]),
            yaxis_title=shared.format_height_axis_title(trace_bundle["z_units"]),
        )
        if trace_bundle["x_range"]:
            fig.update_xaxes(range=trace_bundle["x_range"])
        height_range = global_context.get("height_range") or case_data.get("default_height_range")
        if height_range:
            shared.apply_axis_bounds(fig, "y", shared.padded_range(height_range[0], height_range[1]))
        shared.apply_plot_theme(fig, theme_name)
        return fig

    def _compare_trace_specs(self, files, var_name, case_data, global_context):
        time_indices, time_mode, _time_point = self._time_indices(global_context)
        col_index = int(global_context.get("selected_column") or 0)
        column_mode = global_context.get("column_mode") or "single"
        var_units = ""
        z_units = ""
        var_long_name = ""
        legend_labels = []
        trace_specs = []
        x_values = []
        playback_bounds = []
        source_labels = case_data.get("source_labels") or [f"output {idx + 1}" for idx in range(len(files))]
        enabled_sources = benchmark_overlay.sanitize_enabled_sources(case_data, global_context.get("enabled_benchmarks"))
        for source_name in enabled_sources:
            benchmark_profile = benchmark_overlay.extract_benchmark_profile(
                case_data,
                source_name,
                var_name,
                global_context.get("time_mode") or "range",
                global_context.get("time_range") or case_data.get("default_time_range"),
                global_context.get("time_point"),
            )
            if benchmark_profile is None:
                continue
            trace_specs.append(
                {
                    "x": benchmark_profile["profile"],
                    "y": benchmark_profile["z_values"],
                    "name": benchmark_profile["label"],
                    "line": benchmark_profile["line"],
                    "opacity": None,
                    "showlegend": True,
                }
            )
            x_values.extend([value for value in benchmark_profile["profile"] if value is not None])
            if time_mode == "point":
                playback_bounds.append(shared._finite_bounds(benchmark_profile["profile"]))
        for source_idx, path in enumerate(files):
            extracted = shared.extract_time_avg_profile_for_path(
                path,
                var_name,
                time_indices,
                col_index=col_index,
                column_mode=column_mode,
            )
            if extracted is None:
                continue
            z_units = extracted["z_units"] or z_units
            var_units = extracted["units"] or var_units
            var_long_name = extracted["long_name"] or var_long_name
            label = source_labels[source_idx] if source_idx < len(source_labels) else f"output {source_idx + 1}"
            legend_labels.append(label)
            for idx, profile in enumerate(extracted["profiles"]):
                line = {"width": 1.2 if column_mode == "all" else 2.0}
                if column_mode == "all":
                    line["color"] = shared.column_line_color(idx)
                else:
                    line["color"] = shared.source_line_color(source_idx)
                line["dash"] = shared.source_line_dash(source_idx)
                trace_specs.append(
                    {
                        "x": profile,
                        "y": extracted["z_values"],
                        "name": label,
                        "line": line,
                        "opacity": 0.8 if column_mode == "all" else None,
                        "showlegend": column_mode != "all",
                    }
                )
                x_values.extend([value for value in profile if value is not None])
            if time_mode == "point":
                playback_bounds.append(extracted["bounds"])
        x_range = self._point_mode_x_range(playback_bounds) if time_mode == "point" else shared.padded_data_range(x_values)
        return {
            "trace_specs": trace_specs,
            "x_range": x_range,
            "var_units": var_units,
            "z_units": z_units,
            "var_long_name": var_long_name,
            "legend_labels": legend_labels,
            "column_mode": column_mode,
            "enabled_benchmarks": tuple(enabled_sources),
        }

    def _build_compare_figure(self, files, var_name, case_data, global_context, theme_name):
        trace_bundle = self._compare_trace_specs(files, var_name, case_data, global_context)
        if not trace_bundle["trace_specs"]:
            return shared.make_empty_figure(f"{var_name} is not available across the selected outputs.", theme_name)
        fig = go.Figure()
        for spec in trace_bundle["trace_specs"]:
            fig.add_trace(
                go.Scatter(
                    x=spec["x"],
                    y=spec["y"],
                    mode="lines",
                    name=spec["name"],
                    line=spec["line"],
                    opacity=spec["opacity"],
                    showlegend=spec["showlegend"],
                )
            )
        if trace_bundle["column_mode"] == "all":
            shared.add_directory_legend_traces(fig, trace_bundle["legend_labels"])
        plot_height = shared.figure_height_for_size(global_context.get("size"))
        shared.apply_figure_chrome(
            fig,
            title=shared.format_plot_title(var_name, trace_bundle["var_long_name"]),
            showlegend=bool(case_data.get("compare_mode")),
            height=plot_height,
            uirevision=shared.figure_uirevision(
                self.plot_type_id,
                case_data,
                var_name,
                trace_bundle["column_mode"],
                ",".join(trace_bundle["enabled_benchmarks"]),
            ),
        )
        fig.update_layout(
            xaxis_title=shared.format_axis_title(var_name, trace_bundle["var_units"]),
            yaxis_title=shared.format_height_axis_title(trace_bundle["z_units"]),
        )
        if trace_bundle["x_range"]:
            fig.update_xaxes(range=trace_bundle["x_range"])
        height_range = global_context.get("height_range") or case_data.get("default_height_range")
        if height_range:
            shared.apply_axis_bounds(fig, "y", shared.padded_range(height_range[0], height_range[1]))
        shared.apply_plot_theme(fig, theme_name)
        return fig

    def _build_single_patch(self, files, var_name, case_data, global_context):
        trace_bundle = self._single_trace_specs(files, var_name, case_data, global_context)
        if trace_bundle is None or not trace_bundle["trace_specs"]:
            return None
        patch = Patch()
        for idx, spec in enumerate(trace_bundle["trace_specs"]):
            patch["data"][idx]["x"] = np.asarray(spec["x"], dtype=float).tolist()
        if trace_bundle["x_range"]:
            patch["layout"]["xaxis"]["range"] = list(trace_bundle["x_range"])
        return patch

    def _build_compare_patch(self, files, var_name, case_data, global_context):
        trace_bundle = self._compare_trace_specs(files, var_name, case_data, global_context)
        if not trace_bundle["trace_specs"]:
            return None
        patch = Patch()
        for idx, spec in enumerate(trace_bundle["trace_specs"]):
            patch["data"][idx]["x"] = np.asarray(spec["x"], dtype=float).tolist()
        if trace_bundle["x_range"]:
            patch["layout"]["xaxis"]["range"] = list(trace_bundle["x_range"])
        return patch

    def build_figure(self, state, global_context):
        case_data = global_context.get("case_data") or {}
        theme_name = global_context.get("theme_name")
        files = case_data.get("files") or []
        var_name = state.get("var")
        if not files or not var_name:
            return shared.make_empty_figure("Select a case and variable.", theme_name)
        if case_data.get("compare_mode"):
            return self._build_compare_figure(files, var_name, case_data, global_context, theme_name)
        return self._build_single_figure(files, var_name, case_data, global_context, theme_name)

    def build_patch(self, state, global_context):
        case_data = global_context.get("case_data") or {}
        files = case_data.get("files") or []
        var_name = state.get("var")
        if not files or not var_name or (global_context.get("time_mode") or "range") != "point":
            return None
        if case_data.get("compare_mode"):
            return self._build_compare_patch(files, var_name, case_data, global_context)
        return self._build_single_patch(files, var_name, case_data, global_context)

    def register_callbacks(self, app):
        @app.callback(
            Output(self.graph_id(MATCH), "figure"),
            Output(self.render_signal_id(MATCH), "children"),
            Input(self.var_input_id(MATCH), "value"),
            Input("plots-case-data", "data"),
            Input("plots-enabled-benchmarks", "data"),
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
        def _update_profile_graph(var_name, case_data, enabled_benchmarks, time_mode, time_range, time_point, height_range, selected_column, column_mode, theme_name, size_store_value, graph_id):
            plot_id = int((graph_id or {}).get("index", -1))
            size_value = shared.normalize_plot_size(size_store_value)
            signal = int(time_point) if time_mode == "point" and time_point is not None else ""
            if callback_context.triggered_id == "plots-global-time-point" and plot_id >= 0 and self._has_full_render(plot_id):
                patch = self.build_patch(
                    {"var": var_name, "size": size_value},
                    {
                        "case_data": case_data,
                        "enabled_benchmarks": enabled_benchmarks,
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
                    return patch, signal
            fig = self.build_figure(
                {"var": var_name, "size": size_value},
                {
                    "case_data": case_data,
                    "enabled_benchmarks": enabled_benchmarks,
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
            return fig, signal


PLOT = ProfilePlotType()
