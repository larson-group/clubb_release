import threading

import numpy as np
import plotly.graph_objects as go
from dash import Input, MATCH, Output, Patch, State, callback_context

from . import shared
from .base_plot import BasePlotType
from .budget_groups import BUDGET_GROUPS, DEFAULT_BUDGET_GROUPS


class BudgetPlotType(BasePlotType):
    def __init__(self):
        super().__init__(
            plot_type_id="budget",
            default_vars=DEFAULT_BUDGET_GROUPS,
            var_input_type="budget-var",
            graph_type="budget-graph",
            case_data_var_key="budget_groups",
            subtitle="Uses the global time controls. Budget plots are currently single-column only.",
        )
        self._render_state_lock = threading.Lock()
        self._rendered_plot_ids = set()

    def get_variable_options(self, collection, case_data):
        return case_data.get("budget_groups") or shared.list_budget_groups(collection)

    def _mark_full_render(self, plot_id):
        with self._render_state_lock:
            self._rendered_plot_ids.add(int(plot_id))

    def _has_full_render(self, plot_id):
        with self._render_state_lock:
            return int(plot_id) in self._rendered_plot_ids

    def clear_render_state(self):
        with self._render_state_lock:
            self._rendered_plot_ids.clear()

    def _trace_bundle(self, files, group_name, case_data, global_context):
        group = BUDGET_GROUPS.get(group_name)
        if not group:
            return None
        time_len = max(int(case_data.get("time_len") or 1), 1)
        slider_range = global_context.get("time_range") or case_data.get("default_time_range") or [1, time_len]
        time_mode = global_context.get("time_mode") or "range"
        time_point = global_context.get("time_point") or slider_range[0]
        if time_mode == "point":
            point_idx = shared.slider_value_to_index(time_point, time_len)
            time_indices = [point_idx, point_idx]
        else:
            time_indices = shared.slider_range_to_indices(slider_range, time_len)

        col_index = int(global_context.get("selected_column") or 0)
        trace_specs = []
        x_values = []
        x_units = []
        z_units = ""
        playback_bounds = []
        for term in group["terms"]:
            path, _meta = shared.dataset_info_for_var(files, term)
            if path is None:
                continue
            extracted = shared.extract_time_avg_profile_for_path(
                path,
                term,
                time_indices,
                col_index=col_index,
                column_mode="single",
            )
            if extracted is None:
                continue
            profile = extracted["profiles"][0]
            trace_specs.append(
                {
                    "x": profile,
                    "y": extracted["z_values"],
                    "name": term,
                    "line": {"width": 1.8},
                }
            )
            x_values.extend([value for value in profile if value is not None])
            if extracted["units"]:
                x_units.append(extracted["units"])
            z_units = extracted["z_units"] or z_units
            if time_mode == "point":
                playback_bounds.append(extracted["bounds"])
        if time_mode == "point":
            finite_bounds = [bounds for bounds in playback_bounds if bounds is not None]
            if finite_bounds:
                x_range = shared.padded_range(
                    min(bounds[0] for bounds in finite_bounds),
                    max(bounds[1] for bounds in finite_bounds),
                )
            else:
                x_range = None
        else:
            x_range = shared.padded_data_range(x_values)
        budget_units = x_units[0] if x_units and all(unit == x_units[0] for unit in x_units) else ""
        return {
            "group": group,
            "trace_specs": trace_specs,
            "x_range": x_range,
            "budget_units": budget_units,
            "z_units": z_units,
        }

    def build_figure(self, state, global_context):
        case_data = global_context.get("case_data") or {}
        theme_name = global_context.get("theme_name")
        group_name = state.get("var")
        files = case_data.get("files") or []
        if not files or not group_name:
            return shared.make_empty_figure("Select a case and budget group.", theme_name)
        if case_data.get("compare_mode"):
            return shared.make_empty_figure("Budget plots are disabled in compare mode.", theme_name)
        if (global_context.get("column_mode") or "single") == "all":
            return shared.make_empty_figure("Budget plots are disabled in overplot mode.", theme_name)
        trace_bundle = self._trace_bundle(files, group_name, case_data, global_context)
        if trace_bundle is None:
            return shared.make_empty_figure(f"Unknown budget group: {group_name}", theme_name)
        fig = go.Figure()
        for spec in trace_bundle["trace_specs"]:
            fig.add_trace(
                go.Scatter(
                    x=spec["x"],
                    y=spec["y"],
                    mode="lines",
                    name=spec["name"],
                    line=spec["line"],
                )
            )
        if not fig.data:
            return shared.make_empty_figure(f"No available terms found for the {group_name} budget group.", theme_name)

        plot_height = shared.figure_height_for_size(global_context.get("size"))
        shared.apply_figure_chrome(
            fig,
            title=shared.format_plot_title(group_name, trace_bundle["group"]["label"]),
            showlegend=True,
            height=plot_height,
            uirevision=shared.figure_uirevision(self.plot_type_id, case_data, group_name, "single"),
        )
        fig.update_layout(
            xaxis_title=shared.format_axis_title(group_name, trace_bundle["budget_units"]),
            yaxis_title=shared.format_height_axis_title(trace_bundle["z_units"]),
        )
        if trace_bundle["x_range"]:
            fig.update_xaxes(range=trace_bundle["x_range"])
        height_range = global_context.get("height_range") or case_data.get("default_height_range")
        if height_range:
            shared.apply_axis_bounds(fig, "y", shared.padded_range(height_range[0], height_range[1]))
        shared.apply_plot_theme(fig, theme_name)
        return fig

    def build_patch(self, state, global_context):
        case_data = global_context.get("case_data") or {}
        files = case_data.get("files") or []
        group_name = state.get("var")
        if (
            not files
            or not group_name
            or case_data.get("compare_mode")
            or (global_context.get("column_mode") or "single") == "all"
            or (global_context.get("time_mode") or "range") != "point"
        ):
            return None
        trace_bundle = self._trace_bundle(files, group_name, case_data, global_context)
        if trace_bundle is None or not trace_bundle["trace_specs"]:
            return None
        patch = Patch()
        for idx, spec in enumerate(trace_bundle["trace_specs"]):
            patch["data"][idx]["x"] = np.asarray(spec["x"], dtype=float).tolist()
        if trace_bundle["x_range"]:
            patch["layout"]["xaxis"]["range"] = list(trace_bundle["x_range"])
        return patch

    def register_callbacks(self, app):
        @app.callback(
            Output(self.graph_id(MATCH), "figure"),
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
        def _update_budget_graph(var_name, case_data, time_mode, time_range, time_point, height_range, selected_column, column_mode, theme_name, size_store_value, graph_id):
            plot_id = int((graph_id or {}).get("index", -1))
            size_value = shared.normalize_plot_size(size_store_value)
            signal = int(time_point) if time_mode == "point" and time_point is not None else ""
            if callback_context.triggered_id == "plots-global-time-point" and plot_id >= 0 and self._has_full_render(plot_id):
                patch = self.build_patch(
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
                    return patch, signal
            fig = self.build_figure(
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
            return fig, signal


PLOT = BudgetPlotType()
