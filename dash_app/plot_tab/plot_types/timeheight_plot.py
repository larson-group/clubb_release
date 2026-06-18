import numpy as np
import plotly.graph_objects as go
from dash import Input, Output, MATCH
from plotly.subplots import make_subplots

from ..benchmark_overlay import extract_benchmark_timeheight_panels
from . import shared
from .base_plot import BasePlotType


class TimeHeightPlotType(BasePlotType):
    def __init__(self):
        super().__init__(
            plot_type_id="timeheight",
            default_vars=shared.TIMEHEIGHT_DEFAULT_VARS,
            var_input_type="timeheight-var",
            graph_type="timeheight-graph",
            case_data_var_key="timeheight_vars",
            subtitle="Uses the CLUBB stats time extent and selected average length. In all-column mode, any column dimension is averaged before plotting.",
        )

    def case_data_options(self, case_data):
        return list((case_data or {}).get("timeheight_vars") or [])

    def get_variable_options(self, collection, case_data):
        return case_data.get("timeheight_vars") or shared.list_timeheight_vars(collection)

    def build_figure(self, state, global_context):
        case_data = global_context.get("case_data") or {}
        theme_name = global_context.get("theme_name")
        files = case_data.get("files") or []
        var_name = state.get("var")
        if not files or not var_name:
            return shared.make_empty_figure("Select a case and variable.", theme_name)
        if case_data.get("compare_mode"):
            return shared.make_empty_figure("Time-height plots are disabled in compare mode.", theme_name)
        path, _meta = shared.dataset_info_for_var(files, var_name)
        if path is None:
            return shared.make_empty_figure(f"{var_name} is not available for this case.", theme_name)
        column_mode = global_context.get("column_mode") or "single"
        result = shared.extract_time_height_for_path(
            path,
            var_name,
            col_index=int(global_context.get("selected_column") or 0),
            column_mode=column_mode,
        )
        if result is None:
            return shared.make_empty_figure(f"{var_name} is not compatible with a time-height plot.", theme_name)
        result = self._average_time_height(result, global_context)
        benchmark_panels = extract_benchmark_timeheight_panels(case_data, var_name, global_context.get("enabled_benchmarks"))
        if benchmark_panels:
            global_context = {**global_context, "size": "large"}
            return self._build_panel_figure(
                [
                    {
                        "title": self._clubb_panel_title(global_context),
                        "data": result,
                    },
                    *[
                        {**panel, "data": self._average_time_height(panel["data"], global_context)}
                        for panel in benchmark_panels
                    ],
                ],
                var_name,
                global_context,
            )
        return self._build_single_figure(var_name, result, global_context)

    def _average_time_height(self, result, global_context):
        time_vals, z_vals, data, var_units, var_long_name, z_units = result
        case_data = global_context.get("case_data") or {}
        start_seconds = float(case_data.get("time_slider_start_min_seconds", case_data.get("default_time_start_seconds", 0.0)))
        min_average = float(case_data.get("time_slider_duration_min_minutes") or 0.0)
        average_minutes = max(min_average, float(global_context.get("time_range") if global_context.get("time_range") is not None else case_data.get("default_time_duration_minutes", min_average or 1.0)))
        window_seconds = average_minutes * 60.0
        time_arr = np.asarray(time_vals, dtype=float)
        data_arr = np.asarray(data, dtype=float)
        if time_arr.size == 0 or data_arr.ndim != 2 or window_seconds <= 0.0:
            return result
        final_end = self._clubb_time_final_end(case_data, time_arr)
        if not np.isfinite(final_end) or start_seconds >= final_end:
            return result
        bin_starts = np.arange(start_seconds, final_end, window_seconds)
        averaged = []
        bin_ends = []
        for bin_start in bin_starts:
            bin_end = min(bin_start + window_seconds, final_end)
            mask = (time_arr > bin_start + 1.0e-6) & (time_arr <= bin_end + 1.0e-6)
            if np.any(mask):
                averaged.append(np.nanmean(data_arr[mask, :], axis=0))
            else:
                averaged.append(np.full(data_arr.shape[1], np.nan, dtype=float))
            bin_ends.append(bin_end)
        if not averaged:
            return result
        return (
            np.asarray(bin_ends, dtype=float),
            z_vals,
            np.asarray(averaged, dtype=float),
            var_units,
            var_long_name,
            z_units,
        )

    def _clubb_time_final_end(self, case_data, fallback_time_arr):
        final_end = case_data.get("time_slider_final_end_seconds")
        if final_end is not None:
            return float(final_end)
        time_bounds = np.asarray(case_data.get("time_bounds_seconds") or [], dtype=float)
        if time_bounds.ndim == 2 and time_bounds.shape[1] == 2 and len(time_bounds) > 0:
            return float(np.nanmax(time_bounds[:, 1]))
        return float(np.nanmax(fallback_time_arr))

    def _clubb_panel_title(self, global_context):
        column_mode = global_context.get("column_mode") or "single"
        if column_mode == "all":
            return "CLUBB"
        selected_column = int(global_context.get("selected_column") or 0) + 1
        return f"CLUBB col {selected_column}"

    def _build_single_figure(self, var_name, result, global_context):
        case_data = global_context.get("case_data") or {}
        theme_name = global_context.get("theme_name")
        time_vals, z_vals, data, var_units, var_long_name, z_units = result
        colorbar = {"title": var_units} if var_units else None
        plot_height = shared.figure_height_for_size(global_context.get("size"))
        column_mode = global_context.get("column_mode") or "single"
        fig = go.Figure(
            data=go.Heatmap(
                x=time_vals,
                y=z_vals,
                z=data.T,
                colorscale="Viridis",
                colorbar=colorbar,
            )
        )
        shared.apply_figure_chrome(
            fig,
            title=shared.format_plot_title(var_name, var_long_name),
            showlegend=False,
            height=plot_height,
            uirevision=shared.figure_uirevision(self.plot_type_id, case_data, var_name, column_mode),
        )
        fig.update_layout(
            xaxis_title=shared.format_axis_title("time", "s"),
            yaxis_title=shared.format_height_axis_title(z_units),
        )
        x_range = shared.padded_data_range(time_vals)
        if x_range:
            shared.apply_axis_bounds(fig, "x", x_range)
        height_range = global_context.get("height_range") or case_data.get("default_height_range")
        if height_range:
            shared.apply_axis_bounds(fig, "y", shared.padded_range(height_range[0], height_range[1]))
        shared.apply_plot_theme(fig, theme_name)
        return fig

    def _build_panel_figure(self, panels, field_name, global_context):
        case_data = global_context.get("case_data") or {}
        theme_name = global_context.get("theme_name")
        plot_height = max(shared.figure_height_for_size(global_context.get("size")), 260 * len(panels))
        fig = make_subplots(
            rows=len(panels),
            cols=1,
            shared_xaxes=True,
            shared_yaxes=True,
            vertical_spacing=0.055 if len(panels) > 1 else 0.0,
            subplot_titles=[panel["title"] for panel in panels],
        )
        z_arrays = [np.asarray(panel["data"][2], dtype=float) for panel in panels]
        finite_parts = [arr[np.isfinite(arr)] for arr in z_arrays if np.any(np.isfinite(arr))]
        finite = np.concatenate(finite_parts) if finite_parts else np.array([], dtype=float)
        zmin = float(np.min(finite)) if finite.size else None
        zmax = float(np.max(finite)) if finite.size else None
        var_units = panels[0]["data"][3]
        var_long_name = panels[0]["data"][4]
        z_units = panels[0]["data"][5]
        colorbar = {"title": var_units} if var_units else None
        for row, panel in enumerate(panels, start=1):
            time_vals, z_vals, data, _units, _long_name, _z_units = panel["data"]
            fig.add_trace(
                go.Heatmap(
                    x=time_vals,
                    y=z_vals,
                    z=np.asarray(data, dtype=float).T,
                    coloraxis="coloraxis",
                ),
                row=row,
                col=1,
            )
        coloraxis = {"colorscale": "Viridis", "colorbar": colorbar or {}}
        if zmin is not None and zmax is not None:
            coloraxis["cmin"] = zmin
            coloraxis["cmax"] = zmax
        shared.apply_figure_chrome(
            fig,
            title=shared.format_plot_title(field_name, var_long_name),
            showlegend=False,
            height=plot_height,
            uirevision=shared.figure_uirevision(
                self.plot_type_id,
                case_data,
                field_name,
                global_context.get("column_mode") or "single",
                ",".join(global_context.get("enabled_benchmarks") or []),
            ),
        )
        fig.update_layout(coloraxis=coloraxis)
        fig.update_xaxes(title_text=shared.format_axis_title("time", "s"), row=len(panels), col=1)
        for row in range(1, len(panels) + 1):
            fig.update_yaxes(title_text=shared.format_height_axis_title(z_units), row=row, col=1)
        x_bounds = [shared.padded_data_range(panel["data"][0]) for panel in panels]
        finite_x_bounds = [bounds for bounds in x_bounds if bounds]
        if finite_x_bounds:
            shared.apply_axis_bounds(
                fig,
                "x",
                [min(bounds[0] for bounds in finite_x_bounds), max(bounds[1] for bounds in finite_x_bounds)],
            )
        height_range = global_context.get("height_range") or case_data.get("default_height_range")
        if height_range:
            shared.apply_axis_bounds(fig, "y", shared.padded_range(height_range[0], height_range[1]))
        shared.apply_plot_theme(fig, theme_name)
        return fig

    def register_callbacks(self, app):
        @app.callback(
            Output(self.graph_id(MATCH), "figure"),
            Input(self.var_input_id(MATCH), "value"),
            Input("plots-case-data", "data"),
            Input("plots-global-height-range", "value"),
            Input("plots-global-time-range", "value"),
            Input("plots-selected-column", "data"),
            Input("plots-column-mode", "value"),
            Input("plots-enabled-benchmarks", "data"),
            Input("theme-store", "data"),
            Input(self.size_store_id(MATCH), "data"),
        )
        def _update_timeheight_graph(var_name, case_data, height_range, time_range, selected_column, column_mode, enabled_benchmarks, theme_name, size_store_value):
            size_value = shared.normalize_plot_size(size_store_value)
            return self.build_figure(
                {"var": var_name, "size": size_value},
                {
                    "case_data": case_data,
                    "height_range": height_range,
                    "time_range": time_range,
                    "selected_column": selected_column,
                    "column_mode": column_mode,
                    "enabled_benchmarks": enabled_benchmarks,
                    "size": size_value,
                    "theme_name": theme_name,
                },
            )


PLOT = TimeHeightPlotType()
