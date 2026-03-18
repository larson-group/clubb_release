import plotly.graph_objects as go
from dash import Input, Output, State, MATCH

from . import shared
from .base_plot import BasePlotType


class TimeSeriesPlotType(BasePlotType):
    def __init__(self):
        super().__init__(
            plot_type_id="timeseries",
            default_vars=shared.DEFAULT_TIMESERIES_VARS,
            var_input_type="timeseries-var",
            graph_type="timeseries-graph",
            case_data_var_key="timeseries_vars",
        )

    def get_variable_options(self, collection, case_data):
        return case_data.get("timeseries_vars") or shared.list_timeseries_vars(collection, require_all=bool(case_data.get("compare_mode")))

    def card_subtitle(self, case_data):
        subtitle = "Uses the full file time extent."
        if case_data.get("compare_mode"):
            return subtitle + " In compare mode, it plots one line per output or all columns from all outputs."
        return subtitle + " In all-column mode, it plots one line per column."

    def build_figure(self, state, global_context):
        case_data = global_context.get("case_data") or {}
        theme_name = global_context.get("theme_name")
        files = case_data.get("files") or []
        var_name = state.get("var")
        if not files or not var_name:
            return shared.make_empty_figure("Select a case and variable.", theme_name)
        column_mode = global_context.get("column_mode") or "single"
        col_index = int(global_context.get("selected_column") or 0)
        compare_mode = bool(case_data.get("compare_mode"))
        plot_height = shared.figure_height_for_size(global_context.get("size"))
        source_labels = case_data.get("source_labels") or [f"output {idx + 1}" for idx in range(len(files))]
        fig = go.Figure()
        var_units = ""
        time_units = shared.format_axis_title("time", "s")
        var_long_name = ""
        legend_labels = []
        for source_idx, path in enumerate(files):
            result = shared.extract_timeseries_for_path(path, var_name, col_index=col_index, column_mode=column_mode)
            if result is None:
                continue
            time_vals, lines, labels, units, long_name, time_units_raw = result
            var_units = units or var_units
            var_long_name = long_name or var_long_name
            label = source_labels[source_idx] if source_idx < len(source_labels) else f"output {source_idx + 1}"
            legend_labels.append(label)
            if time_units_raw:
                factor = shared.time_units_factor(time_units_raw)
                unit_name = "s"
                if factor == 60.0:
                    unit_name = "min"
                elif factor == 3600.0:
                    unit_name = "h"
                time_units = shared.format_axis_title("time", unit_name)
            for line_idx in range(lines.shape[1]):
                if compare_mode and column_mode == "all":
                    trace_name = label
                    showlegend = False
                elif compare_mode:
                    trace_name = label
                    showlegend = True
                else:
                    trace_name = labels[line_idx]
                    showlegend = False
                line_style = {"width": 1.5 if column_mode == "all" else 2.0}
                if column_mode == "all":
                    line_style["color"] = shared.column_line_color(line_idx)
                if compare_mode:
                    line_style["dash"] = shared.source_line_dash(source_idx)
                fig.add_trace(
                    go.Scatter(x=time_vals, y=lines[:, line_idx], mode="lines", name=trace_name, line=line_style, opacity=0.9, showlegend=showlegend)
                )
        if not fig.data:
            return shared.make_empty_figure(f"{var_name} is not compatible with a time-series plot.", theme_name)
        if compare_mode and column_mode == "all":
            shared.add_directory_legend_traces(fig, legend_labels)
        shared.apply_figure_chrome(
            fig,
            title=shared.format_plot_title(var_name, var_long_name),
            showlegend=compare_mode,
            height=plot_height,
            uirevision=shared.figure_uirevision(self.plot_type_id, case_data, var_name, column_mode),
        )
        fig.update_layout(
            xaxis_title=time_units,
            yaxis_title=shared.format_axis_title(var_name, var_units),
        )
        x_values = []
        y_values = []
        for trace in fig.data:
            x_values.extend([value for value in trace.x if value is not None])
            y_values.extend([value for value in trace.y if value is not None])
        x_range = shared.padded_data_range(x_values)
        y_range = shared.padded_data_range(y_values)
        if x_range:
            shared.apply_axis_bounds(fig, "x", x_range)
        if y_range:
            fig.update_yaxes(range=y_range)
        shared.apply_plot_theme(fig, theme_name)
        return fig

    def register_callbacks(self, app):
        @app.callback(
            Output(self.graph_id(MATCH), "figure"),
            Input(self.var_input_id(MATCH), "value"),
            Input("plots-case-data", "data"),
            Input("plots-selected-column", "data"),
            Input("plots-column-mode", "value"),
            Input("theme-store", "data"),
            Input(self.size_store_id(MATCH), "data"),
        )
        def _update_timeseries_graph(var_name, case_data, selected_column, column_mode, theme_name, size_store_value):
            size_value = shared.normalize_plot_size(size_store_value)
            return self.build_figure(
                {"var": var_name, "size": size_value},
                {
                    "case_data": case_data,
                    "selected_column": selected_column,
                    "column_mode": column_mode,
                    "size": size_value,
                    "theme_name": theme_name,
                },
            )


PLOT = TimeSeriesPlotType()
