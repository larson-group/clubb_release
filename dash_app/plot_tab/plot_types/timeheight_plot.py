import plotly.graph_objects as go
from dash import Input, Output, State, MATCH

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
            subtitle="Uses the full file time extent. In all-column mode, any column dimension is averaged before plotting.",
        )

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
        time_vals, z_vals, data, var_units, var_long_name, z_units = result
        colorbar = {"title": var_units} if var_units else None
        plot_height = shared.figure_height_for_size(global_context.get("size"))
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

    def register_callbacks(self, app):
        @app.callback(
            Output(self.graph_id(MATCH), "figure"),
            Input(self.var_input_id(MATCH), "value"),
            Input("plots-case-data", "data"),
            Input("plots-global-height-range", "value"),
            Input("plots-selected-column", "data"),
            Input("plots-column-mode", "value"),
            Input("theme-store", "data"),
            Input(self.size_store_id(MATCH), "data"),
        )
        def _update_timeheight_graph(var_name, case_data, height_range, selected_column, column_mode, theme_name, size_store_value):
            size_value = shared.normalize_plot_size(size_store_value)
            return self.build_figure(
                {"var": var_name, "size": size_value},
                {
                    "case_data": case_data,
                    "height_range": height_range,
                    "selected_column": selected_column,
                    "column_mode": column_mode,
                    "size": size_value,
                    "theme_name": theme_name,
                },
            )


PLOT = TimeHeightPlotType()
