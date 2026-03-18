from dash import dcc, html

from . import shared


class BasePlotType:
    def __init__(
        self,
        plot_type_id,
        default_vars,
        var_input_type,
        graph_type,
        case_data_var_key,
        subtitle=None,
    ):
        self.plot_type_id = plot_type_id
        self.default_vars = list(default_vars)
        self.var_input_type = var_input_type
        self.graph_type = graph_type
        self.case_data_var_key = case_data_var_key
        self._subtitle = subtitle

    def get_variable_options(self, collection, case_data):
        raise NotImplementedError

    def supports_case(self, case_data, collection):
        return bool(self.get_variable_options(collection, case_data))

    def var_input_id(self, plot_id):
        return {"type": self.var_input_type, "index": plot_id}

    def graph_id(self, plot_id):
        return {"type": self.graph_type, "index": plot_id}

    def graph_shell_id(self, plot_id):
        return {"type": "plots-graph-shell", "index": plot_id}

    def card_id(self, plot_id):
        return {"type": "plots-card", "index": plot_id}

    def close_button_id(self, plot_id):
        return {"type": "plots-close", "index": plot_id}

    def size_toggle_id(self, plot_id):
        return {"type": "plots-size-toggle", "index": plot_id}

    def size_store_id(self, plot_id):
        return {"type": "plots-size-store", "index": plot_id}

    def render_signal_id(self, plot_id):
        return {"type": f"{self.plot_type_id}-render-signal", "index": plot_id}

    def case_data_options(self, case_data):
        return (case_data or {}).get(self.case_data_var_key, [])

    def default_value(self, options, plot_id):
        return shared.pick_default_var(options, self.default_vars, fallback_index=plot_id)

    def supports_case_data(self, case_data):
        return bool(self.case_data_options(case_data))

    def make_default_state(self, case_data, plot_id):
        options = list(self.case_data_options(case_data))
        return {
            "plot_type": self.plot_type_id,
            "var": self.default_value(options, plot_id),
            "size": "normal",
        }

    def card_subtitle(self, case_data):
        return self._subtitle

    def auxiliary_components(self, plot_id):
        return []

    def render_card(self, plot_id, state, global_context):
        case_data = global_context.get("case_data") or {}
        options = list(self.case_data_options(case_data))
        value = state.get("var") or self.default_value(options, plot_id)
        size_value = shared.normalize_plot_size(state.get("size"))
        option_values = {option["value"] for option in options}
        if value is not None and value not in option_values:
            options = [{"label": f"{value} (missing)", "value": value}] + options
        controls = dcc.Dropdown(
            id=self.var_input_id(plot_id),
            options=options,
            value=value,
            clearable=False,
        )
        size_text, size_class = shared.plot_size_button_props(size_value)
        size_button = html.Button(size_text, id=self.size_toggle_id(plot_id), className=size_class, title="Toggle plot size")
        return shared.make_plot_card(
            subtitle=self.card_subtitle(case_data),
            controls=controls,
            size_button=size_button,
            size_value=size_value,
            size_store_id=self.size_store_id(plot_id),
            graph_id=self.graph_id(plot_id),
            graph_shell_id=self.graph_shell_id(plot_id),
            card_id=self.card_id(plot_id),
            close_button_id=self.close_button_id(plot_id),
            render_signal_id=self.render_signal_id(plot_id),
        )

    def register_callbacks(self, app):
        raise NotImplementedError
