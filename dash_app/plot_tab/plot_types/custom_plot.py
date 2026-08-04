import ast
from math import erf as _scalar_erf

import numpy as np
import plotly.graph_objects as go
from dash import Input, MATCH, Output, State, callback_context, dcc, html

from .. import benchmark_overlay
from . import shared
from .base_plot import BasePlotType


DEFAULT_EXPRESSION = "wp2 / wp3"


def _erf(values):
    values = np.asarray(values, dtype=float)
    if values.size == 0:
        return np.empty_like(values)
    return np.vectorize(_scalar_erf, otypes=[float])(values)


CUSTOM_FUNCTIONS = {
    "abs": np.abs,
    "arccos": np.arccos,
    "arcsin": np.arcsin,
    "arctan": np.arctan,
    "ceil": np.ceil,
    "clip": np.clip,
    "sqrt": np.sqrt,
    "log": np.log,
    "log1p": np.log1p,
    "log10": np.log10,
    "exp": np.exp,
    "erf": _erf,
    "floor": np.floor,
    "sin": np.sin,
    "cos": np.cos,
    "tan": np.tan,
    "minimum": np.minimum,
    "maximum": np.maximum,
    "sign": np.sign,
    "square": np.square,
    "where": np.where,
}
CUSTOM_CONSTANTS = {
    "pi": np.pi,
    "e": np.e,
    "nan": np.nan,
}
NUMPY_ALIASES = {"np", "numpy"}
LES_NAMESPACE = "les"


class ExpressionError(ValueError):
    pass


class _NameCollector(ast.NodeVisitor):
    def __init__(self):
        self.names = set()

    def visit_Call(self, node):
        if _function_name(node.func) not in CUSTOM_FUNCTIONS:
            raise ExpressionError("Unsupported function in custom expression.")
        if node.keywords:
            raise ExpressionError("Keyword arguments are not supported in custom expressions.")
        for arg in node.args:
            self.visit(arg)

    def visit_Attribute(self, node):
        les_name = _les_variable_name(node)
        if les_name is not None:
            self.names.add(les_name)
            return
        if _constant_name(node) is None:
            raise ExpressionError(
                "Only les.<field>, np.<constant>, and np.<function>(...) attribute syntax is supported."
            )

    def visit_Name(self, node):
        if node.id == LES_NAMESPACE:
            raise ExpressionError("Reference an LES field as les.<field>.")
        if node.id not in CUSTOM_CONSTANTS and node.id not in NUMPY_ALIASES:
            self.names.add(node.id)


def _function_name(func_node):
    if isinstance(func_node, ast.Name):
        return func_node.id
    if (
        isinstance(func_node, ast.Attribute)
        and isinstance(func_node.value, ast.Name)
        and func_node.value.id in NUMPY_ALIASES
    ):
        return func_node.attr
    return None


def _constant_name(node):
    if isinstance(node, ast.Name) and node.id in CUSTOM_CONSTANTS:
        return node.id
    if (
        isinstance(node, ast.Attribute)
        and isinstance(node.value, ast.Name)
        and node.value.id in NUMPY_ALIASES
        and node.attr in CUSTOM_CONSTANTS
    ):
        return node.attr
    return None


def _les_variable_name(node):
    if (
        isinstance(node, ast.Attribute)
        and isinstance(node.value, ast.Name)
        and node.value.id == LES_NAMESPACE
        and not node.attr.startswith("_")
    ):
        return f"{LES_NAMESPACE}.{node.attr}"
    return None


def _partition_variable_names(var_names):
    les_prefix = f"{LES_NAMESPACE}."
    clubb_names = [name for name in var_names if not name.startswith(les_prefix)]
    les_names = [name[len(les_prefix) :] for name in var_names if name.startswith(les_prefix)]
    return clubb_names, les_names


def expression_variable_names(expression):
    """Return variable names referenced by a safe custom expression."""
    text = str(expression or "").strip()
    if not text:
        raise ExpressionError("Enter a custom expression.")
    try:
        tree = ast.parse(text, mode="eval")
    except SyntaxError as exc:
        raise ExpressionError(f"Invalid expression: {exc.msg}") from exc
    collector = _NameCollector()
    collector.visit(tree)
    if not collector.names:
        raise ExpressionError("Expression must reference at least one NetCDF variable.")
    return tree, sorted(collector.names)


def _eval_expression_node(node, values):
    if isinstance(node, ast.Expression):
        return _eval_expression_node(node.body, values)
    if isinstance(node, ast.Constant):
        if isinstance(node.value, (int, float)):
            return float(node.value)
        raise ExpressionError("Only numeric constants are allowed.")
    if isinstance(node, ast.Name):
        constant_name = _constant_name(node)
        if constant_name is not None:
            return CUSTOM_CONSTANTS[constant_name]
        if node.id in NUMPY_ALIASES:
            raise ExpressionError("Use np.<function>(...) or np.<constant> in custom expressions.")
        if node.id not in values:
            raise ExpressionError(f"{node.id} is not available.")
        return values[node.id]
    if isinstance(node, ast.Attribute):
        constant_name = _constant_name(node)
        if constant_name is not None:
            return CUSTOM_CONSTANTS[constant_name]
        les_name = _les_variable_name(node)
        if les_name is not None:
            if les_name not in values:
                raise ExpressionError(f"{les_name} is not available.")
            return values[les_name]
        raise ExpressionError(
            "Only les.<field> and np.<constant> attribute syntax is supported outside function calls."
        )
    if isinstance(node, ast.UnaryOp):
        operand = _eval_expression_node(node.operand, values)
        if isinstance(node.op, ast.UAdd):
            return operand
        if isinstance(node.op, ast.USub):
            return -operand
        raise ExpressionError("Unsupported unary operator.")
    if isinstance(node, ast.BinOp):
        left = _eval_expression_node(node.left, values)
        right = _eval_expression_node(node.right, values)
        with np.errstate(all="ignore"):
            if isinstance(node.op, ast.Add):
                return left + right
            if isinstance(node.op, ast.Sub):
                return left - right
            if isinstance(node.op, ast.Mult):
                return left * right
            if isinstance(node.op, ast.Div):
                return left / right
            if isinstance(node.op, ast.Pow):
                return left ** right
            if isinstance(node.op, ast.Mod):
                return left % right
        raise ExpressionError("Unsupported binary operator.")
    if isinstance(node, ast.Call):
        function_name = _function_name(node.func)
        if function_name not in CUSTOM_FUNCTIONS:
            raise ExpressionError("Unsupported function in custom expression.")
        if node.keywords:
            raise ExpressionError("Keyword arguments are not supported in custom expressions.")
        args = [_eval_expression_node(arg, values) for arg in node.args]
        with np.errstate(all="ignore"):
            return CUSTOM_FUNCTIONS[function_name](*args)
    raise ExpressionError("Unsupported expression syntax.")


def _interpolate_profiles(profiles, source_z, target_z):
    profiles = np.asarray(profiles, dtype=float)
    source_z = np.asarray(source_z, dtype=float)
    target_z = np.asarray(target_z, dtype=float)
    if profiles.ndim != 2:
        raise ExpressionError("Custom expression variables must extract to profile arrays.")
    if source_z.shape == target_z.shape and np.allclose(source_z, target_z, equal_nan=True):
        return profiles
    order = np.argsort(source_z)
    source_sorted = source_z[order]
    interpolated = []
    for profile in profiles:
        profile_sorted = np.asarray(profile, dtype=float)[order]
        mask = np.isfinite(source_sorted) & np.isfinite(profile_sorted)
        if np.count_nonzero(mask) < 2:
            interpolated.append(np.full(target_z.shape, np.nan, dtype=float))
            continue
        interpolated.append(
            np.interp(
                target_z,
                source_sorted[mask],
                profile_sorted[mask],
                left=np.nan,
                right=np.nan,
            )
        )
    return np.asarray(interpolated, dtype=float)


def _broadcast_lines(profiles, target_count):
    profiles = np.asarray(profiles, dtype=float)
    if profiles.shape[0] == target_count:
        return profiles
    if profiles.shape[0] == 1 and target_count > 1:
        return np.repeat(profiles, target_count, axis=0)
    raise ExpressionError("Expression variables do not have compatible column selections.")


class CustomPlotType(BasePlotType):
    def __init__(self):
        super().__init__(
            plot_type_id="custom",
            default_vars=[],
            var_input_type="custom-expression",
            graph_type="custom-graph",
            case_data_var_key="profile_vars",
            subtitle=(
                "Type an arithmetic expression over profile variables, for example wp2 / wp3, "
                "np.abs(wp2), or wpthlp + les.wprcp. Prefix a field with les. to use the LES benchmark."
            ),
        )

    def case_data_options(self, case_data):
        return list((case_data or {}).get("profile_vars") or [])

    def get_variable_options(self, collection, case_data):
        require_all = bool(case_data.get("compare_mode"))
        return case_data.get("profile_vars") or shared.list_profile_vars(collection, require_all=require_all)

    def make_default_state(self, case_data, plot_id):
        return {
            "plot_type": self.plot_type_id,
            "expression": DEFAULT_EXPRESSION,
            "size": "normal",
        }

    def render_card(self, plot_id, state, global_context):
        case_data = global_context.get("case_data") or {}
        size_value = shared.normalize_plot_size(state.get("size"))
        expression = state.get("expression") or DEFAULT_EXPRESSION
        controls = dcc.Input(
            id=self.var_input_id(plot_id),
            type="text",
            value=expression,
            debounce=True,
            className="clubb-input custom-expression-input",
            placeholder="wp2 / wp3",
            style={"width": "100%", "boxSizing": "border-box", "fontSize": "14px", "padding": "6px 8px"},
        )
        size_text, size_class = shared.plot_size_button_props(size_value)
        size_button = html.Button(size_text, id=self.size_toggle_id(plot_id), className=size_class, title="Toggle plot size")
        return shared.make_plot_card(
            subtitle=self.card_subtitle(case_data),
            help_button_id=self.help_button_id(plot_id),
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

    def _time_indices(self, global_context):
        case_data = global_context.get("case_data") or {}
        slider_range = global_context.get("time_range")
        time_mode = global_context.get("time_mode") or "range"
        time_point = global_context.get("time_point")
        return shared.selected_time_indices(case_data, slider_range, time_point, time_mode)

    def _path_has_profile_vars(self, path, var_names):
        meta = shared.dataset_metadata_for_path(path)
        for name in var_names:
            info = meta["var_info"].get(name)
            if not info or info.get("t_dim") is None or info.get("z_dim") is None:
                return False
        return True

    def _les_source_name(self, case_data):
        available = benchmark_overlay.available_benchmark_sources(case_data)
        if "sam" in available:
            return "sam"
        if len(available) == 1:
            return available[0]
        if not available:
            raise ExpressionError("les.<field> requires an LES benchmark configured for this case.")
        raise ExpressionError(
            "les.<field> is ambiguous because this case has multiple non-SAM benchmark sources."
        )

    def _extract_les_profiles(self, les_names, global_context):
        if not les_names:
            return {}, None
        case_data = global_context.get("case_data") or {}
        source_name = self._les_source_name(case_data)
        supported = set(benchmark_overlay.normalized_benchmark_supported_fields())
        unsupported = sorted(set(les_names) - supported)
        if unsupported:
            names = ", ".join(f"les.{name}" for name in unsupported)
            raise ExpressionError(f"The LES converter does not support: {names}.")
        profiles = benchmark_overlay.extract_normalized_benchmark_profiles(
            case_data,
            source_name,
            les_names,
            global_context.get("time_mode") or "range",
            global_context.get("time_range"),
            global_context.get("time_point"),
        )
        missing = [name for name in les_names if name not in profiles]
        if missing:
            names = ", ".join(f"les.{name}" for name in missing)
            raise ExpressionError(
                f"Could not create {names} from the configured {source_name.upper()} benchmark."
            )
        return profiles, source_name

    def _extract_expression_profiles(self, path, expression, tree, var_names, global_context):
        time_indices = self._time_indices(global_context)
        col_index = int(global_context.get("selected_column") or 0)
        column_mode = global_context.get("column_mode") or "single"
        extracted = {}
        clubb_names, les_names = _partition_variable_names(var_names)
        for name in clubb_names:
            result = shared.extract_time_avg_profile_for_path(
                path,
                name,
                time_indices,
                col_index=col_index,
                column_mode=column_mode,
                column_filter_indices=global_context.get("column_filter_indices"),
            )
            if result is None:
                raise ExpressionError(f"{name} is not available as a profile variable.")
            extracted[name] = result
        les_profiles, les_source = self._extract_les_profiles(les_names, global_context)
        for name, result in les_profiles.items():
            extracted[f"{LES_NAMESPACE}.{name}"] = {
                "profiles": np.asarray(result["profile"], dtype=float)[np.newaxis, :],
                "z_values": np.asarray(result["z_values"], dtype=float),
                "labels": [f"LES ({les_source.upper()})"],
                "z_units": result.get("z_units", "") or "",
            }
        reference_name = clubb_names[0] if clubb_names else f"{LES_NAMESPACE}.{les_names[0]}"
        reference = extracted[reference_name]
        target_z = np.asarray(reference["z_values"], dtype=float)
        line_count = int(np.asarray(reference["profiles"]).shape[0])
        values = {}
        for name, result in extracted.items():
            profiles = _interpolate_profiles(result["profiles"], result["z_values"], target_z)
            values[name] = _broadcast_lines(profiles, line_count)
        output = np.asarray(_eval_expression_node(tree, values), dtype=float)
        if output.ndim == 1:
            output = output[np.newaxis, :]
        output = _broadcast_lines(output, line_count)
        return {
            "profiles": output,
            "z_values": target_z,
            "labels": reference["labels"],
            "z_units": reference["z_units"],
            "expression": expression,
        }

    def _single_trace_specs(self, files, expression, case_data, global_context):
        tree, var_names = expression_variable_names(expression)
        clubb_names, _les_names = _partition_variable_names(var_names)
        path = next(
            (candidate for candidate in files or [] if self._path_has_profile_vars(candidate, clubb_names)),
            None,
        )
        if path is None:
            raise ExpressionError("No selected stats file contains all CLUBB expression variables as profiles.")
        extracted = self._extract_expression_profiles(path, expression, tree, var_names, global_context)
        column_mode = global_context.get("column_mode") or "single"
        trace_specs = []
        x_values = []
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
                    "name": extracted["labels"][idx] if idx < len(extracted["labels"]) else f"line {idx + 1}",
                    "line": line,
                    "opacity": 0.85 if column_mode == "all" else None,
                    "showlegend": False,
                }
            )
            x_values.extend([value for value in profile if np.isfinite(value)])
        height_range = shared.active_height_range(global_context)
        return {
            "trace_specs": trace_specs,
            "x_range": shared.padded_trace_x_range(trace_specs, height_range, fallback_values=x_values),
            "z_units": extracted["z_units"],
            "column_mode": column_mode,
        }

    def _compare_trace_specs(self, files, expression, case_data, global_context):
        tree, var_names = expression_variable_names(expression)
        clubb_names, _les_names = _partition_variable_names(var_names)
        trace_specs = []
        x_values = []
        legend_labels = []
        source_labels = case_data.get("source_labels") or [f"output {idx + 1}" for idx in range(len(files or []))]
        column_mode = global_context.get("column_mode") or "single"
        z_units = ""
        for source_idx, path in enumerate(files or []):
            if not self._path_has_profile_vars(path, clubb_names):
                continue
            extracted = self._extract_expression_profiles(path, expression, tree, var_names, global_context)
            z_units = extracted["z_units"] or z_units
            label = source_labels[source_idx] if source_idx < len(source_labels) else f"output {source_idx + 1}"
            legend_labels.append(label)
            for idx, profile in enumerate(extracted["profiles"]):
                line = {"width": 1.2 if column_mode == "all" else 2.0, "color": shared.source_line_color(source_idx), "dash": shared.source_line_dash(source_idx)}
                if column_mode == "all":
                    line["color"] = shared.column_line_color(idx)
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
                x_values.extend([value for value in profile if np.isfinite(value)])
        height_range = shared.active_height_range(global_context)
        return {
            "trace_specs": trace_specs,
            "x_range": shared.padded_trace_x_range(trace_specs, height_range, fallback_values=x_values),
            "z_units": z_units,
            "legend_labels": legend_labels,
            "column_mode": column_mode,
        }

    def build_figure(self, state, global_context):
        case_data = global_context.get("case_data") or {}
        theme_name = global_context.get("theme_name")
        files = case_data.get("files") or []
        expression = str(state.get("expression") or "").strip()
        if not files or not expression:
            return shared.make_empty_figure("Select a case and enter a custom expression.", theme_name)
        try:
            if case_data.get("compare_mode"):
                bundle = self._compare_trace_specs(files, expression, case_data, global_context)
            else:
                bundle = self._single_trace_specs(files, expression, case_data, global_context)
        except ExpressionError as exc:
            return shared.make_empty_figure(str(exc), theme_name)
        if not bundle["trace_specs"]:
            return shared.make_empty_figure("No custom-expression data are available for this case.", theme_name)
        fig = go.Figure()
        for spec in bundle["trace_specs"]:
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
        if case_data.get("compare_mode") and bundle["column_mode"] == "all":
            shared.add_directory_legend_traces(fig, bundle["legend_labels"])
        plot_height = shared.figure_height_for_size(global_context.get("size"))
        shared.apply_figure_chrome(
            fig,
            title=f"Custom: {expression}",
            showlegend=bool(case_data.get("compare_mode")),
            height=plot_height,
            uirevision=shared.figure_uirevision(
                self.plot_type_id,
                case_data,
                expression,
                bundle["column_mode"],
            ),
        )
        fig.update_layout(
            xaxis_title=expression,
            yaxis_title=shared.format_height_axis_title(bundle["z_units"]),
        )
        if bundle["x_range"]:
            fig.update_xaxes(range=bundle["x_range"])
        height_range = global_context.get("height_range") or case_data.get("default_height_range")
        if height_range:
            shared.apply_axis_bounds(fig, "y", shared.padded_range(height_range[0], height_range[1]))
        shared.apply_plot_theme(fig, theme_name)
        return fig

    def register_callbacks(self, app):
        @app.callback(
            Output(self.graph_id(MATCH), "figure"),
            Output(self.render_signal_id(MATCH), "children"),
            Input(self.var_input_id(MATCH), "value"),
            Input("plots-case-data", "data"),
            Input("plots-global-time-range", "value"),
            Input("plots-global-time-point", "value"),
            Input("plots-time-override", "data"),
            Input("plots-global-height-range", "value"),
            Input("plots-selected-column", "data"),
            Input("plots-column-mode", "value"),
            Input("plots-column-filters", "data"),
            Input("theme-store", "data"),
            Input(self.size_store_id(MATCH), "data"),
            State(self.graph_id(MATCH), "relayoutData"),
        )
        def _update_custom_graph(
            expression,
            case_data,
            time_range,
            time_point,
            time_override,
            height_range,
            selected_column,
            column_mode,
            column_filters,
            theme_name,
            size_store_value,
            relayout_data,
        ):
            size_value = shared.normalize_plot_size(size_store_value)
            active_time = shared.resolve_active_time_values(case_data, time_range, time_point, time_override)
            signal = int(active_time["start_seconds"]) if active_time["start_seconds"] is not None else ""
            fig = self.build_figure(
                {"expression": expression, "size": size_value},
                {
                    "case_data": case_data,
                    "time_range": active_time["duration_minutes"],
                    "time_point": active_time["start_seconds"],
                    "height_range": height_range,
                    "relayout_data": relayout_data,
                    "selected_column": selected_column,
                    "column_mode": column_mode,
                    "column_filter_indices": shared.column_filter_indices(column_filters),
                    "size": size_value,
                    "theme_name": theme_name,
                },
            )
            if callback_context.triggered_id == "plots-case-data" and (case_data or {}).get("preserve_plot_view"):
                shared.apply_relayout_ranges(fig, relayout_data)
            return fig, signal


PLOT = CustomPlotType()
