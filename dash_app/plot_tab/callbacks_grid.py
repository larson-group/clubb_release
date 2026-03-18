import re

from dash import ALL, MATCH, Input, Output, Patch, State, callback_context, no_update

from .layout import add_plot_controls_card, child_id, render_plot_card, render_plot_grid
from .plot_types.registry import PLOT_TYPES
from .plot_types.specs import PLOT_FAMILY_SPECS
from .plot_types import shared
from .state import default_plot_state


_SEARCH_TOKEN_RE = re.compile(r"[_\-\s]+")


def _normalize_search_text(value):
    return str(value or "").strip().lower()


def _search_tokens(value):
    normalized = _normalize_search_text(value)
    if not normalized:
        return []
    return [token for token in _SEARCH_TOKEN_RE.split(normalized) if token]


def _option_text(option):
    label = option.get("label")
    value = option.get("value")
    return _normalize_search_text(label if label is not None else value)


def _option_rank(option, search_text, search_tokens):
    text = _option_text(option)
    value_text = _normalize_search_text(option.get("value"))
    if not search_text:
        return (99, len(text), text)
    if text == search_text or value_text == search_text:
        return (0, len(text), text)
    tokens = _search_tokens(text)
    if search_text in tokens:
        return (1, len(text), text)
    if search_tokens and all(token in tokens for token in search_tokens):
        return (2, len(text), text)
    if text.startswith(search_text) or value_text.startswith(search_text):
        return (3, len(text), text)
    if any(token.startswith(search_text) for token in tokens):
        return (4, len(text), text)
    if search_text in text or search_text in value_text:
        return (5, len(text), text)
    if search_tokens and all(token in text for token in search_tokens):
        return (6, len(text), text)
    return (99, len(text), text)


def _ensure_selected_option(options, selected_value):
    if selected_value is None:
        return list(options)
    option_values = {option.get("value") for option in options}
    if selected_value in option_values:
        return list(options)
    return [{"label": f"{selected_value} (missing)", "value": selected_value}] + list(options)


def _ranked_dropdown_options(options, search_value, selected_value):
    options = _ensure_selected_option(options, selected_value)
    search_text = _normalize_search_text(search_value)
    if not search_text:
        return options
    search_tokens = _search_tokens(search_text)
    ranked = sorted(
        options,
        key=lambda option: (_option_rank(option, search_text, search_tokens), _option_text(option)),
    )
    matched = [option for option in ranked if _option_rank(option, search_text, search_tokens)[0] < 99]
    if not matched:
        return ranked
    if selected_value is not None and all(option.get("value") != selected_value for option in matched):
        selected = next((option for option in ranked if option.get("value") == selected_value), None)
        if selected is not None:
            matched = [selected] + matched
    return matched


def register_dropdown_search_callback(app, spec):
    """Register one dropdown-options callback that reorders options by search relevance."""
    module = PLOT_TYPES[spec.plot_type_id]

    @app.callback(
        Output({"type": spec.dropdown_type, "index": MATCH}, "options"),
        Input({"type": spec.dropdown_type, "index": MATCH}, "search_value"),
        Input("plots-case-data", "data"),
        State({"type": spec.dropdown_type, "index": MATCH}, "value"),
    )
    def _update_dropdown_options(search_value, case_data, current_value):
        options = list(module.case_data_options(case_data))
        return _ranked_dropdown_options(options, search_value, current_value)

    return _update_dropdown_options


def register_plot_state_sync_callback(app, plot_type_id, dropdown_type):
    """Register one generic dropdown-to-plot-state synchronization callback."""
    @app.callback(
        Output("plots-plot-state", "data", allow_duplicate=True),
        Input({"type": dropdown_type, "index": ALL}, "value"),
        State({"type": dropdown_type, "index": ALL}, "id"),
        State("plots-plot-state", "data"),
        prevent_initial_call=True,
    )
    def _sync_plot_state(values, ids, current_state):
        """Persist plot type and selected variable for one plot-family dropdown set."""
        updated = dict(current_state or {})
        for meta, value in zip(ids or [], values or []):
            idx = str(meta.get("index"))
            if idx not in updated:
                continue
            entry = dict(updated[idx])
            entry["plot_type"] = plot_type_id
            if value is not None:
                entry["var"] = value
            updated[idx] = entry
        return updated

    return _sync_plot_state


def register_grid_callbacks(app):
    """Register plot-grid rendering plus add/remove and dropdown state callbacks."""
    @app.callback(
        Output("plots-plot-container", "children"),
        Input("plots-case-data", "data"),
        State("plots-plot-order", "data"),
        State("plots-plot-state", "data"),
    )
    def render_plot_container(case_data, plot_order, plot_state):
        """Render the plot grid from the current ordered plot state."""
        return render_plot_grid(plot_order or [], plot_state or {}, case_data)

    for spec in PLOT_FAMILY_SPECS.values():
        register_plot_state_sync_callback(app, spec.plot_type_id, spec.dropdown_type)
        register_dropdown_search_callback(app, spec)

    @app.callback(
        Output({"type": "plots-size-store", "index": MATCH}, "data"),
        Input({"type": "plots-size-toggle", "index": MATCH}, "n_clicks_timestamp"),
        State({"type": "plots-size-store", "index": MATCH}, "data"),
        prevent_initial_call=True,
    )
    def sync_plot_sizes(click_ts, current_size):
        """Toggle one plot between normal and large size in its own local store."""
        triggered = callback_context.triggered_id
        if not isinstance(triggered, dict) or triggered.get("type") != "plots-size-toggle":
            return no_update
        if not click_ts or click_ts <= 0:
            return no_update
        size_value = shared.normalize_plot_size(current_size)
        return "large" if size_value != "large" else "normal"

    @app.callback(
        Output({"type": "plots-card", "index": MATCH}, "className"),
        Output({"type": "plots-graph-shell", "index": MATCH}, "style"),
        Output({"type": "plots-size-toggle", "index": MATCH}, "children"),
        Output({"type": "plots-size-toggle", "index": MATCH}, "className"),
        Input({"type": "plots-size-store", "index": MATCH}, "data"),
        prevent_initial_call=False,
    )
    def update_one_plot_size(size_value):
        """Mirror one plot's local size-store state into the card shell and toggle button."""
        current_size = shared.normalize_plot_size(size_value)
        button_text, button_class = shared.plot_size_button_props(current_size)
        return (
            shared.plot_card_class_name(current_size),
            shared.plot_graph_shell_style(current_size),
            button_text,
            button_class,
        )

    @app.callback(
        Output("plots-plot-order", "data", allow_duplicate=True),
        Output("plots-plot-state", "data", allow_duplicate=True),
        Output("plots-next-id", "data", allow_duplicate=True),
        Output("plots-last-add-ts", "data", allow_duplicate=True),
        Output("plots-plot-container", "children", allow_duplicate=True),
        Input("plots-add-budget", "n_clicks_timestamp"),
        Input("plots-add-profile", "n_clicks_timestamp"),
        Input("plots-add-timeseries", "n_clicks_timestamp"),
        Input("plots-add-timeheight", "n_clicks_timestamp"),
        Input("plots-add-subcolumn", "n_clicks_timestamp"),
        State("plots-case-data", "data"),
        State("plots-plot-order", "data"),
        State("plots-plot-state", "data"),
        State("plots-next-id", "data"),
        State("plots-last-add-ts", "data"),
        State("plots-plot-container", "children"),
        prevent_initial_call=True,
    )
    def add_plot(budget_ts, profile_ts, timeseries_ts, timeheight_ts, subcolumn_ts, case_data, plot_order, plot_state, next_id, last_add_ts, current_children):
        """Append one new plot card of the requested family and move the add-card to the end."""
        if not case_data:
            return no_update, no_update, no_update, no_update, no_update
        timestamps = {
            "plots-add-budget": budget_ts or 0,
            "plots-add-profile": profile_ts or 0,
            "plots-add-timeseries": timeseries_ts or 0,
            "plots-add-timeheight": timeheight_ts or 0,
            "plots-add-subcolumn": subcolumn_ts or 0,
        }
        latest_button = max(timestamps, key=timestamps.get)
        latest_ts = timestamps[latest_button]
        if latest_ts <= int(last_add_ts or 0):
            return no_update, no_update, no_update, no_update, no_update
        plot_type = next(
            spec.plot_type_id for spec in PLOT_FAMILY_SPECS.values() if spec.add_button_id == latest_button
        )
        module = PLOT_TYPES[plot_type]
        if not module.supports_case_data(case_data):
            return no_update, no_update, no_update, no_update, no_update
        plot_id = int(next_id or 0)
        state = default_plot_state(case_data, plot_id, plot_type=plot_type, existing_state=plot_state)
        updated_order = list(plot_order or [])
        updated_state = dict(plot_state or {})
        updated_order.append(plot_id)
        updated_state[str(plot_id)] = state
        if current_children:
            patched_children = Patch()
            patched_children[-1] = render_plot_card(plot_id, state, case_data)
            patched_children.append(add_plot_controls_card())
        else:
            patched_children = render_plot_grid(updated_order, updated_state, case_data)
        return updated_order, updated_state, plot_id + 1, latest_ts, patched_children

    @app.callback(
        Output("plots-plot-order", "data", allow_duplicate=True),
        Output("plots-plot-state", "data", allow_duplicate=True),
        Output("plots-plot-container", "children", allow_duplicate=True),
        Input({"type": "plots-close", "index": ALL}, "n_clicks_timestamp"),
        State({"type": "plots-close", "index": ALL}, "id"),
        State("plots-case-data", "data"),
        State("plots-plot-order", "data"),
        State("plots-plot-state", "data"),
        State("plots-plot-container", "children"),
        prevent_initial_call=True,
    )
    def remove_plot(timestamps, ids, case_data, plot_order, plot_state, current_children):
        """Remove one plot card and keep the remaining grid order intact."""
        if not timestamps or not ids:
            return no_update, no_update, no_update
        indexed = [(ts, meta) for ts, meta in zip(timestamps, ids) if ts]
        if not indexed:
            return no_update, no_update, no_update
        _, target = max(indexed, key=lambda pair: pair[0])
        plot_id = int(target.get("index"))
        updated_order = [int(idx) for idx in (plot_order or []) if int(idx) != plot_id]
        updated_state = dict(plot_state or {})
        updated_state.pop(str(plot_id), None)
        if current_children:
            remove_index = None
            for idx, child in enumerate(current_children):
                child_meta = child_id(child)
                if isinstance(child_meta, dict) and child_meta.get("type") == "plots-card" and int(child_meta.get("index")) == plot_id:
                    remove_index = idx
                    break
            if remove_index is not None:
                patched_children = Patch()
                del patched_children[remove_index]
            else:
                patched_children = no_update
        else:
            patched_children = render_plot_grid(updated_order, updated_state, case_data)
        return updated_order, updated_state, patched_children
