import re

from dash import ALL, MATCH, Input, Output, State, callback_context, no_update

from .layout import render_plot_grid
from .plot_types.registry import PLOT_TYPES
from .plot_types.specs import PLOT_FAMILY_SPECS
from .plot_types import shared
from .state import default_plot_state
from dash_app.services import plots as plot_service
from dash_app.shared.activity import acknowledge_ui_request, set_plot_instances


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


def render_plot_help_notecard(trigger, trigger_value, plot_state):
    """Return the requested plot help dialog, close it, or ignore mount events."""
    if not isinstance(trigger_value, (int, float)) or trigger_value <= 0:
        return no_update
    if isinstance(trigger, dict) and trigger.get("type") == "plots-help-close":
        return ""
    if not (
        isinstance(trigger, dict) and trigger.get("type") == "plots-help-open"
    ):
        return no_update
    plot_id = str(trigger.get("index"))
    plot_type = ((plot_state or {}).get(plot_id) or {}).get("plot_type")
    module = PLOT_TYPES.get(plot_type)
    if module is None:
        return no_update
    return module.help_dialog({"type": "plots-help-close", "index": plot_id})


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
        Input("plots-plot-order", "data"),
        State("plots-plot-state", "data"),
    )
    def render_plot_container(case_data, plot_order, plot_state):
        """Render the plot grid from the current ordered plot state."""
        return render_plot_grid(plot_order or [], plot_state or {}, case_data)

    def next_available_plot_id(plot_order, plot_state, next_id):
        """Choose an id not already mounted, even after a stale UI callback."""
        return plot_service.next_plot_id(plot_order, plot_state, next_id)

    @app.callback(
        Output("plots-instance-snapshot", "data"),
        Input("plots-case-data", "data"),
        Input("plots-plot-order", "data"),
        Input("plots-plot-state", "data"),
        State("plots-next-id", "data"),
        State("dashboard-request", "data"),
    )
    def synchronize_plot_instances(case_data, plot_order, plot_state, next_id, dashboard_request):
        """Publish a safe current-card snapshot for typed Plot inspection."""
        snapshot = {
            "case": str((case_data or {}).get("name") or ""),
            "output_dirs": [str(path) for path in (case_data or {}).get("output_dirs") or []],
            "plots": plot_service.list_plot_instances(plot_order, plot_state),
            "next_id": int(next_id or 0),
        }
        set_plot_instances(snapshot)
        request_id = (dashboard_request or {}).get("id") if isinstance(dashboard_request, dict) else None
        if request_id and (dashboard_request or {}).get("operation") in {"set_view", "add_budget", "remove"}:
            acknowledge_ui_request(request_id)
        return snapshot

    @app.callback(
        Output("plots-help-modal", "children"),
        Input({"type": "plots-help-open", "index": ALL}, "n_clicks"),
        Input({"type": "plots-help-close", "index": ALL}, "n_clicks"),
        State("plots-plot-state", "data"),
        prevent_initial_call=True,
    )
    def update_plot_help_notecard(_open_clicks, _close_clicks, plot_state):
        """Open one plot-family help dialog or close the active dialog."""
        trigger = callback_context.triggered_id
        trigger_value = (
            callback_context.triggered[0].get("value")
            if callback_context.triggered
            else None
        )
        return render_plot_help_notecard(trigger, trigger_value, plot_state)
    for spec in PLOT_FAMILY_SPECS.values():
        if spec.dropdown_type is not None:
            register_plot_state_sync_callback(app, spec.plot_type_id, spec.dropdown_type)
            register_dropdown_search_callback(app, spec)

    @app.callback(
        Output("plots-plot-state", "data", allow_duplicate=True),
        Input({"type": "custom-expression", "index": ALL}, "value"),
        State({"type": "custom-expression", "index": ALL}, "id"),
        State("plots-plot-state", "data"),
        prevent_initial_call=True,
    )
    def sync_custom_expression_state(values, ids, current_state):
        """Persist custom expression text for custom plot cards."""
        updated = dict(current_state or {})
        changed = False
        for meta, value in zip(ids or [], values or []):
            idx = str(meta.get("index"))
            if idx not in updated:
                continue
            entry = dict(updated[idx])
            entry["plot_type"] = "custom"
            entry["expression"] = str(value or "")
            updated[idx] = entry
            changed = True
        return updated if changed else no_update

    @app.callback(
        Output("plots-plot-state", "data", allow_duplicate=True),
        Input({"type": "pdf-height", "index": ALL}, "value"),
        State({"type": "pdf-height", "index": ALL}, "id"),
        State("plots-plot-state", "data"),
        prevent_initial_call=True,
    )
    def sync_pdf_height_state(values, ids, current_state):
        """Persist the local height selected on each PDF-contour card."""
        updated = dict(current_state or {})
        changed = False
        for meta, value in zip(ids or [], values or []):
            idx = str(meta.get("index"))
            if idx not in updated or value is None:
                continue
            entry = dict(updated[idx])
            selected_height = float(value)
            if (
                entry.get("plot_type") == "pdf_contour"
                and entry.get("height") == selected_height
            ):
                continue
            entry["plot_type"] = "pdf_contour"
            entry["height"] = selected_height
            updated[idx] = entry
            changed = True
        return updated if changed else no_update

    @app.callback(
        Output("plots-plot-state", "data", allow_duplicate=True),
        Input({"type": "pdf-transport-view", "index": ALL}, "value"),
        State({"type": "pdf-transport-view", "index": ALL}, "id"),
        State("plots-plot-state", "data"),
        prevent_initial_call=True,
    )
    def sync_pdf_transport_view_state(values, ids, current_state):
        """Persist the transport coloring selected on each PDF-contour card."""
        updated = dict(current_state or {})
        changed = False
        for meta, value in zip(ids or [], values or []):
            idx = str(meta.get("index"))
            if idx not in updated or value not in {"upward", "downward", "signed"}:
                continue
            entry = dict(updated[idx])
            entry["plot_type"] = "pdf_contour"
            entry["transport_view"] = value
            updated[idx] = entry
            changed = True
        return updated if changed else no_update

    @app.callback(
        Output("plots-plot-state", "data", allow_duplicate=True),
        Input({"type": "pdf-color-signal", "index": ALL}, "value"),
        State({"type": "pdf-color-signal", "index": ALL}, "id"),
        State("plots-plot-state", "data"),
        prevent_initial_call=True,
    )
    def sync_pdf_color_signal_state(values, ids, current_state):
        """Persist the selected raw-SAM moment signal on each contour card."""
        updated = dict(current_state or {})
        changed = False
        allowed = {"probability", "wprcp", "wprtp", "wpchi", "wpthlp", "wprtp2"}
        for meta, value in zip(ids or [], values or []):
            idx = str(meta.get("index"))
            if idx not in updated or value not in allowed:
                continue
            entry = dict(updated[idx])
            if entry.get("color_signal") == value:
                continue
            entry["plot_type"] = "pdf_contour"
            entry["color_signal"] = value
            updated[idx] = entry
            changed = True
        return updated if changed else no_update

    @app.callback(
        Output("plots-plot-state", "data", allow_duplicate=True),
        Input({"type": "pdf-contour-smoothing", "index": ALL}, "value"),
        State({"type": "pdf-contour-smoothing", "index": ALL}, "id"),
        State("plots-plot-state", "data"),
        prevent_initial_call=True,
    )
    def sync_pdf_contour_smoothing_state(values, ids, current_state):
        """Persist the raw-SAM contour smoothing width for each PDF card."""
        updated = dict(current_state or {})
        changed = False
        for meta, value in zip(ids or [], values or []):
            idx = str(meta.get("index"))
            try:
                smoothing = min(max(float(value), 0.0), 3.0)
            except (TypeError, ValueError):
                continue
            if idx not in updated:
                continue
            entry = dict(updated[idx])
            if entry.get("contour_smoothing_bins") == smoothing:
                continue
            entry["plot_type"] = "pdf_contour"
            entry["contour_smoothing_bins"] = smoothing
            updated[idx] = entry
            changed = True
        return updated if changed else no_update

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
        Output("plots-plot-state", "data", allow_duplicate=True),
        Input({"type": "plots-size-store", "index": ALL}, "data"),
        State({"type": "plots-size-store", "index": ALL}, "id"),
        State("plots-plot-state", "data"),
        prevent_initial_call=True,
    )
    def persist_plot_sizes(size_values, size_ids, plot_state):
        """Persist local plot-card size stores back into the shared plot state."""
        if not size_ids:
            return no_update
        updated_state = dict(plot_state or {})
        changed = False
        for meta, size_value in zip(size_ids or [], size_values or []):
            if not isinstance(meta, dict):
                continue
            plot_id = str(meta.get("index"))
            if plot_id not in updated_state:
                continue
            normalized_size = shared.normalize_plot_size(size_value)
            entry = dict(updated_state[plot_id])
            if entry.get("size") == normalized_size:
                continue
            entry["size"] = normalized_size
            updated_state[plot_id] = entry
            changed = True
        return updated_state if changed else no_update

    @app.callback(
        Output({"type": "plots-card", "index": MATCH}, "className"),
        Output({"type": "plots-graph-shell", "index": MATCH}, "style"),
        Output({"type": "plots-size-toggle", "index": MATCH}, "children"),
        Output({"type": "plots-size-toggle", "index": MATCH}, "className"),
        Output({"type": "plots-size-toggle", "index": MATCH}, "disabled"),
        Input({"type": "plots-size-store", "index": MATCH}, "data"),
        Input("plots-enabled-benchmarks", "data"),
        State({"type": "plots-size-store", "index": MATCH}, "id"),
        State("plots-plot-state", "data"),
        prevent_initial_call=False,
    )
    def update_one_plot_size(size_value, enabled_benchmarks, size_store_id, plot_state):
        """Mirror one plot's local size-store state into the card shell and toggle button."""
        plot_id = str((size_store_id or {}).get("index"))
        plot_entry = (plot_state or {}).get(plot_id, {})
        lock_large = (
            plot_entry.get("plot_type") == "timeheight"
            and bool(enabled_benchmarks)
        )
        current_size = "large" if lock_large else shared.normalize_plot_size(size_value)
        button_text, button_class = shared.plot_size_button_props(current_size)
        return (
            shared.plot_card_class_name(current_size),
            shared.plot_graph_shell_style(current_size),
            button_text,
            button_class,
            lock_large,
        )

    @app.callback(
        Output("plots-plot-order", "data", allow_duplicate=True),
        Output("plots-plot-state", "data", allow_duplicate=True),
        Output("plots-next-id", "data", allow_duplicate=True),
        Output("plots-last-add-ts", "data", allow_duplicate=True),
        Input("plots-add-budget", "n_clicks_timestamp"),
        Input("plots-add-custom", "n_clicks_timestamp"),
        Input("plots-add-pdf-contour", "n_clicks_timestamp"),
        Input("plots-add-profile", "n_clicks_timestamp"),
        Input("plots-add-timeseries", "n_clicks_timestamp"),
        Input("plots-add-timeheight", "n_clicks_timestamp"),
        Input("plots-add-subcolumn", "n_clicks_timestamp"),
        State("plots-case-data", "data"),
        State("plots-plot-order", "data"),
        State("plots-plot-state", "data"),
        State("plots-next-id", "data"),
        State("plots-last-add-ts", "data"),
        prevent_initial_call=True,
    )
    def add_plot(budget_ts, custom_ts, pdf_contour_ts, profile_ts, timeseries_ts, timeheight_ts, subcolumn_ts, case_data, plot_order, plot_state, next_id, last_add_ts):
        """Append one new plot card of the requested family and move the add-card to the end."""
        if not case_data:
            return no_update, no_update, no_update, no_update
        timestamps = {
            "plots-add-budget": budget_ts or 0,
            "plots-add-custom": custom_ts or 0,
            "plots-add-pdf-contour": pdf_contour_ts or 0,
            "plots-add-profile": profile_ts or 0,
            "plots-add-timeseries": timeseries_ts or 0,
            "plots-add-timeheight": timeheight_ts or 0,
            "plots-add-subcolumn": subcolumn_ts or 0,
        }
        latest_button = max(timestamps, key=timestamps.get)
        latest_ts = timestamps[latest_button]
        if latest_ts <= int(last_add_ts or 0):
            return no_update, no_update, no_update, no_update
        plot_type = next(
            spec.plot_type_id for spec in PLOT_FAMILY_SPECS.values() if spec.add_button_id == latest_button
        )
        module = PLOT_TYPES[plot_type]
        if not module.supports_case_data(case_data):
            return no_update, no_update, no_update, no_update
        plot_id = next_available_plot_id(plot_order, plot_state, next_id)
        state = default_plot_state(case_data, plot_id, plot_type=plot_type, existing_state=plot_state)
        transition = plot_service.append_plot_instance(
            case_data,
            plot_order,
            plot_state,
            next_id,
            plot_type=plot_type,
            state=state,
        )
        return transition.plot_order, transition.plot_state, transition.next_id, latest_ts

    @app.callback(
        Output("plots-plot-order", "data", allow_duplicate=True),
        Output("plots-plot-state", "data", allow_duplicate=True),
        Input({"type": "plots-close", "index": ALL}, "n_clicks_timestamp"),
        State({"type": "plots-close", "index": ALL}, "id"),
        Input("dashboard-request", "data"),
        State("plots-case-data", "data"),
        State("plots-plot-order", "data"),
        State("plots-plot-state", "data"),
        State("plots-next-id", "data"),
        prevent_initial_call=True,
    )
    def remove_plot(timestamps, ids, dashboard_request, case_data, plot_order, plot_state, next_id):
        """Remove one native or typed-requested card through the shared service."""
        if (
            isinstance(dashboard_request, dict)
            and dashboard_request.get("tab") == "plots"
            and dashboard_request.get("operation") == "remove"
        ):
            transition = plot_service.remove_plot_instance(
                dashboard_request.get("plot_id"),
                plot_order,
                plot_state,
                next_id,
            )
            return transition.plot_order, transition.plot_state
        if not timestamps or not ids:
            return no_update, no_update
        indexed = [(ts, meta) for ts, meta in zip(timestamps, ids) if ts]
        if not indexed:
            return no_update, no_update
        _, target = max(indexed, key=lambda pair: pair[0])
        transition = plot_service.remove_plot_instance(
            target.get("index"),
            plot_order,
            plot_state,
            next_id,
        )
        return transition.plot_order, transition.plot_state
