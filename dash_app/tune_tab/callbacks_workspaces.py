"""Saved Tune workspace browser and revision-lineage callbacks."""

from __future__ import annotations

from datetime import datetime
import time

from dash import ALL, Input, Output, State, ctx, dcc, html, no_update

from .layout import action_button_style


def _is_positive_click_count(value):
    """Distinguish a real button press from Dash component hydration.

    Workspace-list rows are dynamic.  After a Start/Stop/Reset refresh, Dash
    can report a newly inserted row's ``n_clicks=0`` as the triggering input.
    Treating that as a press navigates away from the selected revision (usually
    to the first/newest workspace).  Only a positive count is a user action.
    """
    try:
        return float(value) > 0.0
    except (TypeError, ValueError):
        return False


def _format_bytes(value):
    value = max(0, int(value or 0))
    for unit in ("B", "KiB", "MiB", "GiB"):
        if value < 1024 or unit == "GiB":
            return f"{value:.0f} {unit}" if unit == "B" else f"{value:.1f} {unit}"
        value /= 1024
    return f"{value:.1f} GiB"


def _format_timestamp(value):
    if not value:
        return "unknown time"
    try:
        return datetime.fromisoformat(str(value).replace("Z", "+00:00")).strftime("%Y-%m-%d %H:%M UTC")
    except ValueError:
        return str(value)


def _workspace_rows(workspaces):
    rows = [
        html.Div(
            html.Button(
                "New workspace",
                id="tune-workspace-new",
                n_clicks=0,
                style=action_button_style("#2563eb"),
            ),
            className="tune-status-card",
            style={"display": "flex", "alignItems": "center", "gap": "10px", "marginBottom": "8px"},
        )
    ]
    if not workspaces:
        rows.append(
            html.Div(
                "No saved Tune workspaces yet. Choose New workspace to configure a run.",
                className="tune-empty-message",
                style={"padding": "12px 20px 20px"},
            )
        )
        return html.Div(rows)
    for item in workspaces:
        workspace_id = str(item["workspace_id"])
        latest_id = str(item.get("latest_revision_id") or "")
        disabled = not latest_id
        rows.append(
            html.Div(
                [
                    html.Button(
                        str(item.get("display_name") or workspace_id),
                        id={"type": "tune-workspace-load", "workspace": workspace_id, "revision": latest_id},
                        n_clicks=0,
                        disabled=disabled,
                        style=action_button_style("#2563eb", disabled=disabled),
                    ),
                    html.Span(
                        f"{item.get('revision_count', 0)} revision(s) · latest {item.get('latest_state', 'draft')} · "
                        f"created {_format_timestamp(item.get('created_at'))} · modified {_format_timestamp(item.get('modified_at'))} · "
                        f"{_format_bytes(item.get('size_bytes'))}",
                        className="tune-info-message",
                        style={"flex": "1 1 360px"},
                    ),
                    html.Div(
                        [
                            html.Button("Rename", id={"type": "tune-workspace-rename", "workspace": workspace_id}, n_clicks=0, style=action_button_style("#475569")),
                            html.Button("×", id={"type": "tune-workspace-delete", "workspace": workspace_id}, n_clicks=0, title="Delete saved workspace", style=action_button_style("#b91c1c")),
                        ],
                        style={"display": "flex", "gap": "6px"},
                    ),
                    dcc.Input(
                        id={"type": "tune-workspace-rename-value", "workspace": workspace_id},
                        value=str(item.get("display_name") or workspace_id),
                        debounce=True,
                        placeholder="Workspace name",
                        style={"width": "190px"},
                    ),
                ],
                className="tune-status-card",
                style={"display": "flex", "alignItems": "center", "gap": "10px", "flexWrap": "wrap", "marginBottom": "8px"},
            )
        )
    return html.Div(rows)


def _selection_for(workspaces, workspace_id, revision_id):
    """Build the UI selection from the durable execution state.

    A revision is editable only before its *first* start.  The selection is
    therefore derived from the persisted status rather than the transient Dash
    active-job store, which keeps the rule correct across a page reload.
    """
    workspace = next(
        (item for item in (workspaces or []) if str(item.get("workspace_id")) == str(workspace_id)),
        {},
    )
    revision = next(
        (item for item in (workspace.get("revisions") or []) if str(item.get("revision_id")) == str(revision_id)),
        {},
    )
    state = str(revision.get("state") or "draft")
    return {
        "workspace_id": str(workspace_id),
        "revision_id": str(revision_id),
        "mode": "draft" if state in {"draft", "idle"} else "readonly",
    }


def _workspace_button_class(selection, activity):
    """Return the cross-workspace activity class without rebuilding controls."""
    selection = dict(selection or {})
    workspace_id = str(selection.get("workspace_id") or "")
    revision_id = str(selection.get("revision_id") or "")
    other_active = any(
        str(entry.get("state") or "") in {"initializing", "running", "stopping"}
        and (
            str(entry.get("workspace_id") or "") != workspace_id
            or str(entry.get("revision_id") or "") != revision_id
        )
        for entry in (activity or [])
    )
    return "tune-workspace-change tune-workspace-change-active" if other_active else "tune-workspace-change"


def _sealed_selection_for_job(active_job, selection):
    """Seal only the draft that launched this job; never navigate the user back.

    ``tune-active-job`` can outlive a browser selection.  In particular, a
    user can open the workspace browser while a previous workspace continues
    in the broker.  That old job must not reassert itself into the persisted
    selection store on a later polling update.
    """
    job = dict(active_job or {})
    workspace_id = str(job.get("workspace_id") or "")
    revision_id = str(job.get("revision_id") or "")
    current = dict(selection or {})
    if not workspace_id or not revision_id:
        return None
    if (
        str(current.get("workspace_id") or "") == workspace_id
        and str(current.get("revision_id") or "") == revision_id
    ):
        if current.get("mode") == "readonly":
            return None
        return {"workspace_id": workspace_id, "revision_id": revision_id, "mode": "readonly"}
    # A newly created workspace has no ID until the Start action returns. A
    # short-lived token binds that one promotion to the exact draft the user
    # started, rather than any older active job still present in the browser.
    if (
        current.get("mode") == "new"
        and current.get("draft_token")
        and str(current.get("draft_token")) == str(job.get("draft_token") or "")
    ):
        return {"workspace_id": workspace_id, "revision_id": revision_id, "mode": "readonly"}
    return None


def register_workspace_callbacks(app):
    """Register browser-only controls; all persistent operations use the broker."""

    @app.callback(
        Output("tune-workspace-list", "data"),
        Output("tune-workspace-picker-list", "children"),
        Input("tune-workspace-refresh", "data"),
    )
    def refresh_workspace_list(_refresh):
        """Load the rich, potentially expensive saved-workspace browser on demand."""
        from dash_app.shared.broker_client import perform_action

        try:
            workspaces = list((perform_action("list_tuning_workspaces", {}, internal=True) or {}).get("workspaces") or [])
            return workspaces, _workspace_rows(workspaces)
        except Exception as exc:
            return [], html.Div(f"Could not load saved Tune workspaces: {exc}", className="tune-validation-message")

    @app.callback(
        Output("tune-workspace-activity", "data"),
        Input("tune-workspace-status-interval", "n_intervals"),
    )
    def refresh_workspace_activity(_status_tick):
        """Poll only compact state records for the cross-workspace indicator."""
        from dash_app.shared.broker_client import perform_action

        try:
            return list((perform_action("list_tuning_workspace_activity", {}, internal=True) or {}).get("activity") or [])
        except Exception:
            # The activity cue is advisory. Do not erase a usable workspace
            # browser or surface a periodic transport failure as a UI error.
            return no_update

    @app.callback(
        Output("tune-workspace-selection", "data"),
        Output("tune-workspace-delete-target", "data"),
        Output("tune-workspace-delete-confirm", "displayed"),
        Output("tune-workspace-refresh", "data", allow_duplicate=True),
        Output("tune-validation-message", "children", allow_duplicate=True),
        Input("tune-workspace-new", "n_clicks"),
        Input("tune-workspace-change", "n_clicks"),
        Input({"type": "tune-workspace-load", "workspace": ALL, "revision": ALL}, "n_clicks"),
        Input("tune-workspace-revision-selector", "value"),
        Input("tune-workspace-new-revision", "n_clicks"),
        Input("tune-workspace-name", "n_submit"),
        Input({"type": "tune-workspace-rename", "workspace": ALL}, "n_clicks"),
        Input({"type": "tune-workspace-delete", "workspace": ALL}, "n_clicks"),
        Input("tune-workspace-delete-confirm", "submit_n_clicks"),
        State("tune-workspace-selection", "data"),
        State("tune-workspace-delete-target", "data"),
        State("tune-workspace-list", "data"),
        State({"type": "tune-workspace-rename-value", "workspace": ALL}, "value"),
        State("tune-workspace-name", "value"),
        prevent_initial_call=True,
    )
    def manage_workspace(
        _new,
        _change,
        _load_clicks,
        revision_value,
        _new_revision,
        _name_submit,
        _rename_clicks,
        _delete_clicks,
        _delete_confirm,
        selection,
        delete_target,
        workspaces,
        rename_values,
        workspace_name,
    ):
        from dash_app.shared.broker_client import perform_action

        selection = dict(selection or {})
        trigger = ctx.triggered_id
        trigger_value = ctx.triggered[0].get("value") if ctx.triggered else None
        refresh = no_update
        if trigger == "tune-workspace-new":
            if not _is_positive_click_count(trigger_value):
                return no_update, no_update, False, no_update, no_update
            from tuner.workspaces import suggested_workspace_display_name

            return (
                {
                    "mode": "new",
                    "display_name": suggested_workspace_display_name(),
                    "draft_token": str(time.time_ns()),
                },
                None,
                False,
                refresh,
                "Configure this new workspace, then Start to save and run its original revision.",
            )
        if trigger == "tune-workspace-change":
            if not _is_positive_click_count(trigger_value):
                return no_update, no_update, False, no_update, no_update
            return {}, None, False, int(time.time() * 1000), ""
        if trigger == "tune-workspace-revision-selector":
            workspace_id = str(selection.get("workspace_id") or "")
            if workspace_id and revision_value:
                next_selection = _selection_for(workspaces, workspace_id, str(revision_value))
                if next_selection != selection:
                    return next_selection, None, False, refresh, ""
            return no_update, no_update, False, refresh, no_update
        if trigger == "tune-workspace-name":
            if not _is_positive_click_count(trigger_value):
                return no_update, no_update, False, no_update, no_update
            label = str(workspace_name or "").strip()
            if not label:
                return no_update, no_update, False, no_update, "Enter a non-empty unique workspace name."
            workspace_id = str(selection.get("workspace_id") or "")
            if not workspace_id:
                if selection.get("mode") == "new":
                    return {**selection, "display_name": label}, None, False, no_update, "New workspace name set."
                return no_update, no_update, False, no_update, "Choose or create a workspace before naming it."
            try:
                perform_action(
                    "rename_tuning_workspace",
                    {"workspace_id": workspace_id, "display_name": label},
                    internal=True,
                )
                return no_update, None, False, int(time.time() * 1000), "Workspace renamed."
            except Exception as exc:
                return no_update, None, False, no_update, f"Could not rename workspace: {exc}"
        if isinstance(trigger, dict):
            workspace_id = str(trigger.get("workspace") or "")
            action = str(trigger.get("type") or "")
            if action == "tune-workspace-load":
                if not _is_positive_click_count(trigger_value):
                    return no_update, no_update, False, no_update, no_update
                return _selection_for(workspaces, workspace_id, str(trigger.get("revision") or "")), None, False, refresh, ""
            if action == "tune-workspace-delete":
                if not _is_positive_click_count(trigger_value):
                    return no_update, no_update, False, no_update, no_update
                return no_update, workspace_id, True, no_update, no_update
            if action == "tune-workspace-rename":
                if not _is_positive_click_count(trigger_value):
                    return no_update, no_update, False, no_update, no_update
                labels = [str(value or "").strip() for value in (rename_values or [])]
                index = next((idx for idx, item in enumerate(workspaces or []) if str(item.get("workspace_id")) == workspace_id), -1)
                label = labels[index] if 0 <= index < len(labels) else ""
                try:
                    if not label:
                        raise ValueError("enter a non-empty workspace name")
                    perform_action("rename_tuning_workspace", {"workspace_id": workspace_id, "display_name": label}, internal=True)
                    return no_update, None, False, int(time.time() * 1000), "Workspace renamed."
                except Exception as exc:
                    return no_update, None, False, no_update, f"Could not rename workspace: {exc}"
        if trigger == "tune-workspace-delete-confirm" and delete_target and _is_positive_click_count(trigger_value):
            try:
                perform_action("delete_tuning_workspace", {"workspace_id": str(delete_target)}, internal=True)
                return {}, None, False, int(time.time() * 1000), "Workspace deleted."
            except Exception as exc:
                return no_update, None, False, no_update, f"Could not delete workspace: {exc}"
        workspace_id = str(selection.get("workspace_id") or "")
        revision_id = str(selection.get("revision_id") or "")
        if (
            trigger == "tune-workspace-new-revision"
            and _is_positive_click_count(trigger_value)
            and workspace_id
            and revision_id
        ):
            try:
                created = perform_action("create_tuning_revision", {"workspace_id": workspace_id, "revision_id": revision_id}, internal=True)
                return {"workspace_id": workspace_id, "revision_id": created["revision_id"], "mode": "draft"}, None, False, int(time.time() * 1000), "New editable revision created."
            except Exception as exc:
                return no_update, None, False, no_update, f"Could not create revision: {exc}"
        return no_update, no_update, False, no_update, no_update

    @app.callback(
        Output("tune-workspace-picker", "style"),
        Output("tune-results-view", "style"),
        Output("tune-right-pane", "style"),
        Output("tune-configure-content", "style"),
        Output("tune-top-controls", "style"),
        Output("tune-workspace-revision-selector", "options"),
        Output("tune-workspace-revision-selector", "placeholder"),
        Output("tune-workspace-new-revision", "disabled"),
        Output("tune-workspace-new-revision", "style"),
        Output("tune-reset-button", "disabled"),
        Output("tune-reset-button", "style"),
        Output("tune-workspace-name", "value"),
        Output("tune-workspace-summary", "children"),
        Output("tune-start-button", "children"),
        Input("tune-workspace-selection", "data"),
        Input("tune-workspace-list", "data"),
        Input("tune-status", "data"),
        Input("tune-active-job", "data"),
    )
    def render_workspace_state(selection, workspaces, status, active_job):
        selection = dict(selection or {})
        workspace_id = str(selection.get("workspace_id") or "")
        revision_id = str(selection.get("revision_id") or "")
        mode = str(selection.get("mode") or "")
        hidden = {"display": "none"}
        visible = {"display": "block"}
        # This callback owns the pane style after first render.  Keep the
        # structural flow properties here too; otherwise returning only the
        # enabled/disabled opacity would erase layout.py's no-nested-scroll
        # policy after a workspace state update.
        right_pane_flow = {
            "paddingLeft": "16px",
            "height": "auto",
            "minHeight": 0,
            "overflow": "visible",
        }
        locked = {**right_pane_flow, "opacity": 0.45, "pointerEvents": "none"}
        editable = {**right_pane_flow, "opacity": 1, "pointerEvents": "auto"}
        workspace = next((item for item in (workspaces or []) if str(item.get("workspace_id")) == workspace_id), None)
        options = [{"label": entry["revision_id"], "value": entry["revision_id"]} for entry in (workspace or {}).get("revisions", [])]
        if not workspace_id:
            creating = mode == "new"
            name = str(selection.get("display_name") or "")
            return (
                hidden if creating else visible, visible if creating else hidden,
                editable if creating else locked, editable if creating else locked,
                editable if creating else locked,
                [], "Revision",
                True, action_button_style("#7c3aed", disabled=True),
                True, action_button_style("#b45309", disabled=True),
                name,
                "Configure a new Tune workspace." if creating else "Choose a saved workspace or create a new one.", "Start",
            )
        is_draft = mode == "draft" or mode == "new"
        selected_revision = next((entry for entry in (workspace or {}).get("revisions", []) if entry.get("revision_id") == revision_id), {})
        listed_state = str(selected_revision.get("state") or "draft")
        status_state = str((status or {}).get("state") or "")
        status_job_dir = str((status or {}).get("job_dir") or "")
        revision_path = str(selected_revision.get("path") or "")
        # The browser list refreshes only when requested.  Prefer the live
        # selected revision status where it identifies that same directory,
        # otherwise a finished run would remain displayed as "running" and
        # keep Reset disabled until a page refresh.
        state = status_state if status_state and revision_path == status_job_dir else listed_state
        summary = f"{workspace.get('display_name', workspace_id) if workspace else workspace_id} · {revision_id or 'new'} · {state}"
        active = bool(active_job)
        can_reset = (
            bool(revision_id)
            and not is_draft
            and not active
            and state not in {"running", "initializing", "stopping"}
        )
        can_branch = can_reset
        start_label = "Continue" if mode == "readonly" and str((status or {}).get("state") or state) == "stopped" else "Start"
        # ``tune-top-controls`` includes Start, Stop, Reset, and the two
        # independently runnable scoreboard actions as well as editable
        # controls.  Locking its whole wrapper makes those actions inert even
        # when their own disabled state correctly says they are available.
        # The display callback disables each setting control individually, so
        # keep this wrapper interactive for every selected workspace.
        return (
            hidden, visible, editable, editable if is_draft else locked, editable,
            options, f"Current: {revision_id}" if revision_id else "Revision",
            not can_branch, action_button_style("#7c3aed", disabled=not can_branch),
            not can_reset, action_button_style("#b45309", disabled=not can_reset),
            str((workspace or {}).get("display_name") or ""),
            summary, start_label,
        )

    @app.callback(
        Output("tune-workspace-change", "className"),
        Input("tune-workspace-selection", "data"),
        Input("tune-workspace-activity", "data"),
    )
    def render_workspace_activity_indicator(selection, activity):
        return _workspace_button_class(selection, activity)

    @app.callback(
        Output("tune-reset-dialog", "style"),
        Input("tune-reset-prompt", "data"),
    )
    def render_reset_prompt(open_prompt):
        if not open_prompt:
            return {"display": "none"}
        return {
            "display": "block", "position": "fixed", "zIndex": 2500, "inset": 0,
            "backgroundColor": "rgba(0, 0, 0, 0.55)", "padding": "24px",
        }

    @app.callback(
        Output("tune-reset-prompt", "data"),
        Input("tune-reset-button", "n_clicks"),
        Input("tune-reset-cancel", "n_clicks"),
        Input("tune-reset-delete-data", "n_clicks"),
        Input("tune-reset-new-revision", "n_clicks"),
        State("tune-workspace-selection", "data"),
        prevent_initial_call=True,
    )
    def toggle_reset_prompt(_reset, _cancel, _delete_data, _new_revision, selection):
        trigger = ctx.triggered_id
        trigger_value = ctx.triggered[0].get("value") if ctx.triggered else None
        if trigger == "tune-reset-button" and _is_positive_click_count(trigger_value):
            current = dict(selection or {})
            if current.get("workspace_id") and current.get("revision_id") and current.get("mode") == "readonly":
                return True
        return False

    @app.callback(
        Output("tune-workspace-selection", "data", allow_duplicate=True),
        Output("tune-workspace-refresh", "data", allow_duplicate=True),
        Output("tune-validation-message", "children", allow_duplicate=True),
        Input("tune-reset-delete-data", "n_clicks"),
        Input("tune-reset-new-revision", "n_clicks"),
        State("tune-workspace-selection", "data"),
        State("tune-active-job", "data"),
        prevent_initial_call=True,
    )
    def apply_reset(_delete_data, _new_revision, selection, active_job):
        trigger = ctx.triggered_id
        if trigger not in {"tune-reset-delete-data", "tune-reset-new-revision"}:
            return no_update, no_update, no_update
        trigger_value = ctx.triggered[0].get("value") if ctx.triggered else None
        if not _is_positive_click_count(trigger_value):
            return no_update, no_update, no_update
        choice = "delete_data" if trigger == "tune-reset-delete-data" else "new_revision"
        current = dict(selection or {})
        workspace_id = str(current.get("workspace_id") or "")
        revision_id = str(current.get("revision_id") or "")
        if active_job:
            return no_update, no_update, "Stop the active Tune revision before resetting it."
        if not workspace_id or not revision_id or current.get("mode") != "readonly":
            return no_update, no_update, "Select a started, inactive revision before resetting it."
        from dash_app.shared.broker_client import perform_action

        try:
            if choice == "delete_data":
                perform_action(
                    "reset_tuning_revision",
                    {"workspace_id": workspace_id, "revision_id": revision_id},
                    internal=True,
                )
                return (
                    {"workspace_id": workspace_id, "revision_id": revision_id, "mode": "draft"},
                    int(time.time() * 1000),
                    "Revision data deleted. Its original settings are editable again.",
                )
            created = perform_action(
                "create_tuning_revision",
                {"workspace_id": workspace_id, "revision_id": revision_id},
                internal=True,
            )
            new_revision = str(created["revision_id"])
            return (
                {"workspace_id": workspace_id, "revision_id": new_revision, "mode": "draft"},
                int(time.time() * 1000),
                f"New editable revision {new_revision} created; the prior data is preserved.",
            )
        except Exception as exc:
            return no_update, no_update, f"Could not reset Tune revision: {exc}"

    @app.callback(
        Output("tune-workspace-selection", "data", allow_duplicate=True),
        Output("tune-workspace-refresh", "data", allow_duplicate=True),
        Input("tune-active-job", "data"),
        State("tune-workspace-selection", "data"),
        prevent_initial_call=True,
    )
    def seal_started_revision(active_job, selection):
        """Seal a draft immediately after launch, including a new workspace."""
        next_selection = _sealed_selection_for_job(active_job, selection)
        if next_selection is None:
            return no_update, no_update
        return (
            next_selection,
            int(time.time() * 1000),
        )

    @app.callback(
        Output("tune-workspace-selection", "data", allow_duplicate=True),
        Input("tune-workspace-list", "data"),
        State("tune-workspace-selection", "data"),
        prevent_initial_call=True,
    )
    def reconcile_persisted_workspace_selection(workspaces, selection):
        """Clear only a browser-restored selection that no longer exists.

        Do not normalize a valid selection's mode here.  The workspace list
        can refresh while a user is editing a draft; rewriting ``draft`` to
        ``readonly`` (or vice versa) re-triggers the large workspace loader
        and overwrites the visible form.  Lifecycle callbacks own those mode
        transitions explicitly.
        """
        current = dict(selection or {})
        workspace_id = str(current.get("workspace_id") or "")
        revision_id = str(current.get("revision_id") or "")
        if not workspace_id or not revision_id:
            return no_update
        workspace = next(
            (item for item in (workspaces or []) if str(item.get("workspace_id")) == workspace_id),
            None,
        )
        if workspace is None or not any(
            str(item.get("revision_id")) == revision_id for item in (workspace.get("revisions") or [])
        ):
            return {}
        return no_update
