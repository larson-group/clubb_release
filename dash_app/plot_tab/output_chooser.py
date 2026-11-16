"""Self-contained Available/Active output-folder chooser for Dash Plots.

build_output_chooser() builds both sides, their controls and private stores.
register_output_chooser() owns discovery, refresh, extra paths, selection and
confirmed deletion. Consumers listen only to plots-output-dirs.data, the
ordered selected-folder list. Catalog refresh never prunes saved selections.

Browser rendering and add/remove transitions live in assets/34_plot_chooser.js;
styles live in assets/34_output_chooser.css. Discovery reuses the profile
service and runs in the existing catalog worker. The component does not call
plot rendering or invalidate scientific selections itself.
"""

from pathlib import Path
import shutil
import time

from dash import ALL, ClientsideFunction, Input, Output, State, callback_context, dcc, html, no_update

from utilities.output_paths import OUTPUT_ROOT
from dash_app.services.profiles import discover_output_directories, normalize_output_directories
from .async_callbacks import task_callback
from .plot_types.shared import normalize_output_directory
from .state import DEFAULT_OUTPUT_DIR


def build_output_chooser(selected_dirs=None, catalog=None):
    """Build one bordered, two-sided chooser, including all of its internal state."""
    selected = [DEFAULT_OUTPUT_DIR] if selected_dirs is None else normalize_output_directories(selected_dirs)
    return html.Section([
        dcc.Store(id='plots-output-dirs', data=selected),
        dcc.Store(id='plots-output-catalog', data=catalog or []),
        dcc.Store(id='plots-path-context', data={
            'repo': str(Path(__file__).resolve().parents[2]),
            'output': str(Path(OUTPUT_ROOT).resolve()), 'home': str(Path.home()),
        }),
        dcc.Store(id='plots-output-delete-confirm', data=None),
        dcc.Store(id='plots-output-expanded', data=False),
        dcc.Interval(id='plots-output-delete-expiry', interval=250, disabled=True, n_intervals=0),
        html.Div([
            html.Div([
                html.Button([
                    html.Span('Output folders'),
                    html.Span('▾', className='plots-output-chooser-chevron', **{'aria-hidden': 'true'}),
                ], id='plots-output-expand', n_clicks=0, n_clicks_timestamp=-1,
                   className='plots-output-chooser-toggle', title='Expand or collapse the output folders',
                   **{'aria-expanded': 'false', 'aria-controls': 'plots-output-chooser-sides'}),
                html.Button('Refresh outputs', id='plots-output-refresh', n_clicks=0),
                html.Button('Add extra folder', id='plots-show-extra-dir', n_clicks=0),
                html.Span(id='plots-catalog-status', role='status'),
            ], className='plots-output-chooser-heading'),
            html.Div([
                html.Div([
                    html.Div(['Available in ', html.Code('output/', className='plots-output-root-path')],
                             className='plots-output-side-heading'),
                    html.Div(id='plots-available-output-list', className='plots-available-output-list'),
                ], className='plots-output-available-side'),
                html.Div(className='plots-output-chooser-divider', role='separator',
                         **{'aria-orientation': 'vertical'}),
                html.Div([
                    html.Div('Active', className='plots-output-side-heading'),
                    html.Div(id='plots-active-output-list', className='plots-output-active-list'),
                ], className='plots-output-active-side'),
            ], id='plots-output-chooser-sides', className='plots-output-chooser-sides'),
            html.Div([
                dcc.Input(id='plots-extra-dir-input', type='text', placeholder='Path outside the discovered list'),
                html.Button('Add', id='plots-add-extra-dir', n_clicks=0, n_clicks_timestamp=-1),
            ], id='plots-extra-dir-control', style={'display': 'none'}),
            html.Small(id='plots-extra-dir-message', className='plots-mutable-output-warning'),
        ], id='plots-output-chooser-surface', className='plots-output-chooser-surface'),
    ], id='plots-output-chooser', className='plots-output-chooser', **{'aria-label': 'Output folders'})


def _output_display(path):
    """Label a directly inspected folder relative to output/, or keep its external path."""
    try:
        relative = Path(path).relative_to(Path(OUTPUT_ROOT).resolve()).as_posix()
        return {'label': 'output' if relative == '.' else 'output/' + relative, 'catalog_origin': 'output'}
    except ValueError:
        return {'label': path, 'catalog_origin': 'external'}


def _triggered_click_is_positive():
    """Ignore mount events when the chooser rebuilds its dynamic action buttons."""
    triggered = callback_context.triggered[0] if callback_context.triggered else None
    value = (triggered or {}).get('value')
    return isinstance(value, (int, float)) and value > 0


def _remove_output_dir(output_dirs, path):
    """Remove a deleted folder from selection while preserving every other path."""
    selected = normalize_output_directories(output_dirs)
    normalized = normalize_output_directory(str(path or ""))
    return [candidate for candidate in selected if candidate != normalized]


def _delete_output_directory(path):
    """Permanently remove one resolved subdirectory beneath output/."""
    target = Path(path).expanduser().resolve()
    output_root = Path(OUTPUT_ROOT).resolve()
    if target == output_root or output_root not in target.parents:
        raise ValueError("Only subdirectories inside output/ can be deleted here")
    shutil.rmtree(target)


def _catalog_tracking_paths(records, selected_dirs):
    """Keep active paths and previously added external paths in the catalog."""
    tracked = normalize_output_directories(selected_dirs)
    for record in records or []:
        if record.get("catalog_origin") != "external":
            continue
        path = str(record.get("path") or "")
        if path and path not in tracked:
            tracked.append(path)
    return tracked


def _sort_output_catalog(records):
    """Keep the catalog ordered by recency after a local insertion or deletion."""
    ordered = sorted(records, key=lambda item: (-float(item.get("modified") or 0), str(item.get("label") or item["path"])))
    found_newest = False
    result = []
    for index, record in enumerate(ordered):
        newest = bool(record.get("available")) and not found_newest
        result.append({**record, "recency_index": index, "is_newest": newest})
        found_newest = found_newest or newest
    return result


def register_output_chooser(app):
    """Register the chooser callbacks; plots consume its selected-folder store."""
    app.clientside_callback(ClientsideFunction("plotsChooser", "toggle"),
        Output('plots-output-expanded', 'data'),
        Output("plots-extra-dir-control", "style"),
        Input('plots-output-expand', 'n_clicks_timestamp'),
        Input("plots-show-extra-dir", "n_clicks"),
        State('plots-output-expanded', 'data'),
        State("plots-extra-dir-control", "style"),
        prevent_initial_call=True,
    )

    @task_callback(app, "catalog",
        Output("plots-output-catalog", "data"),
        Input("plots-output-chooser", "id"),
        Input("plots-output-refresh", "n_clicks"),
        Input("plots-output-dirs", "data"),
        State("plots-output-catalog", "data"),
    )
    def refresh_output_catalog(_mounted, _clicks, selected_dirs, current_catalog):
        """Discover in a worker; inspect only unknown folders on selection changes."""
        trigger = callback_context.triggered_id
        if trigger == "plots-output-dirs" and current_catalog:
            known = {record["path"] for record in (current_catalog or [])}
            missing = [path for path in normalize_output_directories(selected_dirs) if path not in known]
            if not missing:
                return no_update
            catalog = list(current_catalog or [])
            for path in missing:
                records = discover_output_directories(root=path, selected_dirs=[], recursive=False)
                direct = next((record for record in records if record["path"] == path), None)
                if direct is None:
                    direct = {"path": path, "label": path, "case_names": [], "case_count": 0,
                              "available": False, "error": "No readable case output found in this directory."}
                display = _output_display(path)
                catalog.append({**direct, "label": display["label"], "catalog_origin": display["catalog_origin"]})
            return _sort_output_catalog(catalog)
        tracked = _catalog_tracking_paths(current_catalog, selected_dirs)
        catalog = discover_output_directories(selected_dirs=tracked)
        return catalog if catalog != current_catalog else no_update

    app.clientside_callback(ClientsideFunction("plotsChooser", "render"),
        Output("plots-available-output-list", "children"),
        Output("plots-active-output-list", "children"),
        Output('plots-output-chooser-surface', 'className'),
        Output('plots-output-expand', 'aria-expanded'),
        Input("plots-output-catalog", "data"),
        Input("plots-output-dirs", "data"),
        Input("plots-output-delete-confirm", "data"),
        Input("plots-path-context", "data"),
        Input('plots-output-expanded', 'data'),
    )

    @app.callback(
        Output("plots-output-dirs", "data", allow_duplicate=True),
        Output("plots-output-catalog", "data", allow_duplicate=True),
        Output("plots-extra-dir-message", "children", allow_duplicate=True),
        Output("plots-output-delete-confirm", "data"),
        Output("plots-output-delete-expiry", "disabled"),
        Input({"type": "plots-delete-output-dir", "path": ALL}, "n_clicks_timestamp"),
        Input("plots-output-delete-expiry", "n_intervals"),
        State("plots-output-dirs", "data"),
        State("plots-output-catalog", "data"),
        State("plots-output-delete-confirm", "data"),
        prevent_initial_call=True,
    )
    def delete_output_directory(_clicks, _ticks, selected_dirs, catalog, confirmation):
        """Require two clicks within three seconds before deleting an output."""
        trigger = callback_context.triggered_id
        now = time.time()
        if trigger == "plots-output-delete-expiry":
            if confirmation and now >= float(confirmation.get("expires_at") or 0.0):
                return no_update, no_update, no_update, None, True
            return no_update, no_update, no_update, no_update, no_update
        if not isinstance(trigger, dict) or trigger.get("type") != "plots-delete-output-dir":
            return no_update, no_update, no_update, no_update, no_update
        if not _triggered_click_is_positive():
            return no_update, no_update, no_update, no_update, no_update

        target = str(Path(str(trigger.get("path") or "")).expanduser().resolve())
        armed = (
            confirmation
            and confirmation.get("path") == target
            and now < float(confirmation.get("expires_at") or 0.0)
        )
        if not armed:
            return no_update, no_update, "", {"path": target, "expires_at": now + 3.0}, False
        try:
            _delete_output_directory(target)
        except (OSError, ValueError) as exc:
            return no_update, no_update, str(exc), None, True

        selected = _remove_output_dir(selected_dirs, target)
        updated_catalog = _sort_output_catalog([record for record in (catalog or []) if record.get("path") != target])
        return selected, updated_catalog, f"Deleted {target}", None, True

    app.clientside_callback(ClientsideFunction("plotsChooser", "choose"),
        Output("plots-output-dirs", "data"),
        Output("plots-extra-dir-input", "value"),
        Output("plots-extra-dir-message", "children"),
        Input({"type": "plots-add-output-dir", "path": ALL}, "n_clicks_timestamp"),
        Input("plots-add-extra-dir", "n_clicks_timestamp"),
        Input({"type": "plots-remove-output-dir", "path": ALL}, "n_clicks_timestamp"),
        State("plots-output-dirs", "data"),
        State("plots-extra-dir-input", "value"),
        State("plots-output-catalog", "data"),
        State("plots-path-context", "data"),
        prevent_initial_call=True,
    )
