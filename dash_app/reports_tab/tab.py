"""Keep the document browser in sync with its search and files under doc/."""

from dash import Input, Output, State, dcc, no_update

from .catalog import catalog_token, discover_documents
from .layout import build_layout, document_tabs, document_view, selected_value


def build_tab(app, *, lazy=None):
    """Register document discovery, selection and browser-handoff callbacks."""
    app.clientside_callback(
        "function(tab) { return tab !== 'reports'; }",
        Output('reports-catalog-poll', 'disabled'), Input('dashboard-tabs', 'value'),
    )

    @app.callback(
        Output('reports-pages', 'children'),
        Output('reports-pages', 'value', allow_duplicate=True),
        Output('reports-catalog-token', 'data'),
        Input('reports-catalog-poll', 'n_intervals'),
        Input('dashboard-tabs', 'value'), Input('reports-search', 'value'),
        State('reports-catalog-token', 'data'), State('reports-pages', 'value'),
        prevent_initial_call=True,
    )
    def refresh_documents(_tick, tab, query, previous_token, current):
        """Rescan while visible, retaining the selection when it still matches."""
        if tab != 'reports':
            return (no_update,) * 3
        documents = discover_documents(query=query)
        token = catalog_token(documents)
        selection = selected_value(documents, current)
        if token == previous_token and selection == current:
            return (no_update,) * 3
        return document_tabs(documents), selection, token

    @app.callback(
        Output('reports-view', 'children'),
        Input('reports-pages', 'value'), Input('reports-catalog-token', 'data'),
    )
    def display_document(path, _token):
        """Reload a selected file when its version changes, including Markdown edits."""
        document = next((item for item in discover_documents() if item.path == path), None)
        return document_view(document)

    @app.callback(
        Output('reports-pages', 'value', allow_duplicate=True),
        Output('reports-search', 'value'),
        Input('dashboard-request', 'data'), prevent_initial_call=True,
    )
    def select_agent_document(request):
        """Select the exact doc-relative path requested by a browser handoff."""
        request = request or {}
        if request.get('tab') != 'reports' or request.get('operation') != 'open_document':
            return no_update, no_update
        path = request.get('path')
        if not any(item.path == path for item in discover_documents()):
            return no_update, no_update
        return path, ''

    if lazy is not None:
        return lazy.tab(label='Documents', value='reports', build=build_layout)
    return dcc.Tab(label='Documents', value='reports', children=build_layout())
