"""Regression tests for preserving mounted Plot cards during grid edits."""

from dash import Dash, html, no_update

from dash_app.plot_tab import callbacks_grid


def _card(plot_id):
    return html.Div(id={"type": "plots-card", "index": plot_id})


def _current_children(*plot_ids):
    return [*[_card(plot_id) for plot_id in plot_ids], html.Div(id="plots-add-card")]


def test_adding_plot_inserts_only_the_new_card_before_add_controls(monkeypatch):
    monkeypatch.setattr(
        callbacks_grid,
        "render_plot_card",
        lambda plot_id, _state, _case_data: _card(plot_id),
    )

    result = callbacks_grid.update_plot_grid_children(
        [1, 2, 3],
        {"1": {}, "2": {}, "3": {}},
        {"name": "arm"},
        _current_children(1, 2),
    )

    operation = result.to_plotly_json()["operations"]
    assert len(operation) == 1
    assert operation[0]["operation"] == "Insert"
    assert operation[0]["params"]["index"] == 2
    assert operation[0]["params"]["value"].id == {
        "type": "plots-card",
        "index": 3,
    }


def test_removing_plot_deletes_only_the_target_card():
    result = callbacks_grid.update_plot_grid_children(
        [1, 3],
        {"1": {}, "3": {}},
        {"name": "arm"},
        _current_children(1, 2, 3),
    )

    assert result.to_plotly_json()["operations"] == [
        {
            "operation": "Delete",
            "location": [1],
            "params": {},
        }
    ]


def test_unchanged_plot_order_does_not_touch_the_grid():
    assert callbacks_grid.update_plot_grid_children(
        [1, 2],
        {"1": {}, "2": {}},
        {"name": "arm"},
        _current_children(1, 2),
    ) is no_update


def test_case_change_preserves_figures_and_does_not_upload_the_grid():
    app = Dash(__name__)
    callbacks_grid.register_grid_callbacks(app)
    entry = app.callback_map['plots-plot-container.children']
    assert not any(state['id'] == 'plots-plot-container' for state in entry['state'])
    render = entry['callback'].__wrapped__
    assert render({'name': 'rico', 'preserve_plot_view': False}, [1, 2],
                  {'1': {'plot_type': 'profile'}, '2': {'plot_type': 'profile'}},
                  [{'type': 'plots-card', 'index': 1}, {'type': 'plots-card', 'index': 2}]) is no_update


def test_explicit_set_view_can_still_replace_card_variables(monkeypatch):
    from types import SimpleNamespace
    app = Dash(__name__)
    callbacks_grid.register_grid_callbacks(app)
    render = app.callback_map['plots-plot-container.children']['callback'].__wrapped__
    monkeypatch.setattr(callbacks_grid, 'callback_context', SimpleNamespace(
        triggered=[{'prop_id': 'plots-case-data.data'}]))
    replacement = [html.Div('requested cards')]
    monkeypatch.setattr(callbacks_grid, 'render_plot_grid', lambda *args: replacement)
    assert render({'name': 'rico', 'replace_plot_cards': True}, [1], {'1': {}},
                  [{'type': 'plots-card', 'index': 1}]) is replacement
