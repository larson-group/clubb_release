"""Tests for the Plot available-output menu and active-output tray."""

from dash_app.plot_tab.layout import active_output_items, available_output_buttons


def _record(path, label, cases):
    return {
        "path": path,
        "label": label,
        "relative_path": label,
        "case_names": cases,
        "case_count": len(cases),
        "available": True,
    }


def test_collapsed_menu_shows_three_newest_unselected_outputs():
    records = [
        _record("/tmp/newest", "output/newest", ["arm", "bomex"]),
        _record("/tmp/middle", "output/middle", ["rico"]),
        _record("/tmp/older", "output/older", ["atex"]),
        _record("/tmp/oldest", "output/oldest", ["dycoms2_rf01"]),
    ]

    buttons = available_output_buttons(records, [], expanded=False)

    assert len(buttons) == 3
    assert [button.children for button in buttons] == ["newest · 2", "middle · 1", "older · 1"]
    assert buttons[0].id["path"] == "/tmp/newest"
    assert buttons[0].title == "Cases: arm, bomex"


def test_selected_output_leaves_menu_and_promotes_next_newest():
    records = [
        _record("/tmp/newest", "output/newest", ["arm"]),
        _record("/tmp/older", "output/older", ["rico"]),
    ]

    buttons = available_output_buttons(records, ["/tmp/newest"], expanded=False)

    assert [button.id["path"] for button in buttons] == ["/tmp/older"]


def test_expanded_menu_shows_all_and_only_unselected_outputs():
    records = [
        _record("/tmp/newest", "output/newest", ["arm"]),
        _record("/tmp/middle", "output/middle", ["bomex"]),
        _record("/tmp/oldest", "output/oldest", ["rico"]),
    ]

    buttons = available_output_buttons(records, ["/tmp/middle"], expanded=True)

    assert [row.children[0].id["path"] for row in buttons] == ["/tmp/newest", "/tmp/oldest"]


def test_expanded_menu_arms_permanent_delete_with_explicit_confirmation(monkeypatch, tmp_path):
    from dash_app.plot_tab import layout

    monkeypatch.setattr(layout, "OUTPUT_ROOT", tmp_path)
    output = tmp_path / "run"
    record = _record(str(output), "output/run", ["arm"])

    row = available_output_buttons(
        [record],
        [],
        expanded=True,
        delete_confirmation={"path": str(output.resolve())},
    )[0]

    delete_button = row.children[1]
    assert delete_button.children == "Confirm"
    assert delete_button.id == {"type": "plots-delete-output-dir", "path": str(output.resolve())}
    assert "plots-output-dir-delete--confirm" in delete_button.className


def test_active_tray_shows_counts_tooltips_and_remove_controls():
    records = [_record("/tmp/run", "output/run", ["arm", "bomex"])]

    item = active_output_items(records, ["/tmp/run"])[0]

    assert item.children[0].children == "run · 2"
    assert item.title == "Cases: arm, bomex"
    assert item.children[1].id == {"type": "plots-remove-output-dir", "path": "/tmp/run"}
    assert item.children[1].n_clicks_timestamp == -1


def test_unavailable_selected_output_remains_in_active_tray():
    item = active_output_items([], ["/tmp/missing-output"])[0]

    assert "plots-output-active-item--unavailable" in item.className
    assert item.children[0].children == "missing-output · 0"
    assert item.children[1].children == "Unavailable"
