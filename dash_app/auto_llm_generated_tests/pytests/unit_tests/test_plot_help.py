from dash import no_update

from dash_app.plot_tab.callbacks_grid import render_plot_help_notecard
from dash_app.plot_tab.plot_types.help_content import PLOT_HELP
from dash_app.plot_tab.plot_types.registry import PLOT_TYPES
from dash_app.plot_tab.plot_types import shared


def _walk(component):
    yield component
    children = getattr(component, "children", None)
    if isinstance(children, (list, tuple)):
        for child in children:
            if hasattr(child, "_prop_names"):
                yield from _walk(child)
    elif hasattr(children, "_prop_names"):
        yield from _walk(children)


def test_every_plot_family_has_structured_help():
    assert set(PLOT_HELP) == set(PLOT_TYPES)
    for plot_type, spec in PLOT_HELP.items():
        assert spec["title"]
        assert spec["overview"]
        assert spec["sections"]
        dialog = PLOT_TYPES[plot_type].help_dialog(
            {"type": "plots-help-close", "index": plot_type}
        )
        text = " ".join(
            str(getattr(item, "children", "")) for item in _walk(dialog)
        )
        assert "Overview" in text
        assert spec["title"] in text


def test_plot_card_uses_full_size_question_button():
    card = shared.make_plot_card(
        controls="control",
        graph_id="graph",
        graph_shell_id="graph-shell",
        card_id="card",
        help_button_id={"type": "plots-help-open", "index": 4},
        close_button_id={"type": "plots-close", "index": 4},
        subtitle="Legacy concise tooltip",
    )
    help_button = next(
        item for item in _walk(card) if getattr(item, "className", "") == "plots-card-help"
    )

    assert help_button.children == "?"
    assert help_button.id == {"type": "plots-help-open", "index": 4}
    assert help_button.title == "Legacy concise tooltip"


def test_plot_help_open_close_and_mount_events():
    state = {"7": {"plot_type": "budget"}}

    assert render_plot_help_notecard(
        {"type": "plots-help-open", "index": 7}, 0, state
    ) is no_update
    dialog = render_plot_help_notecard(
        {"type": "plots-help-open", "index": 7}, 1, state
    )
    assert getattr(dialog, "className", "") == "shared-notecard-overlay"
    assert render_plot_help_notecard(
        {"type": "plots-help-close", "index": 7}, 1, state
    ) == ""
