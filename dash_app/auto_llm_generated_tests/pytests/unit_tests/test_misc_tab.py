"""Focused contracts for the consolidated Misc tab."""

from dash import Dash, dcc

from dash_app.misc_tab.layout import build_layout, subtab_page_value
from dash_app.misc_tab.registry import discover_subtabs
from dash_app.misc_tab.tab import build_tab


def _walk(component):
    yield component
    children = getattr(component, "children", None)
    if children is None:
        return
    if not isinstance(children, (list, tuple)):
        children = [children]
    for child in children:
        if hasattr(child, "to_plotly_json"):
            yield from _walk(child)


def test_misc_layout_keeps_self_contained_subtabs_in_one_vertical_rail():
    layout = build_layout()
    components = list(_walk(layout))
    directory = next(component for component in components if getattr(component, "id", None) == "misc-pages")

    assert isinstance(directory, dcc.Tabs)
    assert directory.vertical is True
    assert directory.parent_className == "misc-tabs-shell"
    assert directory.className == "misc-directory"
    subtabs = discover_subtabs()
    values = {subtab.slug: subtab_page_value(subtab) for subtab in subtabs}
    assert directory.value == values["sam-w-rt-neighborhood"]
    assert [child.value for child in directory.children] == [
        subtab_page_value(subtab) for subtab in subtabs
    ]
    assert values["mixing-length-trajectories"] in {
        child.value for child in directory.children
    }
    # Dash 3.4 accepts layout components only as tab children, not as labels.
    # Keep labels text-only so the complete app can build during a hot reload.
    assert all(isinstance(child.label, str) for child in directory.children)


def test_misc_tab_builds_as_the_single_top_level_destination():
    app = Dash(__name__, suppress_callback_exceptions=True)
    tab = build_tab(app)
    app.layout = tab
    app._setup_server()

    assert tab.label == "Misc"
    assert tab.value == "misc"


def test_misc_registry_owns_the_self_contained_subtabs():
    subtabs = discover_subtabs()

    assert subtabs
    assert all("misc_tab." in subtab.build_layout.__module__ for subtab in subtabs)
    assert all(".reports." not in subtab.build_layout.__module__ for subtab in subtabs)
