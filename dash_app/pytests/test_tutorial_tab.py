from dash import Dash

from dash_app.tutorial_tab.layout import build_layout
from dash_app.tutorial_tab.tab import build_tab


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


def test_tutorial_layout_has_vertical_welcome_equations_and_adg1_pages():
    layout = build_layout()
    pages = next(component for component in _walk(layout) if getattr(component, "id", None) == "tutorial-pages")
    assert pages.vertical is True
    assert pages.value == "tutorial-welcome"
    assert [page.label for page in pages.children] == [
        "Welcome",
        "CLUBB Equations",
        "ADG1 two-Gaussian explorer",
    ]
    assert any(
        getattr(component, "id", None) == "equation-term-inspector"
        for component in _walk(layout)
    )
    assert any(getattr(component, "id", None) == "notes-adg1-gaussian-figure" for component in _walk(layout))


def test_tutorial_tab_registers_navigation_and_adg1_callbacks():
    app = Dash(__name__, suppress_callback_exceptions=True)
    tutorial = build_tab(app)
    app.layout = tutorial
    app._setup_server()
    assert tutorial.label == "Tutorial"
    assert len(app.callback_map) == 4
