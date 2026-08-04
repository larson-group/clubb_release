from dash import Dash

from dash_app.tutorial_tab.clubb_equations_demo.app import build_layout, register_callbacks
from dash_app.tutorial_tab.clubb_equations_demo.components import DEFAULT_TERM, _term_min_width
from dash_app.tutorial_tab.clubb_equations_demo.content import EQUATION_GROUPS, OWNERSHIP_META, TERM_OCCURRENCES


def _walk(component):
    yield component
    children = getattr(component, "children", None)
    if children is None:
        return
    if not isinstance(children, (list, tuple)):
        children = [children]
    for child in children:
        if hasattr(child, "children"):
            yield from _walk(child)


def test_every_equation_term_has_unique_clickable_metadata():
    layout = build_layout()
    buttons = [
        component
        for component in _walk(layout)
        if isinstance(getattr(component, "id", None), dict)
        and component.id.get("type") == "equation-guide-term"
    ]
    occurrences = [button.id["index"] for button in buttons]
    assert len(occurrences) == len(TERM_OCCURRENCES)
    assert len(occurrences) == len(set(occurrences))
    assert set(occurrences) == set(TERM_OCCURRENCES)
    assert DEFAULT_TERM in TERM_OCCURRENCES
    fraction_buttons = [
        button for button in buttons if "eq-term-fraction" in button.className
    ]
    assert fraction_buttons
    assert all(
        any(
            marker in TERM_OCCURRENCES[button.id["index"]]["latex"]
            for marker in (r"\frac", r"\dfrac", r"\tfrac")
        )
        for button in fraction_buttons
    )
    widths = {
        button.id["index"]: float(button.style["--eq-term-min"].removesuffix("rem"))
        for button in buttons
    }
    assert min(widths.values()) >= 10.0
    assert max(widths.values()) <= 30.0
    assert _term_min_width(r"+\frac{\partial}{\partial z}[K\frac{\partial q}{\partial z}]") > _term_min_width(r"+q")


def test_equation_content_uses_official_ownership_palette_and_sources():
    assert [group["value"] for group in EQUATION_GROUPS] == ["budgets", "pdf", "cloud"]
    for group in EQUATION_GROUPS:
        for equation in group["equations"]:
            assert equation["source"]
            assert equation["lhs_owner"] in OWNERSHIP_META
            for item in equation["terms"]:
                assert item["ownership"] in OWNERSHIP_META
                assert item["meaning"]
                assert item["obtained"]
                assert item["why"]


def test_equation_guide_registers_one_pattern_callback():
    app = Dash(__name__, suppress_callback_exceptions=True)
    app.layout = build_layout()
    register_callbacks(app)
    app._setup_server()
    assert len(app.callback_map) == 1
