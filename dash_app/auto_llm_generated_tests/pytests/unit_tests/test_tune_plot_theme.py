"""Check that Tune diagnostics retain their data when themed."""

import plotly.graph_objects as go

from dash_app.tune_tab.callbacks_display import _theme_diagnostic_figure


def test_dark_theme_preserves_data_and_makes_diverging_colors_visible():
    figure = go.Figure([
        go.Heatmap(z=[[-1, 1]], colorscale="RdBu"),
        go.Bar(x=["a", "b"], y=[-1, 1], marker=dict(color=[-1, 1], colorscale="RdBu")),
    ])
    figure.add_annotation(text="diagnostic", x=0, y=0)
    light = _theme_diagnostic_figure(go.Figure(figure), "light")
    dark = _theme_diagnostic_figure(go.Figure(figure), "dark")

    assert dark.layout.paper_bgcolor != light.layout.paper_bgcolor
    assert list(dark.data[0].z[0]) == [-1, 1]
    assert list(dark.data[1].y) == [-1, 1]
    assert dark.data[0].colorscale == dark.data[1].marker.colorscale
    assert dark.data[0].colorscale != light.data[0].colorscale
    assert dark.layout.annotations[0].font.color != dark.layout.paper_bgcolor
