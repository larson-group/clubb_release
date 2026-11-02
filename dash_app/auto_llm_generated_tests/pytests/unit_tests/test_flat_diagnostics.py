"""Keep undefined diagnostic cells selectable and avoid fabricated tangent lines."""
from dash_app.misc_tab.gradient_response.callbacks import response_figure
from dash_app.misc_tab.gradient_response.multi_callbacks import matrix_figure


def test_undefined_maps_retain_clickable_parameter_and_field_coordinates():
    """An undefined cloud-fraction diagnostic remains a selectable named cell."""
    mapping = {'parameters': ['C8'], 'fields': ['wp2', 'cloud_frac'], 'parameter_ranges': [{'name': 'C8', 'min': .2, 'max': .8}], 'baseline': [.5], 'metrics': {'correlation': {'values': [.9, None], 'gradients': [[.2], [None]]}}}
    figure = matrix_figure(mapping, 'correlation')
    assert figure.data[0].z[1][0] is None
    assert tuple(figure.data[1].x) == ('C8',)
    assert tuple(figure.data[1].y) == ('cloud_frac',)
    assert tuple(figure.data[1].text) == ('undefined',)


def test_undefined_curve_has_explanation_and_no_tangent():
    """Missing values produce gaps and a useful explanation, with no bogus slope."""
    data = {'values': [.2, .5, .8], 'baseline': .5, 'parameter': 'C8', 'metrics': {'correlation': {'values': [None]*3, 'baseline': None, 'gradient': None}}}
    figure = response_figure('correlation', 'Correlation', data)
    assert len(figure.data) == 1
    assert tuple(figure.data[0].y) == (None, None, None)
    assert 'undefined' in figure.layout.title.text.lower()
    assert any('variation' in a.text for a in figure.layout.annotations)


def test_flat_model_std_ratio_zero_is_shown_without_a_gradient():
    """Zero std ratio is a real value; its missing derivative is labeled separately."""
    data = {'values': [.2, .5, .8], 'baseline': .5, 'parameter': 'C8', 'metrics': {'std_ratio': {'values': [0., .2, .4], 'baseline': 0., 'gradient': None}}}
    figure = response_figure('std_ratio', 'Std ratio', data)
    assert tuple(figure.data[0].y) == (0., .2, .4)
    assert len(figure.data) == 2
    assert 'gradient undefined' in figure.layout.title.text.lower()


def test_taylor_diagram_names_fields_with_undefined_coordinates():
    """The tuner explains a missing flat field instead of silently dropping it."""
    from dash_app.tune_tab.callbacks_display import build_taylor_figure
    results = [{'total_loss': 0., 'field_metrics': {'bomex': {'cloud_frac': {'loss': 0., 'rmse': 0., 'correlation': None, 'std_ratio': None, 'bias_norm': 0.}}}}]
    figure = build_taylor_figure(results)
    assert any('Undefined Taylor coordinates: bomex/cloud_frac' in a.text for a in figure.layout.annotations)


def test_average_diagnostic_is_undefined_if_one_included_window_is_undefined():
    """Do not silently discard a flat window when displaying average diagnostics."""
    from dash_app.tune_tab.callbacks_display import _window_display_metrics
    base = {'correlation': 1., 'subwindows': [{'loss': 0., 'correlation': None}, {'loss': 2., 'correlation': 1.}]}
    displayed, _ = _window_display_metrics(base, 'average')
    assert displayed['correlation'] is None
    assert displayed['loss'] == 1.
