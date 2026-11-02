"""Build the field-by-parameter gradient explorer in the Misc tab."""
from dash import html
from dash_app.tune_tab.discovery import available_tuning_fields
from tuner.case_defaults import available_case_names
from tuner.parameter_ranges import tunable_parameter_defaults
from .multi_layout import build_multi_layout


def build_layout():
    """Return gradient maps and selected-cell response controls and plots."""
    return html.Div(build_multi_layout(sorted(tunable_parameter_defaults(), key=str.lower),
                                       available_tuning_fields(), available_case_names()),
                    className='gradient-response-page')
