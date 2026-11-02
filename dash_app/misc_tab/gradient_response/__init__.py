"""Register the field-by-parameter gradient explorer in the Misc tab."""
from dash_app.misc_tab.registry import SubtabSpec, register_subtab
from .multi_callbacks import register_callbacks
from .layout import build_layout

SUBTAB = register_subtab(SubtabSpec(
    slug='gradient-response', title='Gradient response',
    summary='Compare field and parameter sensitivities and explore response curves.',
    category='JAX', tags=('gradient', 'loss'), order=50,
    build_layout=build_layout, register_callbacks=register_callbacks,
))
