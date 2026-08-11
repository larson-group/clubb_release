"""Interactive viewer for CLUBB's dedicated hole-filling diagnostics."""

from ..registry import SubtabSpec, register_subtab

from .callbacks import register_callbacks
from .layout import build_layout


SUBTAB = register_subtab(
    SubtabSpec(
        slug="fill-holes-viewer",
        title="Fill Holes Viewer",
        summary="Inspect CLUBB profiles immediately before and after hole filling.",
        category="Numerical diagnostics",
        updated="Dedicated NetCDF diagnostics",
        tags=("hole filling", "NetCDF", "positive definite"),
        order=40,
        page_value="misc-fill-holes-viewer",
        build_layout=build_layout,
        register_callbacks=register_callbacks,
    )
)
