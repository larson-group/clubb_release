"""Interactive parcel-path diagnostics for CLUBB mixing lengths."""

from ..registry import SubtabSpec, register_subtab

from .callbacks import register_callbacks
from .layout import build_layout


SUBTAB = register_subtab(
    SubtabSpec(
        slug="mixing-length-trajectories",
        title="Mixing Length Trajectories",
        summary="Trace the upward and downward parcel paths used to diagnose Lscale.",
        category="Mixing length",
        updated="Interactive NetCDF diagnostic",
        tags=("Lscale", "parcel trajectories", "mixing length"),
        order=30,
        page_value="misc-mixing-length-trajectories",
        build_layout=build_layout,
        register_callbacks=register_callbacks,
    )
)
