"""Register the Affine Coherence Mixing Length page under Dash's Misc tab.

Importing this package registers SUBTAB with its layout builder and callback
installer. The page recomputes entraining parcel paths and survival-weighted
lengths from saved NetCDF columns; it does not run CLUBB. See analysis.py for
matrix conventions and doc/coherence_mixing_length.md for the method.
"""

from ..registry import SubtabSpec, register_subtab

from .callbacks import register_callbacks
from .layout import build_layout


SUBTAB = register_subtab(
    SubtabSpec(
        slug="affine-mixing-length-trajectories",
        title="Affine Coherence Mixing Length",
        summary=(
            "Explore inverse-square coherence and directional averaging "
            "over one shared launch-target barrier calculation."
        ),
        category="Mixing length",
        updated="Interactive affine NetCDF diagnostic",
        tags=(
            "Lscale",
            "affine parcels",
            "launch-target matrix",
            "survival weight",
            "inverse-square coherence",
        ),
        order=31,
        page_value="misc-affine-mixing-length-trajectories",
        build_layout=build_layout,
        register_callbacks=register_callbacks,
    )
)
