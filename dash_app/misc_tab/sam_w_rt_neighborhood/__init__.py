"""Self-contained raw-SAM w–r_t neighborhood explorer subtab."""

from ..registry import SubtabSpec, register_subtab

from .tab import build_layout, register_callbacks


SUBTAB = register_subtab(
    SubtabSpec(
        slug="sam-w-rt-neighborhood",
        title="SAM w–rₜ time-height grid",
        summary="An instant 5×5 browser over every raw ARM SAM w–rₜ plane.",
        category="SAM explorer",
        updated="Pre-rendered atlas",
        tags=("SAM", "w–rₜ", "cloud transport"),
        order=20,
        page_value="misc-report-sam-w-rt-neighborhood",
        build_layout=build_layout,
        register_callbacks=register_callbacks,
    )
)
