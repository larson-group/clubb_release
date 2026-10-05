"""Structural checks for the supported non-BUGS JAX radiation port."""

from __future__ import annotations
from utilities.output_paths import REPO_ROOT as _REPO_ROOT


ROOT = _REPO_ROOT


def test_compiled_radiation_modules_have_no_tracer_compatibility_path():
    for name in (
        "radiation_module.py",
        "simple_rad_module.py",
        "rad_lwsw_module.py",
        "cos_solar_zen_module.py",
        "soil_vegetation.py",
        "parameters_radiation.py",
    ):
        source = (ROOT / "clubb_jax/src/Radiation" / name).read_text(encoding="utf-8")
        assert "tracer_numpy" not in source
        assert "_is_tracer" not in source
