"""Structural checks for the supported non-BUGS JAX radiation port."""

from __future__ import annotations

import ast
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def _functions(path: Path) -> list[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    return [node.name for node in tree.body if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))]


def _function(path: Path, name: str) -> ast.FunctionDef:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    return next(
        node for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name == name
    )


def _class_fields(path: Path, name: str) -> list[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    class_node = next(
        node for node in tree.body
        if isinstance(node, ast.ClassDef) and node.name == name
    )
    return [
        node.target.id for node in class_node.body
        if isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name)
    ]


def _return_names(function: ast.FunctionDef) -> list[str]:
    returns = [node for node in function.body if isinstance(node, ast.Return)]
    assert len(returns) == 1
    result = returns[0].value
    assert isinstance(result, ast.Tuple)
    return [value.id for value in result.elts if isinstance(value, ast.Name)]


def test_supported_routine_inventory_and_order_match_fortran_modules():
    assert _functions(ROOT / "clubb_jax/src/Radiation/radiation_module.py") == [
        "advance_clubb_radiation", "radiation_driver", "update_radiation_variables",
    ]
    assert _functions(ROOT / "clubb_jax/src/Radiation/simple_rad_module.py") == [
        "simple_rad", "simple_rad_bomex", "simple_rad_lba", "simple_rad_lba_init", "liq_water_path",
        "lba_radhtz",
    ]
    assert _functions(ROOT / "clubb_jax/src/Radiation/rad_lwsw_module.py") == ["sunray_sw"]
    assert _functions(ROOT / "clubb_jax/src/Radiation/cos_solar_zen_module.py") == ["cos_solar_zen"]
    assert _functions(ROOT / "clubb_jax/src/Radiation/soil_vegetation.py") == [
        "advance_soil_veg", "initialize_soil_veg",
    ]


def test_jax_only_helpers_and_boundaries_are_explicitly_documented():
    simple_rad = (ROOT / "clubb_jax/src/Radiation/simple_rad_module.py").read_text(encoding="utf-8")
    shortwave = (ROOT / "clubb_jax/src/Radiation/rad_lwsw_module.py").read_text(encoding="utf-8")
    radiation = (ROOT / "clubb_jax/src/Radiation/radiation_module.py").read_text(encoding="utf-8")
    parameters = (ROOT / "clubb_jax/src/Radiation/parameters_radiation.py").read_text(encoding="utf-8")

    assert "TODO(port-mirror): ``simple_rad_lba`` has this time-interpolation loop inline" in simple_rad
    assert "TODO(port-mirror): ``shortwave_column`` and its scan body" in shortwave
    assert "TODO(port-mirror): The JAX lifecycle calls this driver" in radiation
    assert "BUGSrad-only source fields" in parameters
    assert "including ``sol_const``" in parameters
    assert "excluded because that configuration is rejected" in parameters


def test_functional_returns_follow_fortran_inout_and_out_order():
    radiation = ROOT / "clubb_jax/src/Radiation/radiation_module.py"
    simple_rad = ROOT / "clubb_jax/src/Radiation/simple_rad_module.py"
    soil = ROOT / "clubb_jax/src/Radiation/soil_vegetation.py"

    assert _return_names(_function(simple_rad, "simple_rad")) == [
        "stats", "Frad_LW", "radht_LW",
    ]
    assert _return_names(_function(soil, "advance_soil_veg")) == [
        "stats", "deep_soil_T_in_K", "sfc_soil_T_in_K", "veg_T_in_K",
    ]
    assert _return_names(_function(radiation, "radiation_driver")) == [
        "stats", "err_info", "radht", "Frad", "Frad_SW_up", "Frad_LW_up",
        "Frad_SW_down", "Frad_LW_down", "radht_SW", "radht_LW", "Frad_SW", "Frad_LW",
    ]
    assert _return_names(_function(radiation, "advance_clubb_radiation")) == [
        "stats", "err_info", "deep_soil_T_in_K", "sfc_soil_T_in_K", "veg_T_in_K",
        "radht", "Frad", "Frad_SW_up", "Frad_LW_up", "Frad_SW_down", "Frad_LW_down",
        "radht_SW", "radht_LW", "Frad_SW", "Frad_LW",
    ]


def test_retained_parameter_order_matches_the_supported_source_surface():
    parameters = ROOT / "clubb_jax/src/Radiation/parameters_radiation.py"
    assert _class_fields(parameters, "RadiationParameters") == [
        "rad_scheme", "alvdr", "kappa", "F0", "F1", "eff_drop_radius", "gc", "omega",
        "Fs_values", "cos_solar_zen_times", "cos_solar_zen_values", "l_fix_cos_solar_zen",
        "l_sw_radiation", "l_rad_above_cloud", "nparam", "l_soil_veg", "lba_zrad", "lba_krad",
    ]
    source = parameters.read_text(encoding="utf-8")
    assert "sol_const:" not in source
    assert '"lba_radhtz"' not in (
        ROOT / "clubb_jax/src/Radiation/simple_rad_module.py"
    ).read_text(encoding="utf-8")


def test_active_driver_uses_only_the_canonical_radiation_module():
    driver = (ROOT / "clubb_jax/src/advance_clubb_to_end.py").read_text(encoding="utf-8")
    assert "from clubb_jax.src.Radiation.radiation_module import advance_clubb_radiation" in driver
    assert not (ROOT / "clubb_jax/src/Radiation/radiation.py").exists()


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


def test_conceptual_fortran_comments_are_preserved_verbatim():
    radiation_dir = ROOT / "clubb_jax/src/Radiation"
    expected = {
        "parameters_radiation.py": [
            "Albedo values (alvdr is used in the simplifed schemes as well)",
            "Use DYCOMS II RF02 heaviside step function",
        ],
        "cos_solar_zen_module.py": [
            "The angle  longang  is equivalent to the",
            "June 6, 2006",
        ],
        "soil_vegetation.py": [
            "We consider a semi half infinite medium, initially at the",
            "The equations given above are analogous to those used by",
            "Duynkerke, Peter G.  \"Radiation Fog: A Comparison of Model",
        ],
        "simple_rad_module.py": [
            "SPECIAL METHOD USED TO CALCULATE RADIATION",
            "The following diagram describes the differences in model grids",
            "Computation of radiative fluxes on staggered grid",
            "If you consider Frad to take place on mass levels, then computing",
            "compatible with the use of a stretched",
        ],
        "rad_lwsw_module.py": [
            "The code for sunray_sw was slightly reconstructed in order to more",
            "These variables come from Duynkerke eqn.20, which, with slight",
            "The negative sign for F_diff and F_dir is related to the definition",
            "Computation of shortwave fluxes on staggered grid",
        ],
        "radiation_module.py": [
            "With this call ordering, snow and ice water mixing ratio will be",
            "We update stats here each sample timestep - even if radiation is not advanced",
        ],
    }

    for name, comments in expected.items():
        source = (radiation_dir / name).read_text(encoding="utf-8")
        for comment in comments:
            assert comment in source, f"{name} is missing source comment: {comment}"
