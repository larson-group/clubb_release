"""Structural checks for the Fortran-to-JAX forcing mirrors."""

from __future__ import annotations

import ast
from pathlib import Path


BENCHMARK_CASES = Path(__file__).parents[1] / "src" / "Benchmark_cases"
INPUT_FIELDS = Path(__file__).parents[1] / "src" / "Input_fields"
CLUBB_CORE = Path(__file__).parents[1] / "src" / "CLUBB_core"
ADVANCE_CLUBB_TO_END = Path(__file__).parents[1] / "src" / "advance_clubb_to_end.py"


# Source order for the routines supported by the standalone JAX driver. The
# dycore-only forcing routine is intentionally excluded as documented at the
# top of time_dependent_input.py.
SUPPORTED_SOURCE_ROUTINE_ORDER = {
    "arm.py": ["arm_sfclyr"],
    "arm_0003.py": ["arm_0003_sfclyr"],
    "arm_3year.py": ["arm_3year_sfclyr"],
    "arm_97.py": ["arm_97_sfclyr"],
    "astex_a209.py": ["astex_a209_tndcy", "astex_a209_sfclyr"],
    "atex.py": ["calc_forcings", "atex_tndcy", "atex_sfclyr"],
    "atex_long.py": ["calc_forcings", "atex_long_tndcy", "atex_long_sfclyr"],
    "bomex.py": ["bomex_tndcy", "bomex_sfclyr"],
    "cloud_feedback.py": ["cloud_feedback_sfclyr"],
    "cobra.py": ["cobra_sfclyr"],
    "diag_ustar_module.py": ["diag_ustar"],
    "dycoms2_rf01.py": ["dycoms2_rf01_tndcy", "dycoms2_rf01_sfclyr"],
    "dycoms2_rf02.py": ["dycoms2_rf02_tndcy", "dycoms2_rf02_sfclyr"],
    "ekman.py": ["ekman_sfclyr"],
    "fire.py": ["fire_sfclyr"],
    "gabls2.py": ["gabls2_tndcy", "gabls2_sfclyr"],
    "gabls3.py": ["gabls3_sfclyr"],
    "gabls3_night.py": ["gabls3_night_sfclyr", "psi_h", "gm1", "gh1", "fm1", "fh1", "landflx"],
    "jun25.py": ["jun25_altocu_read_t_dependent"],
    "lba.py": ["lba_tndcy", "lba_sfclyr"],
    "mpace_a.py": ["mpace_a_tndcy", "mpace_a_sfclyr", "mpace_a_init"],
    "mpace_b.py": ["mpace_b_tndcy", "mpace_b_sfclyr"],
    "neutral_case.py": ["neutral_case_sfclyr"],
    "nov11.py": ["nov11_altocu_rtm_adjust", "nov11_altocu_read_t_dependent"],
    "prescribe_forcings.py": ["prescribe_forcings", "read_surface_var_for_bc"],
    "rico.py": ["rico_tndcy", "rico_sfclyr"],
    "sfc_flux.py": [
        "compute_momentum_flux",
        "compute_ubar",
        "compute_ht_mostr_flux",
        "compute_wpthlp_sfc",
        "compute_wprtp_sfc",
        "set_sclr_sfc_rtm_thlm",
        "convert_sens_ht_to_km_s",
        "convert_latent_ht_to_m_s",
    ],
    "spec_hum_to_mixing_ratio.py": [
        "flux_spec_hum_to_mixing_ratio",
        "force_spec_hum_to_mixing_ratio",
    ],
    "time_dependent_input.py": [
        "initialize_t_dependent_input",
        "finalize_t_dependent_input",
        "initialize_t_dependent_sfc",
        "initialize_t_dependent_forcings",
        "finalize_t_dependent_forcings",
        "finalize_t_dependent_sfc",
        "read_to_grid",
        "apply_time_dependent_forcings_from_array",
        "apply_time_dependent_forcings",
        "time_select",
    ],
    "twp_ice.py": ["twp_ice_sfclyr"],
    "wangara.py": ["wangara_tndcy", "wangara_sfclyr"],
}

PRESCRIBE_FORCINGS_ARGUMENTS = [
    "gr",
    "nzm",
    "nzt",
    "ngrdcol",
    "sclr_dim",
    "edsclr_dim",
    "sclr_idx",
    "runtype",
    "sfctype",
    "time_current",
    "time_initial",
    "dt",
    "um",
    "vm",
    "thlm",
    "p_in_Pa",
    "exner",
    "rho",
    "rho_zm",
    "thvm",
    "veg_T_in_K",
    "l_modify_bc_for_cnvg_test",
    "saturation_formula",
    "stats",
    "rtm",
    "wm_zm",
    "wm_zt",
    "ug",
    "vg",
    "um_ref",
    "vm_ref",
    "thlm_forcing",
    "rtm_forcing",
    "um_forcing",
    "vm_forcing",
    "wprtp_forcing",
    "wpthlp_forcing",
    "rtp2_forcing",
    "thlp2_forcing",
    "rtpthlp_forcing",
    "wpsclrp",
    "sclrm_forcing",
    "edsclrm_forcing",
    "wpthlp_sfc",
    "wprtp_sfc",
    "upwp_sfc",
    "vpwp_sfc",
    "T_sfc",
    "p_sfc",
    "sens_ht",
    "latent_ht",
    "wpsclrp_sfc",
    "wpedsclrp_sfc",
    "err_info",
]

PRESCRIBE_FORCINGS_RESULTS = [
    "stats",
    "rtm",
    "wm_zm",
    "wm_zt",
    "ug",
    "vg",
    "um_ref",
    "vm_ref",
    "thlm_forcing",
    "rtm_forcing",
    "um_forcing",
    "vm_forcing",
    "wprtp_forcing",
    "wpthlp_forcing",
    "rtp2_forcing",
    "thlp2_forcing",
    "rtpthlp_forcing",
    "wpsclrp",
    "sclrm_forcing",
    "edsclrm_forcing",
    "wpthlp_sfc",
    "wprtp_sfc",
    "upwp_sfc",
    "vpwp_sfc",
    "T_sfc",
    "p_sfc",
    "sens_ht",
    "latent_ht",
    "wpsclrp_sfc",
    "wpedsclrp_sfc",
    "err_info",
]


HOST_BOUNDARY_FILES = {
    "mpace_a.py",  # mpace_a_init reads the source data files.
    "time_dependent_input.py",  # Source-equivalent input initialization and parsing.
}
DEVICE_PHYSICS_FILES = tuple(
    filename
    for filename in SUPPORTED_SOURCE_ROUTINE_ORDER
    if filename not in HOST_BOUNDARY_FILES
)


def _tree(filename: str) -> ast.Module:
    return ast.parse((BENCHMARK_CASES / filename).read_text())


def _routine_names(filename: str) -> list[str]:
    return [
        node.name
        for node in _tree(filename).body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
    ]


def test_complete_forcing_modules_keep_source_routine_inventory_and_order():
    for filename, expected in SUPPORTED_SOURCE_ROUTINE_ORDER.items():
        assert _routine_names(filename) == expected, filename


def test_prescribe_forcings_keeps_supported_source_argument_order():
    routine = next(
        node
        for node in _tree("prescribe_forcings.py").body
        if isinstance(node, ast.FunctionDef) and node.name == "prescribe_forcings"
    )
    assert [argument.arg for argument in routine.args.args] == PRESCRIBE_FORCINGS_ARGUMENTS
    for return_node in (
        node for node in ast.walk(routine) if isinstance(node, ast.Return)
    ):
        assert isinstance(return_node.value, ast.Tuple)
        assert [value.id for value in return_node.value.elts] == PRESCRIBE_FORCINGS_RESULTS


def test_forcing_driver_keeps_explicit_call_and_result_order():
    tree = ast.parse(ADVANCE_CLUBB_TO_END.read_text())
    driver = next(
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name == "advance_clubb_to_end"
    )
    wrapper = next(
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name == "_prescribe_forcings"
    )
    driver_calls = [
        node
        for node in ast.walk(driver)
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Name)
        and node.func.id == "_prescribe_forcings"
    ]
    assert len(driver_calls) == 1
    assert [argument.id for argument in driver_calls[0].args] == ["state", "time_current"]

    assignment = next(
        node
        for node in ast.walk(wrapper)
        if isinstance(node, ast.Assign)
        and isinstance(node.value, ast.Call)
        and isinstance(node.value.func, ast.Name)
        and node.value.func.id == "prescribe_forcings"
    )

    def state_key(node):
        if (
            isinstance(node, ast.Subscript)
            and isinstance(node.value, ast.Name)
            and node.value.id == "state"
            and isinstance(node.slice, ast.Constant)
        ):
            return node.slice.value
        if (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr == "asarray"
            and len(node.args) == 1
            and isinstance(node.args[0], ast.Name)
        ):
            return node.args[0].id
        if isinstance(node, ast.Attribute) and node.attr == "saturation_formula":
            return "saturation_formula"
        raise AssertionError(ast.dump(node))

    expected_arguments = [
        {
            "dt": "dt_main",
            "stats": "_jax_stats",
        }.get(argument, argument)
        for argument in PRESCRIBE_FORCINGS_ARGUMENTS
    ]
    expected_results = [
        "_jax_stats" if result == "stats" else result
        for result in PRESCRIBE_FORCINGS_RESULTS
    ]

    assert [state_key(argument) for argument in assignment.value.args] == expected_arguments
    assert len(assignment.targets) == 1
    assert isinstance(assignment.targets[0], ast.Tuple)
    assert [state_key(target) for target in assignment.targets[0].elts] == expected_results


def test_input_reader_keeps_source_routine_inventory_and_order():
    tree = ast.parse((INPUT_FIELDS / "input_reader.py").read_text())
    routines = [node.name for node in tree.body if isinstance(node, ast.FunctionDef)]
    assert routines == [
        "read_two_dim_file",
        "read_one_dim_file",
        "fill_blanks_one_dim_vars",
        "fill_blanks_two_dim_vars",
        "linear_fill_blanks",
        "deallocate_one_dim_vars",
        "deallocate_two_dim_vars",
        "read_x_table",
        "read_x_profile",
        "get_target_index",
        "count_columns",
    ]


def test_file_functions_keeps_source_routine_inventory_and_order():
    tree = ast.parse((CLUBB_CORE / "file_functions.py").read_text())
    routines = [node.name for node in tree.body if isinstance(node, ast.FunctionDef)]
    assert routines == ["file_read_1d", "file_read_2d"]


def test_device_physics_does_not_import_numpy_or_raw_f2py():
    for filename in DEVICE_PHYSICS_FILES:
        source = (BENCHMARK_CASES / filename).read_text()
        tree = ast.parse(source)
        imported_modules = {
            alias.name
            for node in tree.body
            if isinstance(node, ast.Import)
            for alias in node.names
        }
        imported_from = {
            node.module or ""
            for node in tree.body
            if isinstance(node, ast.ImportFrom)
        }
        assert "numpy" not in imported_modules, filename
        assert not any("clubb_f2py" in name for name in imported_modules | imported_from), filename


def test_removed_duplicate_forcing_helpers_do_not_return():
    forbidden = {
        "_diag_ustar",
        "_landflx_scalar",
        "prescribe_forcings_arm",
        "prescribe_forcings_generic",
        "_stats_surface_update",
        "_zero_forcings",
        "_time_interp",
        "_read_mpace_dat",
        "_mpace_time_select",
        "lba_diurnal_factor",
        "initialize_surface_bc_metadata",
        "_read_surface_var_for_bc_interp",
        "load_arm_forcings_data",
        "load_generic_forcings_data",
        "_parse_forcings_file",
        "_parse_sfc_file",
        "apply_time_dependent_forcings_from_dycore",
    }
    present = set()
    for path in BENCHMARK_CASES.glob("*.py"):
        present.update(_routine_names(path.name))
    assert forbidden.isdisjoint(present)
