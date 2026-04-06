"""Meta checks for call-tree test coverage."""
import ast
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
PY_API_FILE = REPO_ROOT / "clubb_python_api" / "clubb_python" / "clubb_api.py"
TESTS_DIR = REPO_ROOT / "clubb_python_api" / "tests"


ROUTINE_TO_API_NAMES = {
    "advance_clubb_core_module::advance_clubb_core": ("advance_clubb_core",),
    "advance_helper_module::calc_brunt_vaisala_freq_sqd": ("calc_brunt_vaisala_freq_sqd",),
    "advance_helper_module::calc_ri_zm": ("calc_ri_zm",),
    "advance_helper_module::calc_stability_correction": ("calc_stability_correction",),
    "advance_helper_module::compute_cx_fnc_richardson": ("compute_cx_fnc_richardson",),
    "advance_windm_edsclrm_module::advance_windm_edsclrm": ("advance_windm_edsclrm",),
    "advance_wp2_wp3_module::advance_wp2_wp3": ("advance_wp2_wp3",),
    "advance_xm_wpxp_module::advance_xm_wpxp": ("advance_xm_wpxp",),
    "advance_xp2_xpyp_module::advance_xp2_xpyp": ("advance_xp2_xpyp",),
    "advance_xp3_module::advance_xp3": ("advance_xp3",),
    "adg1_adg2_3d_luhar_pdf::adg1_w_closure": ("adg1_w_closure",),
    "adg1_adg2_3d_luhar_pdf::adg1_pdf_driver": ("adg1_pdf_driver",),
    "adg1_adg2_3d_luhar_pdf::adg2_pdf_driver": ("adg2_pdf_driver",),
    "adg1_adg2_3d_luhar_pdf::calc_luhar_params": ("calc_luhar_params",),
    "adg1_adg2_3d_luhar_pdf::close_luhar_pdf": ("close_luhar_pdf",),
    "adg1_adg2_3d_luhar_pdf::luhar_3d_pdf_driver": ("luhar_3d_pdf_driver",),
    "advance_helper_module::calc_xpwp": ("calc_xpwp",),
    "advance_helper_module::lscale_width_vert_avg": ("lscale_width_vert_avg",),
    "advance_helper_module::pvertinterp": ("pvertinterp",),
    "advance_helper_module::smooth_min": ("smooth_min",),
    "advance_helper_module::smooth_max": ("smooth_max",),
    "advance_helper_module::smooth_heaviside_peskin": ("smooth_heaviside_peskin",),
    "advance_helper_module::vertical_avg": ("vertical_avg",),
    "advance_helper_module::vertical_integral": ("vertical_integral",),
    "advance_helper_module::wp2_term_splat_lhs": ("wp2_term_splat_lhs",),
    "advance_helper_module::wp3_term_splat_lhs": ("wp3_term_splat_lhs",),
    "calc_pressure::calculate_thvm": ("calculate_thvm",),
    "calendar::compute_current_date_api": ("compute_current_date",),
    "calendar::julian2gregorian_date": ("julian2gregorian_date",),
    "clip_explicit::clip_covar": ("clip_covar",),
    "clip_explicit::clip_covars_denom": ("clip_covars_denom",),
    "clip_explicit::clip_rcm": ("clip_rcm",),
    "clip_explicit::clip_variance": ("clip_variance",),
    "diffusion::diffusion_zm_lhs": ("diffusion_zm_lhs",),
    "diffusion::diffusion_zt_lhs": ("diffusion_zt_lhs",),
    "fill_holes::fill_holes_vertical_api": ("fill_holes_vertical",),
    "grid_class::zm2zt2zm": ("zm2zt2zm",),
    "grid_class::zt2zm2zt": ("zt2zm2zt",),
    "grid_class::ddzm": ("ddzm",),
    "grid_class::ddzt": ("ddzt",),
    "interpolation::lin_interpolate_two_points": ("lin_interpolate_two_points",),
    "matrix_solver_wrapper::band_solve": ("band_solve",),
    "matrix_solver_wrapper::tridiag_solve": ("tridiag_solve",),
    "mean_adv::term_ma_zm_lhs": ("term_ma_zm_lhs",),
    "mean_adv::term_ma_zt_lhs": ("term_ma_zt_lhs",),
    "mixing_length::set_lscale_max": ("set_lscale_max",),
    "mono_flux_limiter::calc_turb_adv_range": ("calc_turb_adv_range",),
    "mono_flux_limiter::monotonic_turbulent_flux_limit": ("monotonic_turbulent_flux_limit",),
    "new_hybrid_pdf::calc_coef_wp2xp_implicit": ("calc_coef_wp2xp_implicit",),
    "new_hybrid_pdf::calc_coefs_wpxp2_semiimpl": ("calc_coefs_wpxp2_semiimpl",),
    "new_hybrid_pdf::calculate_coef_wp4_implicit": ("calculate_coef_wp4_implicit",),
    "new_hybrid_pdf::calculate_responder_params": ("calculate_responder_params",),
    "new_hybrid_pdf::calculate_w_params": ("calculate_w_params",),
    "new_hybrid_pdf_main::new_hybrid_pdf_driver": ("new_hybrid_pdf_driver",),
    "ly93_pdf::calc_params_ly93": ("calc_params_ly93",),
    "ly93_pdf::ly93_driver": ("ly93_driver",),
    "new_pdf::calc_coef_wp4_implicit": ("calc_coef_wp4_implicit",),
    "new_pdf::calc_coef_wpxp2_implicit": ("calc_coef_wpxp2_implicit",),
    "new_pdf::calc_coefs_wp2xp_semiimpl": ("calc_coefs_wp2xp_semiimpl",),
    "new_pdf::calc_coefs_wpxpyp_semiimpl": ("calc_coefs_wpxpyp_semiimpl",),
    "new_pdf::calc_limits_F_x_responder": ("calc_limits_f_x_responder",),
    "new_pdf::calc_responder_params": ("calc_responder_params",),
    "new_pdf::calc_setter_var_params": ("calc_setter_var_params",),
    "new_pdf_main::new_pdf_driver": ("new_pdf_driver",),
    "new_tsdadg_pdf::calc_l_x_skx_fnc": ("calc_l_x_skx_fnc",),
    "new_tsdadg_pdf::calc_setter_parameters": ("calc_setter_parameters_tsdadg",),
    "new_tsdadg_pdf::tsdadg_pdf_driver": ("tsdadg_pdf_driver",),
    "mixing_length::calc_lscale_directly": ("calc_lscale_directly",),
    "mixing_length::diagnose_lscale_from_tau": ("diagnose_lscale_from_tau",),
    "numerical_check::calculate_spurious_source": ("calculate_spurious_source",),
    "numerical_check::parameterization_check": ("parameterization_check",),
    "numerical_check::pdf_closure_check": ("pdf_closure_check",),
    "numerical_check::sfc_varnce_check": ("sfc_varnce_check",),
    "pdf_closure_module::calc_w_up_in_cloud": ("calc_w_up_in_cloud",),
    "pdf_closure_module::calc_wp2xp_pdf": ("calc_wp2xp_pdf",),
    "pdf_closure_module::calc_wp4_pdf": ("calc_wp4_pdf",),
    "pdf_closure_module::pdf_closure_driver": ("pdf_closure_driver",),
    "pdf_closure_module::calc_wpxp2_pdf": ("calc_wpxp2_pdf",),
    "pdf_closure_module::calc_wpxpyp_pdf": ("calc_wpxpyp_pdf",),
    "pdf_parameter_module::init_pdf_implicit_coefs_terms_api": ("init_pdf_implicit",),
    "pdf_parameter_module::zero_pdf_implicit_coefs_terms_api": ("zero_pdf_implicit_coefs_terms",),
    "pdf_utilities::calc_comp_corrs_binormal": ("calc_comp_corrs_binormal",),
    "pdf_utilities::smooth_corr_quotient": ("smooth_corr_quotient",),
    "skx_module::lg_2005_ansatz": ("lg_2005_ansatz",),
    "skx_module::skx_func": ("skx_func",),
    "skx_module::xp3_lg_2005_ansatz": ("xp3_lg_2005_ansatz",),
    "sponge_layer_damping::sponge_damp_xm": ("sponge_damp_xm",),
    "sponge_layer_damping::sponge_damp_xp2": ("sponge_damp_xp2",),
    "sponge_layer_damping::sponge_damp_xp3": ("sponge_damp_xp3",),
    "sfc_varnce_module::calc_sfc_varnce": ("calc_sfc_varnce",),
    "remapping_module::remap_vals_to_target": ("remap_vals_to_target",),
    "sigma_sqd_w_module::compute_sigma_sqd_w": ("compute_sigma_sqd_w",),
    "stats_clubb_utilities::stats_accumulate": ("stats_accumulate",),
    "precipitation_fraction::precip_fraction": ("precip_fraction",),
    "saturation::sat_mixrat_ice": ("sat_mixrat_ice",),
}


def _collect_python_api_defs():
    """Collect top-level exported function names in clubb_python/clubb_api.py."""
    tree = ast.parse(PY_API_FILE.read_text())
    defs = {node.name for node in tree.body if isinstance(node, ast.FunctionDef)}
    for node in tree.body:
        if isinstance(node, ast.ImportFrom):
            for alias in node.names:
                if alias.name == "*":
                    continue
                defs.add(alias.asname or alias.name)
    return defs


def _collect_test_api_usage():
    """Collect clubb_python.clubb_api function usages from test files."""
    usage = set()
    for test_file in sorted(TESTS_DIR.glob("test_*.py")):
        if test_file.name == Path(__file__).name:
            continue
        tree = ast.parse(test_file.read_text())

        alias_to_original = {}
        module_aliases = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and node.module == "clubb_python.clubb_api":
                for name in node.names:
                    alias_to_original[name.asname or name.name] = name.name
            elif isinstance(node, ast.ImportFrom) and node.module == "clubb_python":
                for name in node.names:
                    if name.name == "clubb_api":
                        module_aliases.add(name.asname or name.name)

        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            func = node.func
            if (
                isinstance(func, ast.Attribute)
                and isinstance(func.value, ast.Name)
                and func.value.id in module_aliases
            ):
                usage.add(func.attr)
            elif isinstance(func, ast.Name) and func.id in alias_to_original:
                usage.add(alias_to_original[func.id])
    return usage

def test_call_tree_api_bindings_exist():
    """Every covered call-tree routine should have declared Python API binding(s)."""
    api_defs = _collect_python_api_defs()
    required = {name for names in ROUTINE_TO_API_NAMES.values() for name in names}
    missing_defs = sorted(required - api_defs)
    assert not missing_defs, f"Missing API function definitions: {missing_defs}"


def test_call_tree_api_bindings_are_tested():
    """Every covered call-tree routine should be exercised by at least one test."""
    usage = _collect_test_api_usage()
    missing_test_coverage = []
    for routine, api_names in sorted(ROUTINE_TO_API_NAMES.items()):
        if not any(name in usage for name in api_names):
            missing_test_coverage.append((routine, api_names))

    assert not missing_test_coverage, (
        "Covered call-tree routine(s) missing test usage: "
        + ", ".join(f"{routine} -> {api_names}" for routine, api_names in missing_test_coverage)
    )
