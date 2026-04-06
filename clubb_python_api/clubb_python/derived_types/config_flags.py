"""Python representation of model_flags.F90 config flag UDT."""

from typing import NamedTuple

_LOGICAL_FLAG_NAMES = [
    'l_use_precip_frac', 'l_predict_upwp_vpwp',
    'l_ho_nontrad_coriolis', 'l_ho_trad_coriolis',
    'l_min_wp2_from_corr_wx', 'l_min_xp2_from_corr_wx',
    'l_C2_cloud_frac', 'l_diffuse_rtm_and_thlm',
    'l_stability_correct_Kh_N2_zm', 'l_calc_thlp2_rad',
    'l_upwind_xpyp_ta', 'l_upwind_xm_ma',
    'l_uv_nudge', 'l_rtm_nudge', 'l_tke_aniso',
    'l_vert_avg_closure', 'l_trapezoidal_rule_zt',
    'l_trapezoidal_rule_zm', 'l_call_pdf_closure_twice',
    'l_standard_term_ta', 'l_partial_upwind_wp3',
    'l_godunov_upwind_wpxp_ta', 'l_godunov_upwind_xpyp_ta',
    'l_use_cloud_cover', 'l_diagnose_correlations',
    'l_calc_w_corr', 'l_const_Nc_in_cloud',
    'l_fix_w_chi_eta_correlations', 'l_stability_correct_tau_zm',
    'l_damp_wp2_using_em', 'l_do_expldiff_rtm_thlm',
    'l_Lscale_plume_centered', 'l_diag_Lscale_from_tau',
    'l_use_C7_Richardson', 'l_use_C11_Richardson',
    'l_use_shear_Richardson', 'l_brunt_vaisala_freq_moist',
    'l_use_thvm_in_bv_freq', 'l_rcm_supersat_adj',
    'l_damp_wp3_Skw_squared', 'l_prescribed_avg_deltaz',
    'l_lmm_stepping', 'l_e3sm_config', 'l_vary_convect_depth',
    'l_use_tke_in_wp3_pr_turb_term', 'l_use_tke_in_wp2_wp3_K_dfsn',
    'l_use_wp3_lim_with_smth_Heaviside',
    'l_smooth_Heaviside_tau_wpxp',
    'l_modify_limiters_for_cnvg_test',
    'l_enable_relaxed_clipping', 'l_linearize_pbl_winds',
    'l_mono_flux_lim_thlm', 'l_mono_flux_lim_rtm',
    'l_mono_flux_lim_um', 'l_mono_flux_lim_vm',
    'l_mono_flux_lim_spikefix', 'l_host_applies_sfc_fluxes',
    'l_wp2_fill_holes_tke', 'l_add_dycore_grid',
]

_INT_FLAG_NAMES = [
    'iiPDF_type', 'ipdf_call_placement',
    'penta_solve_method', 'tridiag_solve_method',
    'saturation_formula', 'grid_remap_method',
    'grid_adapt_in_time_method', 'fill_holes_type',
]

CONFIG_FLAG_NAMES = _INT_FLAG_NAMES + _LOGICAL_FLAG_NAMES


class ConfigFlags(NamedTuple):
    """CLUBB configuration flags in Fortran argument order."""

    iiPDF_type: int
    ipdf_call_placement: int
    penta_solve_method: int
    tridiag_solve_method: int
    saturation_formula: int
    grid_remap_method: int
    grid_adapt_in_time_method: int
    fill_holes_type: int
    l_use_precip_frac: bool
    l_predict_upwp_vpwp: bool
    l_ho_nontrad_coriolis: bool
    l_ho_trad_coriolis: bool
    l_min_wp2_from_corr_wx: bool
    l_min_xp2_from_corr_wx: bool
    l_C2_cloud_frac: bool
    l_diffuse_rtm_and_thlm: bool
    l_stability_correct_Kh_N2_zm: bool
    l_calc_thlp2_rad: bool
    l_upwind_xpyp_ta: bool
    l_upwind_xm_ma: bool
    l_uv_nudge: bool
    l_rtm_nudge: bool
    l_tke_aniso: bool
    l_vert_avg_closure: bool
    l_trapezoidal_rule_zt: bool
    l_trapezoidal_rule_zm: bool
    l_call_pdf_closure_twice: bool
    l_standard_term_ta: bool
    l_partial_upwind_wp3: bool
    l_godunov_upwind_wpxp_ta: bool
    l_godunov_upwind_xpyp_ta: bool
    l_use_cloud_cover: bool
    l_diagnose_correlations: bool
    l_calc_w_corr: bool
    l_const_Nc_in_cloud: bool
    l_fix_w_chi_eta_correlations: bool
    l_stability_correct_tau_zm: bool
    l_damp_wp2_using_em: bool
    l_do_expldiff_rtm_thlm: bool
    l_Lscale_plume_centered: bool
    l_diag_Lscale_from_tau: bool
    l_use_C7_Richardson: bool
    l_use_C11_Richardson: bool
    l_use_shear_Richardson: bool
    l_brunt_vaisala_freq_moist: bool
    l_use_thvm_in_bv_freq: bool
    l_rcm_supersat_adj: bool
    l_damp_wp3_Skw_squared: bool
    l_prescribed_avg_deltaz: bool
    l_lmm_stepping: bool
    l_e3sm_config: bool
    l_vary_convect_depth: bool
    l_use_tke_in_wp3_pr_turb_term: bool
    l_use_tke_in_wp2_wp3_K_dfsn: bool
    l_use_wp3_lim_with_smth_Heaviside: bool
    l_smooth_Heaviside_tau_wpxp: bool
    l_modify_limiters_for_cnvg_test: bool
    l_enable_relaxed_clipping: bool
    l_linearize_pbl_winds: bool
    l_mono_flux_lim_thlm: bool
    l_mono_flux_lim_rtm: bool
    l_mono_flux_lim_um: bool
    l_mono_flux_lim_vm: bool
    l_mono_flux_lim_spikefix: bool
    l_host_applies_sfc_fluxes: bool
    l_wp2_fill_holes_tke: bool
    l_add_dycore_grid: bool
