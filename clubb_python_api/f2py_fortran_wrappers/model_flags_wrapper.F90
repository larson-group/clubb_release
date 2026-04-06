! model_flags_wrapper.F90 — wrappers for CLUBB_core/model_flags.F90 routines

module model_flags_wrapper

  use clubb_precision, only: core_rknd
  use derived_type_storage

  implicit none
  private

  public :: f2py_set_default_config_flags
  public :: f2py_init_config_flags

contains

  subroutine f2py_set_default_config_flags( &
      iiPDF_type, ipdf_call_placement, penta_solve_method, &
      tridiag_solve_method, saturation_formula, grid_remap_method, &
      grid_adapt_in_time_method, fill_holes_type, &
      l_use_precip_frac, l_predict_upwp_vpwp, &
      l_ho_nontrad_coriolis, l_ho_trad_coriolis, &
      l_min_wp2_from_corr_wx, l_min_xp2_from_corr_wx, &
      l_C2_cloud_frac, l_diffuse_rtm_and_thlm, &
      l_stability_correct_Kh_N2_zm, l_calc_thlp2_rad, &
      l_upwind_xpyp_ta, l_upwind_xm_ma, &
      l_uv_nudge, l_rtm_nudge, l_tke_aniso, &
      l_vert_avg_closure, l_trapezoidal_rule_zt, &
      l_trapezoidal_rule_zm, l_call_pdf_closure_twice, &
      l_standard_term_ta, l_partial_upwind_wp3, &
      l_godunov_upwind_wpxp_ta, l_godunov_upwind_xpyp_ta, &
      l_use_cloud_cover, l_diagnose_correlations, &
      l_calc_w_corr, l_const_Nc_in_cloud, &
      l_fix_w_chi_eta_correlations, l_stability_correct_tau_zm, &
      l_damp_wp2_using_em, l_do_expldiff_rtm_thlm, &
      l_Lscale_plume_centered, l_diag_Lscale_from_tau, &
      l_use_C7_Richardson, l_use_C11_Richardson, &
      l_use_shear_Richardson, l_brunt_vaisala_freq_moist, &
      l_use_thvm_in_bv_freq, l_rcm_supersat_adj, &
      l_damp_wp3_Skw_squared, l_prescribed_avg_deltaz, &
      l_lmm_stepping, l_e3sm_config, l_vary_convect_depth, &
      l_use_tke_in_wp3_pr_turb_term, l_use_tke_in_wp2_wp3_K_dfsn, &
      l_use_wp3_lim_with_smth_Heaviside, &
      l_smooth_Heaviside_tau_wpxp, &
      l_modify_limiters_for_cnvg_test, &
      l_enable_relaxed_clipping, l_linearize_pbl_winds, &
      l_mono_flux_lim_thlm, l_mono_flux_lim_rtm, &
      l_mono_flux_lim_um, l_mono_flux_lim_vm, &
      l_mono_flux_lim_spikefix, l_host_applies_sfc_fluxes, &
      l_wp2_fill_holes_tke, l_add_dycore_grid) &
    bind(C, name="f2py_set_default_config_flags_")

    use model_flags, only: set_default_clubb_config_flags_api

    ! 8 integer flags (intent out)
    integer, intent(out) :: iiPDF_type, ipdf_call_placement, &
      penta_solve_method, tridiag_solve_method, saturation_formula, &
      grid_remap_method, grid_adapt_in_time_method, fill_holes_type

    ! 59 logical flags (intent out)
    logical, intent(out) :: &
      l_use_precip_frac, l_predict_upwp_vpwp, &
      l_ho_nontrad_coriolis, l_ho_trad_coriolis, &
      l_min_wp2_from_corr_wx, l_min_xp2_from_corr_wx, &
      l_C2_cloud_frac, l_diffuse_rtm_and_thlm, &
      l_stability_correct_Kh_N2_zm, l_calc_thlp2_rad, &
      l_upwind_xpyp_ta, l_upwind_xm_ma, &
      l_uv_nudge, l_rtm_nudge, l_tke_aniso, &
      l_vert_avg_closure, l_trapezoidal_rule_zt, &
      l_trapezoidal_rule_zm, l_call_pdf_closure_twice, &
      l_standard_term_ta, l_partial_upwind_wp3, &
      l_godunov_upwind_wpxp_ta, l_godunov_upwind_xpyp_ta, &
      l_use_cloud_cover, l_diagnose_correlations, &
      l_calc_w_corr, l_const_Nc_in_cloud, &
      l_fix_w_chi_eta_correlations, l_stability_correct_tau_zm, &
      l_damp_wp2_using_em, l_do_expldiff_rtm_thlm, &
      l_Lscale_plume_centered, l_diag_Lscale_from_tau, &
      l_use_C7_Richardson, l_use_C11_Richardson, &
      l_use_shear_Richardson, l_brunt_vaisala_freq_moist, &
      l_use_thvm_in_bv_freq, l_rcm_supersat_adj, &
      l_damp_wp3_Skw_squared, l_prescribed_avg_deltaz, &
      l_lmm_stepping, l_e3sm_config, l_vary_convect_depth, &
      l_use_tke_in_wp3_pr_turb_term, l_use_tke_in_wp2_wp3_K_dfsn, &
      l_use_wp3_lim_with_smth_Heaviside, &
      l_smooth_Heaviside_tau_wpxp, &
      l_modify_limiters_for_cnvg_test, &
      l_enable_relaxed_clipping, l_linearize_pbl_winds, &
      l_mono_flux_lim_thlm, l_mono_flux_lim_rtm, &
      l_mono_flux_lim_um, l_mono_flux_lim_vm, &
      l_mono_flux_lim_spikefix, l_host_applies_sfc_fluxes, &
      l_wp2_fill_holes_tke, l_add_dycore_grid


    call set_default_clubb_config_flags_api( &
      iiPDF_type, ipdf_call_placement, penta_solve_method, &
      tridiag_solve_method, saturation_formula, grid_remap_method, &
      grid_adapt_in_time_method, fill_holes_type, &
      l_use_precip_frac, l_predict_upwp_vpwp, &
      l_ho_nontrad_coriolis, l_ho_trad_coriolis, &
      l_min_wp2_from_corr_wx, l_min_xp2_from_corr_wx, &
      l_C2_cloud_frac, l_diffuse_rtm_and_thlm, &
      l_stability_correct_Kh_N2_zm, l_calc_thlp2_rad, &
      l_upwind_xpyp_ta, l_upwind_xm_ma, &
      l_uv_nudge, l_rtm_nudge, l_tke_aniso, &
      l_vert_avg_closure, l_trapezoidal_rule_zt, &
      l_trapezoidal_rule_zm, l_call_pdf_closure_twice, &
      l_standard_term_ta, l_partial_upwind_wp3, &
      l_godunov_upwind_wpxp_ta, l_godunov_upwind_xpyp_ta, &
      l_use_cloud_cover, l_diagnose_correlations, &
      l_calc_w_corr, l_const_Nc_in_cloud, &
      l_fix_w_chi_eta_correlations, l_stability_correct_tau_zm, &
      l_damp_wp2_using_em, l_do_expldiff_rtm_thlm, &
      l_Lscale_plume_centered, l_diag_Lscale_from_tau, &
      l_use_C7_Richardson, l_use_C11_Richardson, &
      l_use_shear_Richardson, l_brunt_vaisala_freq_moist, &
      l_use_thvm_in_bv_freq, l_rcm_supersat_adj, &
      l_damp_wp3_Skw_squared, l_prescribed_avg_deltaz, &
      l_lmm_stepping, l_e3sm_config, l_vary_convect_depth, &
      l_use_tke_in_wp3_pr_turb_term, l_use_tke_in_wp2_wp3_K_dfsn, &
      l_use_wp3_lim_with_smth_Heaviside, &
      l_smooth_Heaviside_tau_wpxp, &
      l_modify_limiters_for_cnvg_test, &
      l_enable_relaxed_clipping, l_linearize_pbl_winds, &
      l_mono_flux_lim_thlm, l_mono_flux_lim_rtm, &
      l_mono_flux_lim_um, l_mono_flux_lim_vm, &
      l_mono_flux_lim_spikefix, l_host_applies_sfc_fluxes, &
      l_wp2_fill_holes_tke, l_add_dycore_grid)


  end subroutine f2py_set_default_config_flags

  subroutine f2py_init_config_flags( &
      iiPDF_type, ipdf_call_placement, penta_solve_method, &
      tridiag_solve_method, saturation_formula, grid_remap_method, &
      grid_adapt_in_time_method, fill_holes_type, &
      l_use_precip_frac, l_predict_upwp_vpwp, &
      l_ho_nontrad_coriolis, l_ho_trad_coriolis, &
      l_min_wp2_from_corr_wx, l_min_xp2_from_corr_wx, &
      l_C2_cloud_frac, l_diffuse_rtm_and_thlm, &
      l_stability_correct_Kh_N2_zm, l_calc_thlp2_rad, &
      l_upwind_xpyp_ta, l_upwind_xm_ma, &
      l_uv_nudge, l_rtm_nudge, l_tke_aniso, &
      l_vert_avg_closure, l_trapezoidal_rule_zt, &
      l_trapezoidal_rule_zm, l_call_pdf_closure_twice, &
      l_standard_term_ta, l_partial_upwind_wp3, &
      l_godunov_upwind_wpxp_ta, l_godunov_upwind_xpyp_ta, &
      l_use_cloud_cover, l_diagnose_correlations, &
      l_calc_w_corr, l_const_Nc_in_cloud, &
      l_fix_w_chi_eta_correlations, l_stability_correct_tau_zm, &
      l_damp_wp2_using_em, l_do_expldiff_rtm_thlm, &
      l_Lscale_plume_centered, l_diag_Lscale_from_tau, &
      l_use_C7_Richardson, l_use_C11_Richardson, &
      l_use_shear_Richardson, l_brunt_vaisala_freq_moist, &
      l_use_thvm_in_bv_freq, l_rcm_supersat_adj, &
      l_damp_wp3_Skw_squared, l_prescribed_avg_deltaz, &
      l_lmm_stepping, l_e3sm_config, l_vary_convect_depth, &
      l_use_tke_in_wp3_pr_turb_term, l_use_tke_in_wp2_wp3_K_dfsn, &
      l_use_wp3_lim_with_smth_Heaviside, &
      l_smooth_Heaviside_tau_wpxp, &
      l_modify_limiters_for_cnvg_test, &
      l_enable_relaxed_clipping, l_linearize_pbl_winds, &
      l_mono_flux_lim_thlm, l_mono_flux_lim_rtm, &
      l_mono_flux_lim_um, l_mono_flux_lim_vm, &
      l_mono_flux_lim_spikefix, l_host_applies_sfc_fluxes, &
      l_wp2_fill_holes_tke, l_add_dycore_grid) &
    bind(C, name="f2py_init_config_flags_")

    use model_flags, only: initialize_clubb_config_flags_type_api

    ! 8 integer flags
    integer, intent(in) :: iiPDF_type, ipdf_call_placement, &
      penta_solve_method, tridiag_solve_method, saturation_formula, &
      grid_remap_method, grid_adapt_in_time_method, fill_holes_type

    ! 59 logical flags
    logical, intent(in) :: &
      l_use_precip_frac, l_predict_upwp_vpwp, &
      l_ho_nontrad_coriolis, l_ho_trad_coriolis, &
      l_min_wp2_from_corr_wx, l_min_xp2_from_corr_wx, &
      l_C2_cloud_frac, l_diffuse_rtm_and_thlm, &
      l_stability_correct_Kh_N2_zm, l_calc_thlp2_rad, &
      l_upwind_xpyp_ta, l_upwind_xm_ma, &
      l_uv_nudge, l_rtm_nudge, l_tke_aniso, &
      l_vert_avg_closure, l_trapezoidal_rule_zt, &
      l_trapezoidal_rule_zm, l_call_pdf_closure_twice, &
      l_standard_term_ta, l_partial_upwind_wp3, &
      l_godunov_upwind_wpxp_ta, l_godunov_upwind_xpyp_ta, &
      l_use_cloud_cover, l_diagnose_correlations, &
      l_calc_w_corr, l_const_Nc_in_cloud, &
      l_fix_w_chi_eta_correlations, l_stability_correct_tau_zm, &
      l_damp_wp2_using_em, l_do_expldiff_rtm_thlm, &
      l_Lscale_plume_centered, l_diag_Lscale_from_tau, &
      l_use_C7_Richardson, l_use_C11_Richardson, &
      l_use_shear_Richardson, l_brunt_vaisala_freq_moist, &
      l_use_thvm_in_bv_freq, l_rcm_supersat_adj, &
      l_damp_wp3_Skw_squared, l_prescribed_avg_deltaz, &
      l_lmm_stepping, l_e3sm_config, l_vary_convect_depth, &
      l_use_tke_in_wp3_pr_turb_term, l_use_tke_in_wp2_wp3_K_dfsn, &
      l_use_wp3_lim_with_smth_Heaviside, &
      l_smooth_Heaviside_tau_wpxp, &
      l_modify_limiters_for_cnvg_test, &
      l_enable_relaxed_clipping, l_linearize_pbl_winds, &
      l_mono_flux_lim_thlm, l_mono_flux_lim_rtm, &
      l_mono_flux_lim_um, l_mono_flux_lim_vm, &
      l_mono_flux_lim_spikefix, l_host_applies_sfc_fluxes, &
      l_wp2_fill_holes_tke, l_add_dycore_grid

    call initialize_clubb_config_flags_type_api( &
      iiPDF_type, ipdf_call_placement, penta_solve_method, &
      tridiag_solve_method, saturation_formula, grid_remap_method, &
      grid_adapt_in_time_method, fill_holes_type, &
      l_use_precip_frac, l_predict_upwp_vpwp, &
      l_ho_nontrad_coriolis, l_ho_trad_coriolis, &
      l_min_wp2_from_corr_wx, l_min_xp2_from_corr_wx, &
      l_C2_cloud_frac, l_diffuse_rtm_and_thlm, &
      l_stability_correct_Kh_N2_zm, l_calc_thlp2_rad, &
      l_upwind_xpyp_ta, l_upwind_xm_ma, &
      l_uv_nudge, l_rtm_nudge, l_tke_aniso, &
      l_vert_avg_closure, l_trapezoidal_rule_zt, &
      l_trapezoidal_rule_zm, l_call_pdf_closure_twice, &
      l_standard_term_ta, l_partial_upwind_wp3, &
      l_godunov_upwind_wpxp_ta, l_godunov_upwind_xpyp_ta, &
      l_use_cloud_cover, l_diagnose_correlations, &
      l_calc_w_corr, l_const_Nc_in_cloud, &
      l_fix_w_chi_eta_correlations, l_stability_correct_tau_zm, &
      l_damp_wp2_using_em, l_do_expldiff_rtm_thlm, &
      l_Lscale_plume_centered, l_diag_Lscale_from_tau, &
      l_use_C7_Richardson, l_use_C11_Richardson, &
      l_use_shear_Richardson, l_brunt_vaisala_freq_moist, &
      l_use_thvm_in_bv_freq, l_rcm_supersat_adj, &
      l_damp_wp3_Skw_squared, l_prescribed_avg_deltaz, &
      l_lmm_stepping, l_e3sm_config, &
      l_vary_convect_depth, &
      l_use_tke_in_wp3_pr_turb_term, &
      l_use_tke_in_wp2_wp3_K_dfsn, &
      l_use_wp3_lim_with_smth_Heaviside, &
      l_smooth_Heaviside_tau_wpxp, &
      l_modify_limiters_for_cnvg_test, &
      l_enable_relaxed_clipping, l_linearize_pbl_winds, &
      l_mono_flux_lim_thlm, l_mono_flux_lim_rtm, &
      l_mono_flux_lim_um, l_mono_flux_lim_vm, &
      l_mono_flux_lim_spikefix, l_host_applies_sfc_fluxes, &
      l_wp2_fill_holes_tke, l_add_dycore_grid, &
      stored_config_flags)

  end subroutine f2py_init_config_flags

end module model_flags_wrapper
