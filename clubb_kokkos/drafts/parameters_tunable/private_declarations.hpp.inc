namespace parameters_tunable_detail {
template <class Matrix>
void check_parameter_hard_bounds(int ngrdcol, Matrix clubb_params, err_info_type& err_info // InOut
);
template <class Column, class Matrix>
void pack_parameters(
    int ngrdcol,                                                                  // In
    Column C1, Column C1b, Column C1c, Column C2rt, Column C2thl, Column C2rtthl, // In
    Column C4, Column C_uu_shr, Column C_uu_buoy, Column C6rt, Column C6rtb, Column C6rtc,
    Column C6thl, Column C6thlb, Column C6thlc,                           // In
    Column C7, Column C7b, Column C7c, Column C8, Column C8b, Column C10, // In
    Column C11, Column C11b, Column C11c, Column C12, Column C13, Column C14, Column C_wp2_pr_dfsn,
    Column C_wp3_pr_tp,                                                                 // In
    Column C_wp3_pr_turb, Column C_wp3_pr_dfsn, Column C_wp2_splat,                     // In
    Column C6rt_Lscale0, Column C6thl_Lscale0, Column C7_Lscale0, Column wpxp_L_thresh, // In
    Column c_K, Column c_K1, Column nu1, Column c_K2, Column nu2, Column c_K6, Column nu6,
    Column c_K8, Column nu8, // In
    Column c_K9, Column nu9, Column nu10, Column c_K_hm, Column c_K_hmb, Column K_hm_min_coef,
    Column nu_hm,                                                             // In
    Column slope_coef_spread_DG_means_w, Column pdf_component_stdev_factor_w, // In
    Column coef_spread_DG_means_rt, Column coef_spread_DG_means_thl,          // In
    Column gamma_coef, Column gamma_coefb, Column gamma_coefc, Column Lscale_mu_entrain,
    Column beta, Column lmin_coef,                                                      // In
    Column omicron, Column zeta_vrnce_rat, Column upsilon_precip_frac_rat,              // In
    Column lambda0_stability_coef, Column mult_coef, Column taumin, Column taumax,      // In
    Column Lscale_mu_coef, Column Lscale_pert_coef, Column alpha_corr,                  // In
    Column Skw_denom_coef, Column c_K10, Column c_K10h, Column thlp2_rad_coef,          // In
    Column thlp2_rad_cloud_frac_thresh, Column up2_sfc_coef,                            // In
    Column Skw_max_mag, Column xp3_coef_base, Column xp3_coef_slope,                    // In
    Column altitude_threshold, Column rtp2_clip_coef, Column C_invrs_tau_bkgnd,         // In
    Column C_invrs_tau_sfc, Column C_invrs_tau_shear, Column C_invrs_tau_N2,            // In
    Column C_invrs_tau_N2_wp2, Column C_invrs_tau_N2_xp2,                               // In
    Column C_invrs_tau_N2_wpxp, Column C_invrs_tau_N2_clear_wp3,                        // In
    Column C_invrs_tau_wpxp_Ri, Column C_invrs_tau_wpxp_N2_thresh,                      // In
    Column Cx_min, Column Cx_max, Column Richardson_num_min, Column Richardson_num_max, // In
    Column wpxp_Ri_exp, Column a3_coef_min, Column a_const, Column bv_efold,
    Column z_displace, // In
    Column Lscale_tke_multiplier, Column Lscale_pdf_launch_frac, Column Lscale_up_weight,
    Matrix clubb_params // Out
);
template <class Vector>
void unpack_parameters(
    Vector clubb_params,                                                                // In
    real_t& C1, real_t& C1b, real_t& C1c, real_t& C2rt, real_t& C2thl, real_t& C2rtthl, // Out
    real_t& C4, real_t& C_uu_shr, real_t& C_uu_buoy, real_t& C6rt, real_t& C6rtb, real_t& C6rtc,
    real_t& C6thl, real_t& C6thlb, real_t& C6thlc,                              // Out
    real_t& C7, real_t& C7b, real_t& C7c, real_t& C8, real_t& C8b, real_t& C10, // Out
    real_t& C11, real_t& C11b, real_t& C11c, real_t& C12, real_t& C13, real_t& C14,
    real_t& C_wp2_pr_dfsn, real_t& C_wp3_pr_tp,                                             // Out
    real_t& C_wp3_pr_turb, real_t& C_wp3_pr_dfsn, real_t& C_wp2_splat,                      // Out
    real_t& C6rt_Lscale0, real_t& C6thl_Lscale0, real_t& C7_Lscale0, real_t& wpxp_L_thresh, // Out
    real_t& c_K, real_t& c_K1, real_t& nu1, real_t& c_K2, real_t& nu2, real_t& c_K6, real_t& nu6,
    real_t& c_K8, real_t& nu8, // Out
    real_t& c_K9, real_t& nu9, real_t& nu10, real_t& c_K_hm, real_t& c_K_hmb, real_t& K_hm_min_coef,
    real_t& nu_hm,                                                              // Out
    real_t& slope_coef_spread_DG_means_w, real_t& pdf_component_stdev_factor_w, // Out
    real_t& coef_spread_DG_means_rt, real_t& coef_spread_DG_means_thl,          // Out
    real_t& gamma_coef, real_t& gamma_coefb, real_t& gamma_coefc, real_t& Lscale_mu_entrain,
    real_t& beta, real_t& lmin_coef,                                                        // Out
    real_t& omicron, real_t& zeta_vrnce_rat, real_t& upsilon_precip_frac_rat,               // Out
    real_t& lambda0_stability_coef, real_t& mult_coef, real_t& taumin, real_t& taumax,      // Out
    real_t& Lscale_mu_coef, real_t& Lscale_pert_coef, real_t& alpha_corr,                   // Out
    real_t& Skw_denom_coef, real_t& c_K10, real_t& c_K10h, real_t& thlp2_rad_coef,          // Out
    real_t& thlp2_rad_cloud_frac_thresh, real_t& up2_sfc_coef,                              // Out
    real_t& Skw_max_mag, real_t& xp3_coef_base, real_t& xp3_coef_slope,                     // Out
    real_t& altitude_threshold, real_t& rtp2_clip_coef, real_t& C_invrs_tau_bkgnd,          // Out
    real_t& C_invrs_tau_sfc, real_t& C_invrs_tau_shear, real_t& C_invrs_tau_N2,             // Out
    real_t& C_invrs_tau_N2_wp2, real_t& C_invrs_tau_N2_xp2,                                 // Out
    real_t& C_invrs_tau_N2_wpxp, real_t& C_invrs_tau_N2_clear_wp3,                          // Out
    real_t& C_invrs_tau_wpxp_Ri, real_t& C_invrs_tau_wpxp_N2_thresh,                        // Out
    real_t& Cx_min, real_t& Cx_max, real_t& Richardson_num_min, real_t& Richardson_num_max, // Out
    real_t& wpxp_Ri_exp, real_t& a3_coef_min, real_t& a_const, real_t& bv_efold,
    real_t& z_displace,                                                                     // Out
    real_t& Lscale_tke_multiplier, real_t& Lscale_pdf_launch_frac, real_t& Lscale_up_weight // Out
);
} // namespace parameters_tunable_detail
