// Merge constants after the generated nu_vertical_res_dep type and before procedure fragments.
// Native parameters are immutable. No global View allocation or mutable cache is introduced.
namespace parameters_tunable_detail {
inline constexpr real_t no_lower_bound = -std::numeric_limits<real_t>::max();
inline constexpr real_t no_upper_bound = std::numeric_limits<real_t>::max();
inline constexpr real_t strictly_positive = std::numeric_limits<real_t>::min();
inline constexpr real_t strictly_less_than_one = one - std::numeric_limits<real_t>::epsilon();
inline constexpr real_t init_value =
    real_t(-999.); // Native missing-value sentinel (inactive check).
} // namespace parameters_tunable_detail

// These are referenced together often enough to keep a list. lmin_coef is the
// input parameter; the actual lmin model constant is computed from it.
// IMPORTANT: Keep this order identical to parameter_indices when adding parameters.
inline constexpr int parameter_name_length = 32;
inline constexpr std::array<std::string_view, nparams> params_list = {
    "C1                              ", "C1b                             ",
    "C1c                             ", "C2rt                            ",
    "C2thl                           ", "C2rtthl                         ",
    "C4                              ", "C_uu_shr                        ",
    "C_uu_buoy                       ", "C6rt                            ",
    "C6rtb                           ", "C6rtc                           ",
    "C6thl                           ", "C6thlb                          ",
    "C6thlc                          ", "C7                              ",
    "C7b                             ", "C7c                             ",
    "C8                              ", "C8b                             ",
    "C10                             ", "C11                             ",
    "C11b                            ", "C11c                            ",
    "C12                             ", "C13                             ",
    "C14                             ", "C_wp2_pr_dfsn                   ",
    "C_wp3_pr_tp                     ", "C_wp3_pr_turb                   ",
    "C_wp3_pr_dfsn                   ", "C_wp2_splat                     ",
    "C6rt_Lscale0                    ", "C6thl_Lscale0                   ",
    "C7_Lscale0                      ", "wpxp_L_thresh                   ",
    "c_K                             ", "c_K1                            ",
    "nu1                             ", "c_K2                            ",
    "nu2                             ", "c_K6                            ",
    "nu6                             ", "c_K8                            ",
    "nu8                             ", "c_K9                            ",
    "nu9                             ", "nu10                            ",
    "c_K_hm                          ", "c_K_hmb                         ",
    "K_hm_min_coef                   ", "nu_hm                           ",
    "slope_coef_spread_DG_means_w    ", "pdf_component_stdev_factor_w    ",
    "coef_spread_DG_means_rt         ", "coef_spread_DG_means_thl        ",
    "gamma_coef                      ", "gamma_coefb                     ",
    "gamma_coefc                     ", "Lscale_mu_entrain               ",
    "beta                            ", "lmin_coef                       ",
    "omicron                         ", "zeta_vrnce_rat                  ",
    "upsilon_precip_frac_rat         ", "lambda0_stability_coef          ",
    "mult_coef                       ", "taumin                          ",
    "taumax                          ", "Lscale_mu_coef                  ",
    "Lscale_pert_coef                ", "alpha_corr                      ",
    "Skw_denom_coef                  ", "c_K10                           ",
    "c_K10h                          ", "thlp2_rad_coef                  ",
    "thlp2_rad_cloud_frac_thresh     ", "up2_sfc_coef                    ",
    "Skw_max_mag                     ", "C_invrs_tau_bkgnd               ",
    "C_invrs_tau_sfc                 ", "C_invrs_tau_shear               ",
    "C_invrs_tau_N2                  ", "C_invrs_tau_N2_wp2              ",
    "C_invrs_tau_N2_xp2              ", "C_invrs_tau_N2_wpxp             ",
    "C_invrs_tau_N2_clear_wp3        ", "C_invrs_tau_wpxp_Ri             ",
    "C_invrs_tau_wpxp_N2_thresh      ", "xp3_coef_base                   ",
    "xp3_coef_slope                  ", "altitude_threshold              ",
    "rtp2_clip_coef                  ", "Cx_min                          ",
    "Cx_max                          ", "Richardson_num_min              ",
    "Richardson_num_max              ", "a3_coef_min                     ",
    "a_const                         ", "bv_efold                        ",
    "wpxp_Ri_exp                     ", "z_displace                      ",
    "Lscale_tke_multiplier           ", "Lscale_pdf_launch_frac          ",
    "Lscale_up_weight                ",
};

// Hard configuration bounds in the same order as params_list.
// Rows are [minimum,maximum]. no_*_bound denotes an open endpoint;
// strictly_positive/strictly_less_than_one encode open physical endpoints.
// These are rejection limits, not recommended tuning intervals. The initializer
// is Fortran-flat: minimum,maximum for each parameter in order.
inline constexpr std::array<real_t, 2 * nparams> parameter_hard_bounds = {
    zero,
    parameters_tunable_detail::no_upper_bound, // C1
    zero,
    parameters_tunable_detail::no_upper_bound, // C1b
    zero,
    parameters_tunable_detail::no_upper_bound, // C1c
    zero,
    parameters_tunable_detail::no_upper_bound, // C2rt
    zero,
    parameters_tunable_detail::no_upper_bound, // C2thl
    zero,
    parameters_tunable_detail::no_upper_bound, // C2rtthl
    zero,
    parameters_tunable_detail::no_upper_bound, // C4
    zero,
    one, // C_uu_shr
    zero,
    one, // C_uu_buoy
    zero,
    parameters_tunable_detail::no_upper_bound, // C6rt
    zero,
    parameters_tunable_detail::no_upper_bound, // C6rtb
    zero,
    parameters_tunable_detail::no_upper_bound, // C6rtc
    zero,
    parameters_tunable_detail::no_upper_bound, // C6thl
    zero,
    parameters_tunable_detail::no_upper_bound, // C6thlb
    zero,
    parameters_tunable_detail::no_upper_bound, // C6thlc
    zero,
    one, // C7
    zero,
    one, // C7b
    zero,
    parameters_tunable_detail::no_upper_bound, // C7c
    zero,
    parameters_tunable_detail::no_upper_bound, // C8
    zero,
    parameters_tunable_detail::no_upper_bound, // C8b
    zero,
    parameters_tunable_detail::no_upper_bound, // C10
    zero,
    one, // C11
    zero,
    one, // C11b
    zero,
    parameters_tunable_detail::no_upper_bound, // C11c
    zero,
    parameters_tunable_detail::no_upper_bound, // C12
    zero,
    parameters_tunable_detail::no_upper_bound, // C13
    zero,
    parameters_tunable_detail::no_upper_bound, // C14
    zero,
    parameters_tunable_detail::no_upper_bound, // C_wp2_pr_dfsn
    zero,
    parameters_tunable_detail::no_upper_bound, // C_wp3_pr_tp
    zero,
    parameters_tunable_detail::no_upper_bound, // C_wp3_pr_turb
    zero,
    parameters_tunable_detail::no_upper_bound, // C_wp3_pr_dfsn
    zero,
    parameters_tunable_detail::no_upper_bound, // C_wp2_splat
    zero,
    parameters_tunable_detail::no_upper_bound, // C6rt_Lscale0
    zero,
    parameters_tunable_detail::no_upper_bound, // C6thl_Lscale0
    zero,
    parameters_tunable_detail::no_upper_bound, // C7_Lscale0
    zero,
    parameters_tunable_detail::no_upper_bound, // wpxp_L_thresh
    zero,
    parameters_tunable_detail::no_upper_bound, // c_K
    zero,
    parameters_tunable_detail::no_upper_bound, // c_K1
    zero,
    parameters_tunable_detail::no_upper_bound, // nu1
    zero,
    parameters_tunable_detail::no_upper_bound, // c_K2
    zero,
    parameters_tunable_detail::no_upper_bound, // nu2
    zero,
    parameters_tunable_detail::no_upper_bound, // c_K6
    zero,
    parameters_tunable_detail::no_upper_bound, // nu6
    zero,
    parameters_tunable_detail::no_upper_bound, // c_K8
    zero,
    parameters_tunable_detail::no_upper_bound, // nu8
    zero,
    parameters_tunable_detail::no_upper_bound, // c_K9
    zero,
    parameters_tunable_detail::no_upper_bound, // nu9
    zero,
    parameters_tunable_detail::no_upper_bound, // nu10
    zero,
    parameters_tunable_detail::no_upper_bound, // c_K_hm
    zero,
    parameters_tunable_detail::no_upper_bound, // c_K_hmb
    zero,
    parameters_tunable_detail::no_upper_bound, // K_hm_min_coef
    zero,
    parameters_tunable_detail::no_upper_bound, // nu_hm
    parameters_tunable_detail::strictly_positive,
    parameters_tunable_detail::no_upper_bound, // slope_coef_spread_DG_means_w
    parameters_tunable_detail::strictly_positive,
    parameters_tunable_detail::no_upper_bound, // pdf_component_stdev_factor_w
    zero,
    parameters_tunable_detail::strictly_less_than_one, // coef_spread_DG_means_rt
    zero,
    parameters_tunable_detail::strictly_less_than_one, // coef_spread_DG_means_thl
    zero,
    parameters_tunable_detail::no_upper_bound, // gamma_coef
    zero,
    parameters_tunable_detail::no_upper_bound, // gamma_coefb
    zero,
    parameters_tunable_detail::no_upper_bound, // gamma_coefc
    parameters_tunable_detail::strictly_positive,
    parameters_tunable_detail::no_upper_bound, // Lscale_mu_entrain
    zero,
    real_t(3.0), // beta
    zero,
    parameters_tunable_detail::no_upper_bound, // lmin_coef
    parameters_tunable_detail::strictly_positive,
    one, // omicron
    -one + std::numeric_limits<real_t>::epsilon(),
    parameters_tunable_detail::no_upper_bound, // zeta_vrnce_rat
    zero,
    one, // upsilon_precip_frac_rat
    zero,
    parameters_tunable_detail::no_upper_bound, // lambda0_stability_coef
    zero,
    parameters_tunable_detail::no_upper_bound, // mult_coef
    zero,
    parameters_tunable_detail::no_upper_bound, // taumin
    zero,
    parameters_tunable_detail::no_upper_bound, // taumax
    zero,
    parameters_tunable_detail::no_upper_bound, // Lscale_mu_coef
    zero,
    parameters_tunable_detail::no_upper_bound, // Lscale_pert_coef
    zero,
    parameters_tunable_detail::no_upper_bound, // alpha_corr
    zero,
    parameters_tunable_detail::no_upper_bound, // Skw_denom_coef
    zero,
    parameters_tunable_detail::no_upper_bound, // c_K10
    zero,
    parameters_tunable_detail::no_upper_bound, // c_K10h
    zero,
    parameters_tunable_detail::no_upper_bound, // thlp2_rad_coef
    zero,
    parameters_tunable_detail::no_upper_bound, // thlp2_rad_cloud_frac_thresh
    zero,
    parameters_tunable_detail::no_upper_bound, // up2_sfc_coef
    zero,
    parameters_tunable_detail::no_upper_bound, // Skw_max_mag
    zero,
    parameters_tunable_detail::no_upper_bound, // C_invrs_tau_bkgnd
    zero,
    parameters_tunable_detail::no_upper_bound, // C_invrs_tau_sfc
    zero,
    parameters_tunable_detail::no_upper_bound, // C_invrs_tau_shear
    zero,
    parameters_tunable_detail::no_upper_bound, // C_invrs_tau_N2
    zero,
    parameters_tunable_detail::no_upper_bound, // C_invrs_tau_N2_wp2
    zero,
    parameters_tunable_detail::no_upper_bound, // C_invrs_tau_N2_xp2
    zero,
    parameters_tunable_detail::no_upper_bound, // C_invrs_tau_N2_wpxp
    zero,
    parameters_tunable_detail::no_upper_bound, // C_invrs_tau_N2_clear_wp3
    zero,
    parameters_tunable_detail::no_upper_bound, // C_invrs_tau_wpxp_Ri
    zero,
    parameters_tunable_detail::no_upper_bound, // C_invrs_tau_wpxp_N2_thresh
    zero,
    parameters_tunable_detail::no_upper_bound, // xp3_coef_base
    zero,
    parameters_tunable_detail::no_upper_bound, // xp3_coef_slope
    zero,
    parameters_tunable_detail::no_upper_bound, // altitude_threshold
    zero,
    parameters_tunable_detail::no_upper_bound, // rtp2_clip_coef
    zero,
    parameters_tunable_detail::no_upper_bound, // Cx_min
    zero,
    parameters_tunable_detail::no_upper_bound, // Cx_max
    zero,
    parameters_tunable_detail::no_upper_bound, // Richardson_num_min
    zero,
    parameters_tunable_detail::no_upper_bound, // Richardson_num_max
    one,
    real_t(3.0), // a3_coef_min
    zero,
    parameters_tunable_detail::no_upper_bound, // a_const
    zero,
    parameters_tunable_detail::no_upper_bound, // bv_efold
    zero,
    parameters_tunable_detail::no_upper_bound, // wpxp_Ri_exp
    zero,
    parameters_tunable_detail::no_upper_bound, // z_displace
    parameters_tunable_detail::strictly_positive,
    parameters_tunable_detail::no_upper_bound, // Lscale_tke_multiplier
    zero,
    one, // Lscale_pdf_launch_frac
    zero,
    one, // Lscale_up_weight
};
