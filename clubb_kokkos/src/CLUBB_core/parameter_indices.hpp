// Tunable-parameter column indices from parameter_indices.F90.
// nparams is the column count; every i-prefixed index is converted once to zero-based C++.
// Keep params_list and the native packing/unpacking order consistent when adding parameters.
#pragma once

namespace clubb_kokkos {

inline constexpr int nparams = 105; // Total tunable parameters

//***************************************************************
//                    ***** IMPORTANT *****
// If you change the order of these parameters, you will need to
// change the order of params_list as well or the tuner will
// break!
//                    ***** IMPORTANT *****
//***************************************************************

inline constexpr int iC1            = 1 - 1;
inline constexpr int iC1b           = 2 - 1;
inline constexpr int iC1c           = 3 - 1;
inline constexpr int iC2rt          = 4 - 1;
inline constexpr int iC2thl         = 5 - 1;
inline constexpr int iC2rtthl       = 6 - 1;
inline constexpr int iC4            = 7 - 1;
inline constexpr int iC_uu_shr      = 8 - 1;
inline constexpr int iC_uu_buoy     = 9 - 1;
inline constexpr int iC6rt          = 10 - 1;
inline constexpr int iC6rtb         = 11 - 1;
inline constexpr int iC6rtc         = 12 - 1;
inline constexpr int iC6thl         = 13 - 1;
inline constexpr int iC6thlb        = 14 - 1;
inline constexpr int iC6thlc        = 15 - 1;
inline constexpr int iC7            = 16 - 1;
inline constexpr int iC7b           = 17 - 1;
inline constexpr int iC7c           = 18 - 1;
inline constexpr int iC8            = 19 - 1;
inline constexpr int iC8b           = 20 - 1;
inline constexpr int iC10           = 21 - 1;
inline constexpr int iC11           = 22 - 1;
inline constexpr int iC11b          = 23 - 1;
inline constexpr int iC11c          = 24 - 1;
inline constexpr int iC12           = 25 - 1;
inline constexpr int iC13           = 26 - 1;
inline constexpr int iC14           = 27 - 1;
inline constexpr int iC_wp2_pr_dfsn = 28 - 1;
inline constexpr int iC_wp3_pr_tp   = 29 - 1;
inline constexpr int iC_wp3_pr_turb = 30 - 1;
inline constexpr int iC_wp3_pr_dfsn = 31 - 1;
inline constexpr int iC_wp2_splat   = 32 - 1;

inline constexpr int iC6rt_Lscale0  = 33 - 1;
inline constexpr int iC6thl_Lscale0 = 34 - 1;
inline constexpr int iC7_Lscale0    = 35 - 1;
inline constexpr int iwpxp_L_thresh = 36 - 1;

inline constexpr int ic_K           = 37 - 1;
inline constexpr int ic_K1          = 38 - 1;
inline constexpr int inu1           = 39 - 1;
inline constexpr int ic_K2          = 40 - 1;
inline constexpr int inu2           = 41 - 1;
inline constexpr int ic_K6          = 42 - 1;
inline constexpr int inu6           = 43 - 1;
inline constexpr int ic_K8          = 44 - 1;
inline constexpr int inu8           = 45 - 1;
inline constexpr int ic_K9          = 46 - 1;
inline constexpr int inu9           = 47 - 1;
inline constexpr int inu10          = 48 - 1;
inline constexpr int ic_K_hm        = 49 - 1;
inline constexpr int ic_K_hmb       = 50 - 1;
inline constexpr int iK_hm_min_coef = 51 - 1;
inline constexpr int inu_hm         = 52 - 1;

inline constexpr int islope_coef_spread_DG_means_w = 53 - 1;
inline constexpr int ipdf_component_stdev_factor_w = 54 - 1;
inline constexpr int icoef_spread_DG_means_rt      = 55 - 1;
inline constexpr int icoef_spread_DG_means_thl     = 56 - 1;
inline constexpr int igamma_coef                   = 57 - 1;
inline constexpr int igamma_coefb                  = 58 - 1;
inline constexpr int igamma_coefc                  = 59 - 1;
inline constexpr int iLscale_mu_entrain            = 60 - 1;
inline constexpr int ibeta                         = 61 - 1;
inline constexpr int ilmin_coef                    = 62 - 1;
inline constexpr int iomicron                      = 63 - 1;
inline constexpr int izeta_vrnce_rat               = 64 - 1;
inline constexpr int iupsilon_precip_frac_rat      = 65 - 1;
inline constexpr int ilambda0_stability_coef       = 66 - 1;
inline constexpr int imult_coef                    = 67 - 1;
inline constexpr int itaumin                       = 68 - 1;
inline constexpr int itaumax                       = 69 - 1;
inline constexpr int iLscale_mu_coef               = 70 - 1;
inline constexpr int iLscale_pert_coef             = 71 - 1;
inline constexpr int ialpha_corr                   = 72 - 1;
inline constexpr int iSkw_denom_coef               = 73 - 1;
inline constexpr int ic_K10                        = 74 - 1;
inline constexpr int ic_K10h                       = 75 - 1;
inline constexpr int ithlp2_rad_coef               = 76 - 1;
inline constexpr int ithlp2_rad_cloud_frac_thresh  = 77 - 1;
inline constexpr int iup2_sfc_coef                 = 78 - 1;
inline constexpr int iSkw_max_mag                  = 79 - 1;
inline constexpr int iC_invrs_tau_bkgnd            = 80 - 1;
inline constexpr int iC_invrs_tau_sfc              = 81 - 1;
inline constexpr int iC_invrs_tau_shear            = 82 - 1;
inline constexpr int iC_invrs_tau_N2               = 83 - 1;
inline constexpr int iC_invrs_tau_N2_wp2           = 84 - 1;
inline constexpr int iC_invrs_tau_N2_xp2           = 85 - 1;
inline constexpr int iC_invrs_tau_N2_wpxp          = 86 - 1;
inline constexpr int iC_invrs_tau_N2_clear_wp3     = 87 - 1;
inline constexpr int iC_invrs_tau_wpxp_Ri          = 88 - 1;
inline constexpr int iC_invrs_tau_wpxp_N2_thresh   = 89 - 1;
inline constexpr int ixp3_coef_base                = 90 - 1;
inline constexpr int ixp3_coef_slope               = 91 - 1;
inline constexpr int ialtitude_threshold           = 92 - 1;
inline constexpr int irtp2_clip_coef               = 93 - 1;
inline constexpr int iCx_min                       = 94 - 1;
inline constexpr int iCx_max                       = 95 - 1;
inline constexpr int iRichardson_num_min           = 96 - 1;
inline constexpr int iRichardson_num_max           = 97 - 1;
inline constexpr int ia3_coef_min                  = 98 - 1;
inline constexpr int ia_const                      = 99 - 1;
inline constexpr int ibv_efold                     = 100 - 1;
inline constexpr int iwpxp_Ri_exp                  = 101 - 1;
inline constexpr int iz_displace                   = 102 - 1;
inline constexpr int iLscale_tke_multiplier        = 103 - 1;
inline constexpr int iLscale_pdf_launch_frac       = 104 - 1;
inline constexpr int iLscale_up_weight             = 105 - 1;

} // namespace clubb_kokkos
