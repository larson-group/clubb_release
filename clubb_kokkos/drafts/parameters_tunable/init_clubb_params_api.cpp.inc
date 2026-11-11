// Merge inside namespace clubb_kokkos in parameters_tunable.cpp.
void init_clubb_params_api(
    int ngrdcol, int iunit, const std::string& filename,                       // In
    Kokkos::View<real_t**, Kokkos::LayoutLeft, Kokkos::HostSpace> clubb_params // Out
) {

  // Description:
  // Read a namelist containing the model parameters

  // References:
  // None
  //-----------------------------------------------------------------------
  // ngrdcol: Number of columns [#]. iunit: Native required logical-unit argument [-].
  // filename: Optional-by-value (empty means defaults) parameter namelist path [-].
  // clubb_params: Output parameter matrix, shape (ngrdcol,nparams) [-].
  // One call-local HostSpace allocation owns all Native local coefficient arrays.
  // Column subviews are contiguous because LayoutLeft stores columns first.
  using Host2D = Kokkos::View<real_t**, Kokkos::LayoutLeft, Kokkos::HostSpace>;
  Host2D local_parameters(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "init_clubb_params_api::local_parameters"),
      ngrdcol, nparams);
  auto C1 = Kokkos::subview(local_parameters, Kokkos::ALL, iC1);
  auto C1b = Kokkos::subview(local_parameters, Kokkos::ALL, iC1b);
  auto C1c = Kokkos::subview(local_parameters, Kokkos::ALL, iC1c);
  auto C2rt = Kokkos::subview(local_parameters, Kokkos::ALL, iC2rt);
  auto C2thl = Kokkos::subview(local_parameters, Kokkos::ALL, iC2thl);
  auto C2rtthl = Kokkos::subview(local_parameters, Kokkos::ALL, iC2rtthl);
  auto C4 = Kokkos::subview(local_parameters, Kokkos::ALL, iC4);
  auto C_uu_shr = Kokkos::subview(local_parameters, Kokkos::ALL, iC_uu_shr);
  auto C_uu_buoy = Kokkos::subview(local_parameters, Kokkos::ALL, iC_uu_buoy);
  auto C6rt = Kokkos::subview(local_parameters, Kokkos::ALL, iC6rt);
  auto C6rtb = Kokkos::subview(local_parameters, Kokkos::ALL, iC6rtb);
  auto C6rtc = Kokkos::subview(local_parameters, Kokkos::ALL, iC6rtc);
  auto C6thl = Kokkos::subview(local_parameters, Kokkos::ALL, iC6thl);
  auto C6thlb = Kokkos::subview(local_parameters, Kokkos::ALL, iC6thlb);
  auto C6thlc = Kokkos::subview(local_parameters, Kokkos::ALL, iC6thlc);
  auto C7 = Kokkos::subview(local_parameters, Kokkos::ALL, iC7);
  auto C7b = Kokkos::subview(local_parameters, Kokkos::ALL, iC7b);
  auto C7c = Kokkos::subview(local_parameters, Kokkos::ALL, iC7c);
  auto C8 = Kokkos::subview(local_parameters, Kokkos::ALL, iC8);
  auto C8b = Kokkos::subview(local_parameters, Kokkos::ALL, iC8b);
  auto C10 = Kokkos::subview(local_parameters, Kokkos::ALL, iC10);
  auto C11 = Kokkos::subview(local_parameters, Kokkos::ALL, iC11);
  auto C11b = Kokkos::subview(local_parameters, Kokkos::ALL, iC11b);
  auto C11c = Kokkos::subview(local_parameters, Kokkos::ALL, iC11c);
  auto C12 = Kokkos::subview(local_parameters, Kokkos::ALL, iC12);
  auto C13 = Kokkos::subview(local_parameters, Kokkos::ALL, iC13);
  auto C14 = Kokkos::subview(local_parameters, Kokkos::ALL, iC14);
  auto C_wp2_pr_dfsn = Kokkos::subview(local_parameters, Kokkos::ALL, iC_wp2_pr_dfsn);
  auto C_wp3_pr_tp = Kokkos::subview(local_parameters, Kokkos::ALL, iC_wp3_pr_tp);
  auto C_wp3_pr_turb = Kokkos::subview(local_parameters, Kokkos::ALL, iC_wp3_pr_turb);
  auto C_wp3_pr_dfsn = Kokkos::subview(local_parameters, Kokkos::ALL, iC_wp3_pr_dfsn);
  auto C_wp2_splat = Kokkos::subview(local_parameters, Kokkos::ALL, iC_wp2_splat);
  auto C6rt_Lscale0 = Kokkos::subview(local_parameters, Kokkos::ALL, iC6rt_Lscale0);
  auto C6thl_Lscale0 = Kokkos::subview(local_parameters, Kokkos::ALL, iC6thl_Lscale0);
  auto C7_Lscale0 = Kokkos::subview(local_parameters, Kokkos::ALL, iC7_Lscale0);
  auto wpxp_L_thresh = Kokkos::subview(local_parameters, Kokkos::ALL, iwpxp_L_thresh);
  auto c_K = Kokkos::subview(local_parameters, Kokkos::ALL, ic_K);
  auto c_K1 = Kokkos::subview(local_parameters, Kokkos::ALL, ic_K1);
  auto nu1 = Kokkos::subview(local_parameters, Kokkos::ALL, inu1);
  auto c_K2 = Kokkos::subview(local_parameters, Kokkos::ALL, ic_K2);
  auto nu2 = Kokkos::subview(local_parameters, Kokkos::ALL, inu2);
  auto c_K6 = Kokkos::subview(local_parameters, Kokkos::ALL, ic_K6);
  auto nu6 = Kokkos::subview(local_parameters, Kokkos::ALL, inu6);
  auto c_K8 = Kokkos::subview(local_parameters, Kokkos::ALL, ic_K8);
  auto nu8 = Kokkos::subview(local_parameters, Kokkos::ALL, inu8);
  auto c_K9 = Kokkos::subview(local_parameters, Kokkos::ALL, ic_K9);
  auto nu9 = Kokkos::subview(local_parameters, Kokkos::ALL, inu9);
  auto nu10 = Kokkos::subview(local_parameters, Kokkos::ALL, inu10);
  auto c_K_hm = Kokkos::subview(local_parameters, Kokkos::ALL, ic_K_hm);
  auto c_K_hmb = Kokkos::subview(local_parameters, Kokkos::ALL, ic_K_hmb);
  auto K_hm_min_coef = Kokkos::subview(local_parameters, Kokkos::ALL, iK_hm_min_coef);
  auto nu_hm = Kokkos::subview(local_parameters, Kokkos::ALL, inu_hm);
  auto slope_coef_spread_DG_means_w =
      Kokkos::subview(local_parameters, Kokkos::ALL, islope_coef_spread_DG_means_w);
  auto pdf_component_stdev_factor_w =
      Kokkos::subview(local_parameters, Kokkos::ALL, ipdf_component_stdev_factor_w);
  auto coef_spread_DG_means_rt =
      Kokkos::subview(local_parameters, Kokkos::ALL, icoef_spread_DG_means_rt);
  auto coef_spread_DG_means_thl =
      Kokkos::subview(local_parameters, Kokkos::ALL, icoef_spread_DG_means_thl);
  auto gamma_coef = Kokkos::subview(local_parameters, Kokkos::ALL, igamma_coef);
  auto gamma_coefb = Kokkos::subview(local_parameters, Kokkos::ALL, igamma_coefb);
  auto gamma_coefc = Kokkos::subview(local_parameters, Kokkos::ALL, igamma_coefc);
  auto Lscale_mu_entrain = Kokkos::subview(local_parameters, Kokkos::ALL, iLscale_mu_entrain);
  auto beta = Kokkos::subview(local_parameters, Kokkos::ALL, ibeta);
  auto lmin_coef = Kokkos::subview(local_parameters, Kokkos::ALL, ilmin_coef);
  auto omicron = Kokkos::subview(local_parameters, Kokkos::ALL, iomicron);
  auto zeta_vrnce_rat = Kokkos::subview(local_parameters, Kokkos::ALL, izeta_vrnce_rat);
  auto upsilon_precip_frac_rat =
      Kokkos::subview(local_parameters, Kokkos::ALL, iupsilon_precip_frac_rat);
  auto lambda0_stability_coef =
      Kokkos::subview(local_parameters, Kokkos::ALL, ilambda0_stability_coef);
  auto mult_coef = Kokkos::subview(local_parameters, Kokkos::ALL, imult_coef);
  auto taumin = Kokkos::subview(local_parameters, Kokkos::ALL, itaumin);
  auto taumax = Kokkos::subview(local_parameters, Kokkos::ALL, itaumax);
  auto Lscale_mu_coef = Kokkos::subview(local_parameters, Kokkos::ALL, iLscale_mu_coef);
  auto Lscale_pert_coef = Kokkos::subview(local_parameters, Kokkos::ALL, iLscale_pert_coef);
  auto alpha_corr = Kokkos::subview(local_parameters, Kokkos::ALL, ialpha_corr);
  auto Skw_denom_coef = Kokkos::subview(local_parameters, Kokkos::ALL, iSkw_denom_coef);
  auto c_K10 = Kokkos::subview(local_parameters, Kokkos::ALL, ic_K10);
  auto c_K10h = Kokkos::subview(local_parameters, Kokkos::ALL, ic_K10h);
  auto thlp2_rad_coef = Kokkos::subview(local_parameters, Kokkos::ALL, ithlp2_rad_coef);
  auto thlp2_rad_cloud_frac_thresh =
      Kokkos::subview(local_parameters, Kokkos::ALL, ithlp2_rad_cloud_frac_thresh);
  auto up2_sfc_coef = Kokkos::subview(local_parameters, Kokkos::ALL, iup2_sfc_coef);
  auto Skw_max_mag = Kokkos::subview(local_parameters, Kokkos::ALL, iSkw_max_mag);
  auto xp3_coef_base = Kokkos::subview(local_parameters, Kokkos::ALL, ixp3_coef_base);
  auto xp3_coef_slope = Kokkos::subview(local_parameters, Kokkos::ALL, ixp3_coef_slope);
  auto altitude_threshold = Kokkos::subview(local_parameters, Kokkos::ALL, ialtitude_threshold);
  auto rtp2_clip_coef = Kokkos::subview(local_parameters, Kokkos::ALL, irtp2_clip_coef);
  auto C_invrs_tau_bkgnd = Kokkos::subview(local_parameters, Kokkos::ALL, iC_invrs_tau_bkgnd);
  auto C_invrs_tau_sfc = Kokkos::subview(local_parameters, Kokkos::ALL, iC_invrs_tau_sfc);
  auto C_invrs_tau_shear = Kokkos::subview(local_parameters, Kokkos::ALL, iC_invrs_tau_shear);
  auto C_invrs_tau_N2 = Kokkos::subview(local_parameters, Kokkos::ALL, iC_invrs_tau_N2);
  auto C_invrs_tau_N2_wp2 = Kokkos::subview(local_parameters, Kokkos::ALL, iC_invrs_tau_N2_wp2);
  auto C_invrs_tau_N2_xp2 = Kokkos::subview(local_parameters, Kokkos::ALL, iC_invrs_tau_N2_xp2);
  auto C_invrs_tau_N2_wpxp = Kokkos::subview(local_parameters, Kokkos::ALL, iC_invrs_tau_N2_wpxp);
  auto C_invrs_tau_N2_clear_wp3 =
      Kokkos::subview(local_parameters, Kokkos::ALL, iC_invrs_tau_N2_clear_wp3);
  auto C_invrs_tau_wpxp_Ri = Kokkos::subview(local_parameters, Kokkos::ALL, iC_invrs_tau_wpxp_Ri);
  auto C_invrs_tau_wpxp_N2_thresh =
      Kokkos::subview(local_parameters, Kokkos::ALL, iC_invrs_tau_wpxp_N2_thresh);
  auto Cx_min = Kokkos::subview(local_parameters, Kokkos::ALL, iCx_min);
  auto Cx_max = Kokkos::subview(local_parameters, Kokkos::ALL, iCx_max);
  auto Richardson_num_min = Kokkos::subview(local_parameters, Kokkos::ALL, iRichardson_num_min);
  auto Richardson_num_max = Kokkos::subview(local_parameters, Kokkos::ALL, iRichardson_num_max);
  auto wpxp_Ri_exp = Kokkos::subview(local_parameters, Kokkos::ALL, iwpxp_Ri_exp);
  auto a3_coef_min = Kokkos::subview(local_parameters, Kokkos::ALL, ia3_coef_min);
  auto a_const = Kokkos::subview(local_parameters, Kokkos::ALL, ia_const);
  auto bv_efold = Kokkos::subview(local_parameters, Kokkos::ALL, ibv_efold);
  auto z_displace = Kokkos::subview(local_parameters, Kokkos::ALL, iz_displace);
  auto Lscale_tke_multiplier =
      Kokkos::subview(local_parameters, Kokkos::ALL, iLscale_tke_multiplier);
  auto Lscale_pdf_launch_frac =
      Kokkos::subview(local_parameters, Kokkos::ALL, iLscale_pdf_launch_frac);
  auto Lscale_up_weight = Kokkos::subview(local_parameters, Kokkos::ALL, iLscale_up_weight);

  // Set the default tunable parameter values.
  set_default_parameters(ngrdcol,                                                    // In
                         C1, C1b, C1c, C2rt, C2thl, C2rtthl,                         // Out
                         C4, C_uu_shr, C_uu_buoy, C6rt, C6rtb, C6rtc,                // Out
                         C6thl, C6thlb, C6thlc, C7, C7b, C7c, C8, C8b, C10,          // Out
                         C11, C11b, C11c, C12, C13, C14, C_wp2_pr_dfsn, C_wp3_pr_tp, // Out
                         C_wp3_pr_turb, C_wp3_pr_dfsn, C_wp2_splat,                  // Out
                         C6rt_Lscale0, C6thl_Lscale0, C7_Lscale0, wpxp_L_thresh,     // Out
                         c_K, c_K1, nu1, c_K2, nu2, c_K6, nu6, c_K8, nu8,            // Out
                         c_K9, nu9, nu10, c_K_hm, c_K_hmb, K_hm_min_coef, nu_hm,     // Out
                         slope_coef_spread_DG_means_w, pdf_component_stdev_factor_w, // Out
                         coef_spread_DG_means_rt, coef_spread_DG_means_thl,          // Out
                         gamma_coef, gamma_coefb, gamma_coefc, Lscale_mu_entrain, beta,
                         lmin_coef,                                                      // Out
                         omicron, zeta_vrnce_rat, upsilon_precip_frac_rat,               // Out
                         lambda0_stability_coef, mult_coef, taumin, taumax,              // Out
                         Lscale_mu_coef, Lscale_pert_coef, alpha_corr,                   // Out
                         Skw_denom_coef, c_K10, c_K10h, thlp2_rad_coef,                  // Out
                         thlp2_rad_cloud_frac_thresh, up2_sfc_coef,                      // Out
                         Skw_max_mag, xp3_coef_base, xp3_coef_slope,                     // Out
                         altitude_threshold, rtp2_clip_coef, C_invrs_tau_bkgnd,          // Out
                         C_invrs_tau_sfc, C_invrs_tau_shear, C_invrs_tau_N2,             // Out
                         C_invrs_tau_N2_wp2, C_invrs_tau_N2_xp2,                         // Out
                         C_invrs_tau_N2_wpxp, C_invrs_tau_N2_clear_wp3,                  // Out
                         C_invrs_tau_wpxp_Ri, C_invrs_tau_wpxp_N2_thresh,                // Out
                         Cx_min, Cx_max, Richardson_num_min, Richardson_num_max,         // Out
                         wpxp_Ri_exp, a3_coef_min, a_const, bv_efold, z_displace,        // Out
                         Lscale_tke_multiplier, Lscale_pdf_launch_frac, Lscale_up_weight // Out
  );

  // If filename is empty, use the statically specified working parameter set.
  // Fortran pads character comparisons: an all-blank filename is also empty.
  if (filename.find_first_not_of(' ') != std::string::npos) {
    // The existing descriptor/namelist owner reads initialized arrays only.
    // Its call-local NEWUNIT replaces OPEN(unit=iunit), avoiding unit collisions;
    // iunit remains required at this exact transient I/O boundary.
    int read_status = 0;
    f2k::read_clubb_params_nl(
        filename, ngrdcol, iunit,                               // In
        C1, C1b, C1c, C2rt,                                     // InOut
        C2thl, C2rtthl, C4, C_uu_shr,                           // InOut
        C_uu_buoy, C6rt, C6rtb, C6rtc,                          // InOut
        C6thl, C6thlb, C6thlc, C7,                              // InOut
        C7b, C7c, C8, C8b,                                      // InOut
        C10, C11, C11b, C11c,                                   // InOut
        C12, C13, C14, C_wp2_pr_dfsn,                           // InOut
        C_wp3_pr_tp, C_wp3_pr_turb, C_wp3_pr_dfsn, C_wp2_splat, // InOut
        C6rt_Lscale0, C6thl_Lscale0, C7_Lscale0, wpxp_L_thresh, // InOut
        c_K, c_K1, nu1, c_K2,                                   // InOut
        nu2, c_K6, nu6, c_K8,                                   // InOut
        nu8, c_K9, nu9, nu10,                                   // InOut
        c_K_hm, c_K_hmb, K_hm_min_coef, nu_hm,                  // InOut
        slope_coef_spread_DG_means_w, pdf_component_stdev_factor_w, coef_spread_DG_means_rt,
        coef_spread_DG_means_thl,                                               // InOut
        beta, gamma_coef, gamma_coefb, gamma_coefc,                             // InOut
        lmin_coef, omicron, zeta_vrnce_rat, upsilon_precip_frac_rat,            // InOut
        lambda0_stability_coef, mult_coef, taumin, taumax,                      // InOut
        Lscale_mu_entrain, Lscale_mu_coef, Lscale_pert_coef, alpha_corr,        // InOut
        Skw_denom_coef, c_K10, c_K10h, thlp2_rad_coef,                          // InOut
        thlp2_rad_cloud_frac_thresh, up2_sfc_coef, Skw_max_mag, xp3_coef_base,  // InOut
        xp3_coef_slope, altitude_threshold, rtp2_clip_coef, C_invrs_tau_bkgnd,  // InOut
        C_invrs_tau_sfc, C_invrs_tau_shear, C_invrs_tau_N2, C_invrs_tau_N2_wp2, // InOut
        C_invrs_tau_N2_xp2, C_invrs_tau_N2_wpxp, C_invrs_tau_N2_clear_wp3,
        C_invrs_tau_wpxp_Ri,                                                 // InOut
        C_invrs_tau_wpxp_N2_thresh, Cx_min, Cx_max, Richardson_num_min,      // InOut
        Richardson_num_max, wpxp_Ri_exp, a3_coef_min, a_const,               // InOut
        bv_efold, z_displace, Lscale_tke_multiplier, Lscale_pdf_launch_frac, // InOut
        Lscale_up_weight,                                                    // InOut
        read_status);                                                        // Out
    if (read_status != 0)
      throw std::runtime_error("Error reading clubb_params_nl from " + filename + ": " +
                               std::to_string(read_status));
  }

  // Put the variables in the output array.
  Host2D packed_parameters(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "init_clubb_params_api::packed_parameters"),
      ngrdcol, nparams);
  parameters_tunable_detail::pack_parameters(
      ngrdcol,                                                                           // In
      C1, C1b, C1c, C2rt, C2thl, C2rtthl,                                                // In
      C4, C_uu_shr, C_uu_buoy, C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc,                // In
      C7, C7b, C7c, C8, C8b, C10,                                                        // In
      C11, C11b, C11c, C12, C13, C14, C_wp2_pr_dfsn, C_wp3_pr_tp,                        // In
      C_wp3_pr_turb, C_wp3_pr_dfsn, C_wp2_splat,                                         // In
      C6rt_Lscale0, C6thl_Lscale0, C7_Lscale0, wpxp_L_thresh,                            // In
      c_K, c_K1, nu1, c_K2, nu2, c_K6, nu6, c_K8, nu8,                                   // In
      c_K9, nu9, nu10, c_K_hm, c_K_hmb, K_hm_min_coef, nu_hm,                            // In
      slope_coef_spread_DG_means_w, pdf_component_stdev_factor_w,                        // In
      coef_spread_DG_means_rt, coef_spread_DG_means_thl,                                 // In
      gamma_coef, gamma_coefb, gamma_coefc, Lscale_mu_entrain, beta, lmin_coef,          // In
      omicron, zeta_vrnce_rat, upsilon_precip_frac_rat,                                  // In
      lambda0_stability_coef, mult_coef, taumin, taumax,                                 // In
      Lscale_mu_coef, Lscale_pert_coef, alpha_corr,                                      // In
      Skw_denom_coef, c_K10, c_K10h, thlp2_rad_coef,                                     // In
      thlp2_rad_cloud_frac_thresh, up2_sfc_coef,                                         // In
      Skw_max_mag, xp3_coef_base, xp3_coef_slope,                                        // In
      altitude_threshold, rtp2_clip_coef, C_invrs_tau_bkgnd,                             // In
      C_invrs_tau_sfc, C_invrs_tau_shear, C_invrs_tau_N2,                                // In
      C_invrs_tau_N2_wp2, C_invrs_tau_N2_xp2,                                            // In
      C_invrs_tau_N2_wpxp, C_invrs_tau_N2_clear_wp3,                                     // In
      C_invrs_tau_wpxp_Ri, C_invrs_tau_wpxp_N2_thresh,                                   // In
      Cx_min, Cx_max, Richardson_num_min, Richardson_num_max,                            // In
      wpxp_Ri_exp, a3_coef_min, a_const, bv_efold, z_displace,                           // In
      Lscale_tke_multiplier, Lscale_pdf_launch_frac, Lscale_up_weight, packed_parameters // Out
  );
  from_execution_space(clubb_params, packed_parameters);

  // Native's missing-parameter check remains disabled: defaults are not
  // initialized to init_value (-999), so that diagnostic is currently broken.
}
