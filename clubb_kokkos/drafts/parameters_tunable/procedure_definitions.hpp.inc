// Definitions follow the Native supported procedure order. Merge after constants/private
// declarations.

template <class Column>
void set_default_parameters(
    int ngrdcol,                                                                           // In
    Column C1, Column C1b, Column C1c, Column C2rt, Column C2thl, Column C2rtthl,          // Out
    Column C4, Column C_uu_shr, Column C_uu_buoy, Column C6rt, Column C6rtb, Column C6rtc, // Out
    Column C6thl, Column C6thlb, Column C6thlc, Column C7, Column C7b, Column C7c, Column C8,
    Column C8b, Column C10, // Out
    Column C11, Column C11b, Column C11c, Column C12, Column C13, Column C14, Column C_wp2_pr_dfsn,
    Column C_wp3_pr_tp,                                                                 // Out
    Column C_wp3_pr_turb, Column C_wp3_pr_dfsn, Column C_wp2_splat,                     // Out
    Column C6rt_Lscale0, Column C6thl_Lscale0, Column C7_Lscale0, Column wpxp_L_thresh, // Out
    Column c_K, Column c_K1, Column nu1, Column c_K2, Column nu2, Column c_K6, Column nu6,
    Column c_K8, Column nu8, // Out
    Column c_K9, Column nu9, Column nu10, Column c_K_hm, Column c_K_hmb, Column K_hm_min_coef,
    Column nu_hm,                                                             // Out
    Column slope_coef_spread_DG_means_w, Column pdf_component_stdev_factor_w, // Out
    Column coef_spread_DG_means_rt, Column coef_spread_DG_means_thl,          // Out
    Column gamma_coef, Column gamma_coefb, Column gamma_coefc, Column Lscale_mu_entrain,
    Column beta, Column lmin_coef,                                                      // Out
    Column omicron, Column zeta_vrnce_rat, Column upsilon_precip_frac_rat,              // Out
    Column lambda0_stability_coef, Column mult_coef, Column taumin, Column taumax,      // Out
    Column Lscale_mu_coef, Column Lscale_pert_coef, Column alpha_corr,                  // Out
    Column Skw_denom_coef, Column c_K10, Column c_K10h, Column thlp2_rad_coef,          // Out
    Column thlp2_rad_cloud_frac_thresh, Column up2_sfc_coef,                            // Out
    Column Skw_max_mag, Column xp3_coef_base, Column xp3_coef_slope,                    // Out
    Column altitude_threshold, Column rtp2_clip_coef, Column C_invrs_tau_bkgnd,         // Out
    Column C_invrs_tau_sfc, Column C_invrs_tau_shear, Column C_invrs_tau_N2,            // Out
    Column C_invrs_tau_N2_wp2, Column C_invrs_tau_N2_xp2,                               // Out
    Column C_invrs_tau_N2_wpxp, Column C_invrs_tau_N2_clear_wp3,                        // Out
    Column C_invrs_tau_wpxp_Ri, Column C_invrs_tau_wpxp_N2_thresh,                      // Out
    Column Cx_min, Column Cx_max, Column Richardson_num_min, Column Richardson_num_max, // Out
    Column wpxp_Ri_exp, Column a3_coef_min, Column a_const, Column bv_efold,
    Column z_displace,                                                                   // Out
    Column Lscale_tke_multiplier, Column Lscale_pdf_launch_frac, Column Lscale_up_weight // Out
) {
  // ngrdcol: Number of columns [#].
  // All remaining outputs are column profiles of Native tunable coefficients;
  // their physical meanings and units accompany the assignments below.
  // Homogeneous column Views select their own execution space, so the
  // call-local HostSpace namelist arrays need no device round trip here.
  using ExecutionSpace = typename Column::execution_space;
  Kokkos::parallel_for(
      "set_default_parameters", Kokkos::RangePolicy<ExecutionSpace>(0, ngrdcol),
      KOKKOS_LAMBDA(const int i) {
        // NOTE: In CLUBB standalone, as well as some host models, the hardcoded
        //       default values of some or all of the parameters below have no
        //       effect, as the values are simply read in using a namelist or set in
        //       host model specific code.

        // Model tunable parameters
        // WARNING: THE DEFAULT VALUES OF THE PARAMETERS BELOW MAY BE OVERWRITTEN
        //    BY NAMELIST VALUES FROM, E.G., clubb_params_nl!!!
        C1(i) = real_t(1.000000);            // Low Skewness in C1 Skw. Function    [-]
        C1b(i) = real_t(1.000000);           // High Skewness in C1 Skw. Function   [-]
        C1c(i) = real_t(1.000000);           // Degree of Slope of C1 Skw. Function [-]
        C2rt(i) = real_t(2.000000);          // C2 coef. for the rtp2_dp1 term      [-]
        C2thl(i) = real_t(2.000000);         // C2 coef. for the thlp2_dp1 term     [-]
        C2rtthl(i) = real_t(2.000000);       // C2 coef. for the rtpthlp_dp1 term   [-]
        C4(i) = real_t(2.000000);            // Pressure term coefficient for w'^2  [-]
        C_uu_shr(i) = real_t(0.400000);      // Coef. in pressure terms (shear): w'^2 eqn   [-]
        C_uu_buoy(i) = real_t(0.300000);     // Coef. in pressure terms (buoyancy): w'^2 eqn [-]
        C6rt(i) = real_t(2.000000);          // Low Skewness in C6rt Skw. Function  [-]
        C6rtb(i) = real_t(2.000000);         // High Skewness in C6rt Skw. Function [-]
        C6rtc(i) = real_t(1.000000);         // Degree of Slope of C6rt Skw. Fnct.  [-]
        C6thl(i) = real_t(2.000000);         // Low Skewness in C6thl Skw. Function [-]
        C6thlb(i) = real_t(2.000000);        // High Skewness in C6thl Skw. Fnct.   [-]
        C6thlc(i) = real_t(1.000000);        // Degree of Slope of C6thl Skw. Fnct. [-]
        C7(i) = real_t(0.500000);            // Low Skewness in C7 Skw. Function    [-]
        C7b(i) = real_t(0.500000);           // High Skewness in C7 Skw. Function   [-]
        C7c(i) = real_t(0.500000);           // Degree of Slope of C7 Skw. Function [-]
        C8(i) = real_t(0.500000);            // Coef. #1 in C8 Skewness Equation    [-]
        C8b(i) = real_t(0.020000);           // Coef. #2 in C8 Skewness Equation    [-]
        C10(i) = real_t(3.300000);           // Currently Not Used in the Model     [-]
        C11(i) = real_t(0.400000);           // Low Skewness in C11 Skw. Function   [-]
        C11b(i) = real_t(0.400000);          // High Skewness in C11 Skw. Function  [-]
        C11c(i) = real_t(0.500000);          // Degree of Slope of C11 Skw. Fnct.   [-]
        C12(i) = real_t(1.000000);           // Constant in w'^3 Crank-Nich. diff.  [-]
        C13(i) = real_t(0.100000);           // Not currently used in model         [-]
        C14(i) = real_t(1.000000);           // Constant for u'^2 and v'^2 terms    [-]
        C_wp2_pr_dfsn(i) = real_t(0.000000); // Coefficient for the wp2_pr_dfsn term [-]
        C_wp3_pr_tp(i) = real_t(0.000000);   // Coefficient for the wp3_pr_tp term [-]
        C_wp3_pr_turb(i) = real_t(0.000000); // Coefficient for the wp3_pr_turb term [-]
        C_wp3_pr_dfsn(i) = real_t(0.000000); // Coefficient for the wp3_pr_dfsn term [-]
        C_wp2_splat(i) = real_t(2.000000);   // Coefficient for gustiness near ground [-]
        C6rt_Lscale0(i) = real_t(14.0);      // Damp C6rt as a fnct. of Lscale  [-]
        C6thl_Lscale0(i) = real_t(14.0);     // Damp C6thl as a fnct. of Lscale [-]
        C7_Lscale0(i) = real_t(0.8500000);   // Damp C7 as a fnct. of Lscale    [-]
        wpxp_L_thresh(i) = real_t(60.0);     // Lscale threshold: damp C6 & C7  [m]
        c_K(i) = real_t(0.200000);           // Constant C_mu^(1/4) in DD 1987 [m^2/s]
        c_K1(i) = real_t(0.200000);          // Coef. of Eddy Diffusion: wp2   [m^2/s]
        c_K2(i) = real_t(0.025000);          // Coef. of Eddy Diffusion: xp2   [m^2/s]
        c_K6(i) = real_t(0.375000);          // Coef. of Eddy Diffusion: wpxp  [m^2/s]
        c_K8(i) = real_t(5.000000);          // Coef. of Eddy Diffusion: wp3   [m^2/s]
        c_K9(i) = real_t(0.100000);          // Coef. of Eddy Diff.: up2/vp2   [m^2/s]
        c_K_hm(i) = real_t(0.750000);        // Coef. of Eddy Diffusion: hmm   [m^2/s]
        c_K_hmb(i) = real_t(0.750000); // Coef. of Non-Local Factor, Eddy Diffusion: hmm   [m^2/s]
        K_hm_min_coef(i) = real_t(0.10000); // Min. of Non-Local Factor, Eddy Diffusion: hmm [m^2/s]
        gamma_coef(i) = real_t(0.250000);   // Low Skw.: gamma coef. Skw. Fnct.   [-]
        gamma_coefb(i) = real_t(0.250000);  // High Skw.: gamma coef. Skw. Fnct.  [-]
        gamma_coefc(i) = real_t(5.000000);  // Deg. Slope: gamma coef. Skw. Fnct. [-]
        Lscale_mu_entrain(i) =
            real_t(1.000E-3);            // Positive fractional parcel entrainment rate [1/m]
        mult_coef(i) = real_t(0.500000); // Coef. applied to log(avg dz/thresh)[-]
        taumin(i) = real_t(90.00000);    // Min. allow. value: time-scale tau  [s]
        taumax(i) = real_t(3600.000);    // Max. allow. value: time-scale tau  [s]
        Lscale_mu_coef(i) =
            real_t(2.0); // Entrainment perturbation factor for averaged method 1 [-]
        Lscale_pert_coef(i) = real_t(0.1); // Coef pert thlm/rtm: av calc Lscale [-]
        alpha_corr(i) = real_t(0.15);      // Coef. for the corr. diagnosis algorithm  [-]
        nu1(i) = real_t(20.00000);         // Bg. Coef. Eddy Diffusion: wp2        [m^2/s]
        nu2(i) = real_t(1.000000);         // Bg. Coef. Eddy Diffusion: xp2        [m^2/s]
        nu6(i) = real_t(5.000000);         // Bg. Coef. Eddy Diffusion: wpxp       [m^2/s]
        nu8(i) = real_t(20.00000);         // Bg. Coef. Eddy Diffusion: wp3        [m^2/s]
        nu9(i) = real_t(10.00000);         // Bg. Coef. Eddy Diffusion: up2/vp2    [m^2/s]
        nu10(i) = real_t(0.000000);        // Bg. Coef. Eddy Diffusion: edsclrm    [m^2/s]
        nu_hm(i) = real_t(1.500000);       // Bg. Coef. Eddy Diffusion: hmm        [m^2/s]

        // Vince Larson added a constant to set plume widths for theta_l and rt
        // beta should vary between 0 and 3.
        beta(i) = real_t(1.000000);      // Beta coefficient     [-]
        lmin_coef(i) = real_t(0.500000); // Coefficient of lmin    [-]
        Skw_max_mag(i) = real_t(10.0);   // Max magnitude of skewness [-]

        C_invrs_tau_bkgnd(i) = real_t(1.1);
        C_invrs_tau_sfc(i) = real_t(0.1);
        C_invrs_tau_shear(i) = real_t(0.15);
        C_invrs_tau_N2(i) = real_t(0.4);
        C_invrs_tau_N2_wp2(i) = real_t(0.2);
        C_invrs_tau_N2_xp2(i) = real_t(0.05);
        C_invrs_tau_N2_wpxp(i) = real_t(0.0);
        C_invrs_tau_N2_clear_wp3(i) = real_t(1.0);
        C_invrs_tau_wpxp_Ri(i) = real_t(0.35);
        C_invrs_tau_wpxp_N2_thresh(i) = real_t(3.3e-4);

        // Parameters for the new PDF (w, rt, and theta-l).
        //
        // Brian Griffin added a tunable parameter for the PDF of w,
        // slope_coef_spread_DG_means_w, to increase or decrease the spread between
        // the two PDF component means of w.  When the value of this slope parameter
        // is larger, F_w is smaller and the PDF component means of w are closer
        // together.
        // Valid values are slope_coef_spread_DG_means_w > 0.
        //
        // A second parameter for the PDF of w, pdf_component_stdev_factor_w, is used
        // to adjust the standard deviations of the 1st PDF component against the 2nd
        // PDF component for w.  This parameter is related to zeta_w, where:
        //
        // 1 + zeta_w = ( mixt_frac * sigma_w_1^2 )
        //              / ( ( 1 - mixt_frac ) * sigma_w_2^2 );
        //
        // The pdf_component_stdev_factor_w is set such that:
        //
        // pdf_component_stdev_factor_w = zeta_w + 1.
        //
        // Valid values are pdf_component_stdev_factor_w > 0.
        //
        // The parameter for the PDF of rt is coef_spread_DG_means_rt.  Valid values
        // are 0 <= coef_spread_DG_means_rt < 1.  When coef_spread_DG_means_rt
        // approaches 0, F_rt approaches min_F_rt, and the two PDF component means
        // become closer together.  When coef_spread_DG_means_rt approaches 1, F_rt
        // approaches max_F_rt, and the two PDF component means are spread farther
        // apart.
        //
        // The parameter for the PDF of theta-l is coef_spread_DG_means_thl.
        // Valid values are 0 <= coef_spread_DG_means_thl < 1.  When
        // coef_spread_DG_means_thl approaches 0, F_thl approaches min_F_thl, and the
        // two PDF component means become closer together.  When
        // coef_spread_DG_means_thl approaches 1, F_thl approaches max_F_thl, and the
        // two PDF component means are spread farther apart.
        // Slope coefficient for the spread between the PDF component means of w.
        slope_coef_spread_DG_means_w(i) = real_t(21.0);
        // Parameter to adjust the PDF component standard deviations of w.
        pdf_component_stdev_factor_w(i) = real_t(1.0);
        // Coefficient for the spread between the PDF component means of rt.
        coef_spread_DG_means_rt(i) = real_t(0.8);
        // Coefficient for the spread between the PDF component means of thl.
        coef_spread_DG_means_thl(i) = real_t(0.8);

        // Parameters for the hydrometeor portion of the PDF.
        //
        // Brian Griffin added a parameter for hydrometeors, omicron, to increase the
        // standard deviation of each component and decrease the spread between the
        // component means as the value of omicron inreases.  Valid value are
        // 0 < omicron <= 1.
        // A second parameter for hydrometeors, zeta, increases the standard
        // deviation of component 1 at the expense of the standard deviation of
        // component 2 when the value of zeta > 0 (and increasingly so as zeta
        // increases).  Valid values are zeta > -1.
        omicron(i) = real_t(0.5);        // Hydromet width/spread-of-means param [-]
        zeta_vrnce_rat(i) = real_t(0.0); // Ratio sigma^2/mu^2 comp. 1 / comp. 2 [-]
        // ratio mixt_frac*precip_frac_1/precip_frac (precip_frac_calc_type=2)    [-]
        upsilon_precip_frac_rat(i) = real_t(0.55);

        // Intensity of stability correction applied to C1 and C6 [-]
        lambda0_stability_coef(i) = real_t(0.03);
        // Factor to decrease sensitivity in the denominator of Skw calculation
        Skw_denom_coef(i) = real_t(4.0);

        // Momentum coefficient of Kh_zm
        c_K10(i) = real_t(1.0);
        // Thermodynamic coefficient of Kh_zm
        c_K10h(i) = real_t(1.0);

        thlp2_rad_coef(i) = real_t(1.0);              // Coefficient of thlp2_rad               [-]
        thlp2_rad_cloud_frac_thresh(i) = real_t(0.1); // Minimum cloud fraction for
        // computation of thlp2_rad  [-]

        up2_sfc_coef(i) = real_t(4.0); // Coefficients of up2 and vp2    [-]

        xp3_coef_base(i) = real_t(0.25);  // "Base" value of xp3_coef in simple eqn
        xp3_coef_slope(i) = real_t(0.01); // Slope in regards to Brunt-Vaisla freq.

        altitude_threshold(i) = real_t(100.0); // Altitude above which damping should occur for wpxp
        rtp2_clip_coef(i) = real_t(0.5);       // Coef. appled the clipping threshold on rtp2

        Cx_min(i) = real_t(0.33);              // Threshold on Cx_fnc_Richardson
        Cx_max(i) = real_t(0.95);              // Threshold on Cx_fnc_Richardson
        Richardson_num_min(i) = real_t(0.25);  // Threshold on Richardson number
        Richardson_num_max(i) = real_t(400.0); // Threshold on Richardson number

        a3_coef_min(i) = real_t(1.0); // Minimum threshold on the a3 coefficient
        a_const(i) = real_t(1.8);     // Used in sfc_varnce_module.F90 to set surface values
        bv_efold(i) = real_t(5.);     // Control parameter for inverse e-folding of
        // cloud fraction in the mixed Brunt Vaisala frequency
        wpxp_Ri_exp(i) = real_t(.5); // Exponent determining the influence of
        // the Richardson number on invrs_tau_wpxp
        z_displace(i) = real_t(25.0);   // displacement of log law profile above ground   [m]
        Lscale_tke_multiplier(i) = one; // Mean launch energy/TKE in method 3 [-]
        Lscale_pdf_launch_frac(i) =
            zero; // Fraction toward directional PDF-component means in method 3 [-]
        Lscale_up_weight(i) = one_half; // Upward exponent in method-3 geometric mixing [-]
      });
}

template <class Matrix>
void check_parameters_api(int ngrdcol, Matrix clubb_params, err_info_type& err_info // InOut
) {

  // Description:
  //   Compatibility parameter-only validation.  New callers should use
  //   check_clubb_settings_api, which invokes this routine before checking
  //   flags and mixed parameter/flag consistency rules.

  // References:
  // None
  //-----------------------------------------------------------------------
  // ngrdcol: Number of columns [#].
  // clubb_params: Column-first tunable parameter matrix [-].
  // err_info: Existing column error codes and diagnostic headers [-].
  parameters_tunable_detail::check_parameter_hard_bounds(ngrdcol, clubb_params, err_info);
  const auto parameters = to_host_space(clubb_params);
  const auto err_header = to_host_space(err_info.err_header);
  auto err_code = to_host_space(err_info.err_code);
  std::ostringstream diagnostic;
  diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
  for (int i = 0; i < ngrdcol; ++i) {
    // The C6rt parameters must equal the C6thl parameters. Otherwise,
    // the wpthlp pr1 term is calculated inconsistently.
    if (std::abs(parameters(i, iC6rt) - parameters(i, iC6thl)) >
        std::abs(parameters(i, iC6rt) + parameters(i, iC6thl)) / 2 * eps) {
      diagnostic.write(err_header.data() + i * err_header.extent(0),
                       static_cast<std::streamsize>(err_header.extent(0)));
      diagnostic << '\n';
      diagnostic << "C6rt = " << parameters(i, iC6rt) << '\n';
      diagnostic << "C6thl = " << parameters(i, iC6thl) << '\n';
      diagnostic << "C6rt and C6thl must be equal.\n";
      // Error in column i: set only this entry to clubb_fatal_error.
      err_code(i) = clubb_fatal_error;
    }
    if (std::abs(parameters(i, iC6rtb) - parameters(i, iC6thlb)) >
        std::abs(parameters(i, iC6rtb) + parameters(i, iC6thlb)) / 2 * eps) {
      diagnostic.write(err_header.data() + i * err_header.extent(0),
                       static_cast<std::streamsize>(err_header.extent(0)));
      diagnostic << '\n';
      diagnostic << "C6rtb = " << parameters(i, iC6rtb) << '\n';
      diagnostic << "C6thlb = " << parameters(i, iC6thlb) << '\n';
      diagnostic << "C6rtb and C6thlb must be equal.\n";
      // Error in column i: set only this entry to clubb_fatal_error.
      err_code(i) = clubb_fatal_error;
    }
    if (std::abs(parameters(i, iC6rtc) - parameters(i, iC6thlc)) >
        std::abs(parameters(i, iC6rtc) + parameters(i, iC6thlc)) / 2 * eps) {
      diagnostic.write(err_header.data() + i * err_header.extent(0),
                       static_cast<std::streamsize>(err_header.extent(0)));
      diagnostic << '\n';
      diagnostic << "C6rtc = " << parameters(i, iC6rtc) << '\n';
      diagnostic << "C6thlc = " << parameters(i, iC6thlc) << '\n';
      diagnostic << "C6rtc and C6thlc must be equal.\n";
      // Error in column i: set only this entry to clubb_fatal_error.
      err_code(i) = clubb_fatal_error;
    }
    if (std::abs(parameters(i, iC6rt_Lscale0) - parameters(i, iC6thl_Lscale0)) >
        std::abs(parameters(i, iC6rt_Lscale0) + parameters(i, iC6thl_Lscale0)) / 2 * eps) {
      diagnostic.write(err_header.data() + i * err_header.extent(0),
                       static_cast<std::streamsize>(err_header.extent(0)));
      diagnostic << '\n';
      diagnostic << "C6rt_Lscale0 = " << parameters(i, iC6rt_Lscale0) << '\n';
      diagnostic << "C6thl_Lscale0 = " << parameters(i, iC6thl_Lscale0) << '\n';
      diagnostic << "C6rt_Lscale0 and C6thl_Lscale0 must be equal.\n";
      // Error in column i: set only this entry to clubb_fatal_error.
      err_code(i) = clubb_fatal_error;
    }
  }
  std::cerr << diagnostic.str();
  from_execution_space(err_info.err_code, err_code);
}

template <class Matrix>
void get_parameter_hard_bounds_api(int nparams_in, Matrix hard_bounds // Out
) {

  // Description:
  //   Return Fortran's authoritative 2 x nparams hard-bound table.  Row 1
  //   is minimum and row 2 is maximum; huge endpoints denote no bound.
  //   These are rejection limits, not recommended tuning intervals.
  //-----------------------------------------------------------------------
  // nparams_in: Caller-supplied parameter count [#].
  // hard_bounds: Parameter bounds, shape (2,nparams_in): minimum,maximum [-].
  if (nparams_in != nparams)
    throw std::runtime_error("get_parameter_hard_bounds_api: nparams mismatch");
  // Immutable host table is exposed only as a call-local const unmanaged View.
  using Bounds = Kokkos::View<const real_t**, Kokkos::LayoutLeft, Kokkos::HostSpace,
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  Bounds bounds(parameter_hard_bounds.data(), 2, nparams);
  from_execution_space(hard_bounds, bounds);
}

namespace parameters_tunable_detail {
template <class Matrix>
void check_parameter_hard_bounds(int ngrdcol, Matrix clubb_params, err_info_type& err_info // InOut
) {
  // Description:
  //   Reject parameters outside the immutable Native hard configuration bounds.
  // ngrdcol: Number of columns [#].
  // clubb_params: Tunable parameter matrix (column,parameter) [-].
  // err_info: Caller-owned column error state and diagnostic headers [-].
  const auto parameters = to_host_space(clubb_params);
  const auto err_header = to_host_space(err_info.err_header);
  auto err_code = to_host_space(err_info.err_code);
  std::ostringstream diagnostic;
  diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
  for (int i = 0; i < ngrdcol; ++i) {
    for (int k = 0; k < nparams; ++k) {
      const bool has_lower = (parameter_hard_bounds[2 * k] > no_lower_bound);
      const bool has_upper = (parameter_hard_bounds[2 * k + 1] < no_upper_bound);
      const bool below_lower = (has_lower && parameters(i, k) < parameter_hard_bounds[2 * k]);
      const bool above_upper = (has_upper && parameters(i, k) > parameter_hard_bounds[2 * k + 1]);
      if (below_lower || above_upper) {
        diagnostic.write(err_header.data() + i * err_header.extent(0),
                         static_cast<std::streamsize>(err_header.extent(0)));
        diagnostic << '\n';
        const auto name = params_list[k].substr(0, params_list[k].find_last_not_of(' ') + 1);
        diagnostic << name << " = " << parameters(i, k) << '\n';
        if (has_lower && has_upper) {
          diagnostic << name
                     << " is outside hard configuration bounds: " << parameter_hard_bounds[2 * k]
                     << ' ' << parameter_hard_bounds[2 * k + 1] << '\n';
        } else if (has_lower) {
          diagnostic << name
                     << " violates lower hard configuration bound: " << parameter_hard_bounds[2 * k]
                     << '\n';
        } else {
          diagnostic << name << " violates upper hard configuration bound: "
                     << parameter_hard_bounds[2 * k + 1] << '\n';
        }
        err_code(i) = clubb_fatal_error;
      }
    }
  }
  std::cerr << diagnostic.str();
  from_execution_space(err_info.err_code, err_code);
}
} // namespace parameters_tunable_detail

template <class Spacing, class Matrix, class MinimumLength, class MixtureLimit>
void calc_derived_params(const grid& gr, int ngrdcol, int grid_type, Spacing deltaz, // In
                         Matrix clubb_params,                                        // In
                         bool l_prescribed_avg_deltaz,                               // In
                         nu_vertical_res_dep& nu_vert_res_dep, MinimumLength lmin,   // InOut
                         MixtureLimit mixt_frac_max_mag                              // InOut
) {
  // Description:
  //   Calculate CLUBB parameters derived from other quantities.
  //   Adjust background eddy diffusivity for vertical grid spacing.
  //   Hosts with grid heights that vary in space or time (for example WRF)
  //   may call this repeatedly. --ldgrant Jul 2010
  // gr: Grid definition, with Native one-based physical boundary metadata [-].
  // ngrdcol: Number of columns [#].
  // grid_type: Uniform(1), stretched thermodynamic(2), or momentum-input(3) grid [-].
  // deltaz: Change per height level [m].
  // clubb_params: Tunable model parameters, shape (ngrdcol,nparams) [-].
  // l_prescribed_avg_deltaz: Use deltaz as the average spacing when true [-].
  // nu_vert_res_dep: Vertical-resolution-dependent background diffusivities [m^2/s].
  // lmin: Minimum length scale [m].
  // mixt_frac_max_mag: Maximum allowable magnitude of mixture fraction [-].

  constexpr real_t lmin_deltaz = real_t(40.0); // Fixed minimum-length reference spacing [m].
  // Turn off to retain the values declared in tunable_parameters.in.
  constexpr bool l_adj_low_res_nu = true;
  // At or below this average grid spacing the background diffusivities remain
  // unchanged. Above it, increase them by the factors below. Traditionally 40m.
  constexpr real_t grid_spacing_thresh = real_t(40.0); // Grid-spacing threshold [m].

  if (!nu_vert_res_dep.nu1.is_allocated()) {
    nu_vert_res_dep.nu1 = Real1D(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "nu_vertical_res_dep::nu1"), ngrdcol);
    nu_vert_res_dep.nu2 = Real1D(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "nu_vertical_res_dep::nu2"), ngrdcol);
    nu_vert_res_dep.nu6 = Real1D(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "nu_vertical_res_dep::nu6"), ngrdcol);
    nu_vert_res_dep.nu8 = Real1D(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "nu_vertical_res_dep::nu8"), ngrdcol);
    nu_vert_res_dep.nu9 = Real1D(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "nu_vertical_res_dep::nu9"), ngrdcol);
    nu_vert_res_dep.nu10 = Real1D(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "nu_vertical_res_dep::nu10"), ngrdcol);
    nu_vert_res_dep.nu_hm = Real1D(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "nu_vertical_res_dep::nu_hm"), ngrdcol);
  }

  const auto parameters = to_execution_space(clubb_params);
  const auto spacing = to_execution_space(deltaz);
  auto minimum_length = to_execution_space(lmin, false);
  auto mixture_limit = to_execution_space(mixt_frac_max_mag, false);

  // Experimentation favored a fixed lmin reference for grid-independent results.
  // -dschanen 21 May 2007. Old: lmin = lmin_coef * deltaz.
  int invalid_lmin = 0;
  Kokkos::parallel_reduce(
      "calc_derived_params::lmin", ngrdcol,
      KOKKOS_LAMBDA(const int i, int& invalid) {
        minimum_length(i) = parameters(i, ilmin_coef) * lmin_deltaz;
        if (minimum_length(i) < real_t(1.0))
          ++invalid;
      },
      invalid_lmin);
  // Preserve the Native write-before-fatal ordering for direct callers.
  from_execution_space(lmin, minimum_length);
  if (invalid_lmin != 0)
    throw std::runtime_error("lmin is < 1.0_core_rknd");

  // Capture component Views and physical index scalars, not the owning UDTs.
  const auto zm = gr.zm;
  const auto zt = gr.zt;
  const int k_ub_zm = gr.k_ub_zm - 1;
  const int k_lb_zm = gr.k_lb_zm - 1;
  const int k_ub_zt = gr.k_ub_zt - 1;
  const int k_lb_zt = gr.k_lb_zt - 1;
  const int nzm = gr.nzm;
  const int nzt = gr.nzt;
  const auto nu1 = nu_vert_res_dep.nu1;
  const auto nu2 = nu_vert_res_dep.nu2;
  const auto nu6 = nu_vert_res_dep.nu6;
  const auto nu8 = nu_vert_res_dep.nu8;
  const auto nu9 = nu_vert_res_dep.nu9;
  const auto nu10 = nu_vert_res_dep.nu10;
  const auto nu_hm = nu_vert_res_dep.nu_hm;
  Kokkos::parallel_for(
      "calc_derived_params", ngrdcol, KOKKOS_LAMBDA(const int i) {
        // Known magic number; preserve the Native parenthesization.
        mixture_limit(i) =
            real_t(1.0) -
            (real_t(0.5) *
             (real_t(1.0) -
              parameters(i, iSkw_max_mag) /
                  Kokkos::sqrt(real_t(4.0) *
                                   ((real_t(1.0) - real_t(0.4)) * (real_t(1.0) - real_t(0.4)) *
                                    (real_t(1.0) - real_t(0.4))) +
                               parameters(i, iSkw_max_mag) * parameters(i, iSkw_max_mag))));

        // Adjust constant background diffusivity coefficients for grid spacing;
        // leave them as declared in parameters.in if the compile flag is off.
        if (l_adj_low_res_nu) {
          real_t avg_deltaz;                     // Average grid-box height [m].
          real_t mult_factor_zt, mult_factor_zm; // Diffusivity multipliers [-].
          if (l_prescribed_avg_deltaz) {
            avg_deltaz = spacing(i);
          } else if (grid_type == 3) {
            // Implemented host or stretched momentum-input grid.
            avg_deltaz = (zm(i, k_ub_zm) - zm(i, k_lb_zm)) / real_t(nzm - 1);
          } else if (grid_type == 1) {
            // Evenly spaced grid.
            avg_deltaz = spacing(i);
          } else if (grid_type == 2) {
            // Stretched thermodynamic-input grid.
            avg_deltaz = (zt(i, k_ub_zt) - zt(i, k_lb_zt)) / real_t(nzt - 1);
          }
          // TODO: Like Native, invalid grid_type without prescribed spacing
          // leaves avg_deltaz undefined. The existing setup-grid gate excludes it.

          // The nu values are chosen for deltaz <= 40m. Adjust larger spacing
          // (Vince Larson), using a constant multiplier over each column.
          if (avg_deltaz > grid_spacing_thresh) {
            mult_factor_zt = real_t(1.0) + parameters(i, imult_coef) *
                                               Kokkos::log(avg_deltaz / grid_spacing_thresh);
            mult_factor_zm = mult_factor_zt;
          } else {
            mult_factor_zt = real_t(1.0);
            mult_factor_zm = real_t(1.0);
          }
          nu1(i) = parameters(i, inu1) * mult_factor_zm;
          nu2(i) = parameters(i, inu2) * mult_factor_zm;
          nu6(i) = parameters(i, inu6) * mult_factor_zm;
          nu8(i) = parameters(i, inu8) * mult_factor_zt;
          nu9(i) = parameters(i, inu9) * mult_factor_zm;
          nu10(i) = parameters(i, inu10) * mult_factor_zt;
          // We're unsure of the grid.
          nu_hm(i) = parameters(i, inu_hm) * mult_factor_zt;
        } else {
          // Nu values are not adjusted.
          nu1(i) = parameters(i, inu1);
          nu2(i) = parameters(i, inu2);
          nu6(i) = parameters(i, inu6);
          nu8(i) = parameters(i, inu8);
          nu9(i) = parameters(i, inu9);
          nu10(i) = parameters(i, inu10);
          nu_hm(i) = parameters(i, inu_hm);
        }
      });
  from_execution_space(mixt_frac_max_mag, mixture_limit);
}

// Native initializer uses HostSpace at its file boundary. Keep its out-of-line
// body separate from generated UDT conversion headers to avoid an include cycle.
void init_clubb_params_api(
    int ngrdcol, int iunit, const std::string& filename,                         // In
    Kokkos::View<real_t**, Kokkos::LayoutLeft, Kokkos::HostSpace> clubb_params); // Out

template <class Matrix>
void init_clubb_params_api(int ngrdcol, int iunit, const std::string& filename, // In
                           Matrix clubb_params) {                               // Out
  // Host/strided/device callers retain their output allocation. The actual
  // HostSpace driver takes the fixed overload directly, without a round trip.
  static_assert(Matrix::rank == 2);
  Kokkos::View<real_t**, Kokkos::LayoutLeft, Kokkos::HostSpace> parameters(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "init_clubb_params_api::output"), ngrdcol,
      nparams);
  init_clubb_params_api(ngrdcol, iunit, filename, // In
                        parameters);              // Out
  from_execution_space(clubb_params, parameters);
}

// TODO: Native tuner-only read_param_minmax and read_param_constraints are
// outside the actual JAX/standalone supported scope; they are not ported here.

namespace parameters_tunable_detail {
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
) {
  // Description:
  //   Put each column's coefficient list into its parameter vector. Keep this
  //   packing code together so it can be updated when new variables are added.
  // References:
  //   None.
  // ngrdcol: Number of columns [#]. Input coefficients retain their Native units.
  // clubb_params: Output matrix (column,parameter) [-].
  using ExecutionSpace = typename Matrix::execution_space;
  Kokkos::parallel_for(
      "pack_parameters", Kokkos::RangePolicy<ExecutionSpace>(0, ngrdcol),
      KOKKOS_LAMBDA(const int i) {
        clubb_params(i, iC1) = C1(i);
        clubb_params(i, iC1b) = C1b(i);
        clubb_params(i, iC1c) = C1c(i);
        clubb_params(i, iC2rt) = C2rt(i);
        clubb_params(i, iC2thl) = C2thl(i);
        clubb_params(i, iC2rtthl) = C2rtthl(i);
        clubb_params(i, iC4) = C4(i);
        clubb_params(i, iC_uu_shr) = C_uu_shr(i);
        clubb_params(i, iC_uu_buoy) = C_uu_buoy(i);
        clubb_params(i, iC6rt) = C6rt(i);
        clubb_params(i, iC6rtb) = C6rtb(i);
        clubb_params(i, iC6rtc) = C6rtc(i);
        clubb_params(i, iC6thl) = C6thl(i);
        clubb_params(i, iC6thlb) = C6thlb(i);
        clubb_params(i, iC6thlc) = C6thlc(i);
        clubb_params(i, iC7) = C7(i);
        clubb_params(i, iC7b) = C7b(i);
        clubb_params(i, iC7c) = C7c(i);
        clubb_params(i, iC8) = C8(i);
        clubb_params(i, iC8b) = C8b(i);
        clubb_params(i, iC10) = C10(i);
        clubb_params(i, iC11) = C11(i);
        clubb_params(i, iC11b) = C11b(i);
        clubb_params(i, iC11c) = C11c(i);
        clubb_params(i, iC12) = C12(i);
        clubb_params(i, iC13) = C13(i);
        clubb_params(i, iC14) = C14(i);
        clubb_params(i, iC_wp2_pr_dfsn) = C_wp2_pr_dfsn(i);
        clubb_params(i, iC_wp3_pr_tp) = C_wp3_pr_tp(i);
        clubb_params(i, iC_wp3_pr_turb) = C_wp3_pr_turb(i);
        clubb_params(i, iC_wp3_pr_dfsn) = C_wp3_pr_dfsn(i);
        clubb_params(i, iC_wp2_splat) = C_wp2_splat(i);
        clubb_params(i, iC6rt_Lscale0) = C6rt_Lscale0(i);
        clubb_params(i, iC6thl_Lscale0) = C6thl_Lscale0(i);
        clubb_params(i, iC7_Lscale0) = C7_Lscale0(i);
        clubb_params(i, iwpxp_L_thresh) = wpxp_L_thresh(i);
        clubb_params(i, ic_K) = c_K(i);
        clubb_params(i, ic_K1) = c_K1(i);
        clubb_params(i, inu1) = nu1(i);
        clubb_params(i, ic_K2) = c_K2(i);
        clubb_params(i, inu2) = nu2(i);
        clubb_params(i, ic_K6) = c_K6(i);
        clubb_params(i, inu6) = nu6(i);
        clubb_params(i, ic_K8) = c_K8(i);
        clubb_params(i, inu8) = nu8(i);
        clubb_params(i, ic_K9) = c_K9(i);
        clubb_params(i, inu9) = nu9(i);
        clubb_params(i, inu10) = nu10(i);
        clubb_params(i, ic_K_hm) = c_K_hm(i);
        clubb_params(i, ic_K_hmb) = c_K_hmb(i);
        clubb_params(i, iK_hm_min_coef) = K_hm_min_coef(i);
        clubb_params(i, inu_hm) = nu_hm(i);
        clubb_params(i, islope_coef_spread_DG_means_w) = slope_coef_spread_DG_means_w(i);
        clubb_params(i, ipdf_component_stdev_factor_w) = pdf_component_stdev_factor_w(i);
        clubb_params(i, icoef_spread_DG_means_rt) = coef_spread_DG_means_rt(i);
        clubb_params(i, icoef_spread_DG_means_thl) = coef_spread_DG_means_thl(i);
        clubb_params(i, igamma_coef) = gamma_coef(i);
        clubb_params(i, igamma_coefb) = gamma_coefb(i);
        clubb_params(i, igamma_coefc) = gamma_coefc(i);
        clubb_params(i, iLscale_mu_entrain) = Lscale_mu_entrain(i);
        clubb_params(i, ibeta) = beta(i);
        clubb_params(i, ilmin_coef) = lmin_coef(i);
        clubb_params(i, iomicron) = omicron(i);
        clubb_params(i, izeta_vrnce_rat) = zeta_vrnce_rat(i);
        clubb_params(i, iupsilon_precip_frac_rat) = upsilon_precip_frac_rat(i);
        clubb_params(i, ilambda0_stability_coef) = lambda0_stability_coef(i);
        clubb_params(i, imult_coef) = mult_coef(i);
        clubb_params(i, itaumin) = taumin(i);
        clubb_params(i, itaumax) = taumax(i);
        clubb_params(i, iLscale_mu_coef) = Lscale_mu_coef(i);
        clubb_params(i, iLscale_pert_coef) = Lscale_pert_coef(i);
        clubb_params(i, ialpha_corr) = alpha_corr(i);
        clubb_params(i, iSkw_denom_coef) = Skw_denom_coef(i);
        clubb_params(i, ic_K10) = c_K10(i);
        clubb_params(i, ic_K10h) = c_K10h(i);
        clubb_params(i, ithlp2_rad_coef) = thlp2_rad_coef(i);
        clubb_params(i, ithlp2_rad_cloud_frac_thresh) = thlp2_rad_cloud_frac_thresh(i);
        clubb_params(i, iup2_sfc_coef) = up2_sfc_coef(i);
        clubb_params(i, iSkw_max_mag) = Skw_max_mag(i);
        clubb_params(i, ixp3_coef_base) = xp3_coef_base(i);
        clubb_params(i, ixp3_coef_slope) = xp3_coef_slope(i);
        clubb_params(i, ialtitude_threshold) = altitude_threshold(i);
        clubb_params(i, irtp2_clip_coef) = rtp2_clip_coef(i);
        clubb_params(i, iC_invrs_tau_bkgnd) = C_invrs_tau_bkgnd(i);
        clubb_params(i, iC_invrs_tau_sfc) = C_invrs_tau_sfc(i);
        clubb_params(i, iC_invrs_tau_shear) = C_invrs_tau_shear(i);
        clubb_params(i, iC_invrs_tau_N2) = C_invrs_tau_N2(i);
        clubb_params(i, iC_invrs_tau_N2_wp2) = C_invrs_tau_N2_wp2(i);
        clubb_params(i, iC_invrs_tau_N2_xp2) = C_invrs_tau_N2_xp2(i);
        clubb_params(i, iC_invrs_tau_N2_wpxp) = C_invrs_tau_N2_wpxp(i);
        clubb_params(i, iC_invrs_tau_N2_clear_wp3) = C_invrs_tau_N2_clear_wp3(i);
        clubb_params(i, iC_invrs_tau_wpxp_Ri) = C_invrs_tau_wpxp_Ri(i);
        clubb_params(i, iC_invrs_tau_wpxp_N2_thresh) = C_invrs_tau_wpxp_N2_thresh(i);
        clubb_params(i, iCx_min) = Cx_min(i);
        clubb_params(i, iCx_max) = Cx_max(i);
        clubb_params(i, iRichardson_num_min) = Richardson_num_min(i);
        clubb_params(i, iRichardson_num_max) = Richardson_num_max(i);
        clubb_params(i, ia3_coef_min) = a3_coef_min(i);
        clubb_params(i, ia_const) = a_const(i);
        clubb_params(i, ibv_efold) = bv_efold(i);
        clubb_params(i, iwpxp_Ri_exp) = wpxp_Ri_exp(i);
        clubb_params(i, iz_displace) = z_displace(i);
        clubb_params(i, iLscale_tke_multiplier) = Lscale_tke_multiplier(i);
        clubb_params(i, iLscale_pdf_launch_frac) = Lscale_pdf_launch_frac(i);
        clubb_params(i, iLscale_up_weight) = Lscale_up_weight(i);
      });
}

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
) {
  // Description:
  //   Take the 1D parameter vector and return the list of scalar variables.
  //   Keep this code together so it can be updated when new variables are added.
  // References:
  //   None.
  // clubb_params: One parameter vector (nparams) [-]; no column argument exists.
  // Output coefficients retain their Native physical meanings and units.
  const auto parameters = to_host_space(clubb_params);
  C1 = parameters(iC1);
  C1b = parameters(iC1b);
  C1c = parameters(iC1c);
  C2rt = parameters(iC2rt);
  C2thl = parameters(iC2thl);
  C2rtthl = parameters(iC2rtthl);
  C4 = parameters(iC4);
  C_uu_shr = parameters(iC_uu_shr);
  C_uu_buoy = parameters(iC_uu_buoy);
  C6rt = parameters(iC6rt);
  C6rtb = parameters(iC6rtb);
  C6rtc = parameters(iC6rtc);
  C6thl = parameters(iC6thl);
  C6thlb = parameters(iC6thlb);
  C6thlc = parameters(iC6thlc);
  C7 = parameters(iC7);
  C7b = parameters(iC7b);
  C7c = parameters(iC7c);
  C8 = parameters(iC8);
  C8b = parameters(iC8b);
  C10 = parameters(iC10);
  C11 = parameters(iC11);
  C11b = parameters(iC11b);
  C11c = parameters(iC11c);
  C12 = parameters(iC12);
  C13 = parameters(iC13);
  C14 = parameters(iC14);
  C_wp2_pr_dfsn = parameters(iC_wp2_pr_dfsn);
  C_wp3_pr_tp = parameters(iC_wp3_pr_tp);
  C_wp3_pr_turb = parameters(iC_wp3_pr_turb);
  C_wp3_pr_dfsn = parameters(iC_wp3_pr_dfsn);
  C_wp2_splat = parameters(iC_wp2_splat);
  C6rt_Lscale0 = parameters(iC6rt_Lscale0);
  C6thl_Lscale0 = parameters(iC6thl_Lscale0);
  C7_Lscale0 = parameters(iC7_Lscale0);
  wpxp_L_thresh = parameters(iwpxp_L_thresh);
  c_K = parameters(ic_K);
  c_K1 = parameters(ic_K1);
  nu1 = parameters(inu1);
  c_K2 = parameters(ic_K2);
  nu2 = parameters(inu2);
  c_K6 = parameters(ic_K6);
  nu6 = parameters(inu6);
  c_K8 = parameters(ic_K8);
  nu8 = parameters(inu8);
  c_K9 = parameters(ic_K9);
  nu9 = parameters(inu9);
  nu10 = parameters(inu10);
  c_K_hm = parameters(ic_K_hm);
  c_K_hmb = parameters(ic_K_hmb);
  K_hm_min_coef = parameters(iK_hm_min_coef);
  nu_hm = parameters(inu_hm);
  slope_coef_spread_DG_means_w = parameters(islope_coef_spread_DG_means_w);
  pdf_component_stdev_factor_w = parameters(ipdf_component_stdev_factor_w);
  coef_spread_DG_means_rt = parameters(icoef_spread_DG_means_rt);
  coef_spread_DG_means_thl = parameters(icoef_spread_DG_means_thl);
  gamma_coef = parameters(igamma_coef);
  gamma_coefb = parameters(igamma_coefb);
  gamma_coefc = parameters(igamma_coefc);
  Lscale_mu_entrain = parameters(iLscale_mu_entrain);
  beta = parameters(ibeta);
  lmin_coef = parameters(ilmin_coef);
  omicron = parameters(iomicron);
  zeta_vrnce_rat = parameters(izeta_vrnce_rat);
  upsilon_precip_frac_rat = parameters(iupsilon_precip_frac_rat);
  lambda0_stability_coef = parameters(ilambda0_stability_coef);
  mult_coef = parameters(imult_coef);
  taumin = parameters(itaumin);
  taumax = parameters(itaumax);
  Lscale_mu_coef = parameters(iLscale_mu_coef);
  Lscale_pert_coef = parameters(iLscale_pert_coef);
  alpha_corr = parameters(ialpha_corr);
  Skw_denom_coef = parameters(iSkw_denom_coef);
  c_K10 = parameters(ic_K10);
  c_K10h = parameters(ic_K10h);
  thlp2_rad_coef = parameters(ithlp2_rad_coef);
  thlp2_rad_cloud_frac_thresh = parameters(ithlp2_rad_cloud_frac_thresh);
  up2_sfc_coef = parameters(iup2_sfc_coef);
  Skw_max_mag = parameters(iSkw_max_mag);
  xp3_coef_base = parameters(ixp3_coef_base);
  xp3_coef_slope = parameters(ixp3_coef_slope);
  altitude_threshold = parameters(ialtitude_threshold);
  rtp2_clip_coef = parameters(irtp2_clip_coef);
  C_invrs_tau_bkgnd = parameters(iC_invrs_tau_bkgnd);
  C_invrs_tau_sfc = parameters(iC_invrs_tau_sfc);
  C_invrs_tau_shear = parameters(iC_invrs_tau_shear);
  C_invrs_tau_N2 = parameters(iC_invrs_tau_N2);
  C_invrs_tau_N2_wp2 = parameters(iC_invrs_tau_N2_wp2);
  C_invrs_tau_N2_xp2 = parameters(iC_invrs_tau_N2_xp2);
  C_invrs_tau_N2_wpxp = parameters(iC_invrs_tau_N2_wpxp);
  C_invrs_tau_N2_clear_wp3 = parameters(iC_invrs_tau_N2_clear_wp3);
  C_invrs_tau_wpxp_Ri = parameters(iC_invrs_tau_wpxp_Ri);
  C_invrs_tau_wpxp_N2_thresh = parameters(iC_invrs_tau_wpxp_N2_thresh);
  Cx_min = parameters(iCx_min);
  Cx_max = parameters(iCx_max);
  Richardson_num_min = parameters(iRichardson_num_min);
  Richardson_num_max = parameters(iRichardson_num_max);
  a3_coef_min = parameters(ia3_coef_min);
  a_const = parameters(ia_const);
  bv_efold = parameters(ibv_efold);
  wpxp_Ri_exp = parameters(iwpxp_Ri_exp);
  z_displace = parameters(iz_displace);
  Lscale_tke_multiplier = parameters(iLscale_tke_multiplier);
  Lscale_pdf_launch_frac = parameters(iLscale_pdf_launch_frac);
  Lscale_up_weight = parameters(iLscale_up_weight);
}
} // namespace parameters_tunable_detail
