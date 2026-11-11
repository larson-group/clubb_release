! Compiled Fortran oracles for the temperature, skewness and one-level PDF kernels.
subroutine f2k_temperature_reference_c(nz, ngrdcol, thlm, exner, rcm, output) bind(C)
  ! Description: Compare the scalar, profile and column temperature formulas and
  ! their elemental inverse on the same inputs, without copying those formulas.
  use, intrinsic :: iso_c_binding, only: c_int
  use clubb_precision, only: core_rknd
  use T_in_K_module, only: thlm2T_in_K_api, T_in_K2thlm_api
  implicit none

  !--------------------------- Input Variables ---------------------------
  integer(c_int), value, intent(in) :: nz, ngrdcol ! Levels and columns [-]
  real( kind = core_rknd ), intent(in) :: thlm(ngrdcol,nz) ! Liquid potential temperature [K]
  real( kind = core_rknd ), intent(in) :: exner(ngrdcol,nz) ! Exner function [-]
  real( kind = core_rknd ), intent(in) :: rcm(ngrdcol,nz) ! Liquid mixing ratio [kg/kg]

  !--------------------------- Output Variables ---------------------------
  real( kind = core_rknd ), intent(out) :: output(ngrdcol,nz,4) ! Temperatures [K]

  !--------------------------- Local Variables ---------------------------
  integer :: i, k ! Column and level indices [-]

  !--------------------------- Begin Code ---------------------------
  do k = 1, nz
    do i = 1, ngrdcol
      output(i,k,1) = thlm2T_in_K_api(thlm(i,k), exner(i,k), rcm(i,k))
    end do
  end do

  do i = 1, ngrdcol
    output(i,:,2) = thlm2T_in_K_api(nz, thlm(i,:), exner(i,:), rcm(i,:))
  end do

  output(:,:,3) = thlm2T_in_K_api(nz, ngrdcol, thlm, exner, rcm)
  output(:,:,4) = T_in_K2thlm_api(output(:,:,3), exner, rcm)
end subroutine f2k_temperature_reference_c

subroutine f2k_skewness_reference_c(nz, ngrdcol, xp2, xp3, wp2, wpxp, sigma_sqd_w, &
                                    Skw, clubb_params, x_tol, output) bind(C)
  ! Description: Evaluate all four native Skx_module routines on supplied moments.
  use, intrinsic :: iso_c_binding, only: c_int
  use clubb_precision, only: core_rknd
  use parameter_indices, only: nparams, ibeta
  use Skx_module, only: Skx_func, compute_gamma_Skw, LG_2005_ansatz, xp3_LG_2005_ansatz
  implicit none

  !--------------------------- Input Variables ---------------------------
  integer(c_int), value, intent(in) :: nz, ngrdcol ! Levels and columns [-]
  real( kind = core_rknd ), intent(in) :: xp2(ngrdcol,nz), xp3(ngrdcol,nz) ! Moments of x
  real( kind = core_rknd ), intent(in) :: wp2(ngrdcol,nz) ! Vertical velocity variance [m^2/s^2]
  real( kind = core_rknd ), intent(in) :: wpxp(ngrdcol,nz) ! Turbulent flux [m/s (x units)]
  real( kind = core_rknd ), intent(in) :: &
    sigma_sqd_w(ngrdcol,nz), Skw(ngrdcol,nz) ! Width/skewness [-]
  real( kind = core_rknd ), intent(in) :: &
    clubb_params(ngrdcol,nparams) ! Tunable parameters [units vary]
  real( kind = core_rknd ), value, intent(in) :: x_tol ! Minimum x tolerance [x units]

  !--------------------------- Output Variables ---------------------------
  real( kind = core_rknd ), intent(out) :: &
    output(ngrdcol,nz,4) ! Skewness, gamma, ansatz, third moment

  !--------------------------- Begin Code ---------------------------
  call Skx_func(nz, ngrdcol, xp2, xp3, & ! In
                 x_tol, clubb_params, & ! In
                 output(:,:,1))        ! Out

  call compute_gamma_Skw(nz, ngrdcol, Skw, clubb_params, & ! In
                          output(:,:,2))                 ! Out

  call LG_2005_ansatz(nz, ngrdcol, Skw, wpxp, wp2, & ! In
                       xp2, clubb_params(:,ibeta), sigma_sqd_w, x_tol, & ! In
                       output(:,:,3))                                    ! Out

  call xp3_LG_2005_ansatz(nz, ngrdcol, Skw, wpxp, wp2, & ! In
                           xp2, sigma_sqd_w, &           ! In
                           clubb_params, x_tol, &        ! In
                           output(:,:,4))                ! Out
end subroutine f2k_skewness_reference_c

subroutine f2k_sigma_one_level_reference_c(ngrdcol, wp3, wp2, thlp2, rtp2, &
                                           up2, vp2, wpthlp, wprtp, upwp, vpwp, &
                                           clubb_params, predict_winds, output) bind(C)
  ! Description: Evaluate the native nzm=1 branch with and without momentum fluxes.
  use, intrinsic :: iso_c_binding, only: c_int
  use clubb_precision, only: core_rknd
  use parameter_indices, only: nparams
  use grid_class, only: grid
  use sigma_sqd_w_module, only: compute_sigma_sqd_w
  implicit none

  !--------------------------- Input Variables ---------------------------
  integer(c_int), value, intent(in) :: ngrdcol ! Columns [-]
  real( kind = core_rknd ), intent(in) :: wp3(ngrdcol,1) ! Third moment [m^3/s^3]
  real( kind = core_rknd ), intent(in) :: wp2(ngrdcol,1), up2(ngrdcol,1), vp2(ngrdcol,1) ! [m^2/s^2]
  real( kind = core_rknd ), intent(in) :: thlp2(ngrdcol,1), rtp2(ngrdcol,1) ! Scalar variances
  real( kind = core_rknd ), intent(in) :: wpthlp(ngrdcol,1), wprtp(ngrdcol,1) ! Scalar fluxes
  real( kind = core_rknd ), intent(in) :: &
    upwp(ngrdcol,1), vpwp(ngrdcol,1) ! Momentum fluxes [m^2/s^2]
  real( kind = core_rknd ), intent(in) :: clubb_params(ngrdcol,nparams) ! Parameters [units vary]
  integer(c_int), value, intent(in) :: predict_winds ! Include momentum fluxes [-]

  !--------------------------- Output Variables ---------------------------
  real( kind = core_rknd ), intent(out) :: output(ngrdcol,1) ! PDF width [-]

  !--------------------------- Local Variables ---------------------------
  type(grid) :: gr ! Unused by the native one-level test branch

  !--------------------------- Begin Code ---------------------------
  call compute_sigma_sqd_w(1, 1, ngrdcol, gr, & ! In
                            wp3, wp2, thlp2, rtp2, & ! In
                            up2, vp2, wpthlp, wprtp, upwp, vpwp, & ! In
                            clubb_params, & ! In
                            predict_winds /= 0, & ! In
                            output) ! Out
end subroutine f2k_sigma_one_level_reference_c

subroutine f2k_all_core_constants_c(values, doubles, ints, indices) bind(C)
  ! Description: Export every standalone constant and tunable index from compiled native modules.
  use, intrinsic :: iso_c_binding, only: c_int, c_double
  use clubb_precision, only: core_rknd
  use constants_clubb
  use parameter_indices
  implicit none

  !--------------------------- Output Variables ---------------------------
  real( kind = core_rknd ), intent(out) :: values(96) ! Core constants [units vary]
  real(c_double), intent(out) :: doubles(9) ! Double constants [units vary]
  integer(c_int), intent(out) :: ints(5) ! I/O and array sizes [-]
  integer(c_int), intent(out) :: indices(106) ! Native indices and count [-]

  !--------------------------- Begin Code ---------------------------
  values = [ &
    parab_cyl_max_input, gamma_over_implicit_ts, pi, invrs_pi, &
    radians_per_deg, sqrt_2pi, sqrt_2, one_hundred, &
    fifty, thirty_six, twenty, eighteen, &
    twelve, ten, nine, eight, &
    seven, six, five, four, &
    three, two, three_halves, four_thirds, &
    one, three_fourths, two_thirds, one_half, &
    one_third, one_fourth, three_sixteenth, one_eighth, &
    zero, Cp, Lv, Ls, &
    Lf, Rd, Rv, stefan_boltzmann, &
    T_freeze_K, ep, ep1, ep2, &
    kappa, omega_planet, grav, p0, &
    vonk, rho_lw, rho_ice, w_tol, &
    thl_tol, rt_tol, chi_tol, eta_tol, &
    thl_tol_mfl, rt_tol_mfl, w_tol_sqd, rc_tol, &
    Nc_tol, Ncn_tol, mvr_cloud_max, Nc_in_cloud_min, &
    rr_tol, ri_tol, rs_tol, rg_tol, &
    mvr_rain_max, mvr_ice_max, mvr_snow_max, mvr_graupel_max, &
    Nr_tol, Ni_tol, Ns_tol, Ng_tol, &
    em_min, eps, max_num_stdevs, zero_threshold, &
    max_mag_correlation, max_mag_correlation_flux, min_max_smth_mag, cloud_frac_min, &
    wp2_max, sec_per_day, sec_per_hr, sec_per_min, &
    min_per_hr, g_per_kg, pascal_per_mb, cm3_per_m3, &
    micron_per_m, cm_per_m, mm_per_m, unused_var]
  doubles = [ &
    pi_dp, radians_per_deg_dp, sqrt_2pi_dp, sqrt_2_dp, &
    two_dp, one_dp, one_half_dp, one_fourth_dp, &
    zero_dp]
  ints = [ &
    fstderr, fstdin, fstdout, var_length, &
    num_hf_draw_points]
  indices = [ &
    nparams, iC1, iC1b, iC1c, &
    iC2rt, iC2thl, iC2rtthl, iC4, &
    iC_uu_shr, iC_uu_buoy, iC6rt, iC6rtb, &
    iC6rtc, iC6thl, iC6thlb, iC6thlc, &
    iC7, iC7b, iC7c, iC8, &
    iC8b, iC10, iC11, iC11b, &
    iC11c, iC12, iC13, iC14, &
    iC_wp2_pr_dfsn, iC_wp3_pr_tp, iC_wp3_pr_turb, iC_wp3_pr_dfsn, &
    iC_wp2_splat, iC6rt_Lscale0, iC6thl_Lscale0, iC7_Lscale0, &
    iwpxp_L_thresh, ic_K, ic_K1, inu1, &
    ic_K2, inu2, ic_K6, inu6, &
    ic_K8, inu8, ic_K9, inu9, &
    inu10, ic_K_hm, ic_K_hmb, iK_hm_min_coef, &
    inu_hm, islope_coef_spread_DG_means_w, ipdf_component_stdev_factor_w, &
    icoef_spread_DG_means_rt, &
    icoef_spread_DG_means_thl, igamma_coef, igamma_coefb, igamma_coefc, &
    iLscale_mu_entrain, ibeta, ilmin_coef, iomicron, &
    izeta_vrnce_rat, iupsilon_precip_frac_rat, ilambda0_stability_coef, imult_coef, &
    itaumin, itaumax, iLscale_mu_coef, iLscale_pert_coef, &
    ialpha_corr, iSkw_denom_coef, ic_K10, ic_K10h, &
    ithlp2_rad_coef, ithlp2_rad_cloud_frac_thresh, iup2_sfc_coef, iSkw_max_mag, &
    iC_invrs_tau_bkgnd, iC_invrs_tau_sfc, iC_invrs_tau_shear, iC_invrs_tau_N2, &
    iC_invrs_tau_N2_wp2, iC_invrs_tau_N2_xp2, iC_invrs_tau_N2_wpxp, iC_invrs_tau_N2_clear_wp3, &
    iC_invrs_tau_wpxp_Ri, iC_invrs_tau_wpxp_N2_thresh, ixp3_coef_base, ixp3_coef_slope, &
    ialtitude_threshold, irtp2_clip_coef, iCx_min, iCx_max, &
    iRichardson_num_min, iRichardson_num_max, ia3_coef_min, ia_const, &
    ibv_efold, iwpxp_Ri_exp, iz_displace, iLscale_tke_multiplier, &
    iLscale_pdf_launch_frac, iLscale_up_weight]
end subroutine f2k_all_core_constants_c
