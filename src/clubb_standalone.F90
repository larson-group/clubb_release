!----------------------------------------------------------------------
! $Id$

program clubb_standalone

! Description:
! This is essentially a minimalist frontend for CLUBB.
!
!----------------------------------------------------------------------

  use clubb_driver, only: run_clubb ! Procedure(s)
  
  use error_code, only: &
        clubb_no_error, &               ! Constants
        clubb_fatal_error, &
        err_code                        ! Error indicator

  use parameter_indices, only: nparams ! Variable(s)

  use parameters_tunable, only: &
      set_default_parameters, & ! Procedure(s)
      read_parameters

  use clubb_precision, only: core_rknd ! Constant

  use constants_clubb, only: fstderr ! Constant

  implicit none

  ! External
  intrinsic :: trim

  ! Constant parameters
  integer, parameter :: iunit = 10

  integer, parameter :: success_code = 6

  character(len=13), parameter :: &
    namelist_filename = "clubb.in"  ! Text file containing namelists
                                    ! concatenated from various input files such as
                                    ! model.in, tunable_parameters.in, error....in.

  logical, parameter :: &
    l_stdout = .true.

  ! Run information
  real( kind = core_rknd ), dimension(nparams) :: & 
    params  ! Array of the model constants

  real( kind = core_rknd ) :: & 
    C1, C1b, C1c, C2rt, C2thl, C2rtthl, & 
    C4, C_uu_shr, C_uu_buoy, C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc, & 
    C7, C7b, C7c, C8, C8b, C10, & 
    C11, C11b, C11c, C12, C13, C14, C_wp2_pr_dfsn, C_wp3_pr_tp, &
    C_wp3_pr_turb, C_wp3_pr_dfsn, C_wp2_splat, & 
    C6rt_Lscale0, C6thl_Lscale0, C7_Lscale0, wpxp_L_thresh, &
    c_K, c_K1, nu1, c_K2, nu2, c_K6, nu6, c_K8, nu8,  & 
    c_K9, nu9, nu10, c_K_hm, c_K_hmb, K_hm_min_coef, nu_hm, &
    slope_coef_spread_DG_means_w, pdf_component_stdev_factor_w, &
    coef_spread_DG_means_rt, coef_spread_DG_means_thl, &
    gamma_coef, gamma_coefb, gamma_coefc, mu, beta, lmin_coef, &
    omicron, zeta_vrnce_rat, upsilon_precip_frac_rat, &
    lambda0_stability_coef, mult_coef, taumin, taumax, Lscale_mu_coef, &
    Lscale_pert_coef, alpha_corr, Skw_denom_coef, c_K10, c_K10h, &
    thlp2_rad_coef, thlp2_rad_cloud_frac_thresh, up2_sfc_coef, &
    Skw_max_mag, xp3_coef_base, xp3_coef_slope, altitude_threshold, &
    rtp2_clip_coef, C_invrs_tau_bkgnd, C_invrs_tau_sfc, &
    C_invrs_tau_shear, C_invrs_tau_N2, C_invrs_tau_N2_wp2, &
    C_invrs_tau_N2_xp2, C_invrs_tau_N2_wpxp, C_invrs_tau_N2_clear_wp3, &
    C_invrs_tau_wpxp_Ri, C_invrs_tau_wpxp_N2_thresh, &
    Cx_min, Cx_max, Richardson_num_min, Richardson_num_max, &
    wpxp_Ri_exp, a3_coef_min, a_const, bv_efold

!-----------------------------------------------------------------------

  ! --- Begin Code ---

  ! Set the default tunable parameter values
  call set_default_parameters( &
               C1, C1b, C1c, C2rt, C2thl, C2rtthl, &
               C4, C_uu_shr, C_uu_buoy, C6rt, C6rtb, C6rtc, &
               C6thl, C6thlb, C6thlc, C7, C7b, C7c, C8, C8b, C10, &
               C11, C11b, C11c, C12, C13, C14, C_wp2_pr_dfsn, C_wp3_pr_tp, &
               C_wp3_pr_turb, C_wp3_pr_dfsn, C_wp2_splat, &
               C6rt_Lscale0, C6thl_Lscale0, C7_Lscale0, wpxp_L_thresh, &
               c_K, c_K1, nu1, c_K2, nu2, c_K6, nu6, c_K8, nu8, &
               c_K9, nu9, nu10, c_K_hm, c_K_hmb, K_hm_min_coef, nu_hm, &
               slope_coef_spread_DG_means_w, pdf_component_stdev_factor_w, &
               coef_spread_DG_means_rt, coef_spread_DG_means_thl, &
               gamma_coef, gamma_coefb, gamma_coefc, mu, beta, lmin_coef, &
               omicron, zeta_vrnce_rat, upsilon_precip_frac_rat, &
               lambda0_stability_coef, mult_coef, taumin, taumax, &
               Lscale_mu_coef, Lscale_pert_coef, alpha_corr, &
               Skw_denom_coef, c_K10, c_K10h, thlp2_rad_coef, &
               thlp2_rad_cloud_frac_thresh, up2_sfc_coef, &
               Skw_max_mag, xp3_coef_base, xp3_coef_slope, &
               altitude_threshold, rtp2_clip_coef, C_invrs_tau_bkgnd, &
               C_invrs_tau_sfc, C_invrs_tau_shear, C_invrs_tau_N2, & 
               C_invrs_tau_N2_wp2, C_invrs_tau_N2_xp2, &
               C_invrs_tau_N2_wpxp, C_invrs_tau_N2_clear_wp3, &
               C_invrs_tau_wpxp_Ri, C_invrs_tau_wpxp_N2_thresh, &
               Cx_min, Cx_max, Richardson_num_min, Richardson_num_max, &
               wpxp_Ri_exp, a3_coef_min, a_const, bv_efold )

  ! Read in model parameter values
  call read_parameters( iunit, namelist_filename, &
                        C1, C1b, C1c, C2rt, C2thl, C2rtthl, &
                        C4, C_uu_shr, C_uu_buoy, C6rt, C6rtb, C6rtc, &
                        C6thl, C6thlb, C6thlc, C7, C7b, C7c, C8, C8b, C10, &
                        C11, C11b, C11c, C12, C13, C14, C_wp2_pr_dfsn, C_wp3_pr_tp, &
                        C_wp3_pr_turb, C_wp3_pr_dfsn, C_wp2_splat, &
                        C6rt_Lscale0, C6thl_Lscale0, C7_Lscale0, wpxp_L_thresh, &
                        c_K, c_K1, nu1, c_K2, nu2, c_K6, nu6, c_K8, nu8, &
                        c_K9, nu9, nu10, c_K_hm, c_K_hmb, K_hm_min_coef, nu_hm, &
                        slope_coef_spread_DG_means_w, pdf_component_stdev_factor_w, &
                        coef_spread_DG_means_rt, coef_spread_DG_means_thl, &
                        gamma_coef, gamma_coefb, gamma_coefc, mu, beta, lmin_coef, &
                        omicron, zeta_vrnce_rat, upsilon_precip_frac_rat, &
                        lambda0_stability_coef, mult_coef, taumin, taumax, &
                        Lscale_mu_coef, Lscale_pert_coef, alpha_corr, &
                        Skw_denom_coef, c_K10, c_K10h, thlp2_rad_coef, &
                        thlp2_rad_cloud_frac_thresh, up2_sfc_coef, &
                        Skw_max_mag, xp3_coef_base, xp3_coef_slope, &
                        altitude_threshold, rtp2_clip_coef, C_invrs_tau_bkgnd, &
                        C_invrs_tau_sfc, C_invrs_tau_shear, C_invrs_tau_N2, & 
                        C_invrs_tau_N2_wp2, C_invrs_tau_N2_xp2, &
                        C_invrs_tau_N2_wpxp, C_invrs_tau_N2_clear_wp3, &
                        C_invrs_tau_wpxp_Ri, C_invrs_tau_wpxp_N2_thresh, &
                        Cx_min, Cx_max, Richardson_num_min, Richardson_num_max, &
                        wpxp_Ri_exp, a3_coef_min, a_const, bv_efold, &
                        params )

  ! Initialize status of run 
  err_code = clubb_no_error

  ! Run the model
  call run_clubb( params, namelist_filename, l_stdout )

  if ( err_code == clubb_fatal_error ) then
    error stop "Fatal error in clubb, check your parameter values and timestep"
  else
    write(fstderr,*) "Program exited normally"
    call exit(success_code)
  end if 

end program clubb_standalone
