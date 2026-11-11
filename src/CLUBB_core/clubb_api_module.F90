!--------------------------------------------------------------------------------------------------
! $Id$
!==================================================================================================
!
!       ########  ###       ###    ### #########  #########           ###     ######### ###########
!     ###    ### ###       ###    ### ###    ### ###    ###        ### ###   ###    ###    ###
!    ###        ###       ###    ### ###    ### ###    ###       ###   ###  ###    ###    ###
!   ###        ###       ###    ### #########  #########       ########### #########     ###
!  ###        ###       ###    ### ###    ### ###    ###      ###     ### ###           ###
! ###    ### ###       ###    ### ###    ### ###    ###      ###     ### ###           ###
! ########  ########## ########  #########  #########       ###     ### ###       ###########
!
! The CLUBB API serves as the doorway through which external models can interact with CLUBB.
!
!               PLEASE REMEMBER, IF ANY CODE IS CHANGED IN THIS DOCUMENT,
!                   THE CHANGES MUST BE PROPOGATED TO ALL HOST MODELS.
!

module clubb_api_module

  use clubb_precision, only : &
    time_precision, &
    core_rknd, &
    stat_nknd, &
    stat_rknd, &
    dp  ! Double Precision

  use constants_clubb, only : &
    cloud_frac_min, & ! Threshold for cloud fractions
    cm3_per_m3, & ! Cubic centimeters per cubic meter
    Cp, & ! Dry air specific heat at constant p [J/kg/K]
    em_min, & ! Minimum value for em (turbulence kinetic energy)
    ep, & ! ep  = 0.622  [-]
    fstderr, & ! Fortran file unit I/O constant
    fstdout, & ! Fortran file unit I/O constant
    grav, & ! Gravitational acceleration     [m/s^2]
    Ls, & ! Latent heat of sublimation          [J/kg]
    Lv, & ! Latent heat of vaporization         [J/kg]
    Lf, & ! Latent heat of fusion               [J/kg]
    pi, & ! The ratio of radii to their circumference
    pi_dp, & ! pi in double precision
    radians_per_deg_dp, &
    Rd, & ! Dry air gas constant                [J/kg/K]
    Rv, & ! Water vapor gas constant            [J/kg/K]
    sec_per_day, & ! Seconds in a day.
    sec_per_hr, &  ! Seconds in an hour.
    sec_per_min, & ! Seconds in a minute.
    T_freeze_K, & ! Freezing point of water [K]
    var_length, & ! Maximum variable name length in CLUBB GrADS or netCDF output
    zero, & ! 0.0_core_rknd
    zero_threshold, & ! Defining a threshold on a physical quantity to be 0.
    ! Tolerances
    Nc_tol, & ! Tolerance value for N_c  [#/kg]
    Ng_tol, & ! Tolerance value for N_s [#/kg]
    Ni_tol, & ! Tolerance value for N_i [#/kg]
    Nr_tol, & ! Tolerance value for N_r [#/kg]
    Ns_tol, & ! Tolerance value for N_s [#/kg]
    rg_tol, & ! Tolerance value for r_g [kg/kg]
    rho_lw, &
    ri_tol, & ! Tolerance value for r_i [kg/kg]
    rr_tol, & ! Tolerance value for r_r [kg/kg]
    rs_tol, & ! Tolerance value for r_s [kg/kg]
    rt_tol, & ! [kg/kg]
    thl_tol, & ! [K]
    w_tol_sqd, & ! [m^2/s^2]
    num_hf_draw_points ! Neighboring points to use in hole filling algorithm

  use corr_varnce_module, only : &
      hmp2_ip_on_hmm2_ip_slope_type,      & ! Types
      hmp2_ip_on_hmm2_ip_intrcpt_type, &
      hm_metadata_type

  use array_index, only: &
      sclr_idx_type

  use stats_netcdf, only: &
      stats_type, &
      stats_init_api, &
      stats_begin_timestep_api, &
      stats_end_timestep_api, &
      stats_finalize_api

  use setup_clubb_pdf_params, only: &
      setup_pdf_parameters_api

  use stats_clubb_utilities, only: &
      stats_accumulate_hydromet_api

  use error_code, only: &
      set_clubb_debug_level_api,   &
      clubb_at_least_debug_level_api, &
      clubb_no_error,              & ! Constants
      clubb_fatal_error

  use hydromet_pdf_parameter_module, only : &
    hydromet_pdf_parameter, &
    precipitation_fractions, &
    init_hydromet_pdf_params, &
    init_precip_fracs_api, &
    zero_precip_fracs_api

  use model_flags, only : &
      clubb_config_flags_type, & ! Type
      iiPDF_ADG1, &
      iiPDF_new, &
      iiPDF_new_hybrid, &
      ipdf_pre_advance_fields, &
      ipdf_post_advance_fields, &
      l_use_boussinesq, &    ! Use Boussinesq form of predictive equations (default is Anelastic).
      print_clubb_config_flags_api, &
      global_fill, & 
      sliding_window, & 
      widening_windows, & 
      smart_window, & 
      smart_window_smooth, & 
      parallel_fill

  use parameters_tunable, only : &
    params_list,            & ! Variable(s)
    check_parameters_api,   & ! Procedure(s)
    nu_vertical_res_dep       ! Type(s)

  use parameter_indices, only:  &
    nparams, & ! Variable(s)
    iC1, iC1b, iC1c, &
    iC2rt, iC2thl, iC2rtthl, iC4, iC_uu_shr, iC_uu_buoy, &
    iC6rt, iC6rtb, iC6rtc, iC6thl, iC6thlb, iC6thlc, &
    iC7, iC7b, iC7c, iC8, iC8b, iC10, iC11, iC11b, iC11c, &
    iC12, iC13, iC14, iC_wp3_pr_turb, iC_wp3_pr_dfsn, iC_wp2_splat, &
    iC6rt_Lscale0, iC6thl_Lscale0, &
    iC7_Lscale0, iwpxp_L_thresh, ic_K, ic_K1, inu1, &
    ic_K2, inu2, ic_K6, inu6, ic_K8, inu8, ic_K9, inu9, &
    inu10, ic_K_hm, ic_K_hmb, iK_hm_min_coef, inu_hm, &
    islope_coef_spread_DG_means_w, ipdf_component_stdev_factor_w, &
    icoef_spread_DG_means_rt, icoef_spread_DG_means_thl, &
    ibeta, igamma_coef, igamma_coefb, igamma_coefc, ilmin_coef, &
    iomicron, izeta_vrnce_rat, iupsilon_precip_frac_rat, &
    ilambda0_stability_coef, imult_coef, itaumin, itaumax, iLscale_mu_entrain, &
    iLscale_mu_coef, iLscale_pert_coef, ialpha_corr, iSkw_denom_coef, &
    ic_K10, ic_K10h, ithlp2_rad_coef, ithlp2_rad_cloud_frac_thresh, &
    iup2_sfc_coef, iSkw_max_mag, iC_invrs_tau_bkgnd, &
    iC_invrs_tau_sfc, iC_invrs_tau_shear, iC_invrs_tau_N2, &
    iC_invrs_tau_N2_wp2, iC_invrs_tau_N2_xp2, iC_invrs_tau_N2_wpxp, &
    iC_invrs_tau_N2_clear_wp3, ialtitude_threshold, irtp2_clip_coef, &
    iRichardson_num_min, iRichardson_num_max, iwpxp_Ri_exp, &
    ia3_coef_min, ia_const, iCx_min, iCx_max, ibv_efold, iz_displace, &
    iLscale_tke_multiplier, iLscale_pdf_launch_frac, iLscale_up_weight


  use pdf_parameter_module, only : &
! The CLUBB_CAM preprocessor directives are being commented out because this
! code is now also used for WRF-CLUBB.
!#ifdef CLUBB_CAM /* Code for storing pdf_parameter structs in pbuf as array */
!#endif
    pdf_parameter, &
    implicit_coefs_terms, &
    zero_pdf_params_api, &
    init_pdf_implicit_coefs_terms_api, &
    zero_pdf_implicit_coefs_terms_api

  use sponge_layer_damping, only : &
    thlm_sponge_damp_settings,    & ! Variable(s)
    rtm_sponge_damp_settings,     &
    uv_sponge_damp_settings,      &
    wp2_sponge_damp_settings,     &
    wp3_sponge_damp_settings,     &
    up2_vp2_sponge_damp_settings, &
    thlm_sponge_damp_profile,     &
    rtm_sponge_damp_profile,      &
    uv_sponge_damp_profile,       &
    wp2_sponge_damp_profile,      &
    wp3_sponge_damp_profile,      &
    up2_vp2_sponge_damp_profile,  &
    finalize_tau_sponge_damp_api, &   ! Procedure(s)
    initialize_tau_sponge_damp_api

  use advance_helper_module, only: vertical_integral_api

  use grid_class, only: grid ! Type

  use grid_class, only: &
    ! Interpolate momentum level variables to thermodynamic levels
    zm2zt_api, &

    ! Interpolate thermodynamic level variables to momentum levels
    zt2zm_api, &

    ! Refresh grid heights and interpolation weights
    setup_grid_heights_api, &

    ! Deallocate grid variable
    cleanup_grid_api

  use saturation, only: &
    ! Used to compute the saturation mixing ratio of liquid water.
    sat_mixrat_liq_api

  use T_in_K_module, only: &
    ! Calculates absolute temperature from liquid water potential
    ! temperature.  (Does not include ice.)
    thlm2T_in_K_api, &

    ! Calculates liquid water potential temperature from absolute temperature
    T_in_K2thlm_api

  use calendar, only: &
    compute_current_date_api

  use interpolation, only: &
    lin_interpolate_on_grid_api

  use numerical_check, only : &
    check_clubb_settings_api

  use fill_holes, only : &
    fill_holes_driver_api, &
    fill_holes_hydromet_api, &
    fill_holes_vertical_api

  use model_flags, only: &
    initialize_clubb_config_flags_type_api, &
    set_default_clubb_config_flags_api, &
    saturation_bolton, & ! Constant for Bolton approximations of saturation
    saturation_gfdl,   & ! Constant for the GFDL approximation of saturation
    saturation_flatau, & ! Constant for Flatau approximations of saturation
    saturation_lookup    ! Use a lookup table for mixing length

  use corr_varnce_module, only : &
    setup_corr_varnce_array_api, &
    init_pdf_hydromet_arrays_api

  use parameters_tunable, only: &
    init_clubb_params_api

  use err_info_type_module, only: &
    err_info_type,              & ! Type
    init_err_info_api,          & ! Procedures
    init_default_err_info_api,  &
    set_err_info_values_api,    &
    cleanup_err_info_api

  implicit none

  private

  public &
    ! To Implement CLUBB:
    check_clubb_settings_api, &
    init_clubb_params_api, &
        ! CLUBB can be set more specifically using these flags:
        iiPDF_ADG1, &
        iiPDF_new, &
        iiPDF_new_hybrid, &
        ipdf_pre_advance_fields, &
        ipdf_post_advance_fields, &
        l_use_boussinesq, &
        saturation_bolton, &
        saturation_gfdl, &
        saturation_flatau, &
        saturation_lookup, &
        ! The parameters of CLUBB can be retrieved and tuned using these indices:
        iC1, iC1b, iC1c, &
        iC2rt, iC2thl, iC2rtthl, iC4, iC_uu_shr, iC_uu_buoy, &
        iC6rt, iC6rtb, iC6rtc, iC6thl, iC6thlb, iC6thlc, &
        iC7, iC7b, iC7c, iC8, iC8b, iC10, iC11, iC11b, iC11c, &
        iC12, iC13, iC14, iC_wp3_pr_turb, iC_wp3_pr_dfsn, iC_wp2_splat, & 
        iC6rt_Lscale0, iC6thl_Lscale0, &
        iC7_Lscale0, iwpxp_L_thresh, ic_K, ic_K1, inu1, &
        ic_K2, inu2, ic_K6, inu6, ic_K8, inu8, ic_K9, inu9, &
        inu10, ic_K_hm, ic_K_hmb, iK_hm_min_coef, inu_hm, &
        islope_coef_spread_DG_means_w, ipdf_component_stdev_factor_w, &
        icoef_spread_DG_means_rt, icoef_spread_DG_means_thl, &
        ibeta, igamma_coef, igamma_coefb, igamma_coefc, ilmin_coef, &
        iomicron, izeta_vrnce_rat, iupsilon_precip_frac_rat, &
        ilambda0_stability_coef, imult_coef, itaumin, itaumax, iLscale_mu_entrain, &
        iLscale_mu_coef, iLscale_pert_coef, ialpha_corr, iSkw_denom_coef, &
        ic_K10, ic_K10h, ithlp2_rad_coef, ithlp2_rad_cloud_frac_thresh, &
        iup2_sfc_coef, iSkw_max_mag, iC_invrs_tau_bkgnd, &
        iC_invrs_tau_sfc, iC_invrs_tau_shear, iC_invrs_tau_N2, &
        iC_invrs_tau_N2_wp2, iC_invrs_tau_N2_xp2, iC_invrs_tau_N2_wpxp, &
        iC_invrs_tau_N2_clear_wp3, ialtitude_threshold, irtp2_clip_coef, &
        iRichardson_num_min, iRichardson_num_max, iwpxp_Ri_exp, &
        ia3_coef_min, ia_const, iCx_min, iCx_max, ibv_efold, iz_displace, &
        iLscale_tke_multiplier, iLscale_pdf_launch_frac, iLscale_up_weight

  public &
        advance_clubb_core_api, &
        pdf_parameter, &
        implicit_coefs_terms

  public &
    ! To Implement SILHS:
    setup_corr_varnce_array_api, &
    hydromet_pdf_parameter, &
    init_pdf_hydromet_arrays_api, &
    ! To use the results, you will need these variables:
    hmp2_ip_on_hmm2_ip_slope_type,      & ! Types
    hmp2_ip_on_hmm2_ip_intrcpt_type, &
    hm_metadata_type, &
    sclr_idx_type, &
    grid

  public &
    ! To Interact With CLUBB's Grid:
    ! For Varying Grids
    setup_grid_heights_api, &
    setup_grid_api

  public &
    ! To Obtain NetCDF Stats Output from CLUBB:
    stats_type, &
    stats_init_api, &
    stats_begin_timestep_api, &
    stats_end_timestep_api, &
    stats_finalize_api, &
    setup_pdf_parameters_api, &
    stats_accumulate_hydromet_api

  public :: &
    calculate_thlp2_rad_api, params_list, &
    update_xp2_mc_api, sat_mixrat_liq_api

  public :: &
    ! To Convert Between Common CLUBB-related quantities:
    lin_interpolate_two_points_api, & ! OR
    lin_interpolate_on_grid_api, &
    T_in_K2thlm_api, &
    thlm2T_in_K_api, &
    zm2zt_api, &
    zt2zm_api, &
    cleanup_grid_api

  public &
    ! To Check For and Handle CLUBB's Errors:
    calculate_spurious_source_api, &
    clubb_at_least_debug_level_api, &
    clubb_fatal_error, &
    clubb_no_error, &
    fill_holes_driver_api, & ! OR
    fill_holes_vertical_api, &
    fill_holes_hydromet_api, &
    set_clubb_debug_level_api, &
    vertical_integral_api, &
    num_hf_draw_points, &
    global_fill, & 
    sliding_window, & 
    widening_windows, & 
    smart_window, & 
    smart_window_smooth, & 
    parallel_fill

  public &
    ! Constants That May be Helpful:
    cloud_frac_min, &
    cm3_per_m3, &
    core_rknd, &
    Cp, &
    dp, &
    em_min, &
    ep, &
    fstderr, &
    fstdout, &
    grav, &
    Lf, &
    Ls, &
    Lv, &
    pi_dp, &
    pi, &
    radians_per_deg_dp, &
    Rd, &
    Rv, &
    sec_per_day, &
    sec_per_hr, &
    sec_per_min, &
    T_freeze_K, &
    time_precision, &
    var_length, &
    zero_threshold, &
    zero, &
    ! Tolerances
    Nc_tol, &
    Ng_tol, &
    Ni_tol, &
    Nr_tol, &
    Ns_tol, &
    rg_tol, &
    rho_lw, &
    ri_tol, &
    rr_tol, &
    rs_tol, &
    rt_tol, &
    thl_tol, &
    w_tol_sqd

  public &
    ! Attempt to Not Use the Following:
! The CLUBB_CAM preprocessor directives are being commented out because this
! code is now also used for WRF-CLUBB.
!#ifdef CLUBB_CAM /* Code for storing pdf_parameter structs in pbuf as array */
!#endif
    init_pdf_params_api, &
    init_precip_fracs_api, &
    init_hydromet_pdf_params, &
    zero_precip_fracs_api, &
    precipitation_fractions, &
    zero_pdf_params_api, &
    init_pdf_implicit_coefs_terms_api, &
    zero_pdf_implicit_coefs_terms_api, &
    calc_derived_params_api, &
    compute_current_date_api, &
    gregorian2julian_day_api, &
    leap_year_api

  public &
    nparams, &
    nu_vertical_res_dep, &
    check_parameters_api, &
    stat_nknd, &
    stat_rknd

  public &
    ! Needed to use the configurable CLUBB flags
    clubb_config_flags_type, &
    set_default_clubb_config_flags_api, &
    initialize_clubb_config_flags_type_api, &
    print_clubb_config_flags_api

  public &
    thlm_sponge_damp_settings,      & ! Variable(s)
    rtm_sponge_damp_settings,       &
    uv_sponge_damp_settings,        &
    wp2_sponge_damp_settings,       &
    wp3_sponge_damp_settings,       &
    up2_vp2_sponge_damp_settings,   &
    thlm_sponge_damp_profile,       &
    rtm_sponge_damp_profile,        &
    uv_sponge_damp_profile,         &
    wp2_sponge_damp_profile,        &
    wp3_sponge_damp_profile,        &
    up2_vp2_sponge_damp_profile,    &
    initialize_tau_sponge_damp_api, & ! Procedure(s)
    finalize_tau_sponge_damp_api
    

  public &
   err_info_type,                   &
   init_err_info_api,               &
   init_default_err_info_api,       &
   set_err_info_values_api,         &
   cleanup_err_info_api


contains


  subroutine advance_clubb_core_api &
            ( gr, nzm, nzt, ngrdcol, & ! In
              l_implemented, dt, fcor, fcor_y, sfc_elevation, & ! In
              hydromet_dim, & ! In
              sclr_dim, sclr_tol, edsclr_dim, sclr_idx, & ! In
              thlm_forcing, rtm_forcing, um_forcing, vm_forcing, & ! In
              sclrm_forcing, edsclrm_forcing, wprtp_forcing, & ! In
              wpthlp_forcing, rtp2_forcing, thlp2_forcing, & ! In
              rtpthlp_forcing, wm_zm, wm_zt, & ! In
              wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, p_sfc, & ! In
              wpsclrp_sfc, wpedsclrp_sfc, & ! In
              upwp_sfc_pert, vpwp_sfc_pert, & ! In
              rtm_ref, thlm_ref, um_ref, vm_ref, ug, vg, & ! In
              p_in_Pa, rho_zm, rho, exner, & ! In
              rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, & ! In
              invrs_rho_ds_zt, thv_ds_zm, thv_ds_zt, & ! In
              l_mix_rat_hm, & ! In
              rfrzm, wphydrometp, & ! In
              wp2hmp, rtphmp_zt, thlphmp_zt, & ! In
              host_dx, host_dy, & ! In
              clubb_params, nu_vert_res_dep, lmin, & ! In
              mixt_frac_max_mag, T0, ts_nudge, & ! In
              rtm_min, rtm_nudge_max_altitude, & ! In
              clubb_config_flags, & ! In
              stats, & ! InOut
              um, vm, upwp, vpwp, up2, vp2, up3, vp3, & ! InOut
              thlm, rtm, wprtp, wpthlp, & ! InOut
              wp2, wp3, rtp2, rtp3, thlp2, thlp3, rtpthlp, & ! InOut
              sclrm, sclrp2, sclrp3, sclrprtp, sclrpthlp, & ! InOut
              wpsclrp, edsclrm, err_info_api, & ! InOut
              rcm, cloud_frac, & ! InOut
              wpthvp, wp2thvp, wp2up, rtpthvp, thlpthvp, & ! InOut
              sclrpthvp, & ! InOut
              wp2rtp, wp2thlp, uprcp, vprcp, rc_coef_zm, wp4, & ! InOut
              wpup2, wpvp2, wp2up2, wp2vp2, ice_supersat_frac, & ! InOut
              um_pert, vm_pert, upwp_pert, vpwp_pert, & ! InOut
              pdf_params, pdf_params_zm, & ! InOut
              pdf_implicit_coefs_terms, & ! InOut
              Kh_zm, Kh_zt, & ! Out
#ifdef CLUBB_CAM
              qclvar, & ! Out
#endif
              thlprcp, wprcp, w_up_in_cloud, w_down_in_cloud, & ! Out
              cloudy_updraft_frac, cloudy_downdraft_frac, & ! Out
              rcm_in_layer, cloud_cover, invrs_tau_zm, & ! Out
              Lscale ) ! Out

    ! Description:
    !   Advance CLUBB for all columns in the batch.

    use advance_clubb_core_module, only : advance_clubb_core

    use pdf_parameter_module, only: &
        implicit_coefs_terms    ! Variable Type(s)

    use model_flags, only: &
        clubb_config_flags_type

    use array_index, only: &
      sclr_idx_type

    implicit none

    !------------------------ Input Variables ------------------------
    integer, intent(in) :: &
      nzm, &
      nzt, &
      ngrdcol   ! Number of grid columns

    type(grid), intent(in) :: &
      gr

      !!! Input Variables
    logical, intent(in) :: &
      l_implemented ! True if CLUBB is being implemented and run in a host model

    real( kind = core_rknd ), intent(in) :: &
      dt  ! Current timestep duration    [s]
      
    real( kind = core_rknd ), intent(in), dimension(ngrdcol) :: &
      fcor, &           ! Coriolis forcing             [s^-1]
      fcor_y, &         ! Nontraditional Coriolis parameter [s^-1]
      sfc_elevation     ! Elevation of ground level    [m AMSL]

    integer, intent(in) :: &
      hydromet_dim,   & ! Total number of hydrometeor species       [#]
      sclr_dim,       & ! Number of passive scalars                 [#]
      edsclr_dim        ! Number of eddy-diff. passive scalars      [#]

    real( kind = core_rknd ), intent(in), dimension(sclr_dim) :: &
      sclr_tol          ! Threshold(s) on the passive scalars  [units vary]

    type (sclr_idx_type), intent(in) :: &
      sclr_idx

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzt) :: &
      thlm_forcing,    & ! theta_l forcing (thermodynamic levels)    [K/s]
      rtm_forcing,     & ! r_t forcing (thermodynamic levels)        [(kg/kg)/s]
      um_forcing,      & ! u wind forcing (thermodynamic levels)     [m/s/s]
      vm_forcing,      & ! v wind forcing (thermodynamic levels)     [m/s/s]
      wm_zt,           & ! w mean wind component on thermo. levels   [m/s]
      rho,             & ! Air density on thermodynamic levels       [kg/m^3]
      rho_ds_zt,       & ! Dry, static density on thermo. levels     [kg/m^3]
      invrs_rho_ds_zt, & ! Inv. dry, static density @ thermo. levs.  [m^3/kg]
      thv_ds_zt,       & ! Dry, base-state theta_v on thermo. levs.  [K]
      rfrzm              ! Total ice-phase water mixing ratio        [kg/kg]

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzm) :: &
      wprtp_forcing,   & ! <w'r_t'> forcing (momentum levels)    [m*K/s^2]
      wpthlp_forcing,  & ! <w'th_l'> forcing (momentum levels)   [m*(kg/kg)/s^2]
      rtp2_forcing,    & ! <r_t'^2> forcing (momentum levels)    [(kg/kg)^2/s]
      thlp2_forcing,   & ! <th_l'^2> forcing (momentum levels)   [K^2/s]
      rtpthlp_forcing, & ! <r_t'th_l'> forcing (momentum levels) [K*(kg/kg)/s]
      wm_zm,           & ! w mean wind component on momentum levels  [m/s]
      rho_zm,          & ! Air density on momentum levels            [kg/m^3]
      rho_ds_zm,       & ! Dry, static density on momentum levels    [kg/m^3]
      invrs_rho_ds_zm, & ! Inv. dry, static density @ momentum levs. [m^3/kg]
      thv_ds_zm          ! Dry, base-state theta_v on momentum levs. [K]

    logical, dimension(hydromet_dim), intent(in) :: &
      l_mix_rat_hm   ! if true, then the quantity is a hydrometeor mixing ratio

    real( kind = core_rknd ), dimension(ngrdcol,nzm, hydromet_dim), intent(in) :: &
      wphydrometp    ! Covariance of w and a hydrometeor   [(m/s) <hm units>]

    real( kind = core_rknd ), dimension(ngrdcol,nzt, hydromet_dim), intent(in) :: &
      wp2hmp,      & ! Third moment: <w'^2> * <hydro.'>    [(m/s)^2 <hm units>]
      rtphmp_zt,   & ! Covariance of rt and a hydrometeor  [(kg/kg) <hm units>]
      thlphmp_zt     ! Covariance of thl and a hydrometeor [K <hm units>]

    real( kind = core_rknd ), intent(in), dimension(ngrdcol) ::  &
      wpthlp_sfc,   & ! w' theta_l' at surface   [(m K)/s]
      wprtp_sfc,    & ! w' r_t' at surface       [(kg m)/( kg s)]
      upwp_sfc,     & ! u'w' at surface          [m^2/s^2]
      vpwp_sfc,     & ! v'w' at surface          [m^2/s^2]
      p_sfc           ! Pressure at surface      [Pa]

    ! Passive scalar variables
    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzt,sclr_dim) :: &
      sclrm_forcing    ! Passive scalar forcing         [{units vary}/s]

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,sclr_dim) ::  &
      wpsclrp_sfc      ! Scalar flux at surface         [{units vary} m/s]

    real( kind = core_rknd ), intent(in), dimension(ngrdcol) :: &
      upwp_sfc_pert, & ! pertubed u'w' at surface    [m^2/s^2]
      vpwp_sfc_pert    ! pertubed v'w' at surface    [m^2/s^2]

    ! Eddy passive scalar variables
    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzt,edsclr_dim) :: &
      edsclrm_forcing  ! Eddy passive scalar forcing    [{units vary}/s]

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,edsclr_dim) ::  &
      wpedsclrp_sfc    ! Eddy-Scalar flux at surface    [{units vary} m/s]

    ! Reference profiles (used for nudging, sponge damping, and Coriolis effect)
    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) ::  &
      rtm_ref,  & ! Initial total water mixing ratio             [kg/kg]
      thlm_ref, & ! Initial liquid water potential temperature   [K]
      um_ref,   & ! Initial u wind; Michael Falk                 [m/s]
      vm_ref,   & ! Initial v wind; Michael Falk                 [m/s]
      ug,       & ! u geostrophic wind                           [m/s]
      vg          ! v geostrophic wind                           [m/s]

    ! Host model horizontal grid spacing, if part of host model.
    real( kind = core_rknd ), intent(in), dimension(ngrdcol) :: &
      host_dx,  & ! East-West horizontal grid spacing     [m]
      host_dy     ! North-South horizontal grid spacing   [m]

    real( kind = core_rknd ), dimension(ngrdcol,nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

    type(nu_vertical_res_dep), intent(in) :: &
      nu_vert_res_dep    ! Vertical resolution dependent nu values

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
      lmin, & ! Min. value for the length scale    [m]
      mixt_frac_max_mag

    real( kind = core_rknd ), intent(in) :: &
      T0, & ! Reference temperature (usually 300)  [K]
      ts_nudge, & ! Timescale of u/v nudging             [s]
      rtm_min, & ! Value below which rtm will be nudged [kg/kg]
      rtm_nudge_max_altitude ! Highest altitude at which to nudge rtm [m]

    type( clubb_config_flags_type ), intent(in) :: &
      clubb_config_flags ! Derived type holding all configurable CLUBB flags

    type(stats_type), intent(inout) :: &
      stats

    ! These are prognostic or are planned to be in the future
    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,nzt) ::  &
      um,      & ! u mean wind component (thermodynamic levels)   [m/s]
      vm,      & ! v mean wind component (thermodynamic levels)   [m/s]
      up3,     & ! u'^3 (thermodynamic levels)                    [m^3/s^3]
      vp3,     & ! v'^3 (thermodynamic levels)                    [m^3/s^3]
      rtm,     & ! total water mixing ratio, r_t (thermo. levels) [kg/kg]
      thlm,    & ! liq. water pot. temp., th_l (thermo. levels)   [K]
      rtp3,    & ! r_t'^3 (thermodynamic levels)                  [(kg/kg)^3]
      thlp3,   & ! th_l'^3 (thermodynamic levels)                 [K^3]
      wp3        ! w'^3 (thermodynamic levels)                    [m^3/s^3]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,nzm) ::  &
      upwp,    & ! u'w' (momentum levels)                         [m^2/s^2]
      vpwp,    & ! v'w' (momentum levels)                         [m^2/s^2]
      up2,     & ! u'^2 (momentum levels)                         [m^2/s^2]
      vp2,     & ! v'^2 (momentum levels)                         [m^2/s^2]
      wprtp,   & ! w' r_t' (momentum levels)                      [(kg/kg) m/s]
      wpthlp,  & ! w' th_l' (momentum levels)                     [(m/s) K]
      rtp2,    & ! r_t'^2 (momentum levels)                       [(kg/kg)^2]
      thlp2,   & ! th_l'^2 (momentum levels)                      [K^2]
      rtpthlp, & ! r_t' th_l' (momentum levels)                   [(kg/kg) K]
      wp2        ! w'^2 (momentum levels)                         [m^2/s^2]

    ! Passive scalar variables
    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,nzt,sclr_dim) :: &
      sclrm,     & ! Passive scalar mean (thermo. levels) [units vary]
      sclrp3       ! sclr'^3 (thermodynamic levels)       [{units vary}^3]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,nzm,sclr_dim) :: &
      wpsclrp,   & ! w'sclr' (momentum levels)            [{units vary} m/s]
      sclrp2,    & ! sclr'^2 (momentum levels)            [{units vary}^2]
      sclrprtp,  & ! sclr'rt' (momentum levels)           [{units vary} (kg/kg)]
      sclrpthlp    ! sclr'thl' (momentum levels)          [{units vary} K]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,nzt) ::  &
      p_in_Pa, & ! Air pressure (thermodynamic levels)       [Pa]
      exner      ! Exner function (thermodynamic levels)     [-]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,nzt) ::  &
      rcm,        & ! cloud water mixing ratio, r_c (thermo. levels) [kg/kg]
      cloud_frac, & ! cloud fraction (thermodynamic levels)          [-]
      wp2thvp,    & ! < w'^2 th_v' > (thermodynamic levels)          [m^2/s^2 K]
      wp2up         ! < w'^2 u' > (thermodynamic levels)             [m^3/s^3]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,nzm) ::  &
      wpthvp,     & ! < w' th_v' > (momentum levels)                 [kg/kg K]
      rtpthvp,    & ! < r_t' th_v' > (momentum levels)               [kg/kg K]
      thlpthvp      ! < th_l' th_v' > (momentum levels)              [K^2]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,nzm,sclr_dim) :: &
      sclrpthvp     ! < sclr' th_v' > (momentum levels)   [units vary]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,nzt) ::  &
      wp2rtp,            & ! w'^2 rt' (thermodynamic levels)      [m^2/s^2 kg/kg]
      wp2thlp,           & ! w'^2 thl' (thermodynamic levels)     [m^2/s^2 K]
      wpup2,             & ! w'u'^2 (thermodynamic levels)        [m^3/s^3]
      wpvp2,             & ! w'v'^2 (thermodynamic levels)        [m^3/s^3]
      ice_supersat_frac    ! ice cloud fraction (thermo. levels)  [-]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,nzm) ::  &
      uprcp,             & ! < u' r_c' > (momentum levels)        [(m/s)(kg/kg)]
      vprcp,             & ! < v' r_c' > (momentum levels)        [(m/s)(kg/kg)]
      rc_coef_zm,        & ! Coef of X'r_c' in Eq. (34) (m-levs.) [K/(kg/kg)]
      wp4,               & ! w'^4 (momentum levels)               [m^4/s^4]
      wp2up2,            & ! w'^2 u'^2 (momentum levels)          [m^4/s^4]
      wp2vp2               ! w'^2 v'^2 (momentum levels)          [m^4/s^4]

    ! Variables used to track perturbed version of winds.
    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,nzt) :: &
      um_pert,   & ! perturbed <u>       [m/s]
      vm_pert      ! perturbed <v>       [m/s]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,nzm) :: &
      upwp_pert, & ! perturbed <u'w'>    [m^2/s^2]
      vpwp_pert    ! perturbed <v'w'>    [m^2/s^2]

    type(pdf_parameter), intent(inout) :: &
      pdf_params,    & ! PDF parameters (thermodynamic levels)    [units vary]
      pdf_params_zm    ! PDF parameters on momentum levels        [units vary]

    type(implicit_coefs_terms), intent(inout) :: &
      pdf_implicit_coefs_terms    ! Implicit coefs / explicit terms [units vary]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,nzt,edsclr_dim) :: &
    edsclrm   ! Eddy passive scalar mean (thermo. levels)   [units vary]

    type(err_info_type), intent(inout) :: &
      err_info_api          ! err_info struct containing err_code and err_header

    !------------------------ Output Variables ------------------------
    real( kind = core_rknd ), intent(out), dimension(ngrdcol,nzt) ::  &
      rcm_in_layer, & ! rcm in cloud layer                              [kg/kg]
      cloud_cover     ! cloud cover                                     [-]

    ! Variables that need to be output for use in host models
    real( kind = core_rknd ), intent(out), dimension(ngrdcol,nzm) ::  &
      wprcp,                 & ! w'r_c' (momentum levels)              [(kg/kg) m/s]
      invrs_tau_zm             ! One divided by tau on zm levels       [1/s]

    real( kind = core_rknd ), intent(out), dimension(ngrdcol,nzt) ::  &
      w_up_in_cloud,         & ! Average cloudy updraft velocity       [m/s]
      w_down_in_cloud,       & ! Average cloudy downdraft velocity     [m/s]
      cloudy_updraft_frac,   & ! cloudy updraft fraction               [-]
      cloudy_downdraft_frac    ! cloudy downdraft fraction             [-]

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) :: &
      Kh_zt    ! Eddy diffusivity coefficient on thermodynamic levels   [m^2/s]

    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(out) :: &
      Kh_zm    ! Eddy diffusivity coefficient on momentum levels        [m^2/s]

#ifdef CLUBB_CAM
    real( kind = core_rknd), intent(out), dimension(ngrdcol,nzt) :: &
      qclvar        ! cloud water variance
#endif

    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(out) :: &
      thlprcp    ! thl'rc'              [K kg/kg]

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) :: &
      Lscale     ! Length scale         [m]

    !$acc data copyin( gr, gr%zm, gr%zt, gr%dzm, gr%dzt, gr%invrs_dzt, gr%invrs_dzm, &
    !$acc              gr%weights_zt2zm, gr%weights_zm2zt, &
    !$acc              nu_vert_res_dep, nu_vert_res_dep%nu2, nu_vert_res_dep%nu9, &
    !$acc              nu_vert_res_dep%nu1, nu_vert_res_dep%nu8, nu_vert_res_dep%nu10, &
    !$acc              nu_vert_res_dep%nu6, &
    !$acc              sclr_idx, clubb_params, lmin, mixt_frac_max_mag, &
    !$acc              fcor, sfc_elevation, thlm_forcing, rtm_forcing, um_forcing, &
    !$acc              vm_forcing, wprtp_forcing, wpthlp_forcing, rtp2_forcing, thlp2_forcing, &
    !$acc              rtpthlp_forcing, wm_zm, wm_zt, rho_zm, rho, rho_ds_zm, rho_ds_zt, &
    !$acc              invrs_rho_ds_zm, invrs_rho_ds_zt, thv_ds_zm, thv_ds_zt, rfrzm, &
    !$acc              wpthlp_sfc, &
    !$acc              wprtp_sfc, upwp_sfc, vpwp_sfc, p_sfc, & 
    !$acc              upwp_sfc_pert, vpwp_sfc_pert, rtm_ref, thlm_ref, um_ref, &
    !$acc              vm_ref, ug, vg, host_dx, host_dy, &
    !$acc              pdf_params, pdf_params_zm, err_info_api, err_info_api%err_header ) &
    !$acc        copy( um, upwp, vm, vpwp, up2, vp2, up3, vp3, rtm, wprtp, thlm, wpthlp, rtp2, &
    !$acc              rtp3, thlp2, thlp3, rtpthlp, wp2, wp3, &
    !$acc              p_in_Pa, exner, rcm, cloud_frac, wpthvp, wp2thvp, wp2up, &
    !$acc              rtpthvp, thlpthvp, wp2rtp, wp2thlp, uprcp, vprcp, rc_coef_zm, &
    !$acc              wp4, wpup2, wpvp2, wp2up2, wp2vp2, ice_supersat_frac, um_pert, &
    !$acc              vm_pert, upwp_pert, vpwp_pert, &
    !$acc              pdf_params%w_1, pdf_params%w_2, &
    !$acc              pdf_params%varnce_w_1, pdf_params%varnce_w_2, &
    !$acc              pdf_params%rt_1, pdf_params%rt_2, &
    !$acc              pdf_params%varnce_rt_1, pdf_params%varnce_rt_2,  &
    !$acc              pdf_params%thl_1, pdf_params%thl_2, &
    !$acc              pdf_params%varnce_thl_1, pdf_params%varnce_thl_2, &
    !$acc              pdf_params%corr_w_rt_1, pdf_params%corr_w_rt_2,  &
    !$acc              pdf_params%corr_w_thl_1, pdf_params%corr_w_thl_2, &
    !$acc              pdf_params%corr_rt_thl_1, pdf_params%corr_rt_thl_2,&
    !$acc              pdf_params%alpha_thl, pdf_params%alpha_rt, &
    !$acc              pdf_params%crt_1, pdf_params%crt_2, pdf_params%cthl_1, &
    !$acc              pdf_params%cthl_2, pdf_params%chi_1, &
    !$acc              pdf_params%chi_2, pdf_params%stdev_chi_1, &
    !$acc              pdf_params%stdev_chi_2, pdf_params%stdev_eta_1, &
    !$acc              pdf_params%stdev_eta_2, pdf_params%covar_chi_eta_1, &
    !$acc              pdf_params%covar_chi_eta_2, pdf_params%corr_w_chi_1, &
    !$acc              pdf_params%corr_w_chi_2, pdf_params%corr_w_eta_1, &
    !$acc              pdf_params%corr_w_eta_2, pdf_params%corr_chi_eta_1, &
    !$acc              pdf_params%corr_chi_eta_2, pdf_params%rsatl_1, &
    !$acc              pdf_params%rsatl_2, pdf_params%rc_1, pdf_params%rc_2, &
    !$acc              pdf_params%cloud_frac_1, pdf_params%cloud_frac_2,  &
    !$acc              pdf_params%mixt_frac, pdf_params%ice_supersat_frac_1, &
    !$acc              pdf_params%ice_supersat_frac_2, &
    !$acc              pdf_params_zm%w_1, pdf_params_zm%w_2, &
    !$acc              pdf_params_zm%varnce_w_1, pdf_params_zm%varnce_w_2, &
    !$acc              pdf_params_zm%rt_1, pdf_params_zm%rt_2, &
    !$acc              pdf_params_zm%varnce_rt_1, pdf_params_zm%varnce_rt_2,  &
    !$acc              pdf_params_zm%thl_1, pdf_params_zm%thl_2, &
    !$acc              pdf_params_zm%varnce_thl_1, pdf_params_zm%varnce_thl_2, &
    !$acc              pdf_params_zm%corr_w_rt_1, pdf_params_zm%corr_w_rt_2,  &
    !$acc              pdf_params_zm%corr_w_thl_1, pdf_params_zm%corr_w_thl_2, &
    !$acc              pdf_params_zm%corr_rt_thl_1, pdf_params_zm%corr_rt_thl_2,&
    !$acc              pdf_params_zm%alpha_thl, pdf_params_zm%alpha_rt, &
    !$acc              pdf_params_zm%crt_1, pdf_params_zm%crt_2, pdf_params_zm%cthl_1, &
    !$acc              pdf_params_zm%cthl_2, pdf_params_zm%chi_1, &
    !$acc              pdf_params_zm%chi_2, pdf_params_zm%stdev_chi_1, &
    !$acc              pdf_params_zm%stdev_chi_2, pdf_params_zm%stdev_eta_1, &
    !$acc              pdf_params_zm%stdev_eta_2, pdf_params_zm%covar_chi_eta_1, &
    !$acc              pdf_params_zm%covar_chi_eta_2, pdf_params_zm%corr_w_chi_1, &
    !$acc              pdf_params_zm%corr_w_chi_2, pdf_params_zm%corr_w_eta_1, &
    !$acc              pdf_params_zm%corr_w_eta_2, pdf_params_zm%corr_chi_eta_1, &
    !$acc              pdf_params_zm%corr_chi_eta_2, pdf_params_zm%rsatl_1, &
    !$acc              pdf_params_zm%rsatl_2, pdf_params_zm%rc_1, pdf_params_zm%rc_2, &
    !$acc              pdf_params_zm%cloud_frac_1, pdf_params_zm%cloud_frac_2,  &
    !$acc              pdf_params_zm%mixt_frac, pdf_params_zm%ice_supersat_frac_1, &
    !$acc              pdf_params_zm%ice_supersat_frac_2 ) &
    !$acc     copyout( rcm_in_layer, cloud_cover, wprcp, w_up_in_cloud, w_down_in_cloud, &
    !$acc              cloudy_updraft_frac, cloudy_downdraft_frac, invrs_tau_zm, Kh_zt, &
    !$acc              Kh_zm, thlprcp, Lscale, err_info_api%err_code )

    !$acc data if( sclr_dim > 0 ) &
    !$acc      copyin( sclr_tol, sclrm_forcing, wpsclrp_sfc ) &
    !$acc        copy( sclrm, wpsclrp, sclrp2, sclrp3, sclrprtp, sclrpthlp, sclrpthvp )

    !$acc data if( edsclr_dim > 0 ) &
    !$acc      copyin( wpedsclrp_sfc, edsclrm_forcing ) &
    !$acc        copy( edsclrm )

    !$acc data if( hydromet_dim > 0 ) &
    !$acc      copyin( wphydrometp, wp2hmp, rtphmp_zt, thlphmp_zt, &
    !$acc              l_mix_rat_hm )

#ifdef CLUBB_CAM
    !$acc data copyout( qclvar )
#endif

    call advance_clubb_core( gr, nzm, nzt, ngrdcol, &         ! intent(in)
      l_implemented, dt, fcor, fcor_y, sfc_elevation, &       ! intent(in)
      hydromet_dim, &                                         ! intent(in)
      sclr_dim, sclr_tol, edsclr_dim, sclr_idx, &             ! intent(in)
      thlm_forcing, rtm_forcing, um_forcing, vm_forcing, &    ! intent(in)
      sclrm_forcing, edsclrm_forcing, wprtp_forcing, &        ! intent(in)
      wpthlp_forcing, rtp2_forcing, thlp2_forcing, &          ! intent(in)
      rtpthlp_forcing, wm_zm, wm_zt, &                        ! intent(in)
      wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, p_sfc, &     ! intent(in)
      wpsclrp_sfc, wpedsclrp_sfc, &                           ! intent(in)
      upwp_sfc_pert, vpwp_sfc_pert, &                         ! intent(in)
      rtm_ref, thlm_ref, um_ref, vm_ref, ug, vg, &            ! Intent(in)
      p_in_Pa, rho_zm, rho, exner, &                          ! intent(in)
      rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &                ! intent(in)
      invrs_rho_ds_zt, thv_ds_zm, thv_ds_zt, &                ! intent(in) 
      l_mix_rat_hm,&                                          ! intent(in)
      rfrzm, &                                                ! intent(in)
      wphydrometp, wp2hmp, rtphmp_zt, thlphmp_zt, &           ! intent(in)
      host_dx, host_dy, &                                     ! intent(in)
      clubb_params, nu_vert_res_dep, lmin, &                  ! intent(in)
      mixt_frac_max_mag, T0, ts_nudge, &                      ! Intent(in)
      rtm_min, rtm_nudge_max_altitude, &                      ! Intent(in)
      clubb_config_flags, &                                   ! intent(in)
      stats,         &                                        ! intent(inout)
      um, vm, upwp, vpwp, up2, vp2, up3, vp3, &               ! intent(inout)
      thlm, rtm, wprtp, wpthlp, &                             ! intent(inout)
      wp2, wp3, rtp2, rtp3, thlp2, thlp3, rtpthlp, &          ! intent(inout)
      sclrm,   &
      sclrp2, sclrp3, sclrprtp, sclrpthlp, &                  ! intent(inout)
      wpsclrp, edsclrm, &                                     ! intent(inout)
      rcm, cloud_frac, &                                      ! intent(inout)
      wpthvp, wp2thvp, wp2up, rtpthvp, thlpthvp, &            ! intent(inout)
      sclrpthvp, &                                            ! intent(inout)
      wp2rtp, wp2thlp, uprcp, vprcp, rc_coef_zm, wp4, &       ! intent(inout)
      wpup2, wpvp2, wp2up2, wp2vp2, ice_supersat_frac, &      ! intent(inout)
      um_pert, vm_pert, upwp_pert, vpwp_pert, &               ! intent(inout)
      pdf_params, pdf_params_zm, &                            ! intent(inout)
      pdf_implicit_coefs_terms, &                             ! intent(inout)
      err_info_api, &                                         ! intent(inout)
      Kh_zm, Kh_zt, &                                         ! intent(out)
#ifdef CLUBB_CAM
               qclvar, &                                      ! intent(out)
#endif
      thlprcp, wprcp, w_up_in_cloud, w_down_in_cloud, &       ! intent(out)
      cloudy_updraft_frac, cloudy_downdraft_frac, &           ! intent(out)
      rcm_in_layer, cloud_cover, invrs_tau_zm, &              ! intent(out)
      Lscale )

    !$acc end data
    !$acc end data
    !$acc end data
    !$acc end data

#ifdef CLUBB_CAM
    !$acc end data
#endif

    if ( clubb_at_least_debug_level_api( 0 ) ) then
      if ( any(err_info_api%err_code == clubb_fatal_error) ) then
        write(fstderr,*) err_info_api%err_header_global
        write(fstderr,*) "Fatal error advance_clubb_core in API module"
      end if
    end if

  end subroutine advance_clubb_core_api

  subroutine calculate_thlp2_rad_api( ngrdcol, nzm, nzt, gr, &
                                      rcm, thlprcp, radht, clubb_params, &
                                      thlp2_forcing )

    use clubb_precision, only: &
      core_rknd                     ! Constant(s)

    use parameter_indices, only: &
      nparams

    use grid_class, only: &
      grid

    use advance_helper_module, only : &
      calculate_thlp2_rad

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      ngrdcol, &
      nzm, &
      nzt

    type( grid ), intent(in) :: &
      gr

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      rcm, &             ! Cloud water mixing ratio on momentum grid      [kg/kg]
      radht              ! SW + LW heating rate (on momentum grid)        [K/s]

    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: &
      thlprcp               ! thl'rc'                                        [K kg/kg]

    real( kind = core_rknd ), dimension(ngrdcol,nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(inout) :: &
      thlp2_forcing         ! <th_l'^2> forcing (momentum levels)            [K^2/s]

    !----------------------------------------------------------------------

    !$acc data copyin( rcm, thlprcp, radht, clubb_params ) &
    !$acc        copy( thlp2_forcing )

    call calculate_thlp2_rad( ngrdcol, nzm, nzt, gr,              & ! intent(in)
                              rcm, thlprcp, radht, clubb_params,  & ! intent(in)
                              thlp2_forcing )                       ! intent(inout)

    !$acc end data

    return

    end subroutine calculate_thlp2_rad_api

  !================================================================================================
  ! gregorian2julian_day - Computes the number of days since 1 January 4713 BC.
  !================================================================================================

  integer function gregorian2julian_day_api( &
    day, month, year )

    use calendar, only : gregorian2julian_day

    implicit none

    ! Input Variables
    integer, intent(in) ::  &
      day,        & ! Gregorian Calendar Day for given Month        [dd]
      month,      & ! Gregorian Calendar Month for given Year       [mm]
      year          ! Gregorian Calendar Year                       [yyyy]

    gregorian2julian_day_api = gregorian2julian_day( &
      day, month, year )
  end function gregorian2julian_day_api

  !================================================================================================
  ! leap_year - Determines if the given year is a leap year.
  !================================================================================================

  logical function leap_year_api( &
    year )

    use calendar, only : leap_year

    implicit none

    ! External
    intrinsic :: mod

    ! Input Variable(s)
    integer, intent(in) :: year ! Gregorian Calendar Year [yyyy]

    leap_year_api = leap_year( &
      year )
  end function leap_year_api


  !================================================================================================
  ! setup_grid - This subroutine sets up the CLUBB vertical grid.
  !================================================================================================
  subroutine setup_grid_api( nzmax, ngrdcol, sfc_elevation, l_implemented, & ! In
                             l_ascending_grid, grid_type, deltaz, zm_init, zm_top, & ! In
                             momentum_heights, & ! In
                             thermodynamic_heights, & ! In
                             gr, err_info_api ) ! InOut

    ! Description:
    !   Set up the vertical grids for all columns in the batch.

    use grid_class, only: &
        grid, & ! Type
        setup_grid

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nzmax,  & ! Number of vertical levels in grid      [#]
      ngrdcol   ! Number of grid columns                 [#]

    type (grid), intent(inout) :: gr

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
      sfc_elevation  ! Elevation of ground level    [m AMSL]

    logical, intent(in) :: l_implemented
    
    logical, intent(in) :: l_ascending_grid

    integer, intent(in) :: grid_type

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
      deltaz,   & ! Vertical grid spacing                  [m]
      zm_init,  & ! Initial grid altitude (momentum level) [m]
      zm_top      ! Maximum grid altitude (momentum level) [m]

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzmax) :: &
      momentum_heights      ! Momentum level altitudes (input)      [m]

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzmax-1) :: &
      thermodynamic_heights ! Thermodynamic level altitudes (input) [m]

    ! ------------------- InOut Variables -------------------
    type(err_info_type), intent(inout) :: &
      err_info_api          ! err_info struct containing err_code and err_header

    ! ------------------- Begin Code -------------------

    call setup_grid( nzmax, ngrdcol, sfc_elevation, l_implemented,         & ! intent(in)
                     l_ascending_grid, grid_type, deltaz, zm_init, zm_top, & ! intent(in)
                     momentum_heights, thermodynamic_heights,              & ! intent(in)
                     gr, err_info_api )                                      ! intent(inout)

    if ( any(err_info_api%err_code == clubb_fatal_error) ) then
      write(fstderr, *) err_info_api%err_header_global
      write(fstderr, *) "Error in CLUBB calling setup_grid"
    endif

  end subroutine setup_grid_api

  !================================================================================================
  ! lin_interpolate_two_points - Computes a linear interpolation of the value of a variable.
  !================================================================================================

  function lin_interpolate_two_points_api( &
    height_int, height_high, height_low, &
    var_high, var_low )

    use interpolation, only : lin_interpolate_two_points

    implicit none

    real( kind = core_rknd ), intent(in) :: &
      height_int,  & ! Height to be interpolated to     [m]
      height_high, & ! Height above the interpolation   [m]
      height_low,  & ! Height below the interpolation   [m]
      var_high,    & ! Variable above the interpolation [units vary]
      var_low        ! Variable below the interpolation [units vary]

    ! Output Variables
    real( kind = core_rknd ) :: lin_interpolate_two_points_api

    lin_interpolate_two_points_api = lin_interpolate_two_points( &
      height_int, height_high, height_low, &
      var_high, var_low )

  end function lin_interpolate_two_points_api


  !================================================================================================
  ! calc_derived_params_api - Adjusts background eddy diffusivity based
  ! on vertical grid spacing
  !                    and calculates 
  !================================================================================================

  subroutine calc_derived_params_api( gr, ngrdcol, grid_type, deltaz,  & ! Intent(in)
                                                 clubb_params,                    & ! Intent(in)
                                                 l_prescribed_avg_deltaz,         & ! Intent(in)
                                                 nu_vert_res_dep, lmin,           & ! Intent(inout)
                                                 mixt_frac_max_mag )                ! Intent(inout)

    ! Description:
    !   Calculate grid-dependent parameters for all columns in the batch.

    use parameters_tunable, only : calc_derived_params

    implicit none

    ! Input Variables

    ! Grid definition
    type(grid), intent(in) :: &
      gr

    integer, intent(in) :: &
      ngrdcol

    ! If CLUBB is running on it's own, this option determines
    ! if it is using:
    ! 1) an evenly-spaced grid,
    ! 2) a stretched (unevenly-spaced) grid entered on the
    !    thermodynamic grid levels (with momentum levels set
    !    halfway between thermodynamic levels), or
    ! 3) a stretched (unevenly-spaced) grid entered on the
    !    momentum grid levels (with thermodynamic levels set
    !    halfway between momentum levels).
    integer, intent(in) :: grid_type

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) ::  &
      deltaz  ! Change per height level        [m]

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nparams) :: & 
      clubb_params  ! Tuneable model parameters      [-]

    logical, intent(in) :: &
      l_prescribed_avg_deltaz ! used in calc_derived_params_api. If .true., avg_deltaz = deltaz

    ! Output Variables
    type(nu_vertical_res_dep), intent(inout) :: &
      nu_vert_res_dep    ! Vertical resolution dependent nu values

    real( kind = core_rknd ), dimension(ngrdcol), intent(inout) :: &
      mixt_frac_max_mag, &
      lmin ! Min. value for the length scale    [m]

    call calc_derived_params( gr, ngrdcol, grid_type, deltaz, & ! Intent(in)
                               clubb_params,                   & ! Intent(in)
                               l_prescribed_avg_deltaz,        & ! Intent(in)
                               nu_vert_res_dep, lmin,          & ! Intent(inout)
                               mixt_frac_max_mag )               ! Intent(inout)

  end subroutine calc_derived_params_api

  !================================================================================================
  ! init_pdf_params - allocates arrays for pdf_params
  !================================================================================================
  subroutine init_pdf_params_api( nz, ngrdcol, pdf_params )
  
    use pdf_parameter_module, only : init_pdf_params
    
    implicit none

    ! Input Variable(s)
    integer, intent(in) :: &
      nz, &   ! Number of vertical grid levels    [-]
      ngrdcol ! Number of grid columns            [-]

    ! Output Variable(s)
    type(pdf_parameter), intent(out) :: &
      pdf_params    ! PDF parameters            [units vary]
    
    call init_pdf_params( nz, ngrdcol, & ! intent(in)
                          pdf_params ) ! intent(out)
    
  end subroutine init_pdf_params_api


  !================================================================================================
  ! calculate_spurious_source - Checks whether there is conservation within the column.
  !================================================================================================
  function calculate_spurious_source_api ( &
    integral_after, integral_before, &
    flux_top, flux_sfc, &
    integral_forcing, dt ) result( spurious_source )

    use numerical_check, only : calculate_spurious_source

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      integral_after, &   ! Vertically-integrated quantity after dt time  [units vary]
      integral_before, &  ! Vertically-integrated quantity before dt time [units vary]
      flux_top, &         ! Total flux at the top of the domain           [units vary]
      flux_sfc, &         ! Total flux at the bottom of the domain        [units vary]
      integral_forcing, & ! Vertically-integrated forcing                 [units vary]
      dt                  ! Timestep size                                 [s]

    ! Return Variable
    real( kind = core_rknd ) :: spurious_source ! [units vary]

    spurious_source = calculate_spurious_source( &
      integral_after, integral_before, &
      flux_top, flux_sfc, &
      integral_forcing, dt )

  end function calculate_spurious_source_api

  !================================================================================================
  ! update_xp2_mc - Calculates the effects of rain evaporation on rtp2 and thlp2
  !================================================================================================
  subroutine update_xp2_mc_api( gr, nzm, nzt, ngrdcol, dt, cloud_frac, rcm, rvm, & ! In
                                thlm, & ! In
                                wm, exner, rrm_evap, pdf_params, & ! In
                                rtp2_mc, thlp2_mc, wprtp_mc, wpthlp_mc, & ! InOut
                                rtpthlp_mc ) ! InOut

    ! Description:
    !   Update second-moment microphysics tendencies for all columns in the batch.

    use advance_xp2_xpyp_module, only: &
        update_xp2_mc

     ! Type

    implicit none

    ! -------------------- Input Variables --------------------
    integer, intent(in) :: &
      nzm,     & ! Momentum points in the Vertical        [-]
      nzt,     & ! Thermodynamic points in the Vertical        [-]
      ngrdcol

    type(grid), intent(in) :: gr

    real( kind = core_rknd ), intent(in) :: dt ! Model timestep        [s]

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      cloud_frac, &       !Cloud fraction                        [-]
      rcm, &              !Cloud water mixing ratio              [kg/kg]
      rvm, &              !Vapor water mixing ratio              [kg/kg]
      thlm, &             !Liquid potential temperature          [K]
      wm, &               !Mean vertical velocity                [m/s]
      exner, &            !Exner function                        [-]
      rrm_evap         !Evaporation of rain                   [kg/kg/s]
                          !It is expected that this variable is negative, as
                          !that is the convention in Morrison microphysics

    type(pdf_parameter), intent(in) :: &
      pdf_params ! PDF parameters

    ! -------------------- Input/Output Variables --------------------
    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(inout) :: &
      rtp2_mc, &    !Tendency of <rt'^2> due to evaporation   [(kg/kg)^2/s]
      thlp2_mc, &   !Tendency of <thl'^2> due to evaporation  [K^2/s]
      wprtp_mc, &   !Tendency of <w'rt'> due to evaporation   [m*(kg/kg)/s^2]
      wpthlp_mc, &  !Tendency of <w'thl'> due to evaporation  [m*K/s^2] 
      rtpthlp_mc    !Tendency of <rt'thl'> due to evaporation [K*(kg/kg)/s]

    call update_xp2_mc( gr, nzm, nzt, ngrdcol, dt, cloud_frac, rcm, rvm, thlm,        & ! intent(in)
                        wm, exner, rrm_evap, pdf_params,        & ! intent(in)
                        rtp2_mc, thlp2_mc, wprtp_mc, wpthlp_mc,    & ! intent(inout)
                        rtpthlp_mc ) ! intent(inout)
    return
  end subroutine update_xp2_mc_api


end module clubb_api_module
