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


  use mt95, only : &
    assignment( = ), &
    genrand_state, & ! Internal representation of the RNG state.
    genrand_srepr, & ! Public representation of the RNG state. Should be used to save the RNG state
    genrand_intg, &
    genrand_init_api => genrand_init

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

  use error_code, only: &
      clubb_at_least_debug_level,  & ! Procedure
      err_code,                    & ! Error Indicator
      clubb_no_error,              & ! Constants
      clubb_fatal_error

  use hydromet_pdf_parameter_module, only : &
    hydromet_pdf_parameter, &
    precipitation_fractions

  use model_flags, only : &
      clubb_config_flags_type, & ! Type
      iiPDF_ADG1, &
      iiPDF_new_hybrid, &
      ipdf_pre_advance_fields, &
      ipdf_post_advance_fields, &
      l_use_boussinesq    ! Use Boussinesq form of predictive equations (default is Anelastic).

  use parameters_tunable, only : &
    params_list,         & ! Variable(s)
    nu_vertical_res_dep    ! Type(s)

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
    ilambda0_stability_coef, imult_coef, itaumin, itaumax, imu, &
    iLscale_mu_coef, iLscale_pert_coef, ialpha_corr, iSkw_denom_coef, &
    ic_K10, ic_K10h, ithlp2_rad_coef, ithlp2_rad_cloud_frac_thresh, &
    iup2_sfc_coef, iSkw_max_mag, iC_invrs_tau_bkgnd, &
    iC_invrs_tau_sfc, iC_invrs_tau_shear, iC_invrs_tau_N2, &
    iC_invrs_tau_N2_wp2, iC_invrs_tau_N2_xp2, iC_invrs_tau_N2_wpxp, &
    iC_invrs_tau_N2_clear_wp3, ialtitude_threshold, irtp2_clip_coef, &
    iRichardson_num_min, iRichardson_num_max, iwpxp_Ri_exp, &
    ia3_coef_min, ia_const, iCx_min, iCx_max, ibv_efold, iz_displace


  use pdf_parameter_module, only : &
! The CLUBB_CAM preprocessor directives are being commented out because this
! code is now also used for WRF-CLUBB.
!#ifdef CLUBB_CAM /* Code for storing pdf_parameter structs in pbuf as array */
    num_pdf_params, &
!#endif
    pdf_parameter, &
    implicit_coefs_terms

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
    up2_vp2_sponge_damp_profile

  use stat_file_module, only : &
    clubb_i, &    ! Used to output multiple columns
    clubb_j       ! The indices must not exceed nlon (for i) or nlat (for j).

  use stats_rad_zm_module, only : &
    nvarmax_rad_zm ! Maximum variables allowed

  use stats_rad_zt_module, only : &
    nvarmax_rad_zt  ! Maximum variables allowed

  use stats_zm_module, only : &
    nvarmax_zm ! Maximum variables allowed

  use stats_zt_module, only : &
    nvarmax_zt ! Maximum variables allowed

  use stats_sfc_module, only : &
    nvarmax_sfc

  use stats_variables, only: & 
    stats_metadata_type

  use grid_class, only: grid ! Type
  
  use grid_class, only: &
    ! Interpolate momentum level variables to thermodynamic levels 
    zm2zt_api => zm2zt, &

    ! Interpolate thermodynamic level variables to momentum levels 
    zt2zm_api => zt2zm

  use saturation, only: &
    ! Used to compute the saturation mixing ratio of liquid water.
    sat_mixrat_liq_api => sat_mixrat_liq

  use T_in_K_module, only : &
    ! Calculates absolute temperature from liquid water potential
    ! temperature.  (Does not include ice.)
    thlm2T_in_K_api => thlm2T_in_K, &

    ! Calculates liquid water potential temperature from absolute temperature 
    T_in_K2thlm_api => T_in_K2thlm

  use advance_clubb_core_module, only: &
    check_clubb_settings_api => check_clubb_settings

  use parameters_model, only: &
    setup_parameters_model_api =>  setup_parameters_model

  use fill_holes, only : &
    fill_holes_driver_api => fill_holes_driver, &
    fill_holes_hydromet_api => fill_holes_hydromet, &
    fill_holes_vertical_api => fill_holes_vertical

  use stats_rad_zm_module, only : &
    stats_init_rad_zm_api => stats_init_rad_zm

  use stats_rad_zt_module, only : &
   stats_init_rad_zt_api => stats_init_rad_zt

  use stats_zm_module, only : &
   stats_init_zm_api => stats_init_zm

  use stats_zt_module, only : &
   stats_init_zt_api => stats_init_zt

  use stats_sfc_module, only : &
   stats_init_sfc_api => stats_init_sfc

  use stats_clubb_utilities, only : &
    stats_begin_timestep_api => stats_begin_timestep

  use model_flags, only: &
    initialize_clubb_config_flags_type_api => initialize_clubb_config_flags_type, &
    set_default_clubb_config_flags_api     => set_default_clubb_config_flags, &
    saturation_bolton, & ! Constant for Bolton approximations of saturation
    saturation_gfdl,   & ! Constant for the GFDL approximation of saturation
    saturation_flatau, & ! Constant for Flatau approximations of saturation
    saturation_lookup    ! Use a lookup table for mixing length

  use corr_varnce_module, only : &
    setup_corr_varnce_array_api => setup_corr_varnce_array, &
    init_pdf_hydromet_arrays_api => init_pdf_hydromet_arrays

  use stats_clubb_utilities, only : &
    stats_accumulate_hydromet_api => stats_accumulate_hydromet, &
    stats_init_api => stats_init

  use stats_type, only: stats ! Type

  implicit none

  private

  public &
    ! To Implement CLUBB:
    set_default_parameters_api, & ! Procedure(s)
    read_parameters_api, &
    check_clubb_settings_api, &
    setup_parameters_model_api, &
        ! CLUBB can be set more specifically using these flags:
        iiPDF_ADG1, &
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
        ilambda0_stability_coef, imult_coef, itaumin, itaumax, imu, &
        iLscale_mu_coef, iLscale_pert_coef, ialpha_corr, iSkw_denom_coef, &
        ic_K10, ic_K10h, ithlp2_rad_coef, ithlp2_rad_cloud_frac_thresh, &
        iup2_sfc_coef, iSkw_max_mag, iC_invrs_tau_bkgnd, &
        iC_invrs_tau_sfc, iC_invrs_tau_shear, iC_invrs_tau_N2, &
        iC_invrs_tau_N2_wp2, iC_invrs_tau_N2_xp2, iC_invrs_tau_N2_wpxp, &
        iC_invrs_tau_N2_clear_wp3, ialtitude_threshold, irtp2_clip_coef, &
        iRichardson_num_min, iRichardson_num_max, iwpxp_Ri_exp, &
        ia3_coef_min, ia_const, iCx_min, iCx_max, ibv_efold, iz_displace



  public &
        advance_clubb_core_api, &
        advance_clubb_core_api_single_col, &
        advance_clubb_core_api_multi_col, &
        pdf_parameter, &
        implicit_coefs_terms, &
    cleanup_clubb_core_api

  public &
    ! To Implement SILHS:
    setup_corr_varnce_array_api, &
    setup_pdf_parameters_api, &
    hydromet_pdf_parameter, &
    init_pdf_hydromet_arrays_api, &
    ! generate_silhs_sample - SILHS API
    genrand_init_api, & ! if you are doing restarts)
    genrand_state, &
    genrand_srepr, &
    genrand_intg, &
    ! To use the results, you will need these variables:
    hmp2_ip_on_hmm2_ip_slope_type,      & ! Types
    hmp2_ip_on_hmm2_ip_intrcpt_type, &
    hm_metadata_type, &
    sclr_idx_type, &
    grid, &
    stats

  public &
    ! To Interact With CLUBB's Grid:
    ! For Varying Grids
    setup_grid_heights_api, &    ! if heights vary with time
    setup_grid_api

  public &
    ! To Obtain More Output from CLUBB for Diagnostics:
    stats_begin_timestep_api, &
    stats_end_timestep_api, &
    stats_finalize_api, &
    stats_init_api

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
    zt2zm_api

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
    num_hf_draw_points

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
    pack_pdf_params_api, &
    unpack_pdf_params_api, &
    num_pdf_params, &
!#endif
    init_pdf_params_api, &
    init_precip_fracs_api, &
    precipitation_fractions, &
    init_pdf_implicit_coefs_terms_api, &
    adj_low_res_nu_api, &
    assignment( = ), &
    clubb_i, &
    clubb_j, &
    compute_current_date_api, &
    gregorian2julian_day_api, &
    leap_year_api, &
    nvarmax_rad_zm, &
    nvarmax_rad_zt, &
    nvarmax_sfc, &
    nvarmax_zm, &
    nvarmax_zt

  public &
    nparams, &
    nu_vertical_res_dep, &
    setup_parameters_api, &
    stat_nknd, &
    stat_rknd, &
    stats_accumulate_hydromet_api, &
    stats_init_rad_zm_api, &
    stats_init_rad_zt_api, &
    stats_init_sfc_api, &
    stats_init_zm_api, &
    stats_init_zt_api, &
    stats_metadata_type

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
   copy_single_pdf_params_to_multi, &
   copy_multi_pdf_params_to_single
  
  interface setup_pdf_parameters_api
    module procedure setup_pdf_parameters_api_single_col
    module procedure setup_pdf_parameters_api_multi_col
  end interface
  
  interface advance_clubb_core_api
    module procedure advance_clubb_core_api_single_col
    module procedure advance_clubb_core_api_multi_col
  end interface
  
  interface setup_grid_api
    module procedure setup_grid_api_single_col
    module procedure setup_grid_api_multi_col
  end interface
  
  interface setup_parameters_api
    module procedure setup_parameters_api_single_col
    module procedure setup_parameters_api_multi_col
  end interface

  interface adj_low_res_nu_api
    module procedure adj_low_res_nu_api_single_col
    module procedure adj_low_res_nu_api_multi_col
  end interface
  
  interface setup_grid_heights_api
    module procedure setup_grid_heights_api_single_col
    module procedure setup_grid_heights_api_multi_col
  end interface
  
  interface update_xp2_mc_api
    module procedure update_xp2_mc_api_single_col
    module procedure update_xp2_mc_api_multi_col
  end interface
  
contains

  !================================================================================================
  ! advance_clubb_core - Advances the model one timestep.
  !================================================================================================

  subroutine advance_clubb_core_api_single_col( gr, &       ! intent(in)
    l_implemented, dt, fcor, sfc_elevation, &               ! intent(in)
    hydromet_dim, &                                         ! intent(in)
    sclr_dim, sclr_tol, edsclr_dim, sclr_idx, &         ! intent(in)
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
    hydromet, l_mix_rat_hm, &                               ! intent(in)
    rfrzm, radf, &                                          ! intent(in)
#ifdef CLUBBND_CAM
    varmu, &                                                ! intent(in)
#endif
    wphydrometp, wp2hmp, rtphmp_zt, thlphmp_zt, &           ! intent(in)
    host_dx, host_dy, &                                     ! intent(in)
    clubb_params, nu_vert_res_dep, lmin, &                  ! intent(in)
    clubb_config_flags, &                                   ! intent(in)
    stats_metadata, &                                       ! intent(in)
    stats_zt, stats_zm, stats_sfc, &                        ! intent(inout)
    um, vm, upwp, vpwp, up2, vp2, up3, vp3, &               ! intent(inout)
    thlm, rtm, wprtp, wpthlp, &                             ! intent(inout)
    wp2, wp3, rtp2, rtp3, thlp2, thlp3, rtpthlp, &          ! intent(inout)
    sclrm,   &
#ifdef GFDL
               sclrm_trsport_only,  &  ! h1g, 2010-06-16    ! intent(inout)
#endif
    sclrp2, sclrp3, sclrprtp, sclrpthlp, &                  ! intent(inout)
    wpsclrp, edsclrm, err_code_api, &                       ! intent(inout)
    rcm, cloud_frac, &                                      ! intent(inout)
    wpthvp, wp2thvp, rtpthvp, thlpthvp, &                   ! intent(inout)
    sclrpthvp, &                                            ! intent(inout)
    wp2rtp, wp2thlp, uprcp, vprcp, rc_coef_zm, wp4, &       ! intent(inout)
    wpup2, wpvp2, wp2up2, wp2vp2, ice_supersat_frac, &      ! intent(inout)
    um_pert, vm_pert, upwp_pert, vpwp_pert, &               ! intent(inout)
    pdf_params, pdf_params_zm, &                            ! intent(inout)
    pdf_implicit_coefs_terms, &                             ! intent(inout)
#ifdef GFDL
               RH_crit, & !h1g, 2010-06-16                  ! intent(inout)
               do_liquid_only_in_clubb, &                   ! intent(in)
#endif
    Kh_zm, Kh_zt, &                                         ! intent(out)
#ifdef CLUBB_CAM
    qclvar, &                                               ! intent(out)
#endif
    thlprcp, wprcp, w_up_in_cloud, w_down_in_cloud, &       ! intent(out)
    cloudy_updraft_frac, cloudy_downdraft_frac, &           ! intent(out)
    rcm_in_layer, cloud_cover, invrs_tau_zm )               ! intent(out)

    use advance_clubb_core_module, only : advance_clubb_core

    use pdf_parameter_module, only: &
        implicit_coefs_terms    ! Variable Type(s)

    use model_flags, only: &
        clubb_config_flags_type
        
    use stats_zm_module, only: &
        stats_init_zm ! Procedure(s)
        
    use stats_zt_module, only: &
        stats_init_zt
        
    use stats_sfc_module, only: &
        stats_init_sfc ! Procedure(s)

    use array_index, only: &
      sclr_idx_type

    implicit none

    !--------------------------- Input Variables ---------------------------
    type (grid), intent(in) :: gr

    logical, intent(in) ::  &
      l_implemented    ! True if CLUBB is being run within a large-scale host model,
                       !   rather than a standalone single-column model.

    real( kind = core_rknd ), intent(in) ::  &
      dt  ! Current timestep duration    [s]

    real( kind = core_rknd ), intent(in) ::  &
      fcor,  &          ! Coriolis forcing             [s^-1]
      sfc_elevation     ! Elevation of ground level    [m above MSL]

    integer, intent(in) :: &
      hydromet_dim,   & ! Total number of hydrometeor species       [#]
      sclr_dim,       & ! Number of passive scalars                 [#]
      edsclr_dim        ! Number of eddy-diff. passive scalars      [#]

    real( kind = core_rknd ), intent(in), dimension(sclr_dim) :: & 
      sclr_tol          ! Threshold(s) on the passive scalars  [units vary]

    type (sclr_idx_type), intent(in) :: &
      sclr_idx

    real( kind = core_rknd ), intent(in), dimension(gr%nzt) ::  &
      thlm_forcing,    & ! liquid potential temp. forcing (thermodynamic levels)    [K/s]
      rtm_forcing,     & ! total water forcing (thermodynamic levels)        [(kg/kg)/s]
      um_forcing,      & ! eastward wind forcing (thermodynamic levels)     [m/s/s]
      vm_forcing,      & ! northward wind forcing (thermodynamic levels)     [m/s/s]
      wm_zt,           & ! vertical mean wind component on thermo. levels   [m/s]
      rho,             & ! Air density on thermodynamic levels       [kg/m^3]
      rho_ds_zt,       & ! Dry, static density on thermo. levels     [kg/m^3]
      invrs_rho_ds_zt, & ! Inverse dry, static density on thermo levs.  [m^3/kg]
      thv_ds_zt,       & ! Dry, base-state theta_v on thermo levs.  [K]
      rfrzm              ! Total ice-phase water mixing ratio        [kg/kg]

    real( kind = core_rknd ), intent(in), dimension(gr%nzm) ::  &
      wprtp_forcing,   & ! total water turbulent flux forcing (momentum levels)    [m*K/s^2]
      wpthlp_forcing,  & ! liq pot temp turb flux forcing (momentum levels)   [m*(kg/kg)/s^2]
      rtp2_forcing,    & ! total water variance forcing (momentum levels)    [(kg/kg)^2/s]
      thlp2_forcing,   & ! liq pot temp variance forcing (momentum levels)   [K^2/s]
      rtpthlp_forcing, & ! <r_t'th_l'> covariance forcing (momentum levels) [K*(kg/kg)/s]
      wm_zm,           & ! vertical mean wind component on momentum levels  [m/s]
      rho_zm,          & ! Air density on momentum levels            [kg/m^3]
      rho_ds_zm,       & ! Dry, static density on momentum levels    [kg/m^3]
      invrs_rho_ds_zm, & ! Inverse dry, static density on momentum levs. [m^3/kg]
      thv_ds_zm          ! Dry, base-state theta_v on momentum levs. [K]

    real( kind = core_rknd ), dimension(gr%nzt,hydromet_dim), intent(in) :: &
      hydromet           ! Array of hydrometeors                [units vary]

    logical, dimension(hydromet_dim), intent(in) :: &
      l_mix_rat_hm   ! if true, then the quantity is a hydrometeor mixing ratio

    real( kind = core_rknd ), dimension(gr%nzm), intent(in) :: &
      radf          ! Buoyancy production at cloud top due to longwave radiative cooling [m^2/s^3]

#ifdef CLUBBND_CAM
    real( kind = core_rknd ), intent(in) :: &
      varmu
#endif

    real( kind = core_rknd ), dimension(gr%nzm,hydromet_dim), intent(in) :: &
      wphydrometp    ! Covariance of w and a hydrometeor      [(m/s) <hm units>]

    real( kind = core_rknd ), dimension(gr%nzt,hydromet_dim), intent(in) :: &
      wp2hmp,      & ! Third-order moment:  < w'^2 hm' > (hm = hydrometeor) [(m/s)^2 <hm units>]
      rtphmp_zt,   & ! Covariance of rt and hm (on thermo levs.) [(kg/kg) <hm units>]
      thlphmp_zt     ! Covariance of thl and hm (on thermo levs.)      [K <hm units>]

    real( kind = core_rknd ), intent(in) ::  &
      wpthlp_sfc,   & ! w' theta_l' at surface   [(m K)/s]
      wprtp_sfc,    & ! w' r_t' at surface       [(kg m)/( kg s)]
      upwp_sfc,     & ! u'w' at surface          [m^2/s^2]
      vpwp_sfc,     & ! v'w' at surface          [m^2/s^2]
      p_sfc           ! Pressure at surface      [Pa]

    ! Passive scalar variables
    real( kind = core_rknd ), intent(in), dimension(gr%nzt,sclr_dim) :: &
      sclrm_forcing    ! Passive scalar forcing         [{units vary}/s]

    real( kind = core_rknd ), intent(in), dimension(sclr_dim) ::  &
      wpsclrp_sfc      ! Passive scalar flux at surface         [{units vary} m/s]

    ! Eddy passive scalar variables
    real( kind = core_rknd ), intent(in), dimension(gr%nzt,edsclr_dim) :: &
      edsclrm_forcing  ! Eddy-diffusion passive scalar forcing    [{units vary}/s]

    real( kind = core_rknd ), intent(in), dimension(edsclr_dim) ::  &
      wpedsclrp_sfc    ! Eddy-diffusion passive scalar flux at surface    [{units vary} m/s

    real( kind = core_rknd ), intent(in) :: &
      upwp_sfc_pert, & ! pertubed u'w' at surface    [m^2/s^2]
      vpwp_sfc_pert    ! pertubed v'w' at surface    [m^2/s^2]

    ! Reference profiles (used for nudging, sponge damping, and Coriolis effect)
    real( kind = core_rknd ), dimension(gr%nzt), intent(in) ::  &
      rtm_ref,  & ! Initial total water mixing ratio             [kg/kg]
      thlm_ref, & ! Initial liquid water potential temperature   [K]
      um_ref,   & ! Initial u wind; Michael Falk                 [m/s]
      vm_ref,   & ! Initial v wind; Michael Falk                 [m/s]
      ug,       & ! u geostrophic wind                           [m/s]
      vg          ! v geostrophic wind                           [m/s]

    ! Host model horizontal grid spacing, if part of host model.
    real( kind = core_rknd ), intent(in) :: &
      host_dx,  & ! East-west horizontal grid spacing     [m]
      host_dy     ! North-south horizontal grid spacing   [m]

    real( kind = core_rknd ), dimension(nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

    type(nu_vertical_res_dep), intent(in) :: &
      nu_vert_res_dep    ! Vertical resolution dependent nu values

    real( kind = core_rknd ), intent(in) :: &
      lmin    ! Min. value for the length scale    [m]

    type( clubb_config_flags_type ), intent(in) :: &
      clubb_config_flags ! Derived type holding all configurable CLUBB flags

    type (stats_metadata_type), intent(in) :: &
      stats_metadata

    !--------------------------- Input/Output Variables ---------------------------
    type (stats), intent(inout) :: &
      stats_zt, &
      stats_zm, &
      stats_sfc

    ! These are prognostic or are planned to be in the future
    real( kind = core_rknd ), intent(inout), dimension(gr%nzt) ::  &
      um,      & ! eastward grid-mean wind component (thermodynamic levels)   [m/s]
      vm,      & ! northward grid-mean wind component (thermodynamic levels)   [m/s]
      up3,     & ! u'^3 (thermodynamic levels)                    [m^3/s^3]
      vp3,     & ! v'^3 (thermodynamic levels)                    [m^3/s^3]
      rtm,     & ! total water mixing ratio, r_t (thermo. levels) [kg/kg]
      thlm,    & ! liq. water pot. temp., th_l (thermo. levels)   [K]
      rtp3,    & ! r_t'^3 (thermodynamic levels)                  [(kg/kg)^3]
      thlp3,   & ! th_l'^3 (thermodynamic levels)                 [K^3]
      wp3        ! w'^3 (thermodynamic levels)                    [m^3/s^3]

    real( kind = core_rknd ), intent(inout), dimension(gr%nzm) ::  &
      upwp,    & ! u'w' (momentum levels)                         [m^2/s^2]
      vpwp,    & ! v'w' (momentum levels)                         [m^2/s^2]
      up2,     & ! u'^2 (momentum levels)                         [m^2/s^2]
      vp2,     & ! v'^2 (momentum levels)                         [m^2/s^2]
      wprtp,   & ! w' r_t' (momentum levels)                      [(kg/kg) m/s]
      wpthlp,  & ! w'th_l' (momentum levels)                      [(m/s) K]
      rtp2,    & ! r_t'^2 (momentum levels)                       [(kg/kg)^2]
      thlp2,   & ! th_l'^2 (momentum levels)                      [K^2]
      rtpthlp, & ! r_t'th_l' (momentum levels)                    [(kg/kg) K]
      wp2        ! w'^2 (momentum levels)                         [m^2/s^2]

    ! Passive scalar variables
    real( kind = core_rknd ), intent(inout), dimension(gr%nzt,sclr_dim) :: &
      sclrm,     & ! Passive scalar mean (thermo. levels) [units vary]
      sclrp3       ! sclr'^3 (thermodynamic levels)       [{units vary}^3]

    real( kind = core_rknd ), intent(inout), dimension(gr%nzm,sclr_dim) :: &
      wpsclrp,   & ! w'sclr' (momentum levels)            [{units vary} m/s]
      sclrp2,    & ! sclr'^2 (momentum levels)            [{units vary}^2]
      sclrprtp,  & ! sclr'rt' (momentum levels)           [{units vary} (kg/kg)]
      sclrpthlp    ! sclr'thl' (momentum levels)          [{units vary} K]

    real( kind = core_rknd ), intent(inout), dimension(gr%nzt) ::  &
      p_in_Pa, & ! Air pressure (thermodynamic levels)       [Pa]
      exner      ! Exner function (thermodynamic levels)     [-]

    real( kind = core_rknd ), intent(inout), dimension(gr%nzt) ::  &
      rcm,        & ! cloud water mixing ratio, r_c (thermo. levels) [kg/kg]
      cloud_frac, & ! cloud fraction (thermodynamic levels)          [-]
      wp2thvp       ! < w'^2 th_v' > (thermodynamic levels)          [m^2/s^2 K]

    real( kind = core_rknd ), intent(inout), dimension(gr%nzm) ::  &
      wpthvp,     & ! < w' th_v' > (momentum levels)                 [kg/kg K]
      rtpthvp,    & ! < r_t' th_v' > (momentum levels)               [kg/kg K]
      thlpthvp      ! < th_l' th_v' > (momentum levels)              [K^2]

    real( kind = core_rknd ), intent(inout), dimension(gr%nzm,sclr_dim) :: &
      sclrpthvp     ! < sclr' th_v' > (momentum levels)   [units vary]

    real( kind = core_rknd ), intent(inout), dimension(gr%nzt) ::  &
      wp2rtp,            & ! w'^2 rt' (thermodynamic levels)      [m^2/s^2 kg/kg]
      wp2thlp,           & ! w'^2 thl' (thermodynamic levels)     [m^2/s^2 K]
      wpup2,             & ! w'u'^2 (thermodynamic levels)        [m^3/s^3]
      wpvp2,             & ! w'v'^2 (thermodynamic levels)        [m^3/s^3]
      ice_supersat_frac    ! ice cloud fraction (thermo. levels)  [-]

    real( kind = core_rknd ), intent(inout), dimension(gr%nzm) ::  &
      uprcp,             & ! < u' r_c' > (momentum levels)        [(m/s)(kg/kg)]
      vprcp,             & ! < v' r_c' > (momentum levels)        [(m/s)(kg/kg)]
      rc_coef_zm,        & ! Coef of X'r_c' in Eq. (34) (m-levs.) [K/(kg/kg)]
      wp4,               & ! w'^4 (momentum levels)               [m^4/s^4]
      wp2up2,            & ! w'^2 u'^2 (momentum levels)          [m^4/s^4]
      wp2vp2               ! w'^2 v'^2 (momentum levels)          [m^4/s^4]

    ! Variables used to track perturbed version of winds.
    real( kind = core_rknd ), intent(inout), dimension(gr%nzt) :: &
      um_pert,   & ! perturbed <u>       [m/s]
      vm_pert      ! perturbed <v>       [m/s]

    real( kind = core_rknd ), intent(inout), dimension(gr%nzm) :: &
      upwp_pert, & ! perturbed <u'w'>    [m^2/s^2]
      vpwp_pert    ! perturbed <v'w'>    [m^2/s^2] 

    type(pdf_parameter), intent(inout) :: &
      pdf_params,    & ! Fortran structure of PDF parameters on thermodynamic levels    [units vary]
      pdf_params_zm    ! Fortran structure of PDF parameters on momentum levels        [units vary]

    type(implicit_coefs_terms), intent(inout) :: &
      pdf_implicit_coefs_terms    ! Implicit coefs / explicit terms [units vary]

#ifdef GFDL
    real( kind = core_rknd ), intent(inout), dimension(gr%nzt,sclr_dim) :: &  ! h1g, 2010-06-16
      sclrm_trsport_only  ! Passive scalar concentration due to pure transport [{units vary}/s]
#endif

    ! Eddy passive scalar variable
    real( kind = core_rknd ), intent(inout), dimension(gr%nzt,edsclr_dim) :: &
      edsclrm   ! Eddy passive scalar grid-mean (thermo. levels)   [units vary]

    ! Variables that need to be output for use in other parts of the CLUBB
    ! code, such as microphysics (rcm, pdf_params), forcings (rcm), and/or
    ! BUGSrad (cloud_cover).
    real( kind = core_rknd ), intent(out), dimension(gr%nzt) ::  &
      rcm_in_layer, & ! rcm within cloud layer                          [kg/kg]
      cloud_cover     ! cloud cover                                     [-]

    ! Variables that need to be output for use in host models
    real( kind = core_rknd ), intent(out), dimension(gr%nzt) ::  &
      w_up_in_cloud,         & ! Average cloudy updraft velocity       [m/s]
      w_down_in_cloud,       & ! Average cloudy downdraft velocity     [m/s]
      cloudy_updraft_frac,   & ! cloudy updraft fraction               [-]
      cloudy_downdraft_frac    ! cloudy downdraft fraction             [-]

    real( kind = core_rknd ), intent(out), dimension(gr%nzm) ::  &
      wprcp,                 & ! w'r_c' (momentum levels)              [(kg/kg) m/s]
      invrs_tau_zm             ! One divided by tau on zm levels       [1/s]

    real( kind = core_rknd ), dimension(gr%nzt), intent(out) :: &
      Kh_zt    ! Eddy diffusivity coefficient on thermodynamic levels   [m^2/s]

    real( kind = core_rknd ), dimension(gr%nzm), intent(out) :: &
      Kh_zm    ! Eddy diffusivity coefficient on momentum levels        [m^2/s]

#ifdef CLUBB_CAM
    real( kind = core_rknd), intent(out), dimension(gr%nzt) :: &
      qclvar        ! cloud water variance
#endif

    real( kind = core_rknd ), dimension(gr%nzm), intent(out) :: &
      thlprcp    ! thl'rc'              [K kg/kg]

    !!! Output Variable 
    integer, intent(inout) :: err_code_api ! Diagnostic, for if some calculation goes amiss.

#ifdef GFDL
    ! hlg, 2010-06-16
    real( kind = core_rknd ), intent(inout), dimension(gr%nzt, min(1,sclr_dim) , 2) :: &
      RH_crit  ! critical relative humidity for droplet and ice nucleation
! ---> h1g, 2012-06-14
    logical, intent(in)                 ::  do_liquid_only_in_clubb
! <--- h1g, 2012-06-14
#endif

    !------------------------- Local Variables -------------------------
    type(stats), dimension(1) :: &
      stats_zt_col, &
      stats_zm_col, &
      stats_sfc_col

    real( kind = core_rknd ), dimension(1) ::  &
      fcor_col,  &          ! Coriolis forcing             [s^-1]
      sfc_elevation_col     ! Elevation of ground level    [m AMSL]

    ! Input Variables
    real( kind = core_rknd ), dimension(1,gr%nzt) ::  &
      thlm_forcing_col,    & ! liquid potential temp. forcing (thermodynamic levels)    [K/s]
      rtm_forcing_col,     & ! total water forcing (thermodynamic levels)        [(kg/kg)/s]
      um_forcing_col,      & ! eastward wind forcing (thermodynamic levels)     [m/s/s]
      vm_forcing_col,      & ! northward wind forcing (thermodynamic levels)     [m/s/s]
      wm_zt_col,           & ! vertical mean wind component on thermo. levels   [m/s]
      rho_col,             & ! Air density on thermodynamic levels       [kg/m^3]
      rho_ds_zt_col,       & ! Dry, static density on thermo. levels     [kg/m^3]
      invrs_rho_ds_zt_col, & ! Inverse dry, static density on thermo levs.  [m^3/kg]
      thv_ds_zt_col,       & ! Dry, base-state theta_v on thermo levs.  [K]
      rfrzm_col              ! Total ice-phase water mixing ratio        [kg/kg]

    real( kind = core_rknd ), dimension(1,gr%nzm) ::  &
      wprtp_forcing_col,   & ! total water turbulent flux forcing (momentum levels)    [m*K/s^2]
      wpthlp_forcing_col,  & ! liq pot temp turb flux forcing (momentum levels)   [m*(kg/kg)/s^2]
      rtp2_forcing_col,    & ! total water variance forcing (momentum levels)    [(kg/kg)^2/s]
      thlp2_forcing_col,   & ! liq pot temp variance forcing (momentum levels)   [K^2/s]
      rtpthlp_forcing_col, & ! <r_t'th_l'> covariance forcing (momentum levels) [K*(kg/kg)/s]
      wm_zm_col,           & ! vertical mean wind component on momentum levels  [m/s]
      rho_zm_col,          & ! Air density on momentum levels            [kg/m^3]
      rho_ds_zm_col,       & ! Dry, static density on momentum levels    [kg/m^3]
      invrs_rho_ds_zm_col, & ! Inverse dry, static density on momentum levs. [m^3/kg]
      thv_ds_zm_col          ! Dry, base-state theta_v on momentum levs. [K]

    real( kind = core_rknd ), dimension(1,gr%nzt,hydromet_dim) :: &
      hydromet_col           ! Array of hydrometeors                [units vary]

    real( kind = core_rknd ), dimension(1,gr%nzm) :: &
      radf_col          ! Buoyancy production at cloud top due to 
                        ! longwave radiative cooling [m^2/s^3]

#ifdef CLUBBND_CAM 
    real( kind = core_rknd ), dimension(1) :: & 
      varmu_col 
#endif 

    real( kind = core_rknd ), dimension(1,gr%nzm,hydromet_dim) :: &
      wphydrometp_col    ! Covariance of w and a hydrometeor      [(m/s) <hm units>]

    real( kind = core_rknd ), dimension(1,gr%nzt,hydromet_dim) :: &
      wp2hmp_col,      & ! Third-order moment:  < w'^2 hm' > (hm = hydrometeor) [(m/s)^2 <hm units>]
      rtphmp_zt_col,   & ! Covariance of rt and hm (on thermo levs.) [(kg/kg) <hm units>]
      thlphmp_zt_col     ! Covariance of thl and hm (on thermo levs.)      [K <hm units>]

    real( kind = core_rknd ), dimension(1) ::  &
      wpthlp_sfc_col,   & ! w' theta_l' at surface   [(m K)/s]
      wprtp_sfc_col,    & ! w' r_t' at surface       [(kg m)/( kg s)]
      upwp_sfc_col,     & ! u'w' at surface          [m^2/s^2]
      vpwp_sfc_col,     & ! v'w' at surface          [m^2/s^2]
      p_sfc_col           ! Pressure at surface      [Pa]

    ! Passive scalar variables
    real( kind = core_rknd ), dimension(1,gr%nzt,sclr_dim) :: &
      sclrm_forcing_col    ! Passive scalar forcing         [{units vary}/s]

    real( kind = core_rknd ), dimension(1,sclr_dim) ::  &
      wpsclrp_sfc_col      ! Scalar flux at surface         [{units vary} m/s]

    ! Eddy passive scalar variables
    real( kind = core_rknd ), dimension(1,gr%nzt,edsclr_dim) :: &
      edsclrm_forcing_col  ! Eddy passive scalar forcing    [{units vary}/s]

    real( kind = core_rknd ), dimension(1) :: &
      upwp_sfc_pert_col, & ! pertubed u'w' at surface    [m^2/s^2]
      vpwp_sfc_pert_col    ! pertubed v'w' at surface    [m^2/s^2]

    real( kind = core_rknd ), dimension(1,edsclr_dim) ::  &
      wpedsclrp_sfc_col    ! Eddy-Scalar flux at surface    [{units vary} m/s]

    ! Reference profiles (used for nudging, sponge damping, and Coriolis effect)
    real( kind = core_rknd ), dimension(1,gr%nzt) ::  &
      rtm_ref_col,  & ! Initial total water mixing ratio             [kg/kg]
      thlm_ref_col, & ! Initial liquid water potential temperature   [K]
      um_ref_col,   & ! Initial u wind; Michael Falk                 [m/s]
      vm_ref_col,   & ! Initial v wind; Michael Falk                 [m/s]
      ug_col,       & ! u geostrophic wind                           [m/s]
      vg_col          ! v geostrophic wind                           [m/s]

    ! Host model horizontal grid spacing, if part of host model.
    real( kind = core_rknd ), dimension(1) :: &
      host_dx_col,  & ! East-West horizontal grid spacing     [m]
      host_dy_col     ! North-South horizontal grid spacing   [m]

    real( kind = core_rknd ), dimension(1,nparams) :: &
      clubb_params_col    ! Array of CLUBB's tunable parameters    [units vary]

    ! These are prognostic or are planned to be in the future
    real( kind = core_rknd ), dimension(1,gr%nzt) ::  &
      um_col,      & ! eastward grid-mean wind component (thermodynamic levels)   [m/s]
      vm_col,      & ! northward grid-mean wind component (thermodynamic levels)   [m/s]
      up3_col,     & ! u'^3 (thermodynamic levels)                    [m^3/s^3]
      vp3_col,     & ! v'^3 (thermodynamic levels)                    [m^3/s^3]
      rtm_col,     & ! total water mixing ratio, r_t (thermo. levels) [kg/kg]
      thlm_col,    & ! liq. water pot. temp., th_l (thermo. levels)   [K]
      rtp3_col,    & ! r_t'^3 (thermodynamic levels)                  [(kg/kg)^3]
      thlp3_col,   & ! th_l'^3 (thermodynamic levels)                 [K^3]
      wp3_col        ! w'^3 (thermodynamic levels)                    [m^3/s^3]

    real( kind = core_rknd ), dimension(1,gr%nzm) ::  &
      upwp_col,    & ! u'w' (momentum levels)                         [m^2/s^2]
      vpwp_col,    & ! v'w' (momentum levels)                         [m^2/s^2]
      up2_col,     & ! u'^2 (momentum levels)                         [m^2/s^2]
      vp2_col,     & ! v'^2 (momentum levels)                         [m^2/s^2]
      wprtp_col,   & ! w' r_t' (momentum levels)                      [(kg/kg) m/s]
      wpthlp_col,  & ! w'th_l' (momentum levels)                      [(m/s) K]
      rtp2_col,    & ! r_t'^2 (momentum levels)                       [(kg/kg)^2]
      thlp2_col,   & ! th_l'^2 (momentum levels)                      [K^2]
      rtpthlp_col, & ! r_t'th_l' (momentum levels)                    [(kg/kg) K]
      wp2_col        ! w'^2 (momentum levels)                         [m^2/s^2]

    ! Passive scalar variables
    real( kind = core_rknd ), dimension(1,gr%nzt,sclr_dim) :: &
      sclrm_col,     & ! Passive scalar mean (thermo. levels) [units vary]
      sclrp3_col       ! sclr'^3 (thermodynamic levels)       [{units vary}^3]

    real( kind = core_rknd ), dimension(1,gr%nzm,sclr_dim) :: &
      wpsclrp_col,   & ! w'sclr' (momentum levels)            [{units vary} m/s]
      sclrp2_col,    & ! sclr'^2 (momentum levels)            [{units vary}^2]
      sclrprtp_col,  & ! sclr'rt' (momentum levels)           [{units vary} (kg/kg)]
      sclrpthlp_col    ! sclr'thl' (momentum levels)          [{units vary} K]

    real( kind = core_rknd ), dimension(1,gr%nzt) ::  &
      p_in_Pa_col, & ! Air pressure (thermodynamic levels)       [Pa]
      exner_col      ! Exner function (thermodynamic levels)     [-]

    real( kind = core_rknd ), dimension(1,gr%nzt) ::  &
      rcm_col,        & ! cloud water mixing ratio, r_c (thermo. levels) [kg/kg]
      cloud_frac_col, & ! cloud fraction (thermodynamic levels)          [-]
      wp2thvp_col       ! < w'^2 th_v' > (thermodynamic levels)          [m^2/s^2 K]

    real( kind = core_rknd ), dimension(1,gr%nzm) ::  &
      wpthvp_col,     & ! < w' th_v' > (momentum levels)                 [kg/kg K]
      rtpthvp_col,    & ! < r_t' th_v' > (momentum levels)               [kg/kg K]
      thlpthvp_col      ! < th_l' th_v' > (momentum levels)              [K^2]

    real( kind = core_rknd ), dimension(1,gr%nzm,sclr_dim) :: &
      sclrpthvp_col     ! < sclr' th_v' > (momentum levels)   [units vary]

    real( kind = core_rknd ), dimension(1,gr%nzt) ::  &
      wp2rtp_col,            & ! w'^2 rt' (thermodynamic levels)      [m^2/s^2 kg/kg]
      wp2thlp_col,           & ! w'^2 thl' (thermodynamic levels)     [m^2/s^2 K]
      wpup2_col,             & ! w'u'^2 (thermodynamic levels)        [m^3/s^3]
      wpvp2_col,             & ! w'v'^2 (thermodynamic levels)        [m^3/s^3]
      ice_supersat_frac_col    ! ice cloud fraction (thermo. levels)  [-]

    real( kind = core_rknd ), dimension(1,gr%nzm) ::  &
      uprcp_col,             & ! < u' r_c' > (momentum levels)        [(m/s)(kg/kg)]
      vprcp_col,             & ! < v' r_c' > (momentum levels)        [(m/s)(kg/kg)]
      rc_coef_zm_col,        & ! Coef of X'r_c' in Eq. (34) (m-levs.) [K/(kg/kg)]
      wp4_col,               & ! w'^4 (momentum levels)               [m^4/s^4]
      wp2up2_col,            & ! w'^2 u'^2 (momentum levels)          [m^4/s^4]
      wp2vp2_col               ! w'^2 v'^2 (momentum levels)          [m^4/s^4]

    ! Variables used to track perturbed version of winds.
    real( kind = core_rknd ), dimension(1,gr%nzt) :: &
      um_pert_col,   & ! perturbed <u>       [m/s]
      vm_pert_col      ! perturbed <v>       [m/s]

    real( kind = core_rknd ), dimension(1,gr%nzm) :: &
      upwp_pert_col, & ! perturbed <u'w'>    [m^2/s^2]
      vpwp_pert_col    ! perturbed <v'w'>    [m^2/s^2]

#ifdef GFDL
    real( kind = core_rknd ), dimension(1,gr%nzt,sclr_dim) :: &  ! h1g, 2010-06-16
      sclrm_trsport_only_col  ! Passive scalar concentration due to pure transport [{units vary}/s]
#endif

      real( kind = core_rknd ), dimension(1,gr%nzt,edsclr_dim) :: &
      edsclrm_col   ! Eddy passive scalar mean (thermo. levels)   [units vary]

    real( kind = core_rknd ), dimension(1,gr%nzt) ::  &
      rcm_in_layer_col, & ! rcm in cloud layer                              [kg/kg]
      cloud_cover_col     ! cloud cover                                     [-]

    ! Variables that need to be output for use in host models
    real( kind = core_rknd ), dimension(1,gr%nzt) ::  &
      w_up_in_cloud_col,         & ! Average cloudy updraft velocity       [m/s]
      w_down_in_cloud_col,       & ! Average cloudy downdraft velocity     [m/s]
      cloudy_updraft_frac_col,   & ! cloudy updraft fraction               [-]
      cloudy_downdraft_frac_col    ! cloudy downdraft fraction             [-]

    real( kind = core_rknd ), dimension(1,gr%nzm) ::  &
      wprcp_col,                 & ! w'r_c' (momentum levels)              [(kg/kg) m/s]
      invrs_tau_zm_col             ! One divided by tau on zm levels       [1/s]

    real( kind = core_rknd ), dimension(1,gr%nzt) :: &
      Kh_zt_col    ! Eddy diffusivity coefficient on thermodynamic levels   [m^2/s]

    real( kind = core_rknd ), dimension(1,gr%nzm) :: &
      Kh_zm_col    ! Eddy diffusivity coefficient on momentum levels        [m^2/s]

#ifdef CLUBB_CAM
    real( kind = core_rknd), dimension(1,gr%nzt) :: &
      qclvar_col        ! cloud water variance
#endif

    real( kind = core_rknd ), dimension(1,gr%nzm) :: &
      thlprcp_col    ! thl'rc'              [K kg/kg]
      
#ifdef GFDL
    ! hlg, 2010-06-16
    real( kind = core_rknd ), dimension(1,gr%nzt, min(1,sclr_dim) , 2) :: &
      RH_crit_col  ! critical relative humidity for droplet and ice nucleation
    logical, intent(in)                 ::  do_liquid_only_in_clubb
#endif

    !------------------------- Begin Code -------------------------

    fcor_col(1) = fcor
    sfc_elevation_col(1) = sfc_elevation
    
    thlm_forcing_col(1,:) = thlm_forcing
    rtm_forcing_col(1,:) = rtm_forcing
    um_forcing_col(1,:) = um_forcing
    vm_forcing_col(1,:) = vm_forcing
    sclrm_forcing_col(1,:,:) = sclrm_forcing
    edsclrm_forcing_col(1,:,:) = edsclrm_forcing
    wprtp_forcing_col(1,:) = wprtp_forcing
    wpthlp_forcing_col(1,:) = wpthlp_forcing
    rtp2_forcing_col(1,:) = rtp2_forcing
    thlp2_forcing_col(1,:) = thlp2_forcing
    rtpthlp_forcing_col(1,:) = rtpthlp_forcing
    wm_zm_col(1,:) = wm_zm
    wm_zt_col(1,:) = wm_zt
    
    wpthlp_sfc_col(1) = wpthlp_sfc
    wprtp_sfc_col(1) = wprtp_sfc
    upwp_sfc_col(1) = upwp_sfc
    vpwp_sfc_col(1) = vpwp_sfc
    p_sfc_col(1) = p_sfc
    
    wpsclrp_sfc_col(1,:) = wpsclrp_sfc
    wpedsclrp_sfc_col(1,:) = wpedsclrp_sfc
    
    upwp_sfc_pert_col(1) = upwp_sfc_pert
    vpwp_sfc_pert_col(1) = vpwp_sfc_pert

    rtm_ref_col(1,:) = rtm_ref
    thlm_ref_col(1,:) = thlm_ref
    um_ref_col(1,:) = um_ref
    vm_ref_col(1,:) = vm_ref
    ug_col(1,:) = ug
    vg_col(1,:) = vg
    
    p_in_Pa_col(1,:) = p_in_Pa
    rho_zm_col(1,:) = rho_zm
    rho_col(1,:) = rho
    exner_col(1,:) = exner
    rho_ds_zm_col(1,:) = rho_ds_zm
    rho_ds_zt_col(1,:) = rho_ds_zt
    invrs_rho_ds_zm_col(1,:) = invrs_rho_ds_zm
    invrs_rho_ds_zt_col(1,:) = invrs_rho_ds_zt
    thv_ds_zm_col(1,:) = thv_ds_zm
    thv_ds_zt_col(1,:) = thv_ds_zt
    rfrzm_col(1,:) = rfrzm
    
    hydromet_col(1,:,:) = hydromet
    
    radf_col(1,:) = radf
#ifdef CLUBBND_CAM
    varmu_col(1) = varmu
#endif
    wphydrometp_col(1,:,:) = wphydrometp
    wp2hmp_col(1,:,:) = wp2hmp
    rtphmp_zt_col(1,:,:) = rtphmp_zt
    thlphmp_zt_col(1,:,:) = thlphmp_zt
    
    host_dx_col(1) = host_dx
    host_dy_col(1) = host_dy

    clubb_params_col(1,:) = clubb_params
    
    stats_zt_col(1) = stats_zt
    stats_zm_col(1) = stats_zm
    stats_sfc_col(1) = stats_sfc
    
    um_col(1,:) = um
    vm_col(1,:) = vm
    upwp_col(1,:) = upwp
    vpwp_col(1,:) = vpwp
    up2_col(1,:) = up2
    vp2_col(1,:) = vp2
    up3_col(1,:) = up3
    vp3_col(1,:) = vp3
    thlm_col(1,:) = thlm
    rtm_col(1,:) = rtm
    wprtp_col(1,:) = wprtp
    wpthlp_col(1,:) = wpthlp
    wp2_col(1,:) = wp2
    wp3_col(1,:) = wp3
    rtp2_col(1,:) = rtp2
    rtp3_col(1,:) = rtp3
    thlp2_col(1,:) = thlp2
    thlp3_col(1,:) = thlp3
    rtpthlp_col(1,:) = rtpthlp
    
    sclrm_col(1,:,:) = sclrm
#ifdef GFDL
    sclrm_trsport_only_col(1,:,:) = sclrm_trsport_only
#endif
    sclrp2_col(1,:,:) = sclrp2
    sclrp3_col(1,:,:) = sclrp3
    sclrprtp_col(1,:,:) = sclrprtp
    sclrpthlp_col(1,:,:) = sclrpthlp
    wpsclrp_col(1,:,:) = wpsclrp
    edsclrm_col(1,:,:) = edsclrm
    
    rcm_col(1,:) = rcm
    cloud_frac_col(1,:) = cloud_frac
    wpthvp_col(1,:) = wpthvp
    wp2thvp_col(1,:) = wp2thvp
    rtpthvp_col(1,:) = rtpthvp
    thlpthvp_col(1,:) = thlpthvp
    sclrpthvp_col(1,:,:) = sclrpthvp
    wp2rtp_col(1,:) = wp2rtp
    wp2thlp_col(1,:) = wp2thlp
    uprcp_col(1,:) = uprcp
    vprcp_col(1,:) = vprcp
    rc_coef_zm_col(1,:) = rc_coef_zm
    wp4_col(1,:) = wp4
    wpup2_col(1,:) = wpup2
    wpvp2_col(1,:) = wpvp2
    wp2up2_col(1,:) = wp2up2
    wp2vp2_col(1,:) = wp2vp2
    ice_supersat_frac_col(1,:) = ice_supersat_frac
    um_pert_col(1,:) = um_pert
    vm_pert_col(1,:) = vm_pert
    upwp_pert_col(1,:) = upwp_pert
    vpwp_pert_col(1,:) = vpwp_pert
#ifdef GFDL
    RH_crit_col(1,:,:,:) = RH_crit
#endif
    Kh_zm_col(1,:) = Kh_zm
    Kh_zt_col(1,:) = Kh_zt
#ifdef CLUBB_CAM
    qclvar_col(1,:) = qclvar
#endif
    thlprcp_col(1,:) = thlprcp
    wprcp_col(1,:) = wprcp
    w_up_in_cloud_col(1,:) = w_up_in_cloud
    w_down_in_cloud_col(1,:) = w_down_in_cloud
    cloudy_updraft_frac_col(1,:) = cloudy_updraft_frac
    cloudy_downdraft_frac_col(1,:) = cloudy_downdraft_frac
    rcm_in_layer_col(1,:) = rcm_in_layer
    cloud_cover_col(1,:) = cloud_cover
    invrs_tau_zm_col(1,:) = invrs_tau_zm

    !$acc data copyin( gr, gr%zm, gr%zt, gr%dzm, gr%dzt, gr%invrs_dzt, gr%invrs_dzm, &
    !$acc              gr%weights_zt2zm, gr%weights_zm2zt, &
    !$acc              nu_vert_res_dep, nu_vert_res_dep%nu2, nu_vert_res_dep%nu9, &
    !$acc              nu_vert_res_dep%nu1, nu_vert_res_dep%nu8, nu_vert_res_dep%nu10, &
    !$acc              nu_vert_res_dep%nu6, &
    !$acc              sclr_idx, clubb_params_col, &
    !$acc              fcor_col, sfc_elevation_col, thlm_forcing_col, rtm_forcing_col, &
    !$acc              um_forcing_col, &
    !$acc              vm_forcing_col, wprtp_forcing_col, wpthlp_forcing_col, rtp2_forcing_col, &
    !$acc              thlp2_forcing_col, &
    !$acc              rtpthlp_forcing_col, wm_zm_col, wm_zt_col, rho_zm_col, rho_col, &
    !$acc              rho_ds_zm_col, rho_ds_zt_col, &
    !$acc              invrs_rho_ds_zm_col, invrs_rho_ds_zt_col, thv_ds_zm_col, thv_ds_zt_col, &
    !$acc              rfrzm_col, &
    !$acc              radf_col, wpthlp_sfc_col, &
    !$acc              wprtp_sfc_col, upwp_sfc_col, vpwp_sfc_col, p_sfc_col, & 
    !$acc              upwp_sfc_pert_col, vpwp_sfc_pert_col, rtm_ref_col, thlm_ref_col, &
    !$acc              um_ref_col, &
    !$acc              vm_ref_col, ug_col, vg_col, host_dx_col, host_dy_col, &
    !$acc              pdf_params, pdf_params_zm ) &
    !$acc        copy( um_col, upwp_col, vm_col, vpwp_col, up2_col, vp2_col, up3_col, vp3_col, &
    !$acc              rtm_col, wprtp_col, thlm_col, wpthlp_col, rtp2_col, &
    !$acc              rtp3_col, thlp2_col, thlp3_col, rtpthlp_col, wp2_col, wp3_col, &
    !$acc              p_in_Pa_col, exner_col, rcm_col, cloud_frac_col, wpthvp_col, wp2thvp_col, &
    !$acc              rtpthvp_col, thlpthvp_col, wp2rtp_col, wp2thlp_col, uprcp_col, vprcp_col, &
    !$acc              rc_coef_zm_col, &
    !$acc              wp4_col, wpup2_col, wpvp2_col, wp2up2_col, wp2vp2_col, &
    !$acc              ice_supersat_frac_col, um_pert_col, &
    !$acc              vm_pert_col, upwp_pert_col, vpwp_pert_col, &
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
    !$acc     copyout( rcm_in_layer_col, cloud_cover_col, wprcp_col, w_up_in_cloud_col, &
    !$acc              w_down_in_cloud_col, &
    !$acc              cloudy_updraft_frac_col, cloudy_downdraft_frac_col, invrs_tau_zm_col, &
    !$acc              Kh_zt_col, &
    !$acc              Kh_zm_col, &
    !$acc              thlprcp_col )

    !$acc data if( sclr_dim > 0 ) &
    !$acc      copyin( sclr_tol, sclrm_forcing_col, wpsclrp_sfc_col ) &
    !$acc        copy( sclrm_col, wpsclrp_col, sclrp2_col, sclrp3_col, sclrprtp_col, &
    !$acc              sclrpthlp_col, sclrpthvp_col )

    !$acc data if( edsclr_dim > 0 ) &
    !$acc      copyin( wpedsclrp_sfc_col, edsclrm_forcing_col ) &
    !$acc        copy( edsclrm_col )

    !$acc data if( hydromet_dim > 0 ) &
    !$acc      copyin( hydromet_col, wphydrometp_col, wp2hmp_col, rtphmp_zt_col, thlphmp_zt_col, &
    !$acc              l_mix_rat_hm )

#ifdef CLUBB_CAM
    !$acc data copyout( qclvar_col )
#endif

#ifdef CLUBBND_CAM
    !$acc data copyin( varmu_col )
#endif

#ifdef GFDL
    !$acc data if( sclr_dim > 0 ) copy( sclrm_trsport_only )
#endif

    call advance_clubb_core( gr, gr%nzm, gr%nzt, 1,             &             ! intent(in)
      l_implemented, dt, fcor_col, sfc_elevation_col,            &            ! intent(in)
      hydromet_dim, &                                                         ! intent(in)
      sclr_dim, sclr_tol, edsclr_dim, sclr_idx, &                         ! intent(in)
      thlm_forcing_col, rtm_forcing_col, um_forcing_col, vm_forcing_col, &    ! intent(in)
      sclrm_forcing_col, edsclrm_forcing_col, wprtp_forcing_col, &            ! intent(in)
      wpthlp_forcing_col, rtp2_forcing_col, thlp2_forcing_col, &              ! intent(in)
      rtpthlp_forcing_col, wm_zm_col, wm_zt_col, &                            ! intent(in)
      wpthlp_sfc_col, wprtp_sfc_col, upwp_sfc_col, vpwp_sfc_col, p_sfc_col, & ! intent(in)
      wpsclrp_sfc_col, wpedsclrp_sfc_col, &                                   ! intent(in)
      upwp_sfc_pert_col, vpwp_sfc_pert_col, &                                 ! intent(in)
      rtm_ref_col, thlm_ref_col, um_ref_col, vm_ref_col, ug_col, vg_col, &    ! Intent(in)
      p_in_Pa_col, rho_zm_col, rho_col, exner_col, &                          ! intent(in)
      rho_ds_zm_col, rho_ds_zt_col, invrs_rho_ds_zm_col, &                    ! intent(in)
      invrs_rho_ds_zt_col, thv_ds_zm_col, thv_ds_zt_col, &                    ! intent(in)
      hydromet_col, l_mix_rat_hm, &                                           ! intent(in)
      rfrzm_col, radf_col, &                                                  ! intent(in)
#ifdef CLUBBND_CAM
      varmu_col, &
#endif
      wphydrometp_col, wp2hmp_col, rtphmp_zt_col, thlphmp_zt_col, &                 ! intent(in)
      host_dx_col, host_dy_col, &                                             ! intent(in)
      clubb_params_col, nu_vert_res_dep, lmin, &                                  ! intent(in)
      clubb_config_flags, &                                                   ! intent(in)
      stats_metadata, &           ! intent(in)
      stats_zt_col, stats_zm_col, stats_sfc_col, &                                ! intent(inout)
      um_col, vm_col, upwp_col, vpwp_col, up2_col, vp2_col, up3_col, vp3_col, &   ! intent(inout)
      thlm_col, rtm_col, wprtp_col, wpthlp_col, &                                 ! intent(inout)
      wp2_col, wp3_col, rtp2_col, rtp3_col, thlp2_col, thlp3_col, rtpthlp_col, &  ! intent(inout)
      sclrm_col,   &
#ifdef GFDL
      sclrm_trsport_only_col,  &  ! h1g, 2010-06-16                               ! intent(inout)
#endif
      sclrp2_col, sclrp3_col, sclrprtp_col, sclrpthlp_col, &                      ! intent(inout)
      wpsclrp_col, edsclrm_col, &                                                 ! intent(inout)
      rcm_col, cloud_frac_col, &                                                  ! intent(inout)
      wpthvp_col, wp2thvp_col, rtpthvp_col, thlpthvp_col, &                       ! intent(inout)
      sclrpthvp_col, &                                                            ! intent(inout)
      wp2rtp_col, wp2thlp_col, uprcp_col, vprcp_col, rc_coef_zm_col, wp4_col, &   ! intent(inout)
      wpup2_col, wpvp2_col, wp2up2_col, wp2vp2_col, ice_supersat_frac_col, &      ! intent(inout)
      um_pert_col, vm_pert_col, upwp_pert_col, vpwp_pert_col, &                   ! intent(inout)
      pdf_params, pdf_params_zm, &                                                ! intent(inout)
      pdf_implicit_coefs_terms, &                                                 ! intent(inout)
#ifdef GFDL
               RH_crit_col, & !h1g, 2010-06-16                                    ! intent(inout)
               do_liquid_only_in_clubb, &                                         ! intent(in)
#endif
      Kh_zm_col, Kh_zt_col, &                                                     ! intent(out)
#ifdef CLUBB_CAM
               qclvar_col, &                                                      ! intent(out)
#endif
      thlprcp_col, wprcp_col, w_up_in_cloud_col, w_down_in_cloud_col, & ! intent(out)
      cloudy_updraft_frac_col, cloudy_downdraft_frac_col, &                       ! intent(out)
      rcm_in_layer_col, cloud_cover_col, invrs_tau_zm_col, &                      ! intent(out)
      err_code_api )                                                              ! intent(out)
    
    !$acc end data
    !$acc end data
    !$acc end data
    !$acc end data
      
#ifdef CLUBB_CAM
    !$acc end data
#endif

#ifdef CLUBBND_CAM
    !$acc end data
#endif

#ifdef GFDL
    !$acc end data
#endif
    
    ! The following does not work for stats 
    !     stats_zt = stats_zt_col(1)
    !     stats_zm = stats_zm_col(1) 
    !     stats_sfc = stats_sfc_col(1)
    ! because of some mysterious pointer issue. However, the only thing that 
    ! updates in stats is the field values, so we can copy only those instead.
    if ( stats_metadata%l_stats ) then 
      stats_zm%accum_field_values = stats_zm_col(1)%accum_field_values
      stats_zm%accum_num_samples = stats_zm_col(1)%accum_num_samples
      
      stats_zt%accum_field_values = stats_zt_col(1)%accum_field_values
      stats_zt%accum_num_samples = stats_zt_col(1)%accum_num_samples
      
      stats_sfc%accum_field_values = stats_sfc_col(1)%accum_field_values
      stats_sfc%accum_num_samples = stats_sfc_col(1)%accum_num_samples
    end if
      
      
    um = um_col(1,:)
    upwp = upwp_col(1,:)
    vm = vm_col(1,:)
    vpwp = vpwp_col(1,:)
    up2 = up2_col(1,:)
    vp2 = vp2_col(1,:)
    up3 = up3_col(1,:)
    vp3 = vp3_col(1,:)
    rtm = rtm_col(1,:)
    wprtp = wprtp_col(1,:)
    thlm = thlm_col(1,:)
    wpthlp = wpthlp_col(1,:)
    rtp2 = rtp2_col(1,:)
    rtp3 = rtp3_col(1,:)
    thlp2 = thlp2_col(1,:)
    thlp3 = thlp3_col(1,:)
    rtpthlp = rtpthlp_col(1,:)
    wp2 = wp2_col(1,:)
    wp3 = wp3_col(1,:)
    sclrm = sclrm_col(1,:,:)
    wpsclrp = wpsclrp_col(1,:,:)
    sclrp2 = sclrp2_col(1,:,:)
    sclrp3 = sclrp3_col(1,:,:)
    sclrprtp = sclrprtp_col(1,:,:)
    sclrpthlp = sclrpthlp_col(1,:,:)
    p_in_Pa = p_in_Pa_col(1,:)
    exner = exner_col(1,:)
    rcm = rcm_col(1,:)
    cloud_frac = cloud_frac_col(1,:)
    wpthvp = wpthvp_col(1,:)
    wp2thvp = wp2thvp_col(1,:)
    rtpthvp = rtpthvp_col(1,:)
    thlpthvp = thlpthvp_col(1,:)
    sclrpthvp = sclrpthvp_col(1,:,:)
    wp2rtp = wp2rtp_col(1,:)
    wp2thlp = wp2thlp_col(1,:)
    uprcp = uprcp_col(1,:)
    vprcp = vprcp_col(1,:)
    rc_coef_zm = rc_coef_zm_col(1,:)
    wp4 = wp4_col(1,:)
    wpup2 = wpup2_col(1,:)
    wpvp2 = wpvp2_col(1,:)
    wp2vp2 = wp2vp2_col(1,:)
    wp2up2 = wp2up2_col(1,:)
    ice_supersat_frac = ice_supersat_frac_col(1,:)
#ifdef GFDL
    sclrm_trsport_only = sclrm_trsport_only_col(1,:,:)
#endif
    edsclrm = edsclrm_col(1,:,:)
    rcm_in_layer = rcm_in_layer_col(1,:)
    cloud_cover  = cloud_cover_col(1,:)
    wprcp = wprcp_col(1,:)
    w_up_in_cloud = w_up_in_cloud_col(1,:)
    w_down_in_cloud = w_down_in_cloud_col(1,:)
    cloudy_updraft_frac = cloudy_updraft_frac_col(1,:)
    cloudy_downdraft_frac = cloudy_downdraft_frac_col(1,:)
    invrs_tau_zm = invrs_tau_zm_col(1,:)
    Kh_zt = Kh_zt_col(1,:)
    Kh_zm = Kh_zm_col(1,:)
#ifdef CLUBB_CAM
    qclvar = qclvar_col(1,:)
#endif
    thlprcp = thlprcp_col(1,:)
#ifdef GFDL
    RH_crit = RH_crit_col(1,:,:,:)
#endif

  end subroutine advance_clubb_core_api_single_col
  
  subroutine advance_clubb_core_api_multi_col( gr, nzm, nzt, ngrdcol, &  ! intent(in)
    l_implemented, dt, fcor, sfc_elevation,            &    ! intent(in)
    hydromet_dim, &                                         ! intent(in)
    sclr_dim, sclr_tol, edsclr_dim, sclr_idx, &         ! intent(in)
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
    hydromet, l_mix_rat_hm, &                               ! intent(in)
    rfrzm, radf, &                                          ! intent(in)
#ifdef CLUBBND_CAM
    varmu, &                                                ! intent(in)
#endif
    wphydrometp, wp2hmp, rtphmp_zt, thlphmp_zt, &           ! intent(in)
    host_dx, host_dy, &                                     ! intent(in)
    clubb_params, nu_vert_res_dep, lmin, &                  ! intent(in)
    clubb_config_flags, &                                   ! intent(in)
    stats_metadata, &                                       ! intent(in)
    stats_zt, stats_zm, stats_sfc, &                        ! intent(inout)
    um, vm, upwp, vpwp, up2, vp2, up3, vp3, &               ! intent(inout)
    thlm, rtm, wprtp, wpthlp, &                             ! intent(inout)
    wp2, wp3, rtp2, rtp3, thlp2, thlp3, rtpthlp, &          ! intent(inout)
    sclrm,   &
#ifdef GFDL
               sclrm_trsport_only,  &  ! h1g, 2010-06-16    ! intent(inout)
#endif
    sclrp2, sclrp3, sclrprtp, sclrpthlp, &                  ! intent(inout)
    wpsclrp, edsclrm, err_code_api, &                       ! intent(inout)
    rcm, cloud_frac, &                                      ! intent(inout)
    wpthvp, wp2thvp, rtpthvp, thlpthvp, &                   ! intent(inout)
    sclrpthvp, &                                            ! intent(inout)
    wp2rtp, wp2thlp, uprcp, vprcp, rc_coef_zm, wp4, &       ! intent(inout)
    wpup2, wpvp2, wp2up2, wp2vp2, ice_supersat_frac, &      ! intent(inout)
    um_pert, vm_pert, upwp_pert, vpwp_pert, &               ! intent(inout)
    pdf_params, pdf_params_zm, &                            ! intent(inout)
    pdf_implicit_coefs_terms, &                             ! intent(inout)
#ifdef GFDL
               RH_crit, & !h1g, 2010-06-16                  ! intent(inout)
               do_liquid_only_in_clubb, &                   ! intent(in)
#endif
    Kh_zm, Kh_zt, &                                         ! intent(out)
#ifdef CLUBB_CAM
    qclvar, &                                               ! intent(out)
#endif
    thlprcp, wprcp, w_up_in_cloud, w_down_in_cloud, &       ! intent(out)
    cloudy_updraft_frac, cloudy_downdraft_frac, &           ! intent(out)
    rcm_in_layer, cloud_cover, invrs_tau_zm )               ! intent(out)

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
    logical, intent(in) ::  &
      l_implemented ! Is this part of a larger host model (T/F) ?

    real( kind = core_rknd ), intent(in) ::  &
      dt  ! Current timestep duration    [s]
      
    real( kind = core_rknd ), intent(in), dimension(ngrdcol) ::  &
      fcor, &           ! Coriolis forcing             [s^-1]
      sfc_elevation     ! Elevation of ground level    [m AMSL]

    integer, intent(in) :: &
      hydromet_dim,   & ! Total number of hydrometeor species       [#]
      sclr_dim,       & ! Number of passive scalars                 [#]
      edsclr_dim        ! Number of eddy-diff. passive scalars      [#]

    real( kind = core_rknd ), intent(in), dimension(sclr_dim) :: & 
      sclr_tol          ! Threshold(s) on the passive scalars  [units vary]

    type (sclr_idx_type), intent(in) :: &
      sclr_idx

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzt) ::  &
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

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzm) ::  &
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

    real( kind = core_rknd ), dimension(ngrdcol,nzt,hydromet_dim), intent(in) :: &
      hydromet           ! Collection of hydrometeors                [units vary]

    logical, dimension(hydromet_dim), intent(in) :: &
      l_mix_rat_hm   ! if true, then the quantity is a hydrometeor mixing ratio

    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: &
      radf          ! Buoyancy production at the CL top due to LW radiative cooling [m^2/s^3]

#ifdef CLUBBND_CAM 
    real( kind = core_rknd ), intent(in), dimension(ngrdcol) :: & 
      varmu 
#endif 

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

    real( kind = core_rknd ), intent(in) :: &
      lmin    ! Min. value for the length scale    [m]

    type( clubb_config_flags_type ), intent(in) :: &
      clubb_config_flags ! Derived type holding all configurable CLUBB flags

    type (stats_metadata_type), intent(in) :: &
      stats_metadata

    !------------------------ Input/Output Variables ------------------------
    type(stats), intent(inout), dimension(ngrdcol) :: &
      stats_zt, &
      stats_zm, &
      stats_sfc

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
      wp2thvp       ! < w'^2 th_v' > (thermodynamic levels)          [m^2/s^2 K]

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

#ifdef GFDL
    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,nzt,sclr_dim) :: &  ! h1g, 2010-06-16
      sclrm_trsport_only  ! Passive scalar concentration due to pure transport [{units vary}/s]
#endif

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,nzt,edsclr_dim) :: &
    edsclrm   ! Eddy passive scalar mean (thermo. levels)   [units vary]


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

    integer, intent(inout) :: err_code_api ! Diagnostic, for if some calculation goes amiss.

#ifdef GFDL
    ! hlg, 2010-06-16
    real( kind = core_rknd ), intent(inOUT), dimension(ngrdcol,nzt, min(1,sclr_dim) , 2) :: &
      RH_crit  ! critical relative humidity for droplet and ice nucleation
    logical, intent(in)                 ::  do_liquid_only_in_clubb
#endif

    !$acc data copyin( gr, gr%zm, gr%zt, gr%dzm, gr%dzt, gr%invrs_dzt, gr%invrs_dzm, &
    !$acc              gr%weights_zt2zm, gr%weights_zm2zt, &
    !$acc              nu_vert_res_dep, nu_vert_res_dep%nu2, nu_vert_res_dep%nu9, &
    !$acc              nu_vert_res_dep%nu1, nu_vert_res_dep%nu8, nu_vert_res_dep%nu10, &
    !$acc              nu_vert_res_dep%nu6, &
    !$acc              sclr_idx, clubb_params, &
    !$acc              fcor, sfc_elevation, thlm_forcing, rtm_forcing, um_forcing, &
    !$acc              vm_forcing, wprtp_forcing, wpthlp_forcing, rtp2_forcing, thlp2_forcing, &
    !$acc              rtpthlp_forcing, wm_zm, wm_zt, rho_zm, rho, rho_ds_zm, rho_ds_zt, &
    !$acc              invrs_rho_ds_zm, invrs_rho_ds_zt, thv_ds_zm, thv_ds_zt, rfrzm, &
    !$acc              radf, wpthlp_sfc, &
    !$acc              wprtp_sfc, upwp_sfc, vpwp_sfc, p_sfc, & 
    !$acc              upwp_sfc_pert, vpwp_sfc_pert, rtm_ref, thlm_ref, um_ref, &
    !$acc              vm_ref, ug, vg, host_dx, host_dy, &
    !$acc              pdf_params, pdf_params_zm ) &
    !$acc        copy( um, upwp, vm, vpwp, up2, vp2, up3, vp3, rtm, wprtp, thlm, wpthlp, rtp2, &
    !$acc              rtp3, thlp2, thlp3, rtpthlp, wp2, wp3, &
    !$acc              p_in_Pa, exner, rcm, cloud_frac, wpthvp, wp2thvp, &
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
    !$acc              Kh_zm, &
    !$acc              thlprcp )

    !$acc data if( sclr_dim > 0 ) &
    !$acc      copyin( sclr_tol, sclrm_forcing, wpsclrp_sfc ) &
    !$acc        copy( sclrm, wpsclrp, sclrp2, sclrp3, sclrprtp, sclrpthlp, sclrpthvp )

    !$acc data if( edsclr_dim > 0 ) &
    !$acc      copyin( wpedsclrp_sfc, edsclrm_forcing ) &
    !$acc        copy( edsclrm )

    !$acc data if( hydromet_dim > 0 ) &
    !$acc      copyin( hydromet, wphydrometp, wp2hmp, rtphmp_zt, thlphmp_zt, &
    !$acc              l_mix_rat_hm )

#ifdef CLUBB_CAM
    !$acc data copyout( qclvar )
#endif

#ifdef CLUBBND_CAM
    !$acc data copyin( varmu )
#endif

#ifdef GFDL
    !$acc data if( sclr_dim > 0 ) copy( sclrm_trsport_only )
#endif

    call advance_clubb_core( gr, nzm, nzt, ngrdcol, &         ! intent(in)
      l_implemented, dt, fcor, sfc_elevation,            &    ! intent(in)
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
      hydromet, l_mix_rat_hm,&                                ! intent(in)
      rfrzm, radf, &                                          ! intent(in)
#ifdef CLUBBND_CAM
      varmu, &
#endif
      wphydrometp, wp2hmp, rtphmp_zt, thlphmp_zt, &           ! intent(in)
      host_dx, host_dy, &                                     ! intent(in)
      clubb_params, nu_vert_res_dep, lmin, &                  ! intent(in)
      clubb_config_flags, &                                   ! intent(in)
      stats_metadata, &                                       ! intent(in)
      stats_zt, stats_zm, stats_sfc, &                        ! intent(inout)
      um, vm, upwp, vpwp, up2, vp2, up3, vp3, &               ! intent(inout)
      thlm, rtm, wprtp, wpthlp, &                             ! intent(inout)
      wp2, wp3, rtp2, rtp3, thlp2, thlp3, rtpthlp, &          ! intent(inout)
      sclrm,   &
#ifdef GFDL
               sclrm_trsport_only,  &  ! h1g, 2010-06-16      ! intent(inout)
#endif
      sclrp2, sclrp3, sclrprtp, sclrpthlp, &                  ! intent(inout)
      wpsclrp, edsclrm, &                                     ! intent(inout)
      rcm, cloud_frac, &                                      ! intent(inout)
      wpthvp, wp2thvp, rtpthvp, thlpthvp, &                   ! intent(inout)
      sclrpthvp, &                                            ! intent(inout)
      wp2rtp, wp2thlp, uprcp, vprcp, rc_coef_zm, wp4, &       ! intent(inout)
      wpup2, wpvp2, wp2up2, wp2vp2, ice_supersat_frac, &      ! intent(inout)
      um_pert, vm_pert, upwp_pert, vpwp_pert, &               ! intent(inout)
      pdf_params, pdf_params_zm, &                            ! intent(inout)
      pdf_implicit_coefs_terms, &                             ! intent(inout)
#ifdef GFDL
               RH_crit, & !h1g, 2010-06-16                    ! intent(inout)
               do_liquid_only_in_clubb, &                     ! intent(in)
#endif
      Kh_zm, Kh_zt, &                                         ! intent(out)
#ifdef CLUBB_CAM
               qclvar, &                                      ! intent(out)
#endif
      thlprcp, wprcp, w_up_in_cloud, w_down_in_cloud, &       ! intent(out)
      cloudy_updraft_frac, cloudy_downdraft_frac, &           ! intent(out)
      rcm_in_layer, cloud_cover, invrs_tau_zm, &              ! intent(out)
      err_code_api )                                          ! intent(out)

    !$acc end data
    !$acc end data
    !$acc end data
    !$acc end data

#ifdef CLUBB_CAM
    !$acc end data
#endif

#ifdef CLUBBND_CAM
    !$acc end data
#endif

#ifdef GFDL
    !$acc end data
#endif


  end subroutine advance_clubb_core_api_multi_col

  !================================================================================================
  ! cleanup_clubb_core_api - Frees memory used by the model.
  !================================================================================================

  subroutine cleanup_clubb_core_api( gr )

    use advance_clubb_core_module, only : cleanup_clubb_core

    implicit none

    type(grid), intent(inout) :: gr
    call cleanup_clubb_core( gr ) ! intent(inout)

  end subroutine cleanup_clubb_core_api

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
  ! compute_current_date - Computes the current date and the seconds since that date.
  !================================================================================================

  subroutine compute_current_date_api( &
    previous_day, previous_month, &
    previous_year,  &
    seconds_since_previous_date, &
    current_day, current_month, &
    current_year, &
    seconds_since_current_date )

    use calendar, only : compute_current_date

    implicit none

    ! Previous date
    integer, intent(in) :: &
      previous_day,    & ! Day of the month      [dd]
      previous_month,  & ! Month of the year     [mm]
      previous_year      ! Year                  [yyyy]

    real(kind=time_precision), intent(in) :: &
      seconds_since_previous_date ! [s]

    ! Output Variable(s)

    ! Current date
    integer, intent(out) :: &
      current_day,     & ! Day of the month      [dd]
      current_month,   & ! Month of the year     [mm]
      current_year       ! Year                  [yyyy]

    real(kind=time_precision), intent(out) :: &
      seconds_since_current_date

    call compute_current_date( &
      previous_day, previous_month, & ! intent(in)
      previous_year,  & ! intent(in)
      seconds_since_previous_date, & ! intent(in)
      current_day, current_month, & ! intent(out)
      current_year, & ! intent(out)
      seconds_since_current_date ) ! intent(out)
  end subroutine compute_current_date_api

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
  ! set_clubb_debug_level - Controls the importance of error messages sent to the console.
  !================================================================================================

  subroutine set_clubb_debug_level_api( &
    level )

    use error_code, only: &
        set_clubb_debug_level ! Procedure

    implicit none

    ! Input variable
    integer, intent(in) :: level ! The debug level being checked against the current setting

    call set_clubb_debug_level( &
      level ) ! intent(in)
  end subroutine set_clubb_debug_level_api

  !================================================================================================
  ! clubb_at_least_debug_level - Checks to see if clubb has been set to a specified debug level.
  !================================================================================================

  logical function clubb_at_least_debug_level_api( &
    level )
    
    use error_code, only: &
        clubb_at_least_debug_level ! Procedure

    implicit none

    ! Input variable
    integer, intent(in) :: level   ! The debug level being checked against the current setting

    clubb_at_least_debug_level_api = clubb_at_least_debug_level( &
      level )
  end function clubb_at_least_debug_level_api

  !================================================================================================
  ! vertical_integral - Computes the vertical integral.
  !================================================================================================

  function vertical_integral_api( &
    total_idx, rho_ds, &
    field, dz )

    use advance_helper_module, only : vertical_integral

    implicit none

    ! Input variables
    integer, intent(in) :: &
      total_idx  ! The total numer of indices within the range of averaging

    real( kind = core_rknd ), dimension(total_idx), intent(in) ::  &
      rho_ds,  & ! Dry, static density                   [kg/m^3]
      field,   & ! The field to be vertically averaged   [Units vary]
      dz         ! Level thickness                       [1/m]
    ! Note:  The rho_ds and field points need to be arranged from
    !        lowest to highest in altitude, with rho_ds(1) and
    !        field(1) actually their respective values at level k = begin_idx.

    real( kind = core_rknd ) :: &
      vertical_integral_api ! Integral in the numerator (see description)

    vertical_integral_api = vertical_integral( &
      total_idx, rho_ds, &
      field, dz )
  end function vertical_integral_api

  !================================================================================================
  ! setup_grid_heights - Sets the heights and interpolation weights of the column.
  !================================================================================================

  subroutine setup_grid_heights_api_single_col( &
    l_implemented, grid_type,  &
    deltaz, zm_init, momentum_heights,  &
    gr, thermodynamic_heights )

    use grid_class, only: & 
        grid, & ! Type
        setup_grid_heights
    
    use error_code, only : &
        clubb_fatal_error       ! Constant

    implicit none

    ! Input Variables
   
    type(grid), intent(inout) :: gr

    ! Flag to see if CLUBB is running on it's own,
    ! or if it's implemented as part of a host model.
    logical, intent(in) :: l_implemented

    ! If CLUBB is running on it's own, this option determines if it is using:
    ! 1) an evenly-spaced grid;
    ! 2) a stretched (unevenly-spaced) grid entered on the thermodynamic grid
    !    levels (with momentum levels set halfway between thermodynamic levels);
    !    or
    ! 3) a stretched (unevenly-spaced) grid entered on the momentum grid levels
    !    (with thermodynamic levels set halfway between momentum levels).
    integer, intent(in) :: grid_type

    ! If the CLUBB model is running by itself, and is using an evenly-spaced
    ! grid (grid_type = 1), it needs the vertical grid spacing and
    ! momentum-level starting altitude as input.
    real( kind = core_rknd ), intent(in) ::  &
      deltaz,   & ! Vertical grid spacing                  [m]
      zm_init     ! Initial grid altitude (momentum level) [m]


    ! If the CLUBB parameterization is implemented in a host model, it needs to
    ! use the host model's momentum level altitudes and thermodynamic level
    ! altitudes.
    ! If the CLUBB model is running by itself, but is using a stretched grid
    ! entered on thermodynamic levels (grid_type = 2), it needs to use the
    ! thermodynamic level altitudes as input.
    ! If the CLUBB model is running by itself, but is using a stretched grid
    ! entered on momentum levels (grid_type = 3), it needs to use the momentum
    ! level altitudes as input.
    real( kind = core_rknd ), intent(in), dimension(gr%nzm) ::  &
      momentum_heights      ! Momentum level altitudes (input)      [m]

    real( kind = core_rknd ), intent(in), dimension(gr%nzt) ::  &
      thermodynamic_heights ! Thermodynamic level altitudes (input) [m]
      
    ! ------------------- Local Variables -------------------
    
    real( kind = core_rknd ), dimension(1) ::  &
      deltaz_col,   & ! Vertical grid spacing                  [m]
      zm_init_col     ! Initial grid altitude (momentum level) [m]
      
    real( kind = core_rknd ), dimension(1,gr%nzm) ::  &
      momentum_heights_col      ! Momentum level altitudes (input)      [m]

    real( kind = core_rknd ), dimension(1,gr%nzt) ::  &
      thermodynamic_heights_col ! Thermodynamic level altitudes (input) [m]
      
    deltaz_col(1) = deltaz
    zm_init_col(1) = zm_init
    momentum_heights_col(1,:) = momentum_heights
    thermodynamic_heights_col(1,:) = thermodynamic_heights

    call setup_grid_heights( &
      gr%nzm, gr%nzt, 1, & ! intent(in)
      l_implemented, grid_type,  & ! intent(in)
      deltaz_col, zm_init_col, momentum_heights_col,  & ! intent(in)
      thermodynamic_heights_col, & ! intent(in)
      gr ) ! intent(inout)

    if ( err_code == clubb_fatal_error ) error stop

  end subroutine setup_grid_heights_api_single_col
  
  !================================================================================================
  ! setup_grid_heights - Sets the heights and interpolation weights of the column.
  !================================================================================================

  subroutine setup_grid_heights_api_multi_col( &
      nzm, nzt, ngrdcol, &
      l_implemented, grid_type,  &
      deltaz, zm_init, momentum_heights,  &
      gr, thermodynamic_heights )

    use grid_class, only: & 
        grid, & ! Type
        setup_grid_heights
    
    use error_code, only : &
        clubb_fatal_error       ! Constant

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nzm, &
      nzt, &
      ngrdcol

    type (grid), intent(inout) :: gr

    ! Flag to see if CLUBB is running on it's own,
    ! or if it's implemented as part of a host model.
    logical, intent(in) :: l_implemented

    ! If CLUBB is running on it's own, this option determines if it is using:
    ! 1) an evenly-spaced grid;
    ! 2) a stretched (unevenly-spaced) grid entered on the thermodynamic grid
    !    levels (with momentum levels set halfway between thermodynamic levels);
    !    or
    ! 3) a stretched (unevenly-spaced) grid entered on the momentum grid levels
    !    (with thermodynamic levels set halfway between momentum levels).
    integer, intent(in) :: grid_type

    ! If the CLUBB model is running by itself, and is using an evenly-spaced
    ! grid (grid_type = 1), it needs the vertical grid spacing and
    ! momentum-level starting altitude as input.
    real( kind = core_rknd ), dimension(ngrdcol), intent(in) ::  &
      deltaz,   & ! Vertical grid spacing                  [m]
      zm_init     ! Initial grid altitude (momentum level) [m]


    ! If the CLUBB parameterization is implemented in a host model, it needs to
    ! use the host model's momentum level altitudes and thermodynamic level
    ! altitudes.
    ! If the CLUBB model is running by itself, but is using a stretched grid
    ! entered on thermodynamic levels (grid_type = 2), it needs to use the
    ! thermodynamic level altitudes as input.
    ! If the CLUBB model is running by itself, but is using a stretched grid
    ! entered on momentum levels (grid_type = 3), it needs to use the momentum
    ! level altitudes as input.
    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzm) ::  &
      momentum_heights      ! Momentum level altitudes (input)      [m]

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzt) ::  &
      thermodynamic_heights ! Thermodynamic level altitudes (input) [m]

    call setup_grid_heights( &
      nzm, nzt, ngrdcol, & ! intent(in)
      l_implemented, grid_type,  & ! intent(in)
      deltaz, zm_init, momentum_heights,  & ! intent(in)
      thermodynamic_heights, & ! intent(in)
      gr ) ! intent(inout)

    if ( err_code == clubb_fatal_error ) error stop

  end subroutine setup_grid_heights_api_multi_col
  
  !================================================================================================
  ! setup_grid - This subroutine sets up the CLUBB vertical grid for a single column
  !================================================================================================
  subroutine setup_grid_api_single_col( nzmax, sfc_elevation, l_implemented, &
                                        grid_type, deltaz, zm_init, zm_top, &
                                        momentum_heights, thermodynamic_heights, &
                                        gr )
                            
    use grid_class, only: & 
        grid, & ! Type
        setup_grid

    implicit none

    ! Input Variables
    integer, intent(in) ::  & 
      nzmax  ! Number of vertical levels in grid      [#]

    type(grid), intent(inout) :: gr

    real( kind = core_rknd ), intent(in) ::  &
      sfc_elevation  ! Elevation of ground level    [m AMSL]
      
    logical, intent(in) :: l_implemented
    
    integer, intent(in) :: grid_type
    
    real( kind = core_rknd ), intent(in) ::  & 
      deltaz,   & ! Vertical grid spacing                  [m]
      zm_init,  & ! Initial grid altitude (momentum level) [m]
      zm_top      ! Maximum grid altitude (momentum level) [m]
      
    real( kind = core_rknd ), intent(in), dimension(nzmax) ::  & 
      momentum_heights      ! Momentum level altitudes (input)      [m]

    real( kind = core_rknd ), intent(in), dimension(nzmax-1) ::  & 
      thermodynamic_heights ! Thermodynamic level altitudes (input) [m]

    real( kind = core_rknd ), dimension(1) ::  &
      sfc_elevation_col  ! Elevation of ground level    [m AMSL]
      
    real( kind = core_rknd ), dimension(1) ::  & 
      deltaz_col,   & ! Vertical grid spacing                  [m]
      zm_init_col,  & ! Initial grid altitude (momentum level) [m]
      zm_top_col      ! Maximum grid altitude (momentum level) [m]
      
    real( kind = core_rknd ), dimension(1,nzmax) ::  & 
      momentum_heights_col      ! Momentum level altitudes (input)      [m]

    real( kind = core_rknd ), dimension(1,nzmax-1) ::  & 
      thermodynamic_heights_col ! Thermodynamic level altitudes (input) [m]

    sfc_elevation_col(1)            = sfc_elevation
    deltaz_col(1)                   = deltaz
    zm_init_col(1)                  = zm_init
    zm_top_col(1)                   = zm_top
    momentum_heights_col(1,:)       = momentum_heights
    thermodynamic_heights_col(1,:)  = thermodynamic_heights

    call setup_grid( nzmax, 1, sfc_elevation_col, l_implemented,      & ! intent(in)
                     grid_type, deltaz_col, zm_init_col, zm_top_col,  & ! intent(in)
                     momentum_heights_col, thermodynamic_heights_col, & ! intent(in)
                     gr )                                               ! intent(out)

  end subroutine setup_grid_api_single_col
  
  !================================================================================================
  ! setup_grid - This subroutine sets up the CLUBB vertical grid.
  !================================================================================================
  subroutine setup_grid_api_multi_col( nzmax, ngrdcol, sfc_elevation, l_implemented, &
                                       grid_type, deltaz, zm_init, zm_top, &
                                       momentum_heights, thermodynamic_heights, &
                                       gr )
                            
    use grid_class, only: & 
        grid, & ! Type
        setup_grid

    implicit none

    ! Input Variables
    integer, intent(in) ::  & 
      nzmax,  & ! Number of vertical levels in grid      [#]
      ngrdcol   ! Number of grid columns                 [#]
      
    type (grid), intent(inout) :: gr

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) ::  &
      sfc_elevation  ! Elevation of ground level    [m AMSL]
      
    logical, intent(in) :: l_implemented
    
    integer, intent(in) :: grid_type
    
    real( kind = core_rknd ), dimension(ngrdcol), intent(in) ::  & 
      deltaz,   & ! Vertical grid spacing                  [m]
      zm_init,  & ! Initial grid altitude (momentum level) [m]
      zm_top      ! Maximum grid altitude (momentum level) [m]
      
    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzmax) ::  & 
      momentum_heights      ! Momentum level altitudes (input)      [m]

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzmax-1) ::  & 
      thermodynamic_heights ! Thermodynamic level altitudes (input) [m]


    call setup_grid( nzmax, ngrdcol, sfc_elevation, l_implemented, & ! intent(in)
                     grid_type, deltaz, zm_init, zm_top,           & ! intent(in)
                     momentum_heights, thermodynamic_heights,      & ! intent(in)
                     gr )                                            ! intent(out)

  end subroutine setup_grid_api_multi_col

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
  ! lin_interpolate_on_grid - Linear interpolation for 25 June 1996 altocumulus case.
  !================================================================================================

  subroutine lin_interpolate_on_grid_api( &
    nparam, xlist, tlist, xvalue, tvalue )

    use interpolation, only : lin_interpolate_on_grid

    implicit none

    ! Input Variables
    integer, intent(in) :: nparam ! Number of parameters in xlist and tlist

    ! Input/Output Variables
    real( kind = core_rknd ), intent(inout), dimension(nparam) ::  &
      xlist,  & ! List of x-values (independent variable)
      tlist     ! List of t-values (dependent variable)

    real( kind = core_rknd ), intent(in) ::  &
      xvalue  ! x-value at which to interpolate

    real( kind = core_rknd ), intent(inout) ::  &
      tvalue  ! t-value solved by interpolation

    call lin_interpolate_on_grid( &
      nparam, xlist, tlist, xvalue, & ! intent(in)
      tvalue ) ! intent(inout)

  end subroutine lin_interpolate_on_grid_api

  !================================================================================================
  ! read_parameters - Read a namelist containing the model parameters.
  !================================================================================================

  subroutine read_parameters_api( ngrdcol, iunit, filename, &
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
                                  wpxp_Ri_exp, a3_coef_min, a_const, bv_efold, z_displace, &
                                  clubb_params )

    use parameters_tunable, only : read_parameters

    use parameter_indices, only:  &
        nparams ! Variable(s)

    implicit none

    ! Input variables
    integer, intent(in) :: &
      ngrdcol

    integer, intent(in) :: iunit

    character(len=*), intent(in) :: filename

    ! Input/Output variables
    real( kind = core_rknd ), intent(inout) :: & 
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
      wpxp_Ri_exp, a3_coef_min, a_const, bv_efold, z_displace

    ! Output variables
    real( kind = core_rknd ), intent(out), dimension(nparams) :: clubb_params

    call read_parameters( ngrdcol, iunit, filename, & ! intent(in)
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
                          wpxp_Ri_exp, a3_coef_min, a_const, bv_efold, z_displace, &
                          clubb_params ) ! intent(out)

  end subroutine read_parameters_api

  !================================================================================================
  ! setup_parameters - Sets up model parameters for a single column
  !================================================================================================
  subroutine setup_parameters_api_single_col( &
             deltaz, clubb_params, gr, grid_type, &
             l_prescribed_avg_deltaz, &
             lmin, nu_vert_res_dep, err_code_api )

    use parameters_tunable, only: &
        setup_parameters

    use parameter_indices, only:  &
        nparams ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) ::  &
      deltaz  ! Change per height level        [m]

    real( kind = core_rknd ), intent(in), dimension(nparams) :: &
      clubb_params  ! Tuneable model parameters      [-]

    ! Grid definition
    type(grid), intent(in) :: &
      gr

    ! If CLUBB is running on its own, this option determines
    ! if it is using:
    ! 1) an evenly-spaced grid,
    ! 2) a stretched (unevenly-spaced) grid entered on the
    !    thermodynamic grid levels (with momentum levels set
    !    halfway between thermodynamic levels), or
    ! 3) a stretched (unevenly-spaced) grid entered on the
    !    momentum grid levels (with thermodynamic levels set
    !    halfway between momentum levels).
    integer, intent(in) :: grid_type

    logical, intent(in) :: &
      l_prescribed_avg_deltaz ! used in adj_low_res_nu. If .true., avg_deltaz = deltaz

    ! Output Variables 
    real( kind = core_rknd ), intent(out) :: &
      lmin    ! Min. value for the length scale    [m]

    type(nu_vertical_res_dep), intent(out) :: &
      nu_vert_res_dep    ! Vertical resolution dependent nu values

    integer, intent(out) ::  & 	 	      
      err_code_api ! Error condition 
      
      
    real( kind = core_rknd ), dimension(1) ::  &
      deltaz_col  ! Change per height level        [m]

    deltaz_col(1) = deltaz

    call setup_parameters( & 
              deltaz_col, clubb_params, gr, 1, grid_type, &
              l_prescribed_avg_deltaz, &
              lmin, nu_vert_res_dep, err_code_api )

  end subroutine setup_parameters_api_single_col
  
  !================================================================================================
  ! setup_parameters - Sets up model parameters.
  !================================================================================================
  
  subroutine setup_parameters_api_multi_col( &
             deltaz, clubb_params, gr, ngrdcol, grid_type, &
             l_prescribed_avg_deltaz, &
             lmin, nu_vert_res_dep, err_code_api )

    use parameters_tunable, only: &
        setup_parameters

    use parameter_indices, only:  &
        nparams ! Variable(s)

    implicit none

    ! Input Variables

    ! Grid definition
    type(grid), intent(in) :: &
      gr

    integer, intent(in) :: &
      ngrdcol   ! Number of grid columns          [#]

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) ::  &
      deltaz  ! Change per height level        [m]

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nparams) :: &
      clubb_params  ! Tuneable model parameters      [-]
      
    ! If CLUBB is running on its own, this option determines
    ! if it is using:
    ! 1) an evenly-spaced grid,
    ! 2) a stretched (unevenly-spaced) grid entered on the
    !    thermodynamic grid levels (with momentum levels set
    !    halfway between thermodynamic levels), or
    ! 3) a stretched (unevenly-spaced) grid entered on the
    !    momentum grid levels (with thermodynamic levels set
    !    halfway between momentum levels).
    integer, intent(in) :: grid_type

    logical, intent(in) :: &
      l_prescribed_avg_deltaz ! used in adj_low_res_nu. If .true., avg_deltaz = deltaz

    ! Output Variables 
    real( kind = core_rknd ), intent(out) :: &
      lmin    ! Min. value for the length scale    [m]

    type(nu_vertical_res_dep), intent(out) :: &
      nu_vert_res_dep    ! Vertical resolution dependent nu values

    integer, intent(out) ::  & 	 	      
      err_code_api ! Error condition 

    call setup_parameters( & 
              deltaz, clubb_params, gr, ngrdcol, grid_type, &
              l_prescribed_avg_deltaz, &
              lmin, nu_vert_res_dep, err_code_api )

  end subroutine setup_parameters_api_multi_col

  !================================================================================================
  ! adj_low_res_nu - Adjusts values of background eddy diffusivity based on vertical grid spacing.
  !================================================================================================

  subroutine adj_low_res_nu_api_single_col( gr, grid_type, deltaz,    & ! Intent(in)
                                            clubb_params,             & ! Intent(in)
                                            l_prescribed_avg_deltaz,  & ! Intent(in)
                                            nu_vert_res_dep )           ! Intent(out)

    use parameters_tunable, only : adj_low_res_nu

    implicit none

    ! Input Variables

    ! Grid definition
    type(grid), intent(in) :: &
      gr

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

    real( kind = core_rknd ), intent(in) ::  &
      deltaz  ! Change per height level        [m]

    real( kind = core_rknd ), intent(in), dimension(nparams) :: & 
      clubb_params  ! Tuneable model parameters      [-]

    logical, intent(in) :: &
      l_prescribed_avg_deltaz ! used in adj_low_res_nu. If .true., avg_deltaz = deltaz

    ! Output Variables
    type(nu_vertical_res_dep), intent(out) :: &
      nu_vert_res_dep    ! Vertical resolution dependent nu values
      
    ! Local variables
    real( kind = core_rknd ), dimension(1) ::  &
      deltaz_col  ! Change per height level        [m]

    real( kind = core_rknd ), dimension(1,nparams) :: & 
      clubb_params_col  ! Tuneable model parameters      [-]
      
    deltaz_col(1) = deltaz
    clubb_params_col(1,:) = clubb_params

    call adj_low_res_nu( gr, 1, grid_type, deltaz_col,  & ! Intent(in)
                         clubb_params,                  & ! Intent(in)
                         l_prescribed_avg_deltaz,       & ! Intent(in)
                         nu_vert_res_dep )                ! Intent(out)

  end subroutine adj_low_res_nu_api_single_col
  
  !================================================================================================
  ! adj_low_res_nu - Adjusts values of background eddy diffusivity based on vertical grid spacing.
  !================================================================================================

  subroutine adj_low_res_nu_api_multi_col( gr, ngrdcol, grid_type, deltaz,  & ! Intent(in)
                                           clubb_params,                    & ! Intent(in)
                                           l_prescribed_avg_deltaz,         & ! Intent(in)
                                           nu_vert_res_dep )                  ! Intent(out)

    use parameters_tunable, only : adj_low_res_nu

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
      l_prescribed_avg_deltaz ! used in adj_low_res_nu. If .true., avg_deltaz = deltaz

    ! Output Variables
    type(nu_vertical_res_dep), intent(out) :: &
      nu_vert_res_dep    ! Vertical resolution dependent nu values

    call adj_low_res_nu( gr, ngrdcol, grid_type, deltaz, & ! Intent(in)
                         clubb_params,                   & ! Intent(in)
                         l_prescribed_avg_deltaz,        & ! Intent(in)
                         nu_vert_res_dep )                 ! Intent(out)

  end subroutine adj_low_res_nu_api_multi_col

! The CLUBB_CAM preprocessor directives are being commented out because this
! code is now also used for WRF-CLUBB.
!#ifdef CLUBB_CAM /* Code for storing pdf_parameter structs in pbuf as array */
  !================================================================================================
  ! pack_pdf_params - Returns a two dimensional real array with all values.
  !================================================================================================

  subroutine pack_pdf_params_api( pdf_params, nz, r_param_array, &
                                  k_start, k_end )

    use pdf_parameter_module, only : pack_pdf_params

    !use statements

    implicit none

    integer, intent(in) :: nz ! Num Vert Model Levs
    
    ! Input a pdf_parameter array with nz instances of pdf_parameter
    type (pdf_parameter), intent(inout) :: pdf_params

    ! Output a two dimensional real array with all values
    real (kind = core_rknd), dimension(nz,num_pdf_params), intent(inout) :: &
      r_param_array
      
    integer, optional, intent(in) :: k_start, k_end
      
    if( present( k_start ) .and. present( k_end ) ) then
        call pack_pdf_params( pdf_params, nz, & ! intent(in)
                              r_param_array, & ! intent(out)
                              k_start, k_end ) ! intent(in/optional)
    else 
        call pack_pdf_params( pdf_params, nz, & ! intent(in)
                              r_param_array ) ! intent(out)
    end if

  end subroutine pack_pdf_params_api

  !================================================================================================
  ! unpack_pdf_params - Returns a pdf_parameter array with nz instances of pdf_parameter.
  !================================================================================================

  subroutine unpack_pdf_params_api( r_param_array, nz, pdf_params, &
                                    k_start, k_end )

    use pdf_parameter_module, only : unpack_pdf_params

    implicit none
    
    integer, intent(in) :: nz ! Num Vert Model Levs
    
    ! Input a two dimensional real array with pdf values
    real (kind = core_rknd), dimension(nz,num_pdf_params), intent(in) :: &
      r_param_array

    ! Output a pdf_parameter array with nz instances of pdf_parameter
    type (pdf_parameter), intent(inout) :: pdf_params
    
    integer, optional, intent(in) :: k_start, k_end

    if( present( k_start ) .and. present( k_end ) ) then
        call unpack_pdf_params( r_param_array, nz, & ! intent(in)
                                pdf_params, & ! intent(inout)
                                k_start, k_end ) ! intent(in/optional)
    else  
        call unpack_pdf_params( r_param_array, nz, & ! intent(in)
                                pdf_params ) ! intent(inout)
    end if
    
    
  end subroutine unpack_pdf_params_api
  
!#endif
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
  ! copy_single_pdf_params_to_multi - copies values of a single column version of pdf_params
  !   to a multiple column version for a specified column.
  !
  ! NOTE: THIS SUBROUTINE IS INTENDED TO BE TEMPORARY AND SHOULD BECOME UNNECESSARY ONCE 
  !       CLUBB IS ABLE TO OPERATE OVER MULTIPLE COLUMNS.
  !       See https://github.com/larson-group/cam/issues/129#issuecomment-827944454
  !================================================================================================
  subroutine copy_single_pdf_params_to_multi( pdf_params_single, icol, &
                                              pdf_params_multi )
    
    implicit none

    ! Input Variable(s)
    integer, intent(in) :: &
      icol   ! Column number to copy to
      
    type(pdf_parameter), intent(in) :: &
      pdf_params_single  ! PDF parameters            [units vary]

    ! Output Variable(s)
    type(pdf_parameter), intent(inout) :: &
      pdf_params_multi  ! PDF parameters            [units vary]
      
    pdf_params_multi%w_1(icol,:)                  = pdf_params_single%w_1(1,:) 
    pdf_params_multi%w_2(icol,:)                  = pdf_params_single%w_2(1,:) 
    pdf_params_multi%varnce_w_1(icol,:)           = pdf_params_single%varnce_w_1(1,:) 
    pdf_params_multi%varnce_w_2(icol,:)           = pdf_params_single%varnce_w_2(1,:) 
    pdf_params_multi%rt_1(icol,:)                 = pdf_params_single%rt_1(1,:) 
    pdf_params_multi%rt_2(icol,:)                 = pdf_params_single%rt_2(1,:) 
    pdf_params_multi%varnce_rt_1(icol,:)          = pdf_params_single%varnce_rt_1(1,:) 
    pdf_params_multi%varnce_rt_2(icol,:)          = pdf_params_single%varnce_rt_2(1,:) 
    pdf_params_multi%thl_1(icol,:)                = pdf_params_single%thl_1(1,:) 
    pdf_params_multi%thl_2(icol,:)                = pdf_params_single%thl_2(1,:) 
    pdf_params_multi%varnce_thl_1(icol,:)         = pdf_params_single%varnce_thl_1(1,:) 
    pdf_params_multi%varnce_thl_2(icol,:)         = pdf_params_single%varnce_thl_2(1,:) 
    pdf_params_multi%corr_w_rt_1(icol,:)          = pdf_params_single%corr_w_rt_1(1,:) 
    pdf_params_multi%corr_w_rt_2(icol,:)          = pdf_params_single%corr_w_rt_2(1,:) 
    pdf_params_multi%corr_w_thl_1(icol,:)         = pdf_params_single%corr_w_thl_1(1,:) 
    pdf_params_multi%corr_w_thl_2(icol,:)         = pdf_params_single%corr_w_thl_2(1,:) 
    pdf_params_multi%corr_rt_thl_1(icol,:)        = pdf_params_single%corr_rt_thl_1(1,:) 
    pdf_params_multi%corr_rt_thl_2(icol,:)        = pdf_params_single%corr_rt_thl_2(1,:) 
    pdf_params_multi%alpha_thl(icol,:)            = pdf_params_single%alpha_thl(1,:) 
    pdf_params_multi%alpha_rt(icol,:)             = pdf_params_single%alpha_rt(1,:) 
    pdf_params_multi%crt_1(icol,:)                = pdf_params_single%crt_1(1,:) 
    pdf_params_multi%crt_2(icol,:)                = pdf_params_single%crt_2(1,:) 
    pdf_params_multi%cthl_1(icol,:)               = pdf_params_single%cthl_1(1,:) 
    pdf_params_multi%cthl_2(icol,:)               = pdf_params_single%cthl_2(1,:) 
    pdf_params_multi%chi_1(icol,:)                = pdf_params_single%chi_1(1,:) 
    pdf_params_multi%chi_2(icol,:)                = pdf_params_single%chi_2(1,:) 
    pdf_params_multi%stdev_chi_1(icol,:)          = pdf_params_single%stdev_chi_1(1,:) 
    pdf_params_multi%stdev_chi_2(icol,:)          = pdf_params_single%stdev_chi_2(1,:) 
    pdf_params_multi%stdev_eta_1(icol,:)          = pdf_params_single%stdev_eta_1(1,:) 
    pdf_params_multi%stdev_eta_2(icol,:)          = pdf_params_single%stdev_eta_2(1,:) 
    pdf_params_multi%covar_chi_eta_1(icol,:)      = pdf_params_single%covar_chi_eta_1(1,:) 
    pdf_params_multi%covar_chi_eta_2(icol,:)      = pdf_params_single%covar_chi_eta_2(1,:) 
    pdf_params_multi%corr_w_chi_1(icol,:)         = pdf_params_single%corr_w_chi_1(1,:) 
    pdf_params_multi%corr_w_chi_2(icol,:)         = pdf_params_single%corr_w_chi_2(1,:) 
    pdf_params_multi%corr_w_eta_1(icol,:)         = pdf_params_single%corr_w_eta_1(1,:) 
    pdf_params_multi%corr_w_eta_2(icol,:)         = pdf_params_single%corr_w_eta_2(1,:) 
    pdf_params_multi%corr_chi_eta_1(icol,:)       = pdf_params_single%corr_chi_eta_1(1,:) 
    pdf_params_multi%corr_chi_eta_2(icol,:)       = pdf_params_single%corr_chi_eta_2(1,:) 
    pdf_params_multi%rsatl_1(icol,:)              = pdf_params_single%rsatl_1(1,:) 
    pdf_params_multi%rsatl_2(icol,:)              = pdf_params_single%rsatl_2(1,:) 
    pdf_params_multi%rc_1(icol,:)                 = pdf_params_single%rc_1(1,:) 
    pdf_params_multi%rc_2(icol,:)                 = pdf_params_single%rc_2(1,:) 
    pdf_params_multi%cloud_frac_1(icol,:)         = pdf_params_single%cloud_frac_1(1,:) 
    pdf_params_multi%cloud_frac_2(icol,:)         = pdf_params_single%cloud_frac_2(1,:) 
    pdf_params_multi%mixt_frac(icol,:)            = pdf_params_single%mixt_frac(1,:) 
    pdf_params_multi%ice_supersat_frac_1(icol,:)  = pdf_params_single%ice_supersat_frac_1(1,:) 
    pdf_params_multi%ice_supersat_frac_2(icol,:)  = pdf_params_single%ice_supersat_frac_2(1,:) 
    
  end subroutine copy_single_pdf_params_to_multi
  
  !================================================================================================
  ! copy_multi_pdf_params_to_single - copies values of a multiple column version of pdf_params
  !   at a specified column to a single column version.
  !
  ! NOTE: THIS SUBROUTINE IS INTENDED TO BE TEMPORARY AND SHOULD BECOME UNNECESSARY ONCE 
  !       CLUBB IS ABLE TO OPERATE OVER MULTIPLE COLUMNS.
  !       See https://github.com/larson-group/cam/issues/129#issuecomment-827944454
  !================================================================================================
  subroutine copy_multi_pdf_params_to_single( pdf_params_multi, icol, &
                                              pdf_params_single )
    
    implicit none

    ! Input Variable(s)
    integer, intent(in) :: &
      icol   ! Column number to copy to
      
    type(pdf_parameter), intent(in) :: &
      pdf_params_multi  ! PDF parameters            [units vary]

    ! Output Variable(s)
    type(pdf_parameter), intent(inout) :: &
      pdf_params_single   ! PDF parameters            [units vary]
      
    pdf_params_single%w_1(1,:)                  = pdf_params_multi%w_1(icol,:) 
    pdf_params_single%w_2(1,:)                  = pdf_params_multi%w_2(icol,:) 
    pdf_params_single%varnce_w_1(1,:)           = pdf_params_multi%varnce_w_1(icol,:) 
    pdf_params_single%varnce_w_2(1,:)           = pdf_params_multi%varnce_w_2(icol,:) 
    pdf_params_single%rt_1(1,:)                 = pdf_params_multi%rt_1(icol,:) 
    pdf_params_single%rt_2(1,:)                 = pdf_params_multi%rt_2(icol,:) 
    pdf_params_single%varnce_rt_1(1,:)          = pdf_params_multi%varnce_rt_1(icol,:) 
    pdf_params_single%varnce_rt_2(1,:)          = pdf_params_multi%varnce_rt_2(icol,:) 
    pdf_params_single%thl_1(1,:)                = pdf_params_multi%thl_1(icol,:) 
    pdf_params_single%thl_2(1,:)                = pdf_params_multi%thl_2(icol,:) 
    pdf_params_single%varnce_thl_1(1,:)         = pdf_params_multi%varnce_thl_1(icol,:) 
    pdf_params_single%varnce_thl_2(1,:)         = pdf_params_multi%varnce_thl_2(icol,:) 
    pdf_params_single%corr_w_rt_1(1,:)          = pdf_params_multi%corr_w_rt_1(icol,:) 
    pdf_params_single%corr_w_rt_2(1,:)          = pdf_params_multi%corr_w_rt_2(icol,:) 
    pdf_params_single%corr_w_thl_1(1,:)         = pdf_params_multi%corr_w_thl_1(icol,:) 
    pdf_params_single%corr_w_thl_2(1,:)         = pdf_params_multi%corr_w_thl_2(icol,:) 
    pdf_params_single%corr_rt_thl_1(1,:)        = pdf_params_multi%corr_rt_thl_1(icol,:) 
    pdf_params_single%corr_rt_thl_2(1,:)        = pdf_params_multi%corr_rt_thl_2(icol,:) 
    pdf_params_single%alpha_thl(1,:)            = pdf_params_multi%alpha_thl(icol,:) 
    pdf_params_single%alpha_rt(1,:)             = pdf_params_multi%alpha_rt(icol,:) 
    pdf_params_single%crt_1(1,:)                = pdf_params_multi%crt_1(icol,:) 
    pdf_params_single%crt_2(1,:)                = pdf_params_multi%crt_2(icol,:) 
    pdf_params_single%cthl_1(1,:)               = pdf_params_multi%cthl_1(icol,:) 
    pdf_params_single%cthl_2(1,:)               = pdf_params_multi%cthl_2(icol,:) 
    pdf_params_single%chi_1(1,:)                = pdf_params_multi%chi_1(icol,:) 
    pdf_params_single%chi_2(1,:)                = pdf_params_multi%chi_2(icol,:) 
    pdf_params_single%stdev_chi_1(1,:)          = pdf_params_multi%stdev_chi_1(icol,:) 
    pdf_params_single%stdev_chi_2(1,:)          = pdf_params_multi%stdev_chi_2(icol,:) 
    pdf_params_single%stdev_eta_1(1,:)          = pdf_params_multi%stdev_eta_1(icol,:) 
    pdf_params_single%stdev_eta_2(1,:)          = pdf_params_multi%stdev_eta_2(icol,:) 
    pdf_params_single%covar_chi_eta_1(1,:)      = pdf_params_multi%covar_chi_eta_1(icol,:) 
    pdf_params_single%covar_chi_eta_2(1,:)      = pdf_params_multi%covar_chi_eta_2(icol,:) 
    pdf_params_single%corr_w_chi_1(1,:)         = pdf_params_multi%corr_w_chi_1(icol,:) 
    pdf_params_single%corr_w_chi_2(1,:)         = pdf_params_multi%corr_w_chi_2(icol,:) 
    pdf_params_single%corr_w_eta_1(1,:)         = pdf_params_multi%corr_w_eta_1(icol,:) 
    pdf_params_single%corr_w_eta_2(1,:)         = pdf_params_multi%corr_w_eta_2(icol,:) 
    pdf_params_single%corr_chi_eta_1(1,:)       = pdf_params_multi%corr_chi_eta_1(icol,:) 
    pdf_params_single%corr_chi_eta_2(1,:)       = pdf_params_multi%corr_chi_eta_2(icol,:) 
    pdf_params_single%rsatl_1(1,:)              = pdf_params_multi%rsatl_1(icol,:) 
    pdf_params_single%rsatl_2(1,:)              = pdf_params_multi%rsatl_2(icol,:) 
    pdf_params_single%rc_1(1,:)                 = pdf_params_multi%rc_1(icol,:) 
    pdf_params_single%rc_2(1,:)                 = pdf_params_multi%rc_2(icol,:) 
    pdf_params_single%cloud_frac_1(1,:)         = pdf_params_multi%cloud_frac_1(icol,:) 
    pdf_params_single%cloud_frac_2(1,:)         = pdf_params_multi%cloud_frac_2(icol,:) 
    pdf_params_single%mixt_frac(1,:)            = pdf_params_multi%mixt_frac(icol,:) 
    pdf_params_single%ice_supersat_frac_1(1,:)  = pdf_params_multi%ice_supersat_frac_1(icol,:) 
    pdf_params_single%ice_supersat_frac_2(1,:)  = pdf_params_multi%ice_supersat_frac_2(icol,:) 
    
  end subroutine copy_multi_pdf_params_to_single

  
  !================================================================================================
  ! init_precip_fracs - allocates arrays for precip_fracs
  !================================================================================================
  subroutine init_precip_fracs_api( nz, ngrdcol, &
                                    precip_fracs )

    use hydromet_pdf_parameter_module, only : init_precip_fracs

    implicit none
    
    ! Input Variable(s)
    integer, intent(in) :: &
      nz,     & ! Number of vertical grid levels    [-]
      ngrdcol   ! Number of grid columns            [-]

    ! Output Variable
    type(precipitation_fractions), intent(out) :: &
      precip_fracs    ! Hydrometeor PDF parameters      [units vary]

    call init_precip_fracs( nz, ngrdcol, & ! intent(in)
                            precip_fracs ) ! intent(out)

    return

  end subroutine init_precip_fracs_api
  
  !================================================================================================
  ! init_pdf_implicit_coefs_terms - allocates arrays for the PDF implicit
  ! coefficient and explicit terms.
  !================================================================================================
  subroutine init_pdf_implicit_coefs_terms_api( nz, ngrdcol, sclr_dim, &
                                                pdf_implicit_coefs_terms )

    use pdf_parameter_module, only: &
        init_pdf_implicit_coefs_terms    ! Procedure(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz,       & ! Number of vertical grid levels    [-]
      ngrdcol,  & ! Number of grid columns            [-]
      sclr_dim    ! Number of scalar variables        [-]

    ! Output Variable
    type(implicit_coefs_terms), intent(out) :: &
      pdf_implicit_coefs_terms    ! Implicit coefs / explicit terms [units vary]

    call init_pdf_implicit_coefs_terms( nz, ngrdcol, sclr_dim, & ! intent(in)
                                        pdf_implicit_coefs_terms ) ! intent(out)

  end subroutine init_pdf_implicit_coefs_terms_api

  !================================================================================================
  ! setup_pdf_parameters
  !================================================================================================

  subroutine setup_pdf_parameters_api_single_col( gr, & ! intent(in)
    nzm, nzt, pdf_dim, hydromet_dim, dt, &      ! Intent(in)
    Nc_in_cloud, cloud_frac, Kh_zm, &           ! Intent(in)
    ice_supersat_frac, hydromet, wphydrometp, & ! Intent(in)
    corr_array_n_cloud, corr_array_n_below, &   ! Intent(in)
    hm_metadata, &                           ! Intent(in)
    pdf_params, &                               ! Intent(in)
    clubb_params, &                             ! Intent(in)
    iiPDF_type, &                               ! Intent(in)
    l_use_precip_frac, &                        ! Intent(in)
    l_predict_upwp_vpwp, &                      ! Intent(in)
    l_diagnose_correlations, &                  ! Intent(in)
    l_calc_w_corr, &                            ! Intent(in)
    l_const_Nc_in_cloud, &                      ! Intent(in)
    l_fix_w_chi_eta_correlations, &             ! Intent(in)
    stats_metadata, &                           ! Intent(in)
    stats_zt, stats_zm, stats_sfc, &            ! intent(inout)
    hydrometp2, &                               ! Intent(inout)
    mu_x_1_n, mu_x_2_n, &                       ! Intent(out)
    sigma_x_1_n, sigma_x_2_n, &                 ! Intent(out)
    corr_array_1_n, corr_array_2_n, &           ! Intent(out)
    corr_cholesky_mtx_1, corr_cholesky_mtx_2, & ! Intent(out)
    precip_fracs,  &                            ! Intent(out)
    hydromet_pdf_params )                       ! Intent(out)

    use setup_clubb_pdf_params, only : &
        setup_pdf_parameters

    use parameter_indices, only: &
        nparams    ! Variable(s)

    use error_code, only : &
        err_code, &         ! Error Indicator
        clubb_fatal_error   ! Constant

    implicit none

    !----------------------- Input Variables -----------------------
    type(grid), intent(in) :: gr

    integer, intent(in) :: &
      nzm,          & ! Number of model momentum vertical grid levels
      nzt,          & ! Number of model thermodynamic vertical grid levels
      pdf_dim,      & ! Number of variables in the correlation array
      hydromet_dim    ! Number of hydrometeor species

    real( kind = core_rknd ), intent(in) ::  &
      dt    ! Model timestep                                           [s]

    real( kind = core_rknd ), dimension(nzt), intent(in) :: &
      Nc_in_cloud,       & ! Mean (in-cloud) cloud droplet conc.       [num/kg]
      cloud_frac,        & ! Cloud fraction                            [-]
      ice_supersat_frac    ! Ice supersaturation fraction              [-]

    real( kind = core_rknd ), dimension(nzm), intent(in) :: &
      Kh_zm                ! Eddy diffusivity coef. on momentum levels [m^2/s]

    real( kind = core_rknd ), dimension(nzt,hydromet_dim), intent(in) :: &
      hydromet       ! Mean of hydrometeor, hm (overall) (t-levs.) [units]

    real( kind = core_rknd ), dimension(nzm,hydromet_dim), intent(in) :: &
      wphydrometp    ! Covariance < w'h_m' > (momentum levels)     [(m/s)units]

    real( kind = core_rknd ), dimension(pdf_dim,pdf_dim), &
      intent(in) :: &
      corr_array_n_cloud, & ! Prescribed norm. space corr. array in cloud    [-]
      corr_array_n_below    ! Prescribed norm. space corr. array below cloud [-]

    type (hm_metadata_type), intent(in) :: &
      hm_metadata

    type(pdf_parameter), intent(in) :: &
      pdf_params    ! PDF parameters                               [units vary]

    real( kind = core_rknd ), dimension(nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

    integer, intent(in) :: &
      iiPDF_type    ! Selected option for the two-component normal (double
                    ! Gaussian) PDF type to use for the w, rt, and theta-l (or
                    ! w, chi, and eta) portion of CLUBB's multivariate,
                    ! two-component PDF.

    logical, intent(in) :: &
      l_use_precip_frac,            & ! Flag to use precipitation fraction in KK microphysics. The
                                      ! precipitation fraction is automatically set to 1 when this
                                      ! flag is turned off.
      l_predict_upwp_vpwp,          & ! Flag to predict <u'w'> and <v'w'> along with <u> and <v>
                                      ! alongside the advancement of <rt>, <w'rt'>, <thl>,
                                      ! <wpthlp>, <sclr>, and <w'sclr'> in subroutine
                                      ! advance_xm_wpxp.  Otherwise, <u'w'> and <v'w'> are still
                                      ! approximated by eddy diffusivity when <u> and <v> are
                                      ! advanced in subroutine advance_windm_edsclrm.
      l_diagnose_correlations,      & ! Diagnose correlations instead of using fixed ones
      l_calc_w_corr,                & ! Calculate the correlations between w and the hydrometeors
      l_const_Nc_in_cloud,          & ! Use a constant cloud droplet conc. within cloud (K&K)
      l_fix_w_chi_eta_correlations    ! Use a fixed correlation for s and t Mellor(chi/eta)

    type (stats_metadata_type), intent(in) :: &
      stats_metadata

    !----------------------- Input/Output Variables -----------------------
    type(stats), intent(inout) :: &
      stats_zt, &
      stats_zm, &
      stats_sfc

    real( kind = core_rknd ), dimension(nzm,hydromet_dim), intent(inout) :: &
      hydrometp2    ! Variance of a hydrometeor (overall) (m-levs.)   [units^2]

    !----------------------- Output Variables -----------------------
    real( kind = core_rknd ), dimension(nzt,pdf_dim), intent(out) :: &
      mu_x_1_n,    & ! Mean array (normal space): PDF vars. (comp. 1) [un. vary]
      mu_x_2_n,    & ! Mean array (normal space): PDF vars. (comp. 2) [un. vary]
      sigma_x_1_n, & ! Std. dev. array (normal space): PDF vars (comp. 1) [u.v.]
      sigma_x_2_n    ! Std. dev. array (normal space): PDF vars (comp. 2) [u.v.]

    real( kind = core_rknd ), dimension(nzt,pdf_dim,pdf_dim), &
      intent(out) :: &
      corr_array_1_n, & ! Corr. array (normal space):  PDF vars. (comp. 1)   [-]
      corr_array_2_n    ! Corr. array (normal space):  PDF vars. (comp. 2)   [-]

    real( kind = core_rknd ), dimension(nzt,pdf_dim,pdf_dim), &
      intent(out) :: &
      corr_cholesky_mtx_1, & ! Transposed corr. cholesky matrix, 1st comp. [-]
      corr_cholesky_mtx_2    ! Transposed corr. cholesky matrix, 2nd comp. [-]

    ! This is only an output, but it contains allocated arrays, so we need to treat it as inout
    type(precipitation_fractions), intent(inout) :: &
      precip_fracs           ! Precipitation fractions      [-]

    type(hydromet_pdf_parameter), dimension(nzt), optional, intent(out) :: &
      hydromet_pdf_params    ! Hydrometeor PDF parameters        [units vary]
      
    !----------------------- Local Variables -----------------------
    real( kind = core_rknd ), dimension(1,nzt) :: &
      Nc_in_cloud_col,       & ! Mean (in-cloud) cloud droplet conc.       [num/kg]
      cloud_frac_col,        & ! Cloud fraction                            [-]
      ice_supersat_frac_col    ! Ice supersaturation fraction              [-]

    real( kind = core_rknd ), dimension(1,nzm) :: &
      Kh_zm_col                ! Eddy diffusivity coef. on momentum levels [m^2/s]

    real( kind = core_rknd ), dimension(1,nzt,hydromet_dim) :: &
      hydromet_col       ! Mean of hydrometeor, hm (overall) (t-levs.) [units]

    real( kind = core_rknd ), dimension(1,nzm,hydromet_dim) :: &
      wphydrometp_col    ! Covariance < w'h_m' > (momentum levels)     [(m/s)units]
      
    ! Input/Output Variables
    real( kind = core_rknd ), dimension(1,nzm,hydromet_dim) :: &
      hydrometp2_col    ! Variance of a hydrometeor (overall) (m-levs.)   [units^2]

    ! Output Variables
    real( kind = core_rknd ), dimension(1,nzt,pdf_dim) :: &
      mu_x_1_n_col,    & ! Mean array (normal space): PDF vars. (comp. 1) [un. vary]
      mu_x_2_n_col,    & ! Mean array (normal space): PDF vars. (comp. 2) [un. vary]
      sigma_x_1_n_col, & ! Std. dev. array (normal space): PDF vars (comp. 1) [u.v.]
      sigma_x_2_n_col    ! Std. dev. array (normal space): PDF vars (comp. 2) [u.v.]

    real( kind = core_rknd ), dimension(1,nzt,pdf_dim,pdf_dim) :: &
      corr_array_1_n_col, & ! Corr. array (normal space):  PDF vars. (comp. 1)   [-]
      corr_array_2_n_col    ! Corr. array (normal space):  PDF vars. (comp. 2)   [-]

    real( kind = core_rknd ), dimension(1,nzt,pdf_dim,pdf_dim) :: &
      corr_cholesky_mtx_1_col, & ! Transposed corr. cholesky matrix, 1st comp. [-]
      corr_cholesky_mtx_2_col    ! Transposed corr. cholesky matrix, 2nd comp. [-]

    type(hydromet_pdf_parameter), dimension(1,nzt) :: &
      hydromet_pdf_params_col    ! Hydrometeor PDF parameters        [units vary]
      
    type(stats), dimension(1) :: &
      stats_zt_col, &
      stats_zm_col, &
      stats_sfc_col

    !----------------------- Begin Code -----------------------

    Nc_in_cloud_col(1,:) = Nc_in_cloud
    cloud_frac_col(1,:) = cloud_frac
    Kh_zm_col(1,:) = Kh_zm
    ice_supersat_frac_col(1,:) = ice_supersat_frac
    
    hydromet_col(1,:,:) = hydromet
    wphydrometp_col(1,:,:) = wphydrometp
    
    stats_zt_col(1) = stats_zt
    stats_zm_col(1) = stats_zm
    stats_sfc_col(1) = stats_sfc

    call setup_pdf_parameters( gr, &                          ! intent(in)
      nzm, nzt, 1, pdf_dim, hydromet_dim, dt, &               ! Intent(in)
      Nc_in_cloud_col, cloud_frac_col, Kh_zm_col, &           ! Intent(in)
      ice_supersat_frac_col, hydromet_col, wphydrometp_col, & ! Intent(in)
      corr_array_n_cloud, corr_array_n_below, &               ! Intent(in)
      hm_metadata, &                                          ! Intent(in)
      pdf_params, &                                           ! Intent(in)
      clubb_params, &                                         ! Intent(in)
      iiPDF_type, &                                           ! Intent(in)
      l_use_precip_frac, &                                    ! Intent(in)
      l_predict_upwp_vpwp, &                                  ! Intent(in)
      l_diagnose_correlations, &                              ! Intent(in)
      l_calc_w_corr, &                                        ! Intent(in)
      l_const_Nc_in_cloud, &                                  ! Intent(in)
      l_fix_w_chi_eta_correlations, &                         ! Intent(in)
      stats_metadata, &                                       ! Intent(in)
      stats_zt_col, stats_zm_col, stats_sfc_col, &            ! Intent(inout)
      hydrometp2_col, &                                       ! Intent(inout)
      mu_x_1_n_col, mu_x_2_n_col, &                           ! Intent(out)
      sigma_x_1_n_col, sigma_x_2_n_col, &                     ! Intent(out)
      corr_array_1_n_col, corr_array_2_n_col, &               ! Intent(out)
      corr_cholesky_mtx_1_col, corr_cholesky_mtx_2_col, &     ! Intent(out)
      precip_fracs, &                                         ! Intent(inout)
      hydromet_pdf_params_col )                               ! Optional(out)

    if ( err_code == clubb_fatal_error ) error stop
    
    ! The following does not work for stats 
    !     stats_zt = stats_zt_col(1)
    !     stats_zm = stats_zm_col(1) 
    !     stats_sfc = stats_sfc_col(1)
    ! because of some mysterious pointer issue. However, the only thing that 
    ! updates in stats is the field values, so we can copy only those instead.
    if ( stats_metadata%l_stats ) then 
      stats_zm%accum_field_values = stats_zm_col(1)%accum_field_values
      stats_zm%accum_num_samples = stats_zm_col(1)%accum_num_samples
      
      stats_zt%accum_field_values = stats_zt_col(1)%accum_field_values
      stats_zt%accum_num_samples = stats_zt_col(1)%accum_num_samples
      
      stats_sfc%accum_field_values = stats_sfc_col(1)%accum_field_values
      stats_sfc%accum_num_samples = stats_sfc_col(1)%accum_num_samples
    end if
    
    hydrometp2 = hydrometp2_col(1,:,:)
    mu_x_1_n = mu_x_1_n_col(1,:,:)
    mu_x_2_n = mu_x_2_n_col(1,:,:)
    sigma_x_1_n = sigma_x_1_n_col(1,:,:)
    sigma_x_2_n = sigma_x_2_n_col(1,:,:)
    corr_array_1_n = corr_array_1_n_col(1,:,:,:)
    corr_array_2_n = corr_array_2_n_col(1,:,:,:)
    corr_cholesky_mtx_1 = corr_cholesky_mtx_1_col(1,:,:,:)
    corr_cholesky_mtx_2 = corr_cholesky_mtx_2_col(1,:,:,:)
    if ( present(hydromet_pdf_params) ) then
      hydromet_pdf_params = hydromet_pdf_params_col(1,:)
    end if

  end subroutine setup_pdf_parameters_api_single_col
!===========================================================================! 
  subroutine setup_pdf_parameters_api_multi_col( gr, &
    nzm, nzt, ngrdcol, pdf_dim, hydromet_dim, dt, & ! Intent(in)
    Nc_in_cloud, cloud_frac, Kh_zm, &           ! Intent(in)
    ice_supersat_frac, hydromet, wphydrometp, & ! Intent(in)
    corr_array_n_cloud, corr_array_n_below, &   ! Intent(in)
    hm_metadata, &                              ! Intent(in)
    pdf_params, &                               ! Intent(in)
    clubb_params, &                             ! Intent(in)
    iiPDF_type, &                               ! Intent(in)
    l_use_precip_frac, &                        ! Intent(in)
    l_predict_upwp_vpwp, &                      ! Intent(in)
    l_diagnose_correlations, &                  ! Intent(in)
    l_calc_w_corr, &                            ! Intent(in)
    l_const_Nc_in_cloud, &                      ! Intent(in)
    l_fix_w_chi_eta_correlations, &             ! Intent(in)
    stats_metadata, &                           ! Intent(in)
    stats_zt, stats_zm, stats_sfc, &            ! Intent(inout)
    hydrometp2, &                               ! Intent(inout)
    mu_x_1_n, mu_x_2_n, &                       ! Intent(out)
    sigma_x_1_n, sigma_x_2_n, &                 ! Intent(out)
    corr_array_1_n, corr_array_2_n, &           ! Intent(out)
    corr_cholesky_mtx_1, corr_cholesky_mtx_2, & ! Intent(out)
    precip_fracs, &                             ! Intent(out)
    hydromet_pdf_params )                       ! Intent(out)

    use setup_clubb_pdf_params, only : setup_pdf_parameters

    use parameter_indices, only: &
        nparams    ! Variable(s)

    use error_code, only : &
        err_code, &         ! Error Indicator
        clubb_fatal_error   ! Constant

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nzm,          & ! Number of model momentum vertical grid levels
      nzt,          & ! Number of model thermodynamic vertical grid levels
      pdf_dim,      & ! Number of variables in the correlation array
      ngrdcol,      & ! Number of grid columns
      hydromet_dim    ! Number of hydrometeor species
      
    type (grid), intent(in) :: gr

    real( kind = core_rknd ), intent(in) ::  &
      dt    ! Model timestep                                           [s]

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      Nc_in_cloud,       & ! Mean (in-cloud) cloud droplet conc.       [num/kg]
      cloud_frac,        & ! Cloud fraction                            [-]
      ice_supersat_frac    ! Ice supersaturation fraction              [-]

    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: &
      Kh_zm                ! Eddy diffusivity coef. on momentum levels [m^2/s]

    real( kind = core_rknd ), dimension(ngrdcol,nzt,hydromet_dim), intent(in) :: &
      hydromet       ! Mean of hydrometeor, hm (overall) (t-levs.) [units]

    real( kind = core_rknd ), dimension(ngrdcol,nzm,hydromet_dim), intent(in) :: &
      wphydrometp    ! Covariance < w'h_m' > (momentum levels)     [(m/s)units]

    real( kind = core_rknd ), dimension(pdf_dim,pdf_dim), &
      intent(in) :: &
      corr_array_n_cloud, & ! Prescribed norm. space corr. array in cloud    [-]
      corr_array_n_below    ! Prescribed norm. space corr. array below cloud [-]

    type (hm_metadata_type), intent(in) :: &
      hm_metadata

    type(pdf_parameter), intent(in) :: &
      pdf_params    ! PDF parameters                               [units vary]

    real( kind = core_rknd ), dimension(nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

    integer, intent(in) :: &
      iiPDF_type    ! Selected option for the two-component normal (double
                    ! Gaussian) PDF type to use for the w, rt, and theta-l (or
                    ! w, chi, and eta) portion of CLUBB's multivariate,
                    ! two-component PDF.

    logical, intent(in) :: &
      l_use_precip_frac,            & ! Flag to use precipitation fraction in KK microphysics. The
                                      ! precipitation fraction is automatically set to 1 when this
                                      ! flag is turned off.
      l_predict_upwp_vpwp,          & ! Flag to predict <u'w'> and <v'w'> along with <u> and <v>
                                      ! alongside the advancement of <rt>, <w'rt'>, <thl>,
                                      ! <wpthlp>, <sclr>, and <w'sclr'> in subroutine
                                      ! advance_xm_wpxp.  Otherwise, <u'w'> and <v'w'> are still
                                      ! approximated by eddy diffusivity when <u> and <v> are
                                      ! advanced in subroutine advance_windm_edsclrm.
      l_diagnose_correlations,      & ! Diagnose correlations instead of using fixed ones
      l_calc_w_corr,                & ! Calculate the correlations between w and the hydrometeors
      l_const_Nc_in_cloud,          & ! Use a constant cloud droplet conc. within cloud (K&K)
      l_fix_w_chi_eta_correlations    ! Use a fixed correlation for s and t Mellor(chi/eta)

    type (stats_metadata_type), intent(in) :: &
      stats_metadata

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(ngrdcol,nzm,hydromet_dim), intent(inout) :: &
      hydrometp2    ! Variance of a hydrometeor (overall) (m-levs.)   [units^2]

    type(stats), dimension(ngrdcol), intent(inout) :: &
      stats_zt, &
      stats_zm, &
      stats_sfc

    ! Output Variables
    real( kind = core_rknd ), dimension(ngrdcol,nzt,pdf_dim), intent(out) :: &
      mu_x_1_n,    & ! Mean array (normal space): PDF vars. (comp. 1) [un. vary]
      mu_x_2_n,    & ! Mean array (normal space): PDF vars. (comp. 2) [un. vary]
      sigma_x_1_n, & ! Std. dev. array (normal space): PDF vars (comp. 1) [u.v.]
      sigma_x_2_n    ! Std. dev. array (normal space): PDF vars (comp. 2) [u.v.]

    real( kind = core_rknd ), dimension(ngrdcol,nzt,pdf_dim,pdf_dim), &
      intent(out) :: &
      corr_array_1_n, & ! Corr. array (normal space):  PDF vars. (comp. 1)   [-]
      corr_array_2_n    ! Corr. array (normal space):  PDF vars. (comp. 2)   [-]

    real( kind = core_rknd ), dimension(ngrdcol,nzt,pdf_dim,pdf_dim), &
      intent(out) :: &
      corr_cholesky_mtx_1, & ! Transposed corr. cholesky matrix, 1st comp. [-]
      corr_cholesky_mtx_2    ! Transposed corr. cholesky matrix, 2nd comp. [-]

    type(hydromet_pdf_parameter), dimension(ngrdcol,nzt), optional, intent(out) :: &
      hydromet_pdf_params    ! Hydrometeor PDF parameters        [units vary]
      
    ! This is only an output, but it contains allocated arrays, so we need to treat it as inout
    type(precipitation_fractions), intent(inout) :: &
      precip_fracs           ! Precipitation fractions      [-]

    call setup_pdf_parameters( gr, &                  ! intent(in)
      nzm, nzt, ngrdcol, pdf_dim, hydromet_dim, dt, & ! Intent(in)
      Nc_in_cloud, cloud_frac, Kh_zm, &               ! Intent(in)
      ice_supersat_frac, hydromet, wphydrometp, &     ! Intent(in)
      corr_array_n_cloud, corr_array_n_below, &       ! Intent(in)
      hm_metadata, &                                  ! Intent(in)
      pdf_params, &                                   ! Intent(in)
      clubb_params, &                                 ! Intent(in)
      iiPDF_type, &                                   ! Intent(in)
      l_use_precip_frac, &                            ! Intent(in)
      l_predict_upwp_vpwp, &                          ! Intent(in)
      l_diagnose_correlations, &                      ! Intent(in)
      l_calc_w_corr, &                                ! Intent(in)
      l_const_Nc_in_cloud, &                          ! Intent(in)
      l_fix_w_chi_eta_correlations, &                 ! Intent(in)
      stats_metadata, &                               ! Intent(in)
      stats_zt, stats_zm, stats_sfc, &                ! intent(inout)
      hydrometp2, &                                   ! Intent(inout)
      mu_x_1_n, mu_x_2_n, &                           ! Intent(out)
      sigma_x_1_n, sigma_x_2_n, &                     ! Intent(out)
      corr_array_1_n, corr_array_2_n, &               ! Intent(out)
      corr_cholesky_mtx_1, corr_cholesky_mtx_2, &     ! Intent(out)
      precip_fracs, &                                 ! Intent(inout)
      hydromet_pdf_params )                           ! Optional(out)

    if ( err_code == clubb_fatal_error ) error stop

  end subroutine setup_pdf_parameters_api_multi_col

  !================================================================================================
  ! stats_end_timestep - Calls statistics to be written to the output format.
  !================================================================================================

  subroutine stats_end_timestep_api( clubb_params, stats_metadata, & 
                                     stats_zt, stats_zm, stats_sfc, &
                                     stats_lh_zt, stats_lh_sfc, &
                                     stats_rad_zt, stats_rad_zm &
#ifdef NETCDF
                                     , l_uv_nudge, & ! Unused
                                     l_tke_aniso, & ! Unused
                                     l_standard_term_ta & ! Unused
#endif
                                      )

    use stats_clubb_utilities, only : stats_end_timestep

    use parameter_indices, only: nparams

    implicit none

    ! Input Variables
    real( kind = core_rknd ), dimension(nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

    type (stats_metadata_type), intent(in) :: &
      stats_metadata

    type(stats), intent(inout) :: &
      stats_zt, &
      stats_zm, &
      stats_sfc, &
      stats_lh_zt, &
      stats_lh_sfc, &
      stats_rad_zt, &
      stats_rad_zm

#ifdef NETCDF
    logical, intent(in) :: &
      l_uv_nudge,         & ! For wind speed nudging.
      l_tke_aniso,        & ! For anisotropic turbulent kinetic energy, i.e.
                            ! TKE = 1/2 (u'^2 + v'^2 + w'^2)
      l_standard_term_ta    ! Use the standard discretization for the turbulent advection terms.
                            ! Setting to .false. means that a_1 and a_3 are pulled outside of the
                            ! derivative in advance_wp2_wp3_module.F90 and in
                            ! advance_xp2_xpyp_module.F90.
#endif

    call stats_end_timestep( stats_metadata,      & ! intent(in)
                             stats_zt, stats_zm, stats_sfc,     & ! intent(inout)
                             stats_lh_zt, stats_lh_sfc,         & ! intent(inout)
                             stats_rad_zt, stats_rad_zm         & ! intent(inout)
                           )

    if ( err_code == clubb_fatal_error ) error stop

  end subroutine stats_end_timestep_api

  !================================================================================================
  ! stats_finalize - Close NetCDF files and deallocate scratch space and stats file structures.
  !================================================================================================

  subroutine stats_finalize_api ( stats_metadata, &
                                  stats_zt, stats_zm, stats_sfc, &
                                  stats_lh_zt, stats_lh_sfc, &
                                  stats_rad_zt, stats_rad_zm )

    use stats_clubb_utilities, only : &
      stats_finalize

    implicit none

    type(stats), intent(inout) :: &
      stats_zt, &
      stats_zm, &
      stats_sfc, &
      stats_lh_zt, &
      stats_lh_sfc, &
      stats_rad_zt, &
      stats_rad_zm

    type (stats_metadata_type), intent(inout) :: &
      stats_metadata

    call stats_finalize ( stats_metadata, & ! intent(inout)
                          stats_zt, stats_zm, stats_sfc, & ! intent(inout)
                          stats_lh_zt, stats_lh_sfc, & ! intent(inout)
                          stats_rad_zt, stats_rad_zm ) ! intent(inout)

  end subroutine stats_finalize_api


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
  ! calculate_thlp2_rad - Computes the contribution of radiative cooling to thlp2
  !================================================================================================
  subroutine calculate_thlp2_rad_api &
                  ( nz, rcm_zm, thlprcp, radht_zm, & ! Intent(in)
                    clubb_params,                  & ! Intent(in)
                    thlp2_forcing )                  ! Intent(inout)

    use clubb_precision, only: &
        core_rknd                     ! Constant(s)

    use advance_clubb_core_module, only: &
        calculate_thlp2_rad

    use parameter_indices, only: &
        nparams

    implicit none

  ! Input Variables
    integer, intent(in) :: &
      nz                    ! Number of vertical levels                      [-]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      rcm_zm, &             ! Cloud water mixing ratio on momentum grid      [kg/kg]
      thlprcp, &            ! thl'rc'                                        [K kg/kg]
      radht_zm              ! SW + LW heating rate (on momentum grid)        [K/s]

    real( kind = core_rknd ), dimension(nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

  ! Input/Output Variables
    real( kind = core_rknd ), dimension(nz), intent(inout) :: &
      thlp2_forcing         ! <th_l'^2> forcing (momentum levels)            [K^2/s]
  !----------------------------------------------------------------------

    call calculate_thlp2_rad( nz, rcm_zm, thlprcp, radht_zm, & ! intent(in)
                              clubb_params,                  & ! intent(in)
                              thlp2_forcing )                  ! intent(inout)

    return
  end subroutine calculate_thlp2_rad_api

  !================================================================================================
  ! update_xp2_mc - Calculates the effects of rain evaporation on rtp2 and thlp2
  !================================================================================================
  subroutine update_xp2_mc_api_multi_col( gr, nzm, nzt, ngrdcol, dt, cloud_frac, rcm, rvm, thlm, &
                                          wm, exner, rrm_evap, pdf_params,        &
                                          rtp2_mc, thlp2_mc, wprtp_mc, wpthlp_mc,    &
                                          rtpthlp_mc )

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
  end subroutine update_xp2_mc_api_multi_col
  
  !================================================================================================
  ! update_xp2_mc - Calculates the effects of rain evaporation on rtp2 and thlp2
  !================================================================================================
  subroutine update_xp2_mc_api_single_col( gr, nzm, nzt, dt, cloud_frac, rcm, rvm, thlm, &
                                           wm, exner, rrm_evap, pdf_params,        &
                                           rtp2_mc, thlp2_mc, wprtp_mc, wpthlp_mc,    &
                                           rtpthlp_mc )

    use advance_xp2_xpyp_module, only: &
        update_xp2_mc

     ! Type

    implicit none

    ! -------------------- Input Variables --------------------
    integer, intent(in) :: &
      nzm, & ! Momentum points in the Vertical        [-]
      nzt    ! Thermodynamic points in the Vertical

    type(grid), intent(in) :: gr

    real( kind = core_rknd ), intent(in) :: dt ! Model timestep        [s]

    real( kind = core_rknd ), dimension(nzt), intent(in) :: &
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
    real( kind = core_rknd ), dimension(nzm), intent(inout) :: &
      rtp2_mc, &    !Tendency of <rt'^2> due to evaporation   [(kg/kg)^2/s]
      thlp2_mc, &   !Tendency of <thl'^2> due to evaporation  [K^2/s]
      wprtp_mc, &   !Tendency of <w'rt'> due to evaporation   [m*(kg/kg)/s^2]
      wpthlp_mc, &  !Tendency of <w'thl'> due to evaporation  [m*K/s^2] 
      rtpthlp_mc    !Tendency of <rt'thl'> due to evaporation [K*(kg/kg)/s]
      
    ! -------------------- Local Variables --------------------
    real( kind = core_rknd ), dimension(1,nzt) :: &
      cloud_frac_col, &       !Cloud fraction                        [-]
      rcm_col, &              !Cloud water mixing ratio              [kg/kg]
      rvm_col, &              !Vapor water mixing ratio              [kg/kg]
      thlm_col, &             !Liquid potential temperature          [K]
      wm_col, &               !Mean vertical velocity                [m/s]
      exner_col, &            !Exner function                        [-]
      rrm_evap_col         !Evaporation of rain                   [kg/kg/s]
                          !It is expected that this variable is negative, as
                          !that is the convention in Morrison microphysics
                          
    real( kind = core_rknd ), dimension(1,nzm) :: &
      rtp2_mc_col, &    !Tendency of <rt'^2> due to evaporation   [(kg/kg)^2/s]
      thlp2_mc_col, &   !Tendency of <thl'^2> due to evaporation  [K^2/s]
      wprtp_mc_col, &   !Tendency of <w'rt'> due to evaporation   [m*(kg/kg)/s^2]
      wpthlp_mc_col, &  !Tendency of <w'thl'> due to evaporation  [m*K/s^2] 
      rtpthlp_mc_col    !Tendency of <rt'thl'> due to evaporation [K*(kg/kg)/s]
      
    cloud_frac_col(1,:) = cloud_frac
    rcm_col(1,:) = rcm
    rvm_col(1,:) = rvm
    thlm_col(1,:) = thlm
    wm_col(1,:) = wm
    exner_col(1,:) = exner
    rrm_evap_col(1,:) = rrm_evap
    rtp2_mc_col(1,:) = rtp2_mc
    thlp2_mc_col(1,:) = thlp2_mc
    wprtp_mc_col(1,:) = wprtp_mc
    wpthlp_mc_col(1,:) = wpthlp_mc
    rtpthlp_mc_col(1,:) = rtpthlp_mc
      
    call update_xp2_mc( gr, nzm, nzt, 1, dt, cloud_frac_col, rcm_col, rvm_col, thlm_col, & ! In
                        wm_col, exner_col, rrm_evap_col, pdf_params,        & ! In
                        rtp2_mc_col, thlp2_mc_col, wprtp_mc_col, wpthlp_mc_col,    & ! In/Out
                        rtpthlp_mc_col ) ! In/Out
                        
    rtp2_mc = rtp2_mc_col(1,:)
    thlp2_mc = thlp2_mc_col(1,:)
    wprtp_mc = wprtp_mc_col(1,:)
    wpthlp_mc = wpthlp_mc_col(1,:)
    rtpthlp_mc = rtpthlp_mc_col(1,:)
    
    return
  end subroutine update_xp2_mc_api_single_col
    
  

  !================================================================================================
  ! set_default_parameters: Sets all CLUBB tunable parameters to a default setting
  !================================================================================================
  subroutine set_default_parameters_api( &
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
               wpxp_Ri_exp, a3_coef_min, a_const, bv_efold, z_displace )

    use parameters_tunable, only: &
        set_default_parameters    ! Procedure(s)

    implicit none

    ! Output variables
    real( kind = core_rknd ), intent(out) :: &
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
      wpxp_Ri_exp, a3_coef_min, a_const, bv_efold, z_displace

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
               wpxp_Ri_exp, a3_coef_min, a_const, bv_efold, z_displace )

    return

  end subroutine set_default_parameters_api

  !================================================================================================
  ! print_clubb_config_flags: Prints the clubb_config_flags
  !================================================================================================
  subroutine print_clubb_config_flags_api( iunit, clubb_config_flags ) ! In

    use model_flags, only: &
        clubb_config_flags_type, &          ! Type
        print_clubb_config_flags            ! Procedure

    implicit none

    ! Input variables
    integer, intent(in) :: &
      iunit ! The file to write to

    type(clubb_config_flags_type), intent(in) :: &
      clubb_config_flags ! Derived type holding all configurable CLUBB flags

    call print_clubb_config_flags( iunit, clubb_config_flags ) ! In

  end subroutine print_clubb_config_flags_api

  !================================================================================================
  ! initialize_tau_sponge_damp
  !================================================================================================
  subroutine initialize_tau_sponge_damp_api( gr, nz, dt, z, settings, damping_profile )

    use sponge_layer_damping, only: &
        sponge_damp_settings,       & ! Variable(s)
        sponge_damp_profile,        &
        initialize_tau_sponge_damp    ! Procedure(s)

    

    implicit none

    type(grid), intent(in) :: gr

    ! Input Variable(s)
    integer, intent(in) :: &
      nz

    real( kind = core_rknd ), intent(in) :: &
      dt    ! Model Timestep    [s]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      z    ! Height of model grid levels    [m]

    type(sponge_damp_settings), intent(in) :: &
      settings

    ! Output Variable(s)
    type(sponge_damp_profile), intent(out) :: &
      damping_profile

    call initialize_tau_sponge_damp( gr, nz, dt, z, settings, & ! intent(in)
                                     damping_profile ) ! intent(inout)

  end subroutine initialize_tau_sponge_damp_api

  !================================================================================================
  ! finalize_tau_sponge_damp
  !================================================================================================
  subroutine finalize_tau_sponge_damp_api( damping_profile )

    use sponge_layer_damping, only: &
        sponge_damp_profile,      & ! Variable(s)
        finalize_tau_sponge_damp    ! Procedure(s)

    implicit none

    ! Input/Output Variable(s)
    type(sponge_damp_profile), intent(inout) :: &
      damping_profile ! Information for damping the profile

    call finalize_tau_sponge_damp( damping_profile ) ! intent(inout)

  end subroutine finalize_tau_sponge_damp_api
    
end module clubb_api_module
