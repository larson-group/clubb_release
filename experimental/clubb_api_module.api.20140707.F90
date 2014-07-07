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

  use grid_class, only : &
    zt2zm_api => zt2zm, & ! The interface implementation of these subroutines
    zm2zt_api => zm2zt    ! requires a use statement "interface" here.

  use array_index, only : &
    iiNgraupelm, & ! Hydrometeor array index for graupel concentration, Ng
    iiNim, & ! Hydrometeor array index for ice concentration, Ni
    iiNrm, & ! Hydrometeor array index for rain drop concentration, Nr
    iiNsnowm, & ! Hydrometeor array index for snow concentration, Ns
    iirgraupelm, & ! Hydrometeor array index for graupel mixing ratio, rg
    iiricem, & ! Hydrometeor array index for ice mixing ratio, ri
    iirrainm, & ! Hydrometeor array index for rain water mixing ratio, rr
    iirsnowm, & ! Hydrometeor array index for snow mixing ratio, rs
    l_frozen_hm, & ! if true, then the hydrometeor is frozen; otherwise liquid
    l_mix_rat_hm, & ! if true, then the quantity is a hydrometeor mixing ratio
    ! Scalars
    iiedsclr_CO2, &
    iiedsclr_rt, &
    iiedsclr_thl, &
    iisclr_CO2, &
    iisclr_rt, &
    iisclr_thl


  use clubb_precision, only : time_precision, &
    ! The precisions below are arbitrary, and could be adjusted as
    ! needed for long simulations or time averaging.  Note that on
    ! most machines 12 digits of precision will use a data type
    ! which is 8 bytes long.
    core_rknd, &
    stat_nknd, &
    stat_rknd, &
    ! This definition of double precision must use a real type that is 64 bits
    ! wide, because (at least) the LAPACK routines depend on this definition being
    ! accurate. Otherwise, LAPACK must be recompiled, or some other trickery must
    ! be done.
    dp

  use constants_clubb, only : &
    cm3_per_m3, & ! Cubic centimeters per cubic meter
    Cp, & ! Dry air specific heat at constant p [J/kg/K]
    em_min, & ! Minimum value for em (turbulence kinetic energy)
    ep, & ! ep  = 0.622  [-]
    fstderr, & ! Fortran file unit I/O constant
    fstdout, & ! Fortran file unit I/O constant
    grav, & ! Gravitational acceleration     [m/s^2]
    Ls, & ! Latent heat of sublimation          [J/kg]
    Lv, & ! Latent heat of vaporization         [J/kg]
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
    zero_threshold ! Defining a threshold on a physical quantity to be 0.
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

  use corr_matrix_module, only : &
    corr_array_cloud, &
    corr_array_below, &
    d_variables, &
    iiPDF_chi, &
    iiPDF_rrain, &
    iiPDF_w, &
    iiPDF_Nr, &
    iiPDF_rice, &
    iiPDF_Ni, &
    iiPDF_Ncn, &
    iiPDF_rsnow, &
    iiPDF_Nsnow, &
    iiPDF_rgraupel, &
    iiPDF_Ngraupel, &
    sigma2_on_mu2_ip_array_cloud, &
    sigma2_on_mu2_ip_array_below

  use error_code, only : &
    clubb_no_error

  use grid_class, only : &
    gr

  use hydromet_pdf_parameter_module, only : &
    hydromet_pdf_parameter

  use model_flags, only : &
    l_use_boussinesq, &
    l_diagnose_correlations, &
    l_calc_w_corr, &
    l_use_cloud_cover, &
    l_use_precip_frac, &
    l_tke_aniso

  use parameters_microphys, only : &
    l_lh_vert_overlap, &
    LH_microphys_type, &
    hydromet_list, &
    LH_microphys_disabled, &
    LH_microphys_interactive, &
    LH_microphys_non_interactive, &
    LH_microphys_calls, &
    LH_sequence_length, &
    LH_seed, &
    l_local_kk, &
    l_fix_s_t_correlations, &
    l_lh_cloud_weighted_sampling, &
    LH_sequence_length, &
    hydromet_tol, &
    Nc0_in_cloud, &
    l_const_Nc_in_cloud, &
    l_silhs_KK_convergence_adj_mean, &
    l_predict_Nc, &
    rr_sigma2_on_mu2_ip_cloud, &
    rr_sigma2_on_mu2_ip_below, &
    ri_sigma2_on_mu2_ip_cloud, &
    ri_sigma2_on_mu2_ip_below, &
    Nr_sigma2_on_mu2_ip_cloud, &
    Nr_sigma2_on_mu2_ip_below, &
    rs_sigma2_on_mu2_ip_cloud, &
    rs_sigma2_on_mu2_ip_below, &
    Ns_sigma2_on_mu2_ip_cloud, &
    Ns_sigma2_on_mu2_ip_below, &
    Ni_sigma2_on_mu2_ip_cloud, &
    Ni_sigma2_on_mu2_ip_below, &
    rg_sigma2_on_mu2_ip_cloud, &
    rg_sigma2_on_mu2_ip_below, &
    Ng_sigma2_on_mu2_ip_cloud, &
    Ng_sigma2_on_mu2_ip_below, &
    Ncnp2_on_Ncnm2

  use parameters_model, only : &
    hydromet_dim, &
    Lscale_max

  use parameters_tunable, only : &
    l_prescribed_avg_deltaz, &
    C7, &
    C8, &
    C11, &
    C11b, &
    gamma_coef, &
    mu, &
    mult_coef

  use pdf_parameter_module, only : &
#ifdef CLUBB_CAM /* Code for storing pdf_parameter structs in pbuf as array */
    num_pdf_params, &
#endif
    pdf_parameter

  use stat_file_module, only : &
    clubb_i, &
    clubb_j

  use stats_rad_zm, only : &
    nvarmax_rad_zm

  use stats_rad_zt, only : &
    nvarmax_rad_zt

  use stats_variables, only : &
    zt, &
    zm, &
    rad_zt, &
    rad_zm, &
    sfc, &
    l_stats_last, &
    stats_tsamp, &
    stats_tout, &
    l_output_rad_files, &
    l_stats, &
    l_stats_samp, &
    l_grads

  use stats_zm, only : &
    nvarmax_zm

  use stats_zt, only : &
    nvarmax_zt

  use variables_diagnostic_module, only : &
    Lscale, &
    wp2_zt, &
    wphydrometp

  implicit none

  private

  public  advance_clubb_core_api, setup_clubb_core_api, set_Lscale_max_api, &
    gregorian2julian_day_api, compute_current_date_api, leap_year_api, &
    setup_corr_varnce_array_api, setup_pdf_indices_api, &
    reportError_api, fatal_error_api, set_clubb_debug_level_api, clubb_at_least_debug_level_api, &
    vertical_avg_api, fill_holes_driver_api, fill_holes_vertical_api, vertical_integral_api, &
    zt2zm_api, zm2zt_api, setup_grid_api, setup_grid_heights_api, cleanup_grid_api, &
    lin_int_api, linear_interpolation_api, &
    read_parameters_api, setup_parameters_api, adj_low_res_nu_api, &
#ifdef CLUBB_CAM /* Code for storing pdf_parameter structs in pbuf as array */
    pack_pdf_params_api, unpack_pdf_params_api, &
#endif
    sat_mixrat_liq_api, &
    setup_pdf_parameters_api, &
    stats_init_api, stats_begin_timestep_api, stats_end_timestep_api, &
    stats_accumulate_hydromet_api, stats_finalize_api, &
    stats_init_rad_zm_api, &
    stats_init_rad_zt_api, &
    stats_init_zm_api, &
    stats_init_zt_api, &
    thlm2T_in_K_api

contains

  !================================================================================================
  ! advance_clubb_core
  !================================================================================================

  subroutine advance_clubb_core_api( &
    l_implemented, dt, fcor, sfc_elevation, hydromet_dim, & ! intent(in)
    thlm_forcing, rtm_forcing, um_forcing, vm_forcing, &    ! intent(in)
    sclrm_forcing, edsclrm_forcing, wprtp_forcing, &        ! intent(in)
    wpthlp_forcing, rtp2_forcing, thlp2_forcing, &          ! intent(in)
    rtpthlp_forcing, wm_zm, wm_zt, &                        ! intent(in)
    wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, &            ! intent(in)
    wpsclrp_sfc, wpedsclrp_sfc, &                           ! intent(in)
    p_in_Pa, rho_zm, rho, exner, &                          ! intent(in)
    rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &                ! intent(in)
    invrs_rho_ds_zt, thv_ds_zm, thv_ds_zt, hydromet, &      ! intent(in)
    rfrzm, radf, wphydrometp, wp2hmp, rtphmp, thlphmp, &    ! intent(in)
    um, vm, upwp, vpwp, up2, vp2, &                         ! intent(inout)
    thlm, rtm, wprtp, wpthlp, &                             ! intent(inout)
    wp2, wp3, rtp2, thlp2, rtpthlp, &                       ! intent(inout)
    sclrm,   &
#ifdef GFDL
               sclrm_trsport_only,  &  ! h1g, 2010-06-16    ! intent(inout)
#endif
    sclrp2, sclrprtp, sclrpthlp, &                          ! intent(inout)
    wpsclrp, edsclrm, err_code, &                           ! intent(inout)
#ifdef GFDL
               RH_crit, & !h1g, 2010-06-16                  ! intent(inout)
               do_liquid_only_in_clubb, &                   ! intent(in)
#endif
    rcm, wprcp, cloud_frac, ice_supersat_frac, &            ! intent(out)
    rcm_in_layer, cloud_cover, &                            ! intent(out)
#if defined(CLUBB_CAM) || defined(GFDL)
    khzm, khzt, &                                           ! intent(out)
#endif
#ifdef CLUBB_CAM
    qclvar, &                                               ! intent(out)
#endif
    pdf_params )                                            ! intent(out)

    ! Description:
    !   Subroutine to advance the model one timestep

    ! References:
    !   ``A PDF-Based Model for Boundary Layer Clouds. Part I:
    !     Method and Model Description'' Golaz, et al. (2002)
    !   JAS, Vol. 59, pp. 3540--3551.
    !-----------------------------------------------------------------------

    use advance_clubb_core_module, only : advance_clubb_core

    use parameters_model, only: &
      sclr_dim, & ! Variable(s)
      edsclr_dim

    implicit none
      !!! Input Variables
    logical, intent(in) ::  &
      l_implemented ! Is this part of a larger host model (T/F) ?

    real(kind=time_precision), intent(in) ::  &
      dt  ! Current timestep duration    [s]

    real( kind = core_rknd ), intent(in) ::  &
      fcor,  &          ! Coriolis forcing             [s^-1]
      sfc_elevation     ! Elevation of ground level    [m AMSL]

    integer, intent(in) :: &
      hydromet_dim      ! Total number of hydrometeors          [#]

    ! Input Variables
    real( kind = core_rknd ), intent(in), dimension(gr%nz) ::  &
      thlm_forcing,    & ! theta_l forcing (thermodynamic levels)    [K/s]
      rtm_forcing,     & ! r_t forcing (thermodynamic levels)        [(kg/kg)/s]
      um_forcing,      & ! u wind forcing (thermodynamic levels)     [m/s/s]
      vm_forcing,      & ! v wind forcing (thermodynamic levels)     [m/s/s]
      wprtp_forcing,   & ! <w'r_t'> forcing (momentum levels)    [m*K/s^2]
      wpthlp_forcing,  & ! <w'th_l'> forcing (momentum levels)   [m*(kg/kg)/s^2]
      rtp2_forcing,    & ! <r_t'^2> forcing (momentum levels)    [(kg/kg)^2/s]
      thlp2_forcing,   & ! <th_l'^2> forcing (momentum levels)   [K^2/s]
      rtpthlp_forcing, & ! <r_t'th_l'> forcing (momentum levels) [K*(kg/kg)/s]
      wm_zm,           & ! w mean wind component on momentum levels  [m/s]
      wm_zt,           & ! w mean wind component on thermo. levels   [m/s]
      p_in_Pa,         & ! Air pressure (thermodynamic levels)       [Pa]
      rho_zm,          & ! Air density on momentum levels            [kg/m^3]
      rho,             & ! Air density on thermodynamic levels       [kg/m^3]
      exner,           & ! Exner function (thermodynamic levels)     [-]
      rho_ds_zm,       & ! Dry, static density on momentum levels    [kg/m^3]
      rho_ds_zt,       & ! Dry, static density on thermo. levels     [kg/m^3]
      invrs_rho_ds_zm, & ! Inv. dry, static density @ momentum levs. [m^3/kg]
      invrs_rho_ds_zt, & ! Inv. dry, static density @ thermo. levs.  [m^3/kg]
      thv_ds_zm,       & ! Dry, base-state theta_v on momentum levs. [K]
      thv_ds_zt,       & ! Dry, base-state theta_v on thermo. levs.  [K]
      rfrzm              ! Total ice-phase water mixing ratio        [kg/kg]

    real( kind = core_rknd ), dimension(gr%nz,hydromet_dim), intent(in) :: &
      hydromet           ! Collection of hydrometeors                [units vary]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      radf          ! Buoyancy production at the CL top due to LW radiative cooling [m^2/s^3]

    real( kind = core_rknd ), dimension(gr%nz, hydromet_dim), intent(in) :: &
      wphydrometp, & ! Covariance of w and a hydrometeor   [(m/s) <hm units>]
      wp2hmp,      & ! Third moment: <w'^2> * <hydro.'>    [(m/s)^2 <hm units>]
      rtphmp,      & ! Covariance of rt and a hydrometeor  [(kg/kg) <hm units>]
      thlphmp        ! Covariance of thl and a hydrometeor [K <hm units>]

    real( kind = core_rknd ), intent(in) ::  &
      wpthlp_sfc,   & ! w' theta_l' at surface   [(m K)/s]
      wprtp_sfc,    & ! w' r_t' at surface       [(kg m)/( kg s)]
      upwp_sfc,     & ! u'w' at surface          [m^2/s^2]
      vpwp_sfc        ! v'w' at surface          [m^2/s^2]

    ! Passive scalar variables
    real( kind = core_rknd ), intent(in), dimension(gr%nz,sclr_dim) :: &
      sclrm_forcing    ! Passive scalar forcing         [{units vary}/s]

    real( kind = core_rknd ), intent(in),  dimension(sclr_dim) ::  &
      wpsclrp_sfc      ! Scalar flux at surface         [{units vary} m/s]

    ! Eddy passive scalar variables
    real( kind = core_rknd ), intent(in), dimension(gr%nz,edsclr_dim) :: &
      edsclrm_forcing  ! Eddy passive scalar forcing    [{units vary}/s]

    real( kind = core_rknd ), intent(in),  dimension(edsclr_dim) ::  &
      wpedsclrp_sfc    ! Eddy-Scalar flux at surface    [{units vary} m/s]

    !!! Input/Output Variables
    ! These are prognostic or are planned to be in the future
    real( kind = core_rknd ), intent(inout), dimension(gr%nz) ::  &
      um,      & ! u mean wind component (thermodynamic levels)   [m/s]
      upwp,    & ! u'w' (momentum levels)                         [m^2/s^2]
      vm,      & ! v mean wind component (thermodynamic levels)   [m/s]
      vpwp,    & ! v'w' (momentum levels)                         [m^2/s^2]
      up2,     & ! u'^2 (momentum levels)                         [m^2/s^2]
      vp2,     & ! v'^2 (momentum levels)                         [m^2/s^2]
      rtm,     & ! total water mixing ratio, r_t (thermo. levels) [kg/kg]
      wprtp,   & ! w' r_t' (momentum levels)                      [(kg/kg) m/s]
      thlm,    & ! liq. water pot. temp., th_l (thermo. levels)   [K]
      wpthlp,  & ! w' th_l' (momentum levels)                     [(m/s) K]
      rtp2,    & ! r_t'^2 (momentum levels)                       [(kg/kg)^2]
      thlp2,   & ! th_l'^2 (momentum levels)                      [K^2]
      rtpthlp, & ! r_t' th_l' (momentum levels)                   [(kg/kg) K]
      wp2,     & ! w'^2 (momentum levels)                         [m^2/s^2]
      wp3        ! w'^3 (thermodynamic levels)                    [m^3/s^3]

    ! Passive scalar variables
    real( kind = core_rknd ), intent(inout), dimension(gr%nz,sclr_dim) :: &
      sclrm,     & ! Passive scalar mean (thermo. levels) [units vary]
      wpsclrp,   & ! w'sclr' (momentum levels)            [{units vary} m/s]
      sclrp2,    & ! sclr'^2 (momentum levels)            [{units vary}^2]
      sclrprtp,  & ! sclr'rt' (momentum levels)           [{units vary} (kg/kg)]
      sclrpthlp    ! sclr'thl' (momentum levels)          [{units vary} K]

#ifdef GFDL
    real( kind = core_rknd ), intent(inout), dimension(gr%nz,sclr_dim) :: &  ! h1g, 2010-06-16
      sclrm_trsport_only  ! Passive scalar concentration due to pure transport [{units vary}/s]
#endif

      ! Eddy passive scalar variable
      real( kind = core_rknd ), intent(inout), dimension(gr%nz,edsclr_dim) :: &
      edsclrm   ! Eddy passive scalar mean (thermo. levels)   [units vary]

    ! Variables that need to be output for use in other parts of the CLUBB
    ! code, such as microphysics (rcm, pdf_params), forcings (rcm), and/or
    ! BUGSrad (cloud_cover).
    real( kind = core_rknd ), intent(out), dimension(gr%nz) ::  &
      rcm,          & ! cloud water mixing ratio, r_c (thermo. levels)  [kg/kg]
      rcm_in_layer, & ! rcm in cloud layer                              [kg/kg]
      cloud_cover     ! cloud cover                                     [-]

    type(pdf_parameter), dimension(gr%nz), intent(out) :: &
      pdf_params      ! PDF parameters   [units vary]

    ! Variables that need to be output for use in host models
    real( kind = core_rknd ), intent(out), dimension(gr%nz) ::  &
      wprcp,            & ! w'r_c' (momentum levels)                  [(kg/kg) m/s]
      cloud_frac,       & ! cloud fraction (thermodynamic levels)     [-]
      ice_supersat_frac   ! ice cloud fraction (thermodynamic levels) [-]

#if defined(CLUBB_CAM) || defined(GFDL)
    real( kind = core_rknd ), intent(out), dimension(gr%nz) :: &
      khzt, &       ! eddy diffusivity on thermo levels
      khzm          ! eddy diffusivity on momentum levels
#endif

#ifdef CLUBB_CAM
    real( kind = core_rknd), intent(out), dimension(gr%nz) :: &
      qclvar        ! cloud water variance
#endif

      !!! Output Variable
      ! Diagnostic, for if some calculation goes amiss.
      integer, intent(inout) :: err_code

#ifdef GFDL
    ! hlg, 2010-06-16
    real( kind = core_rknd ), intent(inOUT), dimension(gr%nz, min(1,sclr_dim) , 2) :: &
      RH_crit  ! critical relative humidity for droplet and ice nucleation
    ! ---> h1g, 2012-06-14
    logical, intent(in)                 ::  do_liquid_only_in_clubb
    ! <--- h1g, 2012-06-14
#endif
    call advance_clubb_core( &
      l_implemented, dt, fcor, sfc_elevation, hydromet_dim, & ! intent(in)
      thlm_forcing, rtm_forcing, um_forcing, vm_forcing, &    ! intent(in)
      sclrm_forcing, edsclrm_forcing, wprtp_forcing, &        ! intent(in)
      wpthlp_forcing, rtp2_forcing, thlp2_forcing, &          ! intent(in)
      rtpthlp_forcing, wm_zm, wm_zt, &                        ! intent(in)
      wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, &            ! intent(in)
      wpsclrp_sfc, wpedsclrp_sfc, &                           ! intent(in)
      p_in_Pa, rho_zm, rho, exner, &                          ! intent(in)
      rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &                ! intent(in)
      invrs_rho_ds_zt, thv_ds_zm, thv_ds_zt, hydromet, &      ! intent(in)
      rfrzm, radf, wphydrometp, wp2hmp, rtphmp, thlphmp, &    ! intent(in)
      um, vm, upwp, vpwp, up2, vp2, &                         ! intent(inout)
      thlm, rtm, wprtp, wpthlp, &                             ! intent(inout)
      wp2, wp3, rtp2, thlp2, rtpthlp, &                       ! intent(inout)
      sclrm,   &
#ifdef GFDL
               sclrm_trsport_only,  &  ! h1g, 2010-06-16               ! intent(inout)
#endif
      sclrp2, sclrprtp, sclrpthlp, &                          ! intent(inout)
      wpsclrp, edsclrm, err_code, &                           ! intent(inout)
#ifdef GFDL
               RH_crit, & !h1g, 2010-06-16                             ! intent(inout)
               do_liquid_only_in_clubb, &                              ! intent(in)
#endif
      rcm, wprcp, cloud_frac, ice_supersat_frac, &            ! intent(out)
      rcm_in_layer, cloud_cover, &                            ! intent(out)
#if defined(CLUBB_CAM) || defined(GFDL)
               khzm, khzt, &                                           ! intent(out)
#endif
#ifdef CLUBB_CAM
               qclvar, &                                               ! intent(out)
#endif
      pdf_params )                                            ! intent(out)
  end subroutine advance_clubb_core_api

  !================================================================================================
  ! setup_clubb_core
  !================================================================================================

  subroutine setup_clubb_core_api( &
    nzmax, T0_in, ts_nudge_in,              & ! intent(in)
    hydromet_dim_in, sclr_dim_in,           & ! intent(in)
    sclr_tol_in, edsclr_dim_in, params,     & ! intent(in)
    l_host_applies_sfc_fluxes,              & ! intent(in)
    l_uv_nudge, saturation_formula,         & ! intent(in)
#ifdef GFDL
      I_sat_sphum,                                       & ! intent(in)  h1g, 2010-06-16
#endif
    l_implemented, grid_type, deltaz, zm_init, zm_top, & ! intent(in)
    momentum_heights, thermodynamic_heights,           & ! intent(in)
    host_dx, host_dy, sfc_elevation,                   & ! intent(in)
#ifdef GFDL
      cloud_frac_min ,                                   & ! intent(in)  h1g, 2010-06-16
#endif
    err_code )                                           ! intent(out)

    !
    ! Description:
    !   Subroutine to set up the model for execution.
    !
    ! References:
    !   None
    !-------------------------------------------------------------------------

    use advance_clubb_core_module, only : setup_clubb_core

    use grid_class, only: &
      setup_grid ! Procedure

    use parameter_indices, only:  &
      nparams ! Variable(s)

    use parameters_model, only: &
      setup_parameters_model ! Procedure

    use variables_diagnostic_module, only: &
      setup_diagnostic_variables ! Procedure

    use variables_prognostic_module, only: &
      setup_prognostic_variables ! Procedure

    use model_flags, only: &
      setup_model_flags    ! Subroutine

#ifdef MKL
      use csr_matrix_class, only: &
        initialize_csr_class, & ! Subroutine
        intlc_5d_5d_ja_size     ! Variable

      use gmres_wrap, only: &
        gmres_init              ! Subroutine

      use gmres_cache, only: &
        gmres_cache_temp_init, &! Subroutine
        gmres_idx_wp2wp3        ! Variable
#endif

      implicit none

    ! Input Variables

    ! Grid definition
    integer, intent(in) :: nzmax  ! Vertical grid levels            [#]
    !                      Only true when used in a host model
    !                      CLUBB determines what nzmax should be
    !                      given zm_init and zm_top when
    !                      running in standalone mode.

    real( kind = core_rknd ), intent(in) ::  &
      sfc_elevation  ! Elevation of ground level    [m AMSL]

    ! Flag to see if CLUBB is running on it's own,
    ! or if it's implemented as part of a host model.
    logical, intent(in) :: l_implemented   ! (T/F)

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

    ! If the CLUBB model is running by itself, and is using an
    ! evenly-spaced grid (grid_type = 1), it needs the vertical
    ! grid spacing, momentum-level starting altitude, and maximum
    ! altitude as input.
    real( kind = core_rknd ), intent(in) :: &
      deltaz,   & ! Change in altitude per level           [m]
      zm_init,  & ! Initial grid altitude (momentum level) [m]
      zm_top      ! Maximum grid altitude (momentum level) [m]

    ! If the CLUBB parameterization is implemented in a host model,
    ! it needs to use the host model's momentum level altitudes
    ! and thermodynamic level altitudes.
    ! If the CLUBB model is running by itself, but is using a
    ! stretched grid entered on thermodynamic levels (grid_type = 2),
    ! it needs to use the thermodynamic level altitudes as input.
    ! If the CLUBB model is running by itself, but is using a
    ! stretched grid entered on momentum levels (grid_type = 3),
    ! it needs to use the momentum level altitudes as input.
    real( kind = core_rknd ), intent(in), dimension(nzmax) :: &
      momentum_heights,      & ! Momentum level altitudes (input)      [m]
      thermodynamic_heights    ! Thermodynamic level altitudes (input) [m]

    ! Host model horizontal grid spacing, if part of host model.
    real( kind = core_rknd ), intent(in) :: &
      host_dx,  & ! East-West horizontal grid spacing     [m]
      host_dy     ! North-South horizontal grid spacing   [m]

    ! Model parameters
    real( kind = core_rknd ), intent(in) ::  &
      T0_in, ts_nudge_in

    integer, intent(in) :: &
      hydromet_dim_in,  & ! Number of hydrometeor species
      sclr_dim_in,      & ! Number of passive scalars
      edsclr_dim_in       ! Number of eddy-diff. passive scalars

    real( kind = core_rknd ), intent(in), dimension(sclr_dim_in) :: &
      sclr_tol_in    ! Thresholds for passive scalars

    real( kind = core_rknd ), intent(in), dimension(nparams) :: &
      params  ! Including C1, nu1, nu2, etc.

    ! Flags
    logical, intent(in) ::  &
      l_uv_nudge,             & ! Wind nudging
      l_host_applies_sfc_fluxes ! Whether to apply for the surface flux

    character(len=*), intent(in) :: &
      saturation_formula ! Approximation for saturation vapor pressure

#ifdef GFDL
      logical, intent(in) :: &  ! h1g, 2010-06-16 begin mod
         I_sat_sphum

      real( kind = core_rknd ), intent(in) :: &
         cloud_frac_min         ! h1g, 2010-06-16 end mod
#endif

      ! Output variables
      integer, intent(out) :: &
      err_code   ! Diagnostic for a problem with the setup

    call setup_clubb_core &
      ( nzmax, T0_in, ts_nudge_in,              & ! intent(in)
      hydromet_dim_in, sclr_dim_in,           & ! intent(in)
      sclr_tol_in, edsclr_dim_in, params,     & ! intent(in)
      l_host_applies_sfc_fluxes,              & ! intent(in)
      l_uv_nudge, saturation_formula,         & ! intent(in)
#ifdef GFDL
      I_sat_sphum,                                       & ! intent(in)  h1g, 2010-06-16
#endif
      l_implemented, grid_type, deltaz, zm_init, zm_top, & ! intent(in)
      momentum_heights, thermodynamic_heights,           & ! intent(in)
      host_dx, host_dy, sfc_elevation,                   & ! intent(in)
#ifdef GFDL
      cloud_frac_min ,                                   & ! intent(in)  h1g, 2010-06-16
#endif
      err_code )                                           ! intent(out)

  end subroutine setup_clubb_core_api

  !================================================================================================
  ! set_Lscale_max
  !================================================================================================

  subroutine set_Lscale_max_api( &
    l_implemented, host_dx, host_dy, &
    Lscale_max )

    ! Description:
    !   This subroutine sets the value of Lscale_max, which is the maximum
    !   allowable value of Lscale.  For standard CLUBB, it is set to a very large
    !   value so that Lscale will not be limited.  However, when CLUBB is running
    !   as part of a host model, the value of Lscale_max is dependent on the size
    !   of the host model's horizontal grid spacing.  The smaller the host model's
    !   horizontal grid spacing, the smaller the value of Lscale_max.  When Lscale
    !   is limited to a small value, the value of time-scale Tau is reduced, which
    !   in turn produces greater damping on CLUBB's turbulent parameters.  This
    !   is the desired effect on turbulent parameters for a host model with small
    !   horizontal grid spacing, for small areas usually contain much less
    !   variation in meteorological quantities than large areas.

    ! References:
    !   None
    !-----------------------------------------------------------------------

    use advance_clubb_core_module, only : set_Lscale_max

    implicit none

    ! Input Variables
    logical, intent(in) :: &
      l_implemented    ! Flag to see if CLUBB is running on it's own,
    !                    or if it's implemented as part of a host model.

    real( kind = core_rknd ), intent(in) :: &
      host_dx, & ! Host model's east-west horizontal grid spacing     [m]
      host_dy    ! Host model's north-south horizontal grid spacing   [m]

    ! Output Variable
    real( kind = core_rknd ), intent(out) :: &
      Lscale_max    ! Maximum allowable value for Lscale   [m]

    call set_Lscale_max( &
      l_implemented, host_dx, host_dy, &
      Lscale_max )
  end subroutine set_Lscale_max_api

  !================================================================================================
  ! gregorian2julian_day
  !================================================================================================

  integer function gregorian2julian_day_api( &
    day, month, year )

    !
    ! Description:
    !   Computes the Julian Date (gregorian2julian), or the number of days since
    !   1 January 4713 BC, given a Gregorian Calender date (day, month, year).
    !
    ! Reference:
    !   Fliegel, H. F. and van Flandern, T. C.,
    !   Communications of the ACM, Vol. 11, No. 10 (October, 1968)
    !----------------------------------------------------------------------

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
  ! compute_current_date
  !================================================================================================

  subroutine compute_current_date_api( &
    previous_day, previous_month, &
    previous_year,  &
    seconds_since_previous_date, &
    current_day, current_month, &
    current_year, &
    seconds_since_current_date )

    !
    ! Description:
    !   Computes the current Gregorian date from a previous date and
    !   the seconds that have transpired since that date.
    !
    ! References:
    !   None
    !----------------------------------------------------------------------------

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
      previous_day, previous_month, &
      previous_year,  &
      seconds_since_previous_date, &
      current_day, current_month, &
      current_year, &
      seconds_since_current_date )
  end subroutine compute_current_date_api

  !================================================================================================
  ! leap_year
  !================================================================================================

  logical function leap_year_api( &
    year )


    !
    ! Description:
    !   Determines if the given year is a leap year.
    !
    ! References:
    !   None
    !-----------------------------------------------------------------------------

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
  ! setup_corr_varnce_array
  !================================================================================================

  subroutine setup_corr_varnce_array_api( &
    input_file_cloud, input_file_below, iunit )

    ! Description:
    !   Setup an array with the x'^2/xm^2 variables on the diagonal and the other
    !   elements to be correlations between various variables.

    ! References:
    !   None.
    !-------------------------------------------------------------------------------

    use corr_matrix_module, only : setup_corr_varnce_array

    use matrix_operations, only: mirror_lower_triangular_matrix ! Procedure

    use error_code, only: &
      clubb_debug ! Procedure

    implicit none

    ! External
    intrinsic :: max, epsilon, trim

    ! Input Variables
    integer, intent(in) :: &
      iunit ! The file unit

    character(len=*), intent(in) :: &
      input_file_cloud, & ! Path to the in cloud correlation file
      input_file_below    ! Path to the out of cloud correlation file

    call setup_corr_varnce_array( &
      input_file_cloud, input_file_below, iunit )

  end subroutine setup_corr_varnce_array_api

  !================================================================================================
  ! setup_pdf_indices
  !================================================================================================

  subroutine setup_pdf_indices_api( &
    hydromet_dim, iirrainm, iiNrm, &
    iiricem, iiNim, iirsnowm, iiNsnowm, &
    iirgraupelm, iiNgraupelm )

    ! Description:
    !
    ! Setup for the iiPDF indices. These indices are used to address s, t, w
    ! and the hydrometeors in the mean/stdev/corr arrays
    !
    ! References:
    !-----------------------------------------------------------------------

    use corr_matrix_module, only : setup_pdf_indices

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      hydromet_dim    ! Total number of hydrometeor species.

    integer, intent(in) :: &
      iirrainm,    & ! Index of rain water mixing ratio
      iiNrm,       & ! Index of rain drop concentration
      iiricem,     & ! Index of ice mixing ratio
      iiNim,       & ! Index of ice crystal concentration
      iirsnowm,    & ! Index of snow mixing ratio
      iiNsnowm,    & ! Index of snow flake concentration
      iirgraupelm, & ! Index of graupel mixing ratio
      iiNgraupelm    ! Index of graupel concentration

    call setup_pdf_indices( &
      hydromet_dim, iirrainm, iiNrm, &
      iiricem, iiNim, iirsnowm, iiNsnowm, &
      iirgraupelm, iiNgraupelm )
  end subroutine setup_pdf_indices_api

  !================================================================================================
  ! reportError
  !================================================================================================

  subroutine reportError_api( &
    err_code)

    !
    ! Description:
    !   Reports meaning of error code to console.
    !
    !-------------------------------------------------------------------------------

    use error_code, only : reportError

    implicit none

    ! Input Variable
    integer, intent(in) :: err_code ! Error Code being examined

    call reportError( &
      err_code)
  end subroutine reportError_api

  !================================================================================================
  ! fatal_error
  !================================================================================================

  elemental function fatal_error_api( &
    err_code )

    !
    ! Description: Checks to see if the err_code is one that usually
    !   causes an exit in other parts of CLUBB.
    ! References:
    !   None
    !-------------------------------------------------------------------------------

    use error_code, only : fatal_error

    implicit none

    ! Input Variable
    integer, intent(in) :: err_code ! Error Code being examined

    ! Output variable
    logical :: fatal_error_api

    fatal_error_api = fatal_error( &
      err_code )
  end function fatal_error_api

  !================================================================================================
  ! set_clubb_debug_level
  !================================================================================================

  subroutine set_clubb_debug_level_api( &
    level )

    !
    !  Description:
    !    Accessor for clubb_debug_level
    !
    !   0 => Print no debug messages to the screen
    !   1 => Print lightweight debug messages, e.g. print statements
    !   2 => Print debug messages that require extra testing,
    !        e.g. checks for NaNs and spurious negative values.
    !  References:
    !    None
    !-------------------------------------------------------------------------------

    use error_code, only : set_clubb_debug_level

    implicit none

    ! Input variable
    integer, intent(in) :: level ! The debug level being checked against the current setting

    call set_clubb_debug_level( &
      level )
  end subroutine set_clubb_debug_level_api

  !================================================================================================
  ! clubb_at_least_debug_level
  !================================================================================================

  logical function clubb_at_least_debug_level_api( &
    level )

    !
    ! Description:
    !   Checks to see if clubb has been set to a specified debug level
    !------------------------------------------------------------------

    use error_code, only : clubb_at_least_debug_level

    implicit none

    ! Input variable
    integer, intent(in) :: level   ! The debug level being checked against the current setting

    clubb_at_least_debug_level_api = clubb_at_least_debug_level( &
      level )
  end function clubb_at_least_debug_level_api

  !================================================================================================
  ! vertical_avg
  !================================================================================================

  function vertical_avg_api( &
    total_idx, rho_ds, &
    field, invrs_dz )

    ! Description:
    ! Computes the density-weighted vertical average of a field.
    !
    ! The average value of a function, f, over a set domain, [a,b], is
    ! calculated by the equation:
    !
    ! f_avg = ( INT(a:b) f*g ) / ( INT(a:b) g );
    !
    ! as long as f is continous and g is nonnegative and integrable.  Therefore,
    ! the density-weighted (by dry, static, base-static density) vertical
    ! average value of any model field, x, is calculated by the equation:
    !
    ! x_avg|_z = ( INT(z_bot:z_top) x rho_ds dz )
    !            / ( INT(z_bot:z_top) rho_ds dz );
    !
    ! where z_bot is the bottom of the vertical domain, and z_top is the top of
    ! the vertical domain.
    !
    ! This calculation is done slightly differently depending on whether x is a
    ! thermodynamic-level or a momentum-level variable.
    !
    ! Thermodynamic-level computation:

    !
    ! For numerical purposes, INT(z_bot:z_top) x rho_ds dz, which is the
    ! numerator integral, is calculated as:
    !
    ! SUM(k_bot:k_top) x(k) rho_ds(k) delta_z(k);
    !
    ! where k is the index of the given thermodynamic level, x and rho_ds are
    ! both thermodynamic-level variables, and delta_z(k) = zm(k) - zm(k-1).  The
    ! indices k_bot and k_top are the indices of the respective lower and upper
    ! thermodynamic levels involved in the integration.
    !
    ! Likewise, INT(z_bot:z_top) rho_ds dz, which is the denominator integral,
    ! is calculated as:
    !
    ! SUM(k_bot:k_top) rho_ds(k) delta_z(k).
    !
    ! The first (k=1) thermodynamic level is below ground (or below the
    ! official lower boundary at the first momentum level), so it should not
    ! count in a vertical average, whether that vertical average is used for
    ! the hole-filling scheme or for statistical purposes. Begin no lower
    ! than level k=2, which is the first thermodynamic level above ground (or
    ! above the model lower boundary).
    !
    ! For cases where hole-filling over the entire (global) vertical domain
    ! is desired, or where statistics over the entire (global) vertical
    ! domain are desired, the lower (thermodynamic-level) index of k = 2 and
    ! the upper (thermodynamic-level) index of k = gr%nz, means that the
    ! overall vertical domain will be gr%zm(gr%nz) - gr%zm(1).
    !
    !
    ! Momentum-level computation:
    !
    ! For numerical purposes, INT(z_bot:z_top) x rho_ds dz, which is the
    ! numerator integral, is calculated as:
    !
    ! SUM(k_bot:k_top) x(k) rho_ds(k) delta_z(k);
    !
    ! where k is the index of the given momentum level, x and rho_ds are both
    ! momentum-level variables, and delta_z(k) = zt(k+1) - zt(k).  The indices
    ! k_bot and k_top are the indices of the respective lower and upper momentum
    ! levels involved in the integration.
    !
    ! Likewise, INT(z_bot:z_top) rho_ds dz, which is the denominator integral,
    ! is calculated as:
    !
    ! SUM(k_bot:k_top) rho_ds(k) delta_z(k).
    !
    ! The first (k=1) momentum level is right at ground level (or right at
    ! the official lower boundary).  The momentum level variables that call
    ! the hole-filling scheme have set values at the surface (or lower
    ! boundary), and those set values should not be changed.  Therefore, the
    ! vertical average (for purposes of hole-filling) should not include the
    ! surface level (or lower boundary level).  For hole-filling purposes,
    ! begin no lower than level k=2, which is the second momentum level above
    ! ground (or above the model lower boundary).  Likewise, the value at the
    ! model upper boundary (k=gr%nz) is also set for momentum level
    ! variables.  That value should also not be changed.
    !
    ! However, this function is also used to keep track (for statistical
    ! purposes) of the vertical average of certain variables.  In that case,
    ! the vertical average needs to be taken over the entire vertical domain
    ! (level 1 to level gr%nz).
    !
    !
    ! In both the thermodynamic-level computation and the momentum-level
    ! computation, the numerator integral is divided by the denominator integral
    ! in order to find the average value (over the vertical domain) of x.

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use fill_holes, only : vertical_avg

    implicit none

    ! Input variables
    integer, intent(in) :: &
      total_idx ! The total numer of indices within the range of averaging

    real( kind = core_rknd ), dimension(total_idx), intent(in) ::  &
      rho_ds, & ! Dry, static density on either thermodynamic or momentum levels    [kg/m^3]
      field,  & ! The field (e.g. wp2) to be vertically averaged                    [Units vary]
      invrs_dz  ! Reciprocal of thermodynamic or momentum level thickness           [1/m]
                ! depending on whether we're on zt or zm grid.
    ! Note:  The rho_ds and field points need to be arranged from
    !        lowest to highest in altitude, with rho_ds(1) and
    !        field(1) actually their respective values at level k = 1.

    ! Output variable
    real( kind = core_rknd ) :: &
      vertical_avg_api  ! Vertical average of field    [Units of field]

    vertical_avg_api = vertical_avg( &
      total_idx, rho_ds, &
      field, invrs_dz )
  end function vertical_avg_api

  !================================================================================================
  ! fill_holes_driver
  !================================================================================================

  subroutine fill_holes_driver_api( &
    nz, dt, hydromet_dim,        & ! Intent(in)
    l_fill_holes_hm,             & ! Intent(in)
    rho_ds_zm, rho_ds_zt, exner, & ! Intent(in)
    thlm_mc, rvm_mc, hydromet )    ! Intent(inout)

    ! Description:
    ! Fills holes between same-phase hydrometeors(i.e. for frozen hydrometeors).
    ! The hole filling conserves water substance between all same-phase (frozen or liquid)
    ! hydrometeors at each height level.
    !
    ! Attention: The hole filling for the liquid phase hydrometeors is not yet implemented
    !
    ! Attention: l_frozen_hm and l_mix_rat_hm need to be set up before this subroutine is called!
    !
    ! References:
    !
    ! None
    !-----------------------------------------------------------------------

    use fill_holes, only : fill_holes_driver

    use constants_clubb, only: &
      four_thirds,     &
      rho_ice

    use array_index, only: &
      l_mix_rat_hm ! Variable(s)

    use index_mapping, only: &
      Nx2rx_hm_idx, & ! Procedure(s)
      mvr_hm_max

    use error_code, only: &
      clubb_at_least_debug_level ! Procedure(s)

    use stats_type_utilities, only: &
      stat_begin_update, & ! Subroutines
      stat_end_update

    implicit none

    intrinsic :: trim

    ! Input Variables
    integer, intent(in) :: hydromet_dim, nz

    logical, intent(in) :: l_fill_holes_hm

    real(kind=time_precision), intent(in) ::  &
      dt           ! Timestep         [s]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      rho_ds_zm, & ! Dry, static density on momentum levels   [kg/m^3]
      rho_ds_zt    ! Dry, static density on thermo. levels    [kg/m^3]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      exner  ! Exner function                                       [-]

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(nz, hydromet_dim), intent(inout) :: &
      hydromet

    real( kind = core_rknd ), dimension(nz), intent(inout) :: &
      rvm_mc,  & ! Microphysics contributions to vapor water            [kg/kg/s]
      thlm_mc    ! Microphysics contributions to liquid potential temp. [K/s]

    call fill_holes_driver( &
      nz, dt, hydromet_dim,        & ! Intent(in)
      l_fill_holes_hm,             & ! Intent(in)
      rho_ds_zm, rho_ds_zt, exner, & ! Intent(in)
      thlm_mc, rvm_mc, hydromet )    ! Intent(inout)
  end subroutine fill_holes_driver_api

  !================================================================================================
  ! fill_holes_vertical
  !================================================================================================

  subroutine fill_holes_vertical_api( &
    num_pts, threshold, field_grid, &
    rho_ds, rho_ds_zm, &
    field )

    ! Description:
    ! This subroutine clips values of 'field' that are below 'threshold' as much
    ! as possible (i.e. "fills holes"), but conserves the total integrated mass
    ! of 'field'.  This prevents clipping from acting as a spurious source.
    !
    ! Mass is conserved by reducing the clipped field everywhere by a constant
    ! multiplicative coefficient.
    !
    ! This subroutine does not guarantee that the clipped field will exceed
    ! threshold everywhere; blunt clipping is needed for that.

    ! References:
    !   ``Numerical Methods for Wave Equations in Geophysical Fluid
    !     Dynamics'', Durran (1999), p. 292.
    !-----------------------------------------------------------------------

    use fill_holes, only : fill_holes_vertical

    implicit none

    ! Input variables
    integer, intent(in) :: &
      num_pts  ! The number of points on either side of the hole;
               ! Mass is drawn from these points to fill the hole.  []

    real( kind = core_rknd ), intent(in) :: &
      threshold  ! A threshold (e.g. w_tol*w_tol) below which field must not
                 ! fall                           [Units vary; same as field]

    character(len=2), intent(in) :: &
      field_grid ! The grid of the field, either zt or zm

    real( kind = core_rknd ), dimension(gr%nz), intent(in) ::  &
      rho_ds,    & ! Dry, static density on thermodynamic levels    [kg/m^3]
      rho_ds_zm    ! Dry, static density on momentum levels         [kg/m^3]

    ! Input/Output variable
    real( kind = core_rknd ), dimension(gr%nz), intent(inout) :: &
      field  ! The field (e.g. wp2) that contains holes [Units same as threshold]

    call fill_holes_vertical( &
      num_pts, threshold, field_grid, &
      rho_ds, rho_ds_zm, &
      field )
  end subroutine fill_holes_vertical_api

  !================================================================================================
  ! vertical_integral
  !================================================================================================

  function vertical_integral_api( &
    total_idx, rho_ds, &
    field, invrs_dz )

    ! Description:
    ! Computes the vertical integral. rho_ds, field, and invrs_dz must all be
    ! of size total_idx and should all start at the same index.
    !

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use fill_holes, only : vertical_integral

    implicit none

    ! Input variables
    integer, intent(in) :: &
      total_idx  ! The total numer of indices within the range of averaging

    real( kind = core_rknd ), dimension(total_idx), intent(in) ::  &
      rho_ds,  & ! Dry, static density                   [kg/m^3]
      field,   & ! The field to be vertically averaged   [Units vary]
      invrs_dz   ! Level thickness                       [1/m]
    ! Note:  The rho_ds and field points need to be arranged from
    !        lowest to highest in altitude, with rho_ds(1) and
    !        field(1) actually their respective values at level k = begin_idx.

    real( kind = core_rknd ) :: &
      vertical_integral_api ! Integral in the numerator (see description)

    vertical_integral_api = vertical_integral( &
      total_idx, rho_ds, &
      field, invrs_dz )
  end function vertical_integral_api

  !================================================================================================
  ! setup_grid
  !================================================================================================

  subroutine setup_grid_api( &
    nzmax, sfc_elevation, l_implemented,      &
    grid_type, deltaz, zm_init, zm_top,      &
    momentum_heights, thermodynamic_heights, &
    begin_height, end_height )

    ! Description:
    !   Grid Constructor
    !
    !   This subroutine sets up the CLUBB vertical grid.
    !
    ! References:
    !   ``Equations for CLUBB'',  Sec. 8,  Grid Configuration.
    !-----------------------------------------------------------------------

    use grid_class, only : setup_grid

    use error_code, only:  &
      clubb_at_least_debug_level ! Procedure(s)

    implicit none

    ! Input Variables
    integer, intent(in) ::  &
      nzmax  ! Number of vertical levels in grid      [#]

    real( kind = core_rknd ), intent(in) ::  &
      sfc_elevation  ! Elevation of ground level    [m AMSL]

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
      zm_init,  & ! Initial grid altitude (momentum level) [m]
      zm_top      ! Maximum grid altitude (momentum level) [m]

    ! If the CLUBB parameterization is implemented in a host model, it needs to
    ! use the host model's momentum level altitudes and thermodynamic level
    ! altitudes.
    ! If the CLUBB model is running by itself, but is using a stretched grid
    ! entered on thermodynamic levels (grid_type = 2), it needs to use the
    ! thermodynamic level altitudes as input.
    ! If the CLUBB model is running by itself, but is using a stretched grid
    ! entered on momentum levels (grid_type = 3), it needs to use the momentum
    ! level altitudes as input.
    real( kind = core_rknd ), intent(in), dimension(nzmax) ::  &
      momentum_heights,   & ! Momentum level altitudes (input)      [m]
      thermodynamic_heights ! Thermodynamic level altitudes (input) [m]

    integer, intent(out) :: &
      begin_height, &  ! Lower bound for *_heights arrays [-]
      end_height       ! Upper bound for *_heights arrays [-]

    call setup_grid( &
      nzmax, sfc_elevation, l_implemented,      &
      grid_type, deltaz, zm_init, zm_top,      &
      momentum_heights, thermodynamic_heights, &
      begin_height, end_height )

  end subroutine setup_grid_api

  !================================================================================================
  ! setup_grid_heights
  !================================================================================================

  subroutine setup_grid_heights_api( &
    l_implemented, grid_type,  &
    deltaz, zm_init, momentum_heights,  &
    thermodynamic_heights )

    ! Description:
    !   Sets the heights and interpolation weights of the column.
    !   This is seperated from setup_grid for those host models that have heights
    !   that vary with time.
    ! References:
    !   None
    !------------------------------------------------------------------------------

    use grid_class, only : setup_grid_heights, gr

    implicit none

    ! Input Variables

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
    real( kind = core_rknd ), intent(in), dimension(gr%nz) ::  &
      momentum_heights,   & ! Momentum level altitudes (input)      [m]
      thermodynamic_heights ! Thermodynamic level altitudes (input) [m]

    call setup_grid_heights( &
      l_implemented, grid_type,  &
      deltaz, zm_init, momentum_heights,  &
      thermodynamic_heights )

  end subroutine setup_grid_heights_api

  !================================================================================================
  ! cleanup_grid
  !================================================================================================

  subroutine cleanup_grid_api

    ! Description:
    !   De-allocates the memory for the grid
    !
    ! References:
    !   None
    !------------------------------------------------------------------------------

    use grid_class, only : cleanup_grid

    implicit none

    call cleanup_grid

  end subroutine cleanup_grid_api

  !================================================================================================
  ! lin_int
  !================================================================================================

  function lin_int_api( &
    height_int, height_high, height_low, &
    var_high, var_low )

    ! Description:
    ! This function computes a linear interpolation of the value of variable.
    ! Given two known values of a variable at two height values, the value
    ! of that variable at a height between those two height levels (rather
    ! than a height outside of those two height levels) is computed.
    !
    ! Here is a diagram:
    !
    !  ################################ Height high, know variable value
    !
    !
    !
    !  -------------------------------- Height to be interpolated to; linear interpolation
    !
    !
    !
    !
    !
    !  ################################ Height low, know variable value
    !
    !
    ! FORMULA:
    !
    ! variable(@ Height interpolation) =
    !
    ! [ (variable(@ Height high) - variable(@ Height low)) / (Height high - Height low) ]
    ! * (Height interpolation - Height low)  +  variable(@ Height low)

    ! Comments from WRF-HOC, Brian Griffin.

    ! References:
    ! None
    !-------------------------------------------------------------------------------

    use interpolation, only : lin_int

    implicit none

    real( kind = core_rknd ), intent(in) :: &
      height_int,  & ! Height to be interpolated to     [m]
      height_high, & ! Height above the interpolation   [m]
      height_low,  & ! Height below the interpolation   [m]
      var_high,    & ! Variable above the interpolation [units vary]
      var_low        ! Variable below the interpolation [units vary]

    ! Output Variables
    real( kind = core_rknd ) :: lin_int_api

    lin_int_api = lin_int( &
      height_int, height_high, height_low, &
      var_high, var_low )

  end function lin_int_api

  !================================================================================================
  ! linear_interpolation
  !================================================================================================

  subroutine linear_interpolation_api( &
    nparam, xlist, tlist, xvalue, tvalue )

    ! Description:
    !   Linear interpolation for 25 June 1996 altocumulus case.

    !   For example, to interpolate between two temperatures in space, put
    !   your spatial coordinates in x-list and your temperature values in
    !   tlist.  The point in question should have its spatial value stored
    !   in xvalue, and tvalue will be the temperature at that point.

    ! Author: Michael Falk for COAMPS.
    !-------------------------------------------------------------------------------

    use interpolation, only : linear_interpolation

    use error_code, only: clubb_debug ! Procedure

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

    call linear_interpolation( &
      nparam, xlist, tlist, xvalue, tvalue )

  end subroutine linear_interpolation_api

  !================================================================================================
  ! read_parameters
  !================================================================================================

  subroutine read_parameters_api( &
    iunit, filename, params )

    ! Description:
    ! Read a namelist containing the model parameters

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use parameters_tunable, only : read_parameters

    use parameter_indices, only:  &
      nparams ! Variable(s)

    implicit none

    ! Input variables
    integer, intent(in) :: iunit

    character(len=*), intent(in) :: filename

    ! Output variables
    real( kind = core_rknd ), intent(out), dimension(nparams) :: params

    call read_parameters( &
      iunit, filename, params )

  end subroutine read_parameters_api

  !================================================================================================
  ! setup_parameters
  !================================================================================================

  subroutine setup_parameters_api( &
    deltaz, params, nzmax, &
    grid_type, momentum_heights, thermodynamic_heights, &
    err_code )

    ! Description:
    ! Subroutine to setup model parameters

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
      fstderr ! Variable(s)

    use error_code, only:  &
      clubb_var_out_of_bounds ! Variable(s)

    use parameter_indices, only:  &
      nparams ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) ::  &
      deltaz  ! Change per height level        [m]

    real( kind = core_rknd ), intent(in), dimension(nparams) :: &
      params  ! Tuneable model parameters      [-]

    ! Grid definition
    integer, intent(in) :: nzmax  ! Vertical grid levels            [#]

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

    ! If the CLUBB parameterization is implemented in a host model,
    ! it needs to use the host model's momentum level altitudes
    ! and thermodynamic level altitudes.
    ! If the CLUBB model is running by itself, but is using a
    ! stretched grid entered on thermodynamic levels (grid_type = 2),
    ! it needs to use the thermodynamic level altitudes as input.
    ! If the CLUBB model is running by itself, but is using a
    ! stretched grid entered on momentum levels (grid_type = 3),
    ! it needs to use the momentum level altitudes as input.
    real( kind = core_rknd ), intent(in), dimension(nzmax) :: &
      momentum_heights,      & ! Momentum level altitudes (input)      [m]
      thermodynamic_heights    ! Thermodynamic level altitudes (input) [m]

    ! Output Variables
    integer, intent(out) ::  &
      err_code ! Error condition

    call setup_parameters( &
      deltaz, params, nzmax, &
      grid_type, momentum_heights, thermodynamic_heights, &
      err_code )

  end subroutine setup_parameters_api

  !================================================================================================
  ! adj_low_res_nu
  !================================================================================================

  subroutine adj_low_res_nu_api( &
    nzmax, grid_type, deltaz, & ! Intent(in)
    momentum_heights, thermodynamic_heights )  ! Intent(in)

    ! Description:
    !   Adjust the values of background eddy diffusivity based on
    !   vertical grid spacing.
    !   This code was made into a public subroutine so that it may be
    !   called multiple times per model run in scenarios where grid
    !   altitudes, and hence average grid spacing, change through space
    !   and/or time.  This occurs, for example, when CLUBB is
    !   implemented in WRF.  --ldgrant Jul 2010
    !----------------------------------------------------------------------

    use parameters_tunable, only : adj_low_res_nu

    implicit none

    ! Input Variables

    ! Grid definition
    integer, intent(in) :: nzmax  ! Vertical grid levels            [#]

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

    ! If the CLUBB parameterization is implemented in a host model,
    ! it needs to use the host model's momentum level altitudes
    ! and thermodynamic level altitudes.
    ! If the CLUBB model is running by itself, but is using a
    ! stretched grid entered on thermodynamic levels (grid_type = 2),
    ! it needs to use the thermodynamic level altitudes as input.
    ! If the CLUBB model is running by itself, but is using a
    ! stretched grid entered on momentum levels (grid_type = 3),
    ! it needs to use the momentum level altitudes as input.
    real( kind = core_rknd ), intent(in), dimension(nzmax) :: &
      momentum_heights,      & ! Momentum level altitudes (input)      [m]
      thermodynamic_heights    ! Thermodynamic level altitudes (input) [m]

    call adj_low_res_nu( &
      nzmax, grid_type, deltaz, & ! Intent(in)
      momentum_heights, thermodynamic_heights )  ! Intent(in)
  end subroutine adj_low_res_nu_api

#ifdef CLUBB_CAM /* Code for storing pdf_parameter structs in pbuf as array */
  !================================================================================================
  ! pack_pdf_params
  !================================================================================================

  subroutine pack_pdf_params_api( &
    pdf_params, nz, r_param_array)



    use pdf_parameter_module, only : pack_pdf_params

    !use statements

    implicit none

    ! Input a pdf_parameter array with nz instances of pdf_parameter
    integer, intent(in) :: nz ! Num Vert Model Levs
    type (pdf_parameter), dimension(nz), intent(in) :: pdf_params

    ! Output a two dimensional real array with all values
    real (kind = core_rknd), dimension(nz,num_pdf_params), intent(out) :: &
      r_param_array

    call pack_pdf_params( &
      pdf_params, nz, r_param_array)

  end subroutine pack_pdf_params_api

  !================================================================================================
  ! unpack_pdf_params
  !================================================================================================

  subroutine unpack_pdf_params_api( &
    r_param_array, nz, pdf_params)

    use pdf_parameter_module, only : unpack_pdf_params

    implicit none

    ! Input a two dimensional real array with pdf values
    integer, intent(in) :: nz ! Num Vert Model Levs
    real (kind = core_rknd), dimension(nz,num_pdf_params), intent(in) :: &
      r_param_array

    ! Output a pdf_parameter array with nz instances of pdf_parameter
    type (pdf_parameter), dimension(nz), intent(out) :: pdf_params

    call unpack_pdf_params( &
      r_param_array, nz, pdf_params)
  end subroutine unpack_pdf_params_api
#endif

  !================================================================================================
  ! sat_mixrat_liq
  !================================================================================================

  elemental real( kind = core_rknd ) function  sat_mixrat_liq_api( &
    p_in_Pa, T_in_K )

    ! Description:
    !   Used to compute the saturation mixing ratio of liquid water.

    ! References:
    !   Formula from Emanuel 1994, 4.4.14
    !-------------------------------------------------------------------------

    use saturation, only : sat_mixrat_liq

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) ::  &
      p_in_Pa,  & ! Pressure    [Pa]
      T_in_K      ! Temperature [K]

    sat_mixrat_liq_api =  sat_mixrat_liq( &
      p_in_Pa, T_in_K )
  end function sat_mixrat_liq_api

  !================================================================================================
  ! setup_pdf_parameters
  !================================================================================================

  subroutine setup_pdf_parameters_api( &
    nz, d_variables, dt, rho, &                  ! Intent(in)
    wp2_zt, Nc_in_cloud, rcm, cloud_frac, &      ! Intent(in)
    ice_supersat_frac, hydromet, wphydrometp, &  ! Intent(in)
    corr_array_cloud, corr_array_below, &        ! Intent(in)
    pdf_params, l_stats_samp, &                  ! Intent(in)
    mu_x_1_n, mu_x_2_n, &                        ! Intent(out)
    sigma_x_1_n, sigma_x_2_n, &                  ! Intent(out)
    corr_array_1_n, corr_array_2_n, &            ! Intent(out)
    corr_cholesky_mtx_1, corr_cholesky_mtx_2, &  ! Intent(out)
    hydromet_pdf_params, hydrometp2 )            ! Intent(out)

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use setup_clubb_pdf_params, only : setup_pdf_parameters

    use constants_clubb, only: &
      one,            & ! Constant(s)
      rc_tol,         &
      Ncn_tol,        &
      cloud_frac_min

    use model_flags, only: &
      l_use_modified_corr

    use Nc_Ncn_eqns, only: &
      Nc_in_cloud_to_Ncnm  ! Procedure(s)

    use advance_windm_edsclrm_module, only: &
      xpwp_fnc

    use variables_diagnostic_module, only: &
      Kh_zm

    use parameters_tunable, only: &
      c_K_hm

    use pdf_utilities, only: &
      calc_xp2  ! Procedure(s)

    use clip_explicit, only: &
      clip_covar_level, & ! Procedure(s)
      clip_wphydrometp    ! Variables(s)

    use matrix_operations, only: &
      Cholesky_factor, & ! Procedure(s)
      mirror_lower_triangular_matrix

    use stats_type_utilities, only: &
      stat_update_var,    & ! Procedure(s)
      stat_update_var_pt

    use stats_variables, only: &
      ihm1,           & ! Variable(s)
      ihm2,           &
      iprecip_frac,   &
      iprecip_frac_1, &
      iprecip_frac_2, &
      iNcnm,          &
      ihmp2_zt,       &
      zt

    use model_flags, only: &
      l_diagnose_correlations ! Variable(s)

    use diagnose_correlations_module, only: &
      diagnose_correlations, & ! Procedure(s)
      calc_cholesky_corr_mtx_approx

    use corr_matrix_module, only: &
      assert_corr_symmetric ! Procedure(s)

    use index_mapping, only: &
      hydromet2pdf_idx    ! Procedure(s)

    use error_code, only : &
      clubb_at_least_debug_level   ! Procedure(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz,          & ! Number of model vertical grid levels
      d_variables    ! Number of variables in the correlation array

    real( kind = time_precision ), intent(in) ::  &
      dt    ! Model timestep                                           [s]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      rho,         & ! Density                                         [kg/m^3]
      wp2_zt,      & ! Variance of w, <w'^2> (interp. to t-levs.)      [m^2/s^2]
      Nc_in_cloud    ! Mean (in-cloud) cloud droplet concentration     [num/kg]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      rcm,               & ! Mean cloud water mixing ratio, < r_c >    [kg/kg]
      cloud_frac,        & ! Cloud fraction                            [-]
      ice_supersat_frac    ! Ice supersaturation fraction              [-]

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(in) :: &
      hydromet,    & ! Mean of hydrometeor, hm (overall) (t-levs.) [units]
      wphydrometp    ! Covariance < w'h_m' > (momentum levels)     [(m/s)units]

    real( kind = core_rknd ), dimension(d_variables,d_variables), &
      intent(in) :: &
      corr_array_cloud, & ! Prescribed correlation array in cloud      [-]
      corr_array_below    ! Prescribed correlation array below cloud   [-]

    type(pdf_parameter), dimension(nz), intent(in) :: &
      pdf_params    ! PDF parameters                               [units vary]

    logical, intent(in) :: &
      l_stats_samp    ! Flag to sample statistics

    ! Output Variables
    real( kind = core_rknd ), dimension(d_variables,d_variables,nz), &
      intent(out) :: &
      corr_array_1_n, & ! Corr. array (normalized) of PDF vars. (comp. 1)    [-]
      corr_array_2_n    ! Corr. array (normalized) of PDF vars. (comp. 2)    [-]

    real( kind = core_rknd ), dimension(d_variables, nz), intent(out) :: &
      mu_x_1_n,    & ! Mean array (normalized) of PDF vars. (comp. 1) [un. vary]
      mu_x_2_n,    & ! Mean array (normalized) of PDF vars. (comp. 2) [un. vary]
      sigma_x_1_n, & ! Std. dev. array (normalized) of PDF vars (comp. 1) [u.v.]
      sigma_x_2_n    ! Std. dev. array (normalized) of PDF vars (comp. 2) [u.v.]

    type(hydromet_pdf_parameter), dimension(nz), intent(out) :: &
      hydromet_pdf_params    ! Hydrometeor PDF parameters        [units vary]

    real( kind = core_rknd ), dimension(d_variables,d_variables,nz), &
      intent(out) :: &
      corr_cholesky_mtx_1, & ! Transposed corr. cholesky matrix, 1st comp. [-]
      corr_cholesky_mtx_2    ! Transposed corr. cholesky matrix, 2nd comp. [-]

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(out) :: &
      hydrometp2    ! Variance of a hydrometeor (overall) (m-levs.)   [units^2]

    call setup_pdf_parameters( &
      nz, d_variables, dt, rho, &                  ! Intent(in)
      wp2_zt, Nc_in_cloud, rcm, cloud_frac, &      ! Intent(in)
      ice_supersat_frac, hydromet, wphydrometp, &  ! Intent(in)
      corr_array_cloud, corr_array_below, &        ! Intent(in)
      pdf_params, l_stats_samp, &                  ! Intent(in)
      mu_x_1_n, mu_x_2_n, &                        ! Intent(out)
      sigma_x_1_n, sigma_x_2_n, &                  ! Intent(out)
      corr_array_1_n, corr_array_2_n, &            ! Intent(out)
      corr_cholesky_mtx_1, corr_cholesky_mtx_2, &  ! Intent(out)
      hydromet_pdf_params, hydrometp2 )            ! Intent(out)

  end subroutine setup_pdf_parameters_api

  !================================================================================================
  ! stats_init
  !================================================================================================

  subroutine stats_init_api( &
    iunit, fname_prefix, fdir, l_stats_in, &
    stats_fmt_in, stats_tsamp_in, stats_tout_in, fnamelist, &
    nzmax, nlon, nlat, gzt, gzm, nnrad_zt, &
    grad_zt, nnrad_zm, grad_zm, day, month, year, &
    rlon, rlat, time_current, delt )

    !
    ! Description:
    !   Initializes the statistics saving functionality of the CLUBB model.
    !
    ! References:
    !   None
    !-----------------------------------------------------------------------

    use stats_clubb_utilities, only : stats_init

    implicit none

    ! Input Variables
    integer, intent(in) :: iunit  ! File unit for fnamelist

    character(len=*), intent(in) ::  &
      fname_prefix, & ! Start of the stats filenames
      fdir            ! Directory to output to

    logical, intent(in) :: &
      l_stats_in      ! Stats on? T/F

    character(len=*), intent(in) :: &
      stats_fmt_in    ! Format of the stats file output

    real(kind=time_precision), intent(in) ::  &
      stats_tsamp_in,  & ! Sampling interval   [s]
      stats_tout_in      ! Output interval     [s]

    character(len=*), intent(in) :: &
      fnamelist          ! Filename holding the &statsnl

    integer, intent(in) :: &
      nlon, & ! Number of points in the X direction [-]
      nlat, & ! Number of points in the Y direction [-]
      nzmax   ! Grid points in the vertical         [-]

    real( kind = core_rknd ), intent(in), dimension(nzmax) ::  &
      gzt, gzm  ! Thermodynamic and momentum levels           [m]

    integer, intent(in) :: nnrad_zt ! Grid points in the radiation grid [count]

    real( kind = core_rknd ), intent(in), dimension(nnrad_zt) :: grad_zt ! Radiation levels [m]

    integer, intent(in) :: nnrad_zm ! Grid points in the radiation grid [count]

    real( kind = core_rknd ), intent(in), dimension(nnrad_zm) :: grad_zm ! Radiation levels [m]

    integer, intent(in) :: day, month, year  ! Time of year

    real( kind = core_rknd ), dimension(nlon), intent(in) ::  &
      rlon  ! Longitude(s) [Degrees E]

    real( kind = core_rknd ), dimension(nlat), intent(in) ::  &
      rlat  ! Latitude(s)  [Degrees N]

    real(kind=time_precision), intent(in) ::  &
      time_current ! Model time                         [s]

    real(kind=time_precision), intent(in) ::  &
      delt         ! Timestep (dt_main in CLUBB)         [s]

    call stats_init( &
      iunit, fname_prefix, fdir, l_stats_in, &
      stats_fmt_in, stats_tsamp_in, stats_tout_in, fnamelist, &
      nzmax, nlon, nlat, gzt, gzm, nnrad_zt, &
      grad_zt, nnrad_zm, grad_zm, day, month, year, &
      rlon, rlat, time_current, delt )
  end subroutine stats_init_api

  !================================================================================================
  ! stats_begin_timestep
  !================================================================================================

  subroutine stats_begin_timestep_api( &
    time_elapsed )

    !     Description:
    !       Given the elapsed time, set flags determining specifics such as
    !       if this time set should be sampled or if this is the first or
    !       last time step.
    !-----------------------------------------------------------------------

    use stats_clubb_utilities, only : stats_begin_timestep

    implicit none

    ! External
    intrinsic :: mod

    ! Input Variable(s)
    real(kind=time_precision), intent(in) ::  &
      time_elapsed ! Elapsed model time       [s]

    call stats_begin_timestep( &
      time_elapsed )
  end subroutine stats_begin_timestep_api

  !================================================================================================
  ! stats_end_timestep
  !================================================================================================

  subroutine stats_end_timestep_api

    ! Description:
    !   Called when the stats timestep has ended. This subroutine
    !   is responsible for calling statistics to be written to the output
    !   format.
    !
    ! References:
    !   None
    !-----------------------------------------------------------------------

    use stats_clubb_utilities, only : stats_end_timestep

    implicit none

    call stats_end_timestep

  end subroutine stats_end_timestep_api

  !================================================================================================
  ! stats_accumulate_hydromet
  !================================================================================================

  subroutine stats_accumulate_hydromet_api( &
    hydromet, rho_ds_zt )

    ! Description:
    !   Compute stats related the hydrometeors

    ! References:
    !   None
    !------------------------------------------------------------------------------

    use stats_clubb_utilities, only : stats_accumulate_hydromet

    implicit none

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz,hydromet_dim), intent(in) :: &
      hydromet ! All hydrometeors except for rcm        [units vary]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      rho_ds_zt ! Dry, static density (thermo. levs.)      [kg/m^3]

    call stats_accumulate_hydromet( &
      hydromet, rho_ds_zt )
  end subroutine stats_accumulate_hydromet_api

  !================================================================================================
  ! stats_finalize
  !================================================================================================

  subroutine stats_finalize_api

    !     Description:
    !     Close NetCDF files and deallocate scratch space and
    !     stats file structures.
    !-----------------------------------------------------------------------

    use stats_clubb_utilities, only : stats_finalize

    implicit none

    call stats_finalize

  end subroutine stats_finalize_api

  !================================================================================================
  ! stats_init_rad_zm
  !================================================================================================

  subroutine stats_init_rad_zm_api( &
    vars_rad_zm, l_error )
    use stats_rad_zm, only : stats_init_rad_zm, nvarmax_rad_zm

    !     Description:
    !     Initializes array indices for rad_zm variables
    !-----------------------------------------------------------------------

    implicit none

    ! Input Variable
    character(len= * ), dimension(nvarmax_rad_zm), intent(in) :: vars_rad_zm

    ! Input/Output Variable
    logical, intent(inout) :: l_error

    call stats_init_rad_zm( &
      vars_rad_zm, l_error )
  end subroutine stats_init_rad_zm_api

  !================================================================================================
  ! stats_init_rad_zt
  !================================================================================================

  subroutine stats_init_rad_zt_api( &
    vars_rad_zt, l_error )

    ! Description:
    !   Initializes array indices for zt
    !
    ! References:
    !   None
    !-----------------------------------------------------------------------

    use stats_rad_zt, only : stats_init_rad_zt, nvarmax_rad_zt

    implicit none

    ! Input Variable
    character(len= * ), dimension(nvarmax_rad_zt), intent(in) :: vars_rad_zt

    ! Input/Output Variable
    logical, intent(inout) :: l_error

    call stats_init_rad_zt( &
      vars_rad_zt, l_error )
  end subroutine stats_init_rad_zt_api

  !================================================================================================
  ! stats_init_zm
  !================================================================================================

  subroutine stats_init_zm_api( &
    vars_zm, l_error )

    ! Description:
    !   Initializes array indices for zm

    ! Note:
    !   All code that is within subroutine stats_init_zm, including variable
    !   allocation code, is not called if l_stats is false.  This subroutine is
    !   called only when l_stats is true.

    !-----------------------------------------------------------------------

    use stats_zm, only : stats_init_zm, nvarmax_zm

    implicit none

    ! Input Variable
    character(len= * ), dimension(nvarmax_zm), intent(in) :: vars_zm ! zm variable names

    ! Input / Output Variable
    logical, intent(inout) :: l_error

    call stats_init_zm( &
      vars_zm, l_error )

  end subroutine stats_init_zm_api

  !================================================================================================
  ! stats_init_zt
  !================================================================================================

  subroutine stats_init_zt_api( &
    vars_zt, l_error )

    ! Description:
    ! Initializes array indices for zt

    ! Note:
    ! All code that is within subroutine stats_init_zt, including variable
    ! allocation code, is not called if l_stats is false.  This subroutine is
    ! called only when l_stats is true.

    !-----------------------------------------------------------------------

    use stats_zt, only : stats_init_zt, nvarmax_zt

    implicit none

    ! Input Variable
    character(len= * ), dimension(nvarmax_zt), intent(in) :: vars_zt

    ! Input / Output Variable
    logical, intent(inout) :: l_error

    call stats_init_zt( &
      vars_zt, l_error )

  end subroutine stats_init_zt_api

  !================================================================================================
  ! thlm2T_in_K
  !================================================================================================

  elemental function thlm2T_in_K_api( &
    thlm, exner, rcm )  &
    result( T_in_K )

    ! Description:
    !   Calculates absolute temperature from liquid water potential
    !   temperature.  (Does not include ice.)

    ! References:
    !   Cotton and Anthes (1989), "Storm and Cloud Dynamics", Eqn. (2.51).
    !-------------------------------------------------------------------------------

    use T_in_K_module, only : thlm2T_in_K

    implicit none

    ! Input
    real( kind = core_rknd ), intent(in) :: &
      thlm,   & ! Liquid potential temperature  [K]
      exner,  & ! Exner function                [-]
      rcm       ! Liquid water mixing ratio     [kg/kg]

    real( kind = core_rknd ) :: &
      T_in_K ! Result temperature [K]

    T_in_K = thlm2T_in_K( &
      thlm, exner, rcm )

  end function thlm2T_in_K_api

end module clubb_api_module
