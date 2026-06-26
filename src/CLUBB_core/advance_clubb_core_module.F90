
module advance_clubb_core_module


! Description:
!   The module containing the `core' of the CLUBB parameterization.
!   It advances CLUBB's equations one model time step.
!
! References:
! https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:overview_clubb
!
!  ``A PDF-Based Model for Boundary Layer Clouds. Part I:
!    Method and Model Description'' Golaz, et al. (2002)
!    JAS, Vol. 59, pp. 3540--3551.
!
!                         Copyright Notice:
!
!   This code and the source code it references are (C) 2006-2020.
!
!   The distribution of this code and derived works thereof
!                   should include this notice.
!
!   Portions of this code derived from other sources (Hugh Morrison,
!   ACM TOMS, et cetera) are the intellectual
!   property of their respective authors as noted and are also subject
!   to copyright.
!
!
!
! Cloud Layers Unified By Binormals (CLUBB) user license
! agreement.
!
! Thank you for your interest in CLUBB. We work hard to create a
! code that implements the best software engineering practices,
! is supported to the extent allowed by our limited resources,
! and is available without cost to non-commercial users. You may
! use CLUBB if, in return, you abide by these conditions:
!
! 1. Please cite CLUBB in presentations and publications that
!  contain results obtained using CLUBB.
!
! 2. You may not use any part of CLUBB to create or modify
!  another single-column (1D) model that is not called CLUBB.
!  However, you may modify or augment CLUBB or parts of CLUBB if
!  you include "CLUBB" in the name of the resulting single-column
!  model. For example, a user at MIT might modify CLUBB and call
!  the modified version "CLUBB-MIT." Or, for example, a user of
!  the CLM land-surface model might interface CLM to CLUBB and
!  call it "CLM-CLUBB." This naming convention recognizes the
!  contributions of both sets of developers.
!
! 3. You may implement CLUBB as a parameterization in a large-
!  scale host model that has 2 or 3 spatial dimensions without
!  including "CLUBB" in the combined model name, but please
!  acknowledge in presentations and publications that CLUBB has
!  been included as a parameterization.
!
! 4. You may not provide all or part of CLUBB to anyone without
!  prior permission from Vincent Larson (vlarson@uwm.edu). If
!  you wish to share CLUBB with your collaborators without
!  seeking permission, please ask your collaborators to register
!  as CLUBB users at https://carson.math.uwm.edu/larson-group/clubb_site/ and to
!  download CLUBB from there.
!
! 5. You may not use CLUBB for commercial purposes unless you
!  receive permission from Vincent Larson.
!
! 6. You may not re-license all or any part of CLUBB.
!
! 7. CLUBB is provided "as is" and without warranty.
!
! We hope that CLUBB will develop into a community resource. We
! encourage users to contribute their CLUBB modifications or
! extensions to the CLUBB development group. We will then
! consider them for inclusion in CLUBB. Such contributions will
! benefit all CLUBB users. We would be pleased to acknowledge
! contributors and list their CLUBB-related papers on our "About
! CLUBB" webpage (https://carson.math.uwm.edu/larson-group/clubb_site/about.html) for
! those contributors who so desire.
!
! Thanks so much and best wishes for your research!
!
! The CLUBB Development Group
! (Present and past contributors to the source code include
! Vincent Larson, Chris Golaz, David Schanen, Brian Griffin,
! Joshua Fasching, Adam Smith, and Michael Falk).
!-----------------------------------------------------------------------

  ! Options for the placement of the call to CLUBB's PDF.
  use model_flags, only: &
      ipdf_pre_advance_fields, &      ! Call before advancing predictive fields
      ipdf_post_advance_fields, &     ! Call after advancing predictive fields
      ipdf_pre_post_advance_fields    ! Call both before and after advancing
                                      ! predictive fields


  implicit none

  public ::  &
    advance_clubb_core

  private ! Default Scope

  contains

  !-----------------------------------------------------------------------

  !#######################################################################
  !#######################################################################
  ! If you change the argument list of advance_clubb_core you also have to
  ! change the calls to this function in the host models CAM, WRF and SAM.
  !#######################################################################
  !#######################################################################

  subroutine advance_clubb_core( gr, nzm, nzt, ngrdcol, &           ! In
               l_implemented, dt, fcor, fcor_y, sfc_elevation, &    ! In
               hydromet_dim,                                      & ! In
               sclr_dim, sclr_tol, edsclr_dim, sclr_idx,      &     ! In
               thlm_forcing, rtm_forcing, um_forcing, vm_forcing, & ! In
               sclrm_forcing, edsclrm_forcing, wprtp_forcing, &     ! In
               wpthlp_forcing, rtp2_forcing, thlp2_forcing, &       ! In
               rtpthlp_forcing, wm_zm, wm_zt, &                     ! In
               wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, p_sfc, &  ! In
               wpsclrp_sfc, wpedsclrp_sfc, &                        ! In
               upwp_sfc_pert, vpwp_sfc_pert, &                      ! In
               rtm_ref, thlm_ref, um_ref, vm_ref, ug, vg, &         ! In
               p_in_Pa, rho_zm, rho, exner, &                       ! In
               rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &             ! In
               invrs_rho_ds_zt, thv_ds_zm, thv_ds_zt, &             ! In
               l_mix_rat_hm, &                                      ! In
               rfrzm, &                                             ! In
               wphydrometp, wp2hmp, rtphmp_zt, thlphmp_zt, &        ! In
               host_dx, host_dy, &                                  ! In
               clubb_params, nu_vert_res_dep, lmin, &               ! In
               mixt_frac_max_mag, T0, ts_nudge, &                   ! In
               rtm_min, rtm_nudge_max_altitude, &                   ! In
               clubb_config_flags, &                                ! In
               stats,         &                                     ! InOut
               um, vm, upwp, vpwp, up2, vp2, up3, vp3, &            ! InOut
               thlm, rtm, wprtp, wpthlp, &                          ! InOut
               wp2, wp3, rtp2, rtp3, thlp2, thlp3, rtpthlp, &       ! InOut
               sclrm,   &                                           ! InOut
               sclrp2, sclrp3, sclrprtp, sclrpthlp, &               ! InOut
               wpsclrp, edsclrm, &                                  ! InOut
               rcm, cloud_frac, &                                   ! InOut
               wpthvp, wp2thvp, wp2up, rtpthvp, thlpthvp, &         ! InOut
               sclrpthvp, &                                         ! InOut
               wp2rtp, wp2thlp, uprcp, vprcp, rc_coef_zm, wp4, &    ! InOut
               wpup2, wpvp2, wp2up2, wp2vp2, ice_supersat_frac, &   ! InOut
               um_pert, vm_pert, upwp_pert, vpwp_pert, &            ! InOut
               pdf_params, pdf_params_zm, &                         ! InOut
               pdf_implicit_coefs_terms, &                          ! InOut
               err_info, &                                          ! InOut
               Kh_zm, Kh_zt, &                                      ! Out
#ifdef CLUBB_CAM
               qclvar, &                                            ! Out
#endif
               thlprcp, wprcp, w_up_in_cloud, w_down_in_cloud, &    ! Out
               cloudy_updraft_frac, cloudy_downdraft_frac, &        ! Out
               rcm_in_layer, cloud_cover, invrs_tau_zm, &           ! Out
               Lscale )                                             ! Out

    ! Description:
    !   Subroutine to advance CLUBB one timestep

    ! References:
    !   https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:overview_clubb
    !
    !   ``A PDF-Based Model for Boundary Layer Clouds. Part I:
    !     Method and Model Description'' Golaz, et al. (2002)
    !   JAS, Vol. 59, pp. 3540--3551.
    !-----------------------------------------------------------------------

    ! Modules to be included

    use constants_clubb, only: &
        em_min, &
        fstderr, &
        three_halves, &
        zero_threshold, &
        zero, &
        two

    use parameter_indices, only: &
        nparams,                 & ! Variable(s)
        ic_K,                    &
        iC_wp2_splat,            &
        iup2_sfc_coef,           &
        ia_const,                &
        ibv_efold,               &
        ilambda0_stability_coef

    use parameters_tunable, only: &
        nu_vertical_res_dep    ! Type(s)

    use model_flags, only: &
        clubb_config_flags_type, & ! Type
        l_advance_xp3, &
        iiPDF_ADG1, &
        order_xm_wpxp, &
        order_xp2_xpyp, &
        order_wp2_wp3, &
        order_windm

    use grid_class, only: &
        grid,       & ! Type
        zm2zt_api,  &
        zt2zm_api   ! Procedure(s)

    use numerical_check, only: &
        parameterization_check ! Procedure(s)

    use pdf_parameter_module, only: &
        pdf_parameter, &
        implicit_coefs_terms

    use advance_xm_wpxp_module, only: &
        advance_xm_wpxp          ! Compute mean/flux terms

    use advance_xp2_xpyp_module, only: &
        advance_xp2_xpyp     ! Computes variance terms

    use sfc_varnce_module, only:  &
        calc_sfc_varnce ! Procedure

    use mixing_length, only: &
        calc_Lscale

    use advance_windm_edsclrm_module, only:  &
        advance_windm_edsclrm  ! Procedure(s)

    use saturation, only:  &
        ! Procedure
        sat_mixrat_liq_api ! Saturation mixing ratio

    use advance_wp2_wp3_module, only:  &
        advance_wp2_wp3 ! Procedure

    use advance_xp3_module, only: &
        advance_xp3, & ! Procedure(s)
        diagnose_xp3

    use clubb_precision, only:  &
        core_rknd ! Variable(s)

    use error_code, only: &
        clubb_at_least_debug_level_api,  & ! Procedure
        clubb_fatal_error              ! Constant

    use clip_explicit, only: &
        clip_covars_denom, & ! Procedure(s)
        clip_rcm

    use calc_pressure, only: &
        calculate_thvm

    use T_in_K_module, only: &
        ! Read values from namelist
        thlm2T_in_K_api ! Procedure

    use stats_clubb_utilities, only: &
        stats_accumulate ! Procedure

    use sigma_sqd_w_module, only: &
        compute_sigma_sqd_w ! Procedure

    use advance_helper_module, only: &
        compute_Cx_fnc_Richardson, &
        calc_brunt_vaisala_freq_sqd, &
        calc_ddzt_umvm_sqd, &
        calc_stability_correction, &
        wp23_term_splat_lhs

    use stats_netcdf, only: &
        stats_type, &
        stats_update, &
        stats_begin_budget, &
        stats_finalize_budget, &
        var_on_stats_list

    use array_index, only: &
      sclr_idx_type

    use err_info_type_module, only: &
      err_info_type        ! Type

    use pdf_closure_module, only: &
      pdf_closure_driver    ! Procedure

    implicit none

    !!! External
    intrinsic :: sqrt, min, max, exp, mod, real

    ! Constant Parameters

    real( kind = core_rknd ), parameter :: &
      tau_const = 1000._core_rknd

    !--------------------------- Input Variables ---------------------------
    integer, intent(in) :: &
      nzm, &   ! Number of momentum vertical levels
      nzt, &   ! Number of thermodynamic vertical levels
      ngrdcol  ! Number of grid columns

    type (grid), intent(in) :: gr

    logical, intent(in) ::  &
      l_implemented    ! True if CLUBB is being implemented and run in a host model

    real( kind = core_rknd ), intent(in) ::  &
      dt  ! Current timestep duration    [s]

    real( kind = core_rknd ) ::  &
      dt_advance  ! General timestep duration for advance_wp2_wp3,
                  ! advance_xm_xpwp, and advance_xp2_xpyp.
                  ! Only differs from dt if l_lmm_stepping is used    [s]

    real( kind = core_rknd ), intent(in), dimension(ngrdcol) ::  &
      fcor,  &          ! Traditional Coriolis parameter    [s^-1]
                        ! Vertical planetary vorticity.   Proportional to sin(latitude)
      fcor_y, &         ! Nontraditional Coriolis parameter [s^-1]
                        ! Meridional planetary vorticity. Proportional to cos(latitude)
      sfc_elevation     ! Elevation of ground level    [m above MSL]

    integer, intent(in) :: &
      hydromet_dim,   & ! Total number of hydrometeor species       [#]
      sclr_dim,       & ! Number of passive scalars                 [#]
      edsclr_dim        ! Number of eddy-diff. passive scalars      [#]

    real( kind = core_rknd ), intent(in), dimension(sclr_dim) :: &
      sclr_tol          ! Threshold(s) on the passive scalars  [units vary]

    type (sclr_idx_type), intent(in) :: &
      sclr_idx

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzt) ::  &
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

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzm) ::  &
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

    logical, dimension(hydromet_dim), intent(in) :: &
      l_mix_rat_hm   ! if true, then the quantity is a hydrometeor mixing ratio

    real( kind = core_rknd ), dimension(ngrdcol,nzm,hydromet_dim), intent(in) :: &
      wphydrometp    ! Covariance of w and a hydrometeor      [(m/s) <hm units>]

    real( kind = core_rknd ), dimension(ngrdcol,nzt,hydromet_dim), intent(in) :: &
      wp2hmp,      & ! Third-order moment:  < w'^2 hm' > (hm = hydrometeor) [(m/s)^2 <hm units>]
      rtphmp_zt,   & ! Covariance of rt and hm (on thermo levs.) [(kg/kg) <hm units>]
      thlphmp_zt     ! Covariance of thl and hm (on thermo levs.)      [K <hm units>]

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
      wpsclrp_sfc      ! Passive scalar flux at surface         [{units vary} m/s]

    ! Eddy passive scalar variables
    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzt,edsclr_dim) :: &
      edsclrm_forcing  ! Eddy-diffusion passive scalar forcing    [{units vary}/s]

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,edsclr_dim) ::  &
      wpedsclrp_sfc    ! Eddy-diffusion passive scalar flux at surface    [{units vary} m/s

    real( kind = core_rknd ), intent(in), dimension(ngrdcol) :: &
      upwp_sfc_pert, & ! pertubed u'w' at surface    [m^2/s^2]
      vpwp_sfc_pert    ! pertubed v'w' at surface    [m^2/s^2]

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
      host_dx,  & ! East-west horizontal grid spacing     [m]
      host_dy     ! North-south horizontal grid spacing   [m]

    real( kind = core_rknd ), dimension(ngrdcol,nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

    type(nu_vertical_res_dep), intent(in) :: &
      nu_vert_res_dep    ! Vertical resolution dependent nu values

    real( kind = core_rknd ), intent(in) :: &
      lmin, &                 ! Min. value for the length scale    [m]
      mixt_frac_max_mag, &
      T0, &                   ! Reference temperature (usually 300)  [K]
      ts_nudge, &             ! Timescale of u/v nudging             [s]
      rtm_min, &              ! Value below which rtm will be nudged [kg/kg]
      rtm_nudge_max_altitude  ! Highest altitude at which to nudge rtm [m]

    type( clubb_config_flags_type ), intent(in) :: &
      clubb_config_flags ! Derived type holding all configurable CLUBB flags

    !--------------------------- Input/Output Variables ---------------------------
    type(stats_type), intent(inout) :: &
      stats

    ! These are prognostic or are planned to be in the future
    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,nzt) ::  &
      um,      & ! eastward grid-mean wind component (thermodynamic levels)   [m/s]
      vm,      & ! northward grid-mean wind component (thermodynamic levels)   [m/s]
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
      wpthlp,  & ! w'th_l' (momentum levels)                      [(m/s) K]
      rtp2,    & ! r_t'^2 (momentum levels)                       [(kg/kg)^2]
      thlp2,   & ! th_l'^2 (momentum levels)                      [K^2]
      rtpthlp, & ! r_t'th_l' (momentum levels)                    [(kg/kg) K]
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
      pdf_params,    & ! Fortran structure of PDF parameters on thermodynamic levels    [units vary]
      pdf_params_zm    ! Fortran structure of PDF parameters on momentum levels        [units vary]

    type(implicit_coefs_terms), intent(inout) :: &
      pdf_implicit_coefs_terms    ! Implicit coefs / explicit terms [units vary]

    type(err_info_type), intent(inout) :: &
      err_info        ! err_info struct containing err_code and err_header

    ! Eddy passive scalar variable
    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,nzt,edsclr_dim) :: &
      edsclrm   ! Eddy passive scalar grid-mean (thermo. levels)   [units vary]

    ! Variables that need to be output for use in other parts of the CLUBB
    ! code, such as microphysics (rcm, pdf_params), forcings (rcm), and/or
    ! BUGSrad (cloud_cover).
    real( kind = core_rknd ), intent(out), dimension(ngrdcol,nzt) ::  &
      rcm_in_layer, & ! rcm within cloud layer                          [kg/kg]
      cloud_cover     ! cloud cover                                     [-]

    ! Variables that need to be output for use in host models
    real( kind = core_rknd ), intent(out), dimension(ngrdcol,nzt) ::  &
      w_up_in_cloud,         & ! Average cloudy updraft velocity       [m/s]
      w_down_in_cloud,       & ! Average cloudy downdraft velocity     [m/s]
      cloudy_updraft_frac,   & ! cloudy updraft fraction               [-]
      cloudy_downdraft_frac    ! cloudy downdraft fraction             [-]

    real( kind = core_rknd ), intent(out), dimension(ngrdcol,nzm) ::  &
      wprcp,                 & ! w'r_c' (momentum levels)              [(kg/kg) m/s]
      invrs_tau_zm             ! One divided by tau on zm levels       [1/s]

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
      Lscale          ! Length scale                          [m]

    !--------------------------- Local Variables ---------------------------
    integer :: i, k

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      thvm            ! Virtual potential temperature                    [K]

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      rsat   ! Saturation mixing ratio  ! Brian

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      sigma_sqd_w_zt   ! PDF width parameter on thermodynamic levels [-]

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      wpthlp2,   & ! w'thl'^2    [m K^2/s]
      wprtp2,    & ! w'rt'^2     [m kg^2/kg^2]
      wprtpthlp    ! w'rt'thl'   [m kg K/kg s]

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      Lscale_up,   & ! Length scale (upwards component)      [m]
      Lscale_down    ! Length scale (downwards component)    [m]

    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
      Lscale_zm      ! Length scale on momentum levels       [m]

    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
      em,     & ! Turbulent Kinetic Energy (TKE)                      [m^2/s^2]
      tau_zm    ! Eddy dissipation time scale on momentum levels      [s]

    real( kind = core_rknd ), dimension(ngrdcol,nzm,edsclr_dim) :: &
      wpedsclrp   ! w'edsclr'

    real( kind = core_rknd ), dimension(ngrdcol,nzt,sclr_dim) :: &
      wp2sclrp,    & ! w'^2 sclr'
      wpsclrp2,    & ! w'sclr'^2
      wpsclrprtp,  & ! w'sclr'rt'
      wpsclrpthlp    ! w'sclr'thl'

    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
      sigma_sqd_w        ! PDF width parameter (momentum levels)         [-]

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      em_zt       ! TKE interpolated to thermodynamic levels [m^2/s^2]

    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
      w_1_zm,        & ! Mean w (1st PDF component)                   [m/s]
      w_2_zm,        & ! Mean w (2nd PDF component)                   [m/s]
      varnce_w_1_zm, & ! Variance of w (1st PDF component)            [m^2/s^2]
      varnce_w_2_zm, & ! Variance of w (2nd PDF component)            [m^2/s^2]
      mixt_frac_zm     ! Weight of 1st PDF component (Sk_w dependent) [-]

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      rcp2_zt              ! r_c'^2 (on thermo. grid)             [kg^2/kg^2]

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      rtm_before, &
      thlm_before, &
      rel_humidity        ! Relative humidity after PDF closure [-]

    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
       stability_correction,         & ! Stability correction factor
       invrs_tau_C6_zm,              & ! Inverse tau values used for C6 (pr1) term in wpxp [1/s]
       invrs_tau_C1_zm,              & ! Inverse tau values used for C1 (dp1) term in wp2 [1/s]
       invrs_tau_xp2_zm,             & ! Inverse tau values used for advance_xp2_wpxp [s^-1]
       invrs_tau_C4_zm,              & ! Inverse tau values used for C4 terms         [s^-1]
       invrs_tau_C14_zm,             & ! Inverse tau valuse used for C14 terms        [s^-1]
       Cx_fnc_Richardson,            & ! Cx_fnc computed from Richardson_num          [-]
       ddzt_umvm_sqd,                & ! Squared vertical norm of derivative of
                                       ! mean horizontal wind speed                   [s^-2]
       brunt_vaisala_freq_sqd,       & ! Buoyancy frequency squared, N^2              [s^-2]
       brunt_vaisala_freq_sqd_mixed, & ! A mixture of dry and moist N^2               [s^-2]
       brunt_vaisala_freq_sqd_smth     ! Mix between dry and moist N^2 that is
                                       ! smoothed in the vertical                     [s^-2]

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
       invrs_tau_zt,                 & ! Inverse time-scale tau on thermodynamics levels [1/s]
       invrs_tau_wp3_zt                ! Inverse tau wp3 at zt levels

    real( kind = core_rknd ), parameter :: &
       ufmin = 0.01_core_rknd           ! minimum value of friction velocity     [m/s]

    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
      tau_max_zm    ! Max. allowable eddy dissipation time scale on m-levs  [s]

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      tau_max_zt    ! Max. allowable eddy dissipation time scale on t-levs  [s]

    ! Flag to sample stats in a particular call to subroutine
    ! pdf_closure_driver.
    logical :: l_samp_stats_in_pdf_call

    ! Flag to determine whether invrs_tau_N2_iso is used in C4 terms.
    ! Important! This flag is only in use when l_diag_Lscale_from_tau = true
    ! Setting l_use_invrs_tau_N2_iso = true will not change anything unless
    ! l_diag_Lscale_from_tau is also true
    logical, parameter :: l_use_invrs_tau_N2_iso = .false.

    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
       lhs_splat_wp2    ! LHS coefficient of wp2 splatting term  [1/s]

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
       lhs_splat_wp3    ! LHS coefficient of wp3 splatting term  [1/s]

    ! Variables associated with upgradient momentum contributions due to cumuli
    !real( kind = core_rknd ), dimension(nz) :: &
    !  Km_Skw_factor ! Factor, with value < 1, that reduces eddy diffusivity,
    !                                          Km_zm, in skewed layers
    !real( kind = core_rknd ),parameter :: &
    !  Km_Skw_thresh = zero_threshold, &  ! Value of Skw at which Skw correction kicks in
    !  Km_Skw_factor_efold = 0.5_core_rknd, & ! E-folding rate of exponential Skw correction
    !  Km_Skw_factor_min   = 0.2_core_rknd    ! Minimum value of Km_Skw_factor

    integer :: &
      advance_order_loop_iter, &
      wprtp_cl_num, &
      wpthlp_cl_num, &
      upwp_cl_num, &
      vpwp_cl_num

    !----- Begin Code -----

    !$acc enter data create( thvm, &
    !$acc              ddzt_umvm_sqd, rtm_before, thlm_before, &
    !$acc              wpthlp2, wprtp2, wprtpthlp, Lscale_up, Lscale_zm, &
    !$acc              Lscale_down, em, tau_zm, &
    !$acc              sigma_sqd_w, &
    !$acc              sigma_sqd_w_zt, em_zt, w_1_zm, w_2_zm, &
    !$acc              varnce_w_1_zm, varnce_w_2_zm, &
    !$acc              mixt_frac_zm, rcp2_zt, stability_correction, &
    !$acc              invrs_tau_C6_zm, invrs_tau_C1_zm, invrs_tau_xp2_zm, &
    !$acc              invrs_tau_C4_zm, invrs_tau_C14_zm, &
    !$acc              invrs_tau_zt, invrs_tau_wp3_zt, Cx_fnc_Richardson, &
    !$acc              brunt_vaisala_freq_sqd, brunt_vaisala_freq_sqd_mixed, &
    !$acc              brunt_vaisala_freq_sqd_smth, &
    !$acc              tau_max_zm, tau_max_zt, lhs_splat_wp2, lhs_splat_wp3 )

    !$acc enter data if( sclr_dim > 0 ) &
    !$acc            create( wp2sclrp, &
    !$acc                    wpsclrp2, wpsclrprtp, wpsclrpthlp )

    !$acc enter data if( edsclr_dim > 0 ) &
    !$acc            create( wpedsclrp )

    if ( clubb_config_flags%l_lmm_stepping ) then
      dt_advance = two * dt
    else
      dt_advance = dt
    end if

    !----------------------------------------------------------------
    ! Test input variables
    !----------------------------------------------------------------
    if ( clubb_at_least_debug_level_api( 2 ) ) then

      call parameterization_check(                                        &
             nzm, nzt, ngrdcol, sclr_dim, edsclr_dim,                     & ! intent(in)
             thlm_forcing, rtm_forcing, um_forcing,                       & ! intent(in)
             vm_forcing, wm_zm, wm_zt, p_in_Pa,                           & ! intent(in)
             rho_zm, rho, exner, rho_ds_zm,                               & ! intent(in)
             rho_ds_zt, invrs_rho_ds_zm, invrs_rho_ds_zt,                 & ! intent(in)
             thv_ds_zm, thv_ds_zt, wpthlp_sfc, wprtp_sfc, upwp_sfc,       & ! intent(in)
             vpwp_sfc, p_sfc, um, upwp, vm, vpwp, up2, vp2,               & ! intent(in)
             rtm, wprtp, thlm, wpthlp, wp2, wp3,                          & ! intent(in)
             rtp2, thlp2, rtpthlp,                                        & ! intent(in)
             "beginning of ",                                             & ! intent(in)
             wpsclrp_sfc, wpedsclrp_sfc, sclrm, wpsclrp,                  & ! intent(in)
             sclrp2,                                                      & ! intent(in)
             sclrprtp, sclrpthlp, sclrm_forcing, edsclrm,                 & ! intent(in)
             edsclrm_forcing,                                             & ! intent(in)
             err_info )                                                     ! intent(inout)

      if ( any(err_info%err_code == clubb_fatal_error) ) then
        write(fstderr,*) "Fatal error detected in parameterization_check when testing input"
        return
      end if

    end if
    !-----------------------------------------------------------------------

    if ( stats%l_sample ) then

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzt
        do i = 1, ngrdcol
          rtm_before(i,k) = rtm(i,k)
          thlm_before(i,k) = thlm(i,k)
        end do
      end do
      !$acc end parallel loop

      !$acc update host( rfrzm, wp2, vp2, up2, wprtp, wpthlp, upwp, vpwp, &
      !$acc              rtp2, thlp2, rtpthlp, rtm, thlm, um, vm, wp3 )

      call stats_update( "rfrzm", rfrzm, stats )

      ! Set up budget stats variables.
      call stats_begin_budget( "wp2_bt", wp2 / dt, stats )
      call stats_begin_budget( "vp2_bt", vp2 / dt, stats )
      call stats_begin_budget( "up2_bt", up2 / dt, stats )
      call stats_begin_budget( "wprtp_bt", wprtp / dt, stats )
      call stats_begin_budget( "wpthlp_bt", wpthlp / dt, stats )

      if ( clubb_config_flags%l_predict_upwp_vpwp ) then
        call stats_begin_budget( "upwp_bt", upwp / dt, stats )
        call stats_begin_budget( "vpwp_bt", vpwp / dt, stats )
      end if

      call stats_begin_budget( "rtp2_bt", rtp2 / dt, stats )
      call stats_begin_budget( "thlp2_bt", thlp2 / dt, stats )
      call stats_begin_budget( "rtpthlp_bt", rtpthlp / dt, stats )

      call stats_begin_budget( "rtm_bt", rtm / dt, stats )
      call stats_begin_budget( "thlm_bt", thlm / dt, stats )
      call stats_begin_budget( "um_bt", um / dt, stats )
      call stats_begin_budget( "vm_bt", vm / dt, stats )
      call stats_begin_budget( "wp3_bt", wp3 / dt, stats )

      ! Set clipping counts to zero
      wprtp_cl_num = 0
      wpthlp_cl_num = 0
      upwp_cl_num = 0
      vpwp_cl_num = 0

    end if

    call set_sfc_value_of_flux_profiles( nzm, ngrdcol, sclr_dim, edsclr_dim, gr,  & ! In
                                          clubb_config_flags%l_host_applies_sfc_fluxes, & ! In
                                          clubb_config_flags%l_linearize_pbl_winds, & ! In
                                          wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, & ! In
                                          upwp_sfc_pert, vpwp_sfc_pert,            & ! In
                                          wpsclrp_sfc, wpedsclrp_sfc,              & ! In
                                          wpthlp, wprtp, upwp, vpwp,               & ! InOut
                                          upwp_pert, vpwp_pert,                    & ! InOut
                                          wpsclrp, wpedsclrp )                       ! InOut

    call compute_sigma_sqd_w( nzm, nzt, ngrdcol, gr,                        & ! In
                              wp3, wp2, thlp2, rtp2,                        & ! In
                              up2, vp2, wpthlp, wprtp, upwp, vpwp,          & ! In
                              clubb_params,                                 & ! In
                              clubb_config_flags%l_predict_upwp_vpwp,       & ! In
                              sigma_sqd_w )                                   ! Out

    if (      clubb_config_flags%ipdf_call_placement == ipdf_pre_advance_fields &
         .or. clubb_config_flags%ipdf_call_placement == ipdf_pre_post_advance_fields ) then

      ! Sample stats in this call to subroutine pdf_closure_driver for
      ! both of these options (ipdf_pre_advance_fields and
      ! ipdf_pre_post_advance_fields).
      if ( clubb_config_flags%ipdf_call_placement == ipdf_pre_advance_fields ) then
        l_samp_stats_in_pdf_call = .true.
      elseif ( clubb_config_flags%ipdf_call_placement == ipdf_pre_post_advance_fields ) then
        l_samp_stats_in_pdf_call = .true.
      end if

      !########################################################################
      !#######                     CALL CLUBB's PDF                     #######
      !#######   AND OUTPUT PDF PARAMETERS AND INTEGRATED QUANTITITES   #######
      !########################################################################
      call pdf_closure_driver( gr, nzm, nzt, ngrdcol,                       & ! In
                               dt, hydromet_dim, sclr_dim, sclr_tol,        & ! In
                               wprtp, thlm, wpthlp, rtp2, rtp3,             & ! In
                               thlp2, thlp3, rtpthlp, wp2,                  & ! In
                               wp3, wm_zm, wm_zt,                           & ! In
                               um, up2, upwp, up3,                          & ! In
                               vm, vp2, vpwp, vp3,                          & ! In
                               p_in_Pa, exner,                              & ! In
                               thv_ds_zm, thv_ds_zt, rtm_ref,               & ! In
                               wphydrometp,                                 & ! In
                               wp2hmp, rtphmp_zt, thlphmp_zt,               & ! In
                               sclrm, wpsclrp, sclrp2,                      & ! In
                               sclrprtp, sclrpthlp, sclrp3,                 & ! In
                               p_sfc, l_samp_stats_in_pdf_call,             & ! In
                               mixt_frac_max_mag, ts_nudge,                 & ! In
                               rtm_min, rtm_nudge_max_altitude,             & ! In
                               clubb_params,                                & ! In
                               clubb_config_flags%iiPDF_type,               & ! In
                               clubb_config_flags%saturation_formula,       & ! In
                               clubb_config_flags%l_predict_upwp_vpwp,      & ! In
                               clubb_config_flags%l_rtm_nudge,              & ! In
                               clubb_config_flags%l_trapezoidal_rule_zt,    & ! In
                               clubb_config_flags%l_trapezoidal_rule_zm,    & ! In
                               clubb_config_flags%l_call_pdf_closure_twice, & ! In
                               clubb_config_flags%l_use_cloud_cover,        & ! In
                               clubb_config_flags%l_rcm_supersat_adj,       & ! In
                               l_mix_rat_hm,                                & ! In
                               stats,                                       & ! InOut
                               rtm, sigma_sqd_w,                            & ! InOut
                               pdf_implicit_coefs_terms,                    & ! InOut
                               pdf_params, pdf_params_zm, err_info,         & ! InOut
                               rcm, cloud_frac,                             & ! Out
                               ice_supersat_frac, wprcp,                    & ! Out
                               wpthvp, wp2thvp, wp2up,                      & ! Out
                               rtpthvp, thlpthvp,                           & ! Out
                               rcm_in_layer, cloud_cover,                   & ! Out
                               rcp2_zt, thlprcp,                            & ! Out
                               rc_coef_zm, sclrpthvp,                       & ! Out
                               wpup2, wpvp2,                                & ! Out
                               wp2up2, wp2vp2, wp4,                         & ! Out
                               wp2rtp, wprtp2, wp2thlp,                     & ! Out
                               wpthlp2, wprtpthlp,                          & ! Out
                               uprcp, vprcp,                                & ! Out
                               w_up_in_cloud, w_down_in_cloud,              & ! Out
                               cloudy_updraft_frac,                         & ! Out
                               cloudy_downdraft_frac,                       & ! Out
                               wp2sclrp, wpsclrp2,                          & ! Out
                               wpsclrprtp, wpsclrpthlp )                      ! Out

      if ( clubb_at_least_debug_level_api( 0 ) ) then
        if ( any(err_info%err_code == clubb_fatal_error) ) then
          write(fstderr,*) err_info%err_header_global
          write(fstderr,*) "Error calling pdf_closure_driver in advance_clubb_core"
          return
        end if
      end if

    end if

    ! This feels like an awkward place to have this, but this keeps things BFB
    ! with where it was before
    if ( stats%l_sample ) then

      if ( clubb_config_flags%iiPDF_type == iiPDF_ADG1 ) then

        sigma_sqd_w_zt(:,:) = zm2zt_api( nzm, nzt, ngrdcol, gr, sigma_sqd_w(:,:), &
                                        zero_threshold )

      else

        !$acc parallel loop gang vector collapse(2) default(present)
        do k = 1, nzt
          do i = 1, ngrdcol
            sigma_sqd_w_zt(i,k) = zero
          end do
        end do

      end if

      !$acc update host( sigma_sqd_w_zt )
      call stats_update( "sigma_sqd_w_zt", sigma_sqd_w_zt, stats )

    end if

    call calculate_thvm( nzt, ngrdcol, &
                         thlm, rtm, rcm, exner, thv_ds_zt, &
                         thvm )

    ! Compute tke (turbulent kinetic energy).
    if ( .not. clubb_config_flags%l_tke_aniso ) then
      ! tke is assumed to be 3/2 of wp2.
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm
        do i = 1, ngrdcol
          em(i,k) = three_halves * wp2(i,k)
        end do
      end do
      !$acc end parallel loop
    else
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm
        do i = 1, ngrdcol
          em(i,k) = 0.5_core_rknd * ( wp2(i,k) + vp2(i,k) + up2(i,k) )
        end do
      end do
      !$acc end parallel loop
    end if

    if ( clubb_config_flags%l_call_pdf_closure_twice ) then
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm
        do i = 1, ngrdcol
          w_1_zm(i,k)        = pdf_params_zm%w_1(i,k)
          w_2_zm(i,k)        = pdf_params_zm%w_2(i,k)
          varnce_w_1_zm(i,k) = pdf_params_zm%varnce_w_1(i,k)
          varnce_w_2_zm(i,k) = pdf_params_zm%varnce_w_2(i,k)
          mixt_frac_zm(i,k)  = pdf_params_zm%mixt_frac(i,k)
        end do
      end do
      !$acc end parallel loop
    else
      w_1_zm(:,:)        = zt2zm_api( nzm, nzt, ngrdcol, gr, pdf_params%w_1(:,:) )
      w_2_zm(:,:)        = zt2zm_api( nzm, nzt, ngrdcol, gr, pdf_params%w_2(:,:) )
      varnce_w_1_zm(:,:) = zt2zm_api( nzm, nzt, ngrdcol, gr, pdf_params%varnce_w_1(:,:) )
      varnce_w_2_zm(:,:) = zt2zm_api( nzm, nzt, ngrdcol, gr, pdf_params%varnce_w_2(:,:) )
      mixt_frac_zm(:,:)  = zt2zm_api( nzm, nzt, ngrdcol, gr, pdf_params%mixt_frac(:,:) )
    end if

    if ( stats%l_sample ) then
      !$acc update host( rtm, rcm, thlm, exner, p_in_Pa )
      call stats_update( "rvm", rtm - rcm, stats )
      ! Output relative humidity (q/q* where q* is the saturation mixing ratio
      ! over liquid).
      ! Added an extra check for rel_humidity or rsat output; otherwise,
      ! if both are omitted, rsat may not be computed, leading to a
      ! floating-point exception when rel_humidity is written.
      if (      var_on_stats_list( stats, "rel_humidity" ) &
           .or. var_on_stats_list( stats, "rsat" )         ) then
        rsat = sat_mixrat_liq_api( nzt, ngrdcol, gr, p_in_Pa, &
                                   thlm2T_in_K_api( nzt, ngrdcol, thlm, exner, rcm ), &
                                   clubb_config_flags%saturation_formula )

        ! Recompute rsat and rel_humidity. They might have changed.
        do i = 1, ngrdcol
          rel_humidity(i,:) = (rtm(i,:) - rcm(i,:)) / rsat(i,:)
        end do
        ! Write the var to stats
        call stats_update( "rel_humidity", rel_humidity, stats )
      end if
    end if

    call calc_brunt_vaisala_freq_sqd( nzm, nzt, ngrdcol, gr, thlm,                        & ! In
                                      exner, rtm, rcm, p_in_Pa, thvm,                     & ! In
                                      ice_supersat_frac,                                  & ! In
                                      clubb_config_flags%saturation_formula,              & ! In
                                      clubb_config_flags%l_brunt_vaisala_freq_moist,      & ! In
                                      clubb_config_flags%l_use_thvm_in_bv_freq,           & ! In
                                      clubb_config_flags%l_modify_limiters_for_cnvg_test, & ! In
                                      clubb_params(:,ibv_efold), T0,                      & ! In
                                      brunt_vaisala_freq_sqd,                             & ! Out
                                      brunt_vaisala_freq_sqd_mixed,                       & ! Out
                                      brunt_vaisala_freq_sqd_smth,                        & ! Out
                                      stats )                                    ! Optional InOut

    !----------------------------------------------------------------
    ! Compute mixing length and dissipation time
    !----------------------------------------------------------------

    call calc_ddzt_umvm_sqd( nzm, nzt, ngrdcol, gr, um, vm, & ! In
                             ddzt_umvm_sqd )                  ! Out

    call calc_Lscale( nzm, nzt, ngrdcol, gr, l_implemented, host_dx, host_dy,       & ! In
                      p_in_Pa, exner, rtm, thlm, thvm,                              & ! In
                      thlp2, rtp2, rtpthlp,                                         & ! In
                      pdf_params, em, thv_ds_zt, lmin,                              & ! In
                      upwp_sfc, vpwp_sfc, ddzt_umvm_sqd, ice_supersat_frac,         & ! In
                      ufmin, tau_const, sfc_elevation, clubb_params,                & ! In
                      clubb_config_flags%saturation_formula,                        & ! In
                      clubb_config_flags%l_Lscale_plume_centered,                   & ! In
                      clubb_config_flags%l_diag_Lscale_from_tau,                    & ! In
                      clubb_config_flags%l_e3sm_config,                             & ! In
                      clubb_config_flags%l_smooth_Heaviside_tau_wpxp,               & ! In
                      clubb_config_flags%l_modify_limiters_for_cnvg_test,           & ! In
                      l_use_invrs_tau_N2_iso,                                       & ! In
                      brunt_vaisala_freq_sqd_smth,                                  & ! In
                      stats, err_info,                                              & ! InOut
                      invrs_tau_zt, invrs_tau_zm, invrs_tau_xp2_zm,                 & ! Out
                      invrs_tau_wp3_zt,                                             & ! Out
                      invrs_tau_C1_zm, invrs_tau_C4_zm,                             & ! Out
                      invrs_tau_C6_zm, invrs_tau_C14_zm,                            & ! Out
                      tau_max_zm, tau_max_zt, tau_zm,                               & ! Out
                      Lscale, Lscale_zm, Lscale_up, Lscale_down )                     ! Out

    if ( clubb_at_least_debug_level_api( 0 ) ) then
      if ( any(err_info%err_code == clubb_fatal_error) ) then
        write(fstderr,*) err_info%err_header_global
        write(fstderr, *) "Error calling calc_Lscale in advance_clubb_core"
        return
      end if
    end if

    ! Here we determine if we're using tau_zm or tau_N2_zm, which is tau
    ! that has been stability corrected for stably stratified regions.
    ! -dschanen 7 Nov 2014
    if ( clubb_config_flags%l_stability_correct_tau_zm ) then

      call calc_stability_correction( nzm, ngrdcol, brunt_vaisala_freq_sqd,    & ! In
                                      Lscale_zm, em,                           & ! In
                                      clubb_params(:,ilambda0_stability_coef), & ! In
                                      stability_correction )                     ! Out

      ! Create a damping time scale that is more strongly damped at the
      ! altitudes where the Brunt-Vaisala frequency (N^2) is large.
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm
        do i = 1, ngrdcol
          invrs_tau_C6_zm(i,k) = invrs_tau_zm(i,k) * stability_correction(i,k)
          invrs_tau_C1_zm(i,k) = invrs_tau_C6_zm(i,k)
        end do
      end do
      !$acc end parallel loop

    end if ! clubb_config_flags%l_stability_correct_tau_zm

    !----------------------------------------------------------------
    ! Eddy diffusivity coefficient
    !----------------------------------------------------------------
    ! c_K is 0.548 usually (Duynkerke and Driedonks 1987)
    ! CLUBB uses a smaller value to better fit empirical data.

    ! Calculate CLUBB's eddy diffusivity as
    !   CLUBB's length scale times a velocity scale.
    em_zt(:,:) = zm2zt_api( nzm, nzt, ngrdcol, gr, em(:,:), em_min )

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzt
      do i = 1, ngrdcol
        Kh_zt(i,k) = clubb_params(i,ic_K) * Lscale(i,k) * sqrt( em_zt(i,k) )
      end do
    end do
    !$acc end parallel loop

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzm
      do i = 1, ngrdcol
        Kh_zm(i,k) = clubb_params(i,ic_K) * max( Lscale_zm(i,k), zero_threshold )  &
                     * sqrt( max( em(i,k), em_min ) )
      end do
    end do
    !$acc end parallel loop

    call wp23_term_splat_lhs( nzm, nzt, ngrdcol, gr, clubb_params(:,iC_wp2_splat), & ! In
                              brunt_vaisala_freq_sqd_mixed, Lscale_zm, rho_ds_zm,  & ! In
                              stats,                                               & ! InOut
                              lhs_splat_wp2, lhs_splat_wp3 )                         ! Out

    !----------------------------------------------------------------
    ! Set Surface variances
    !----------------------------------------------------------------
    ! Surface variances should be set here, before the call to either
    ! advance_xp2_xpyp or advance_wp2_wp3.
    ! Surface effects should not be included with any case where the lowest
    ! level is not the ground level.  Brian Griffin.  December 22, 2005.

    ! Diagnose surface variances based on surface fluxes.
    call calc_sfc_varnce( nzm, nzt, ngrdcol, sclr_dim, sclr_idx,          & ! In
                          gr, dt, sfc_elevation,                          & ! In
                          upwp_sfc, vpwp_sfc, wpthlp, wprtp_sfc,          & ! In
                          um, vm, Lscale_up, wpsclrp_sfc,                 & ! In
                          lhs_splat_wp2, tau_zm,                          & ! In
                          clubb_config_flags%l_vary_convect_depth,        & ! In
                          T0,                                             & ! In
                          clubb_params(:,iup2_sfc_coef),                  & ! In
                          clubb_params(:,ia_const),                       & ! In
                          stats,                                          & ! InOut
                          wp2, up2, vp2,                                  & ! InOut
                          thlp2, rtp2, rtpthlp,                           & ! InOut
                          sclrp2, sclrprtp, sclrpthlp,                    & ! InOut
                          err_info )                                        ! InOut

    if ( clubb_at_least_debug_level_api( 0 ) ) then
      if ( any(err_info%err_code == clubb_fatal_error) ) then
        write(fstderr,*) err_info%err_header_global
        write(fstderr, *) "Error calling calc_sfc_varnce in advance_clubb_core"
        return
      end if
    end if

    !#######################################################################
    !############## ADVANCE PROGNOSTIC VARIABLES ONE TIMESTEP ##############
    !#######################################################################

    ! Cx_fnc_Richardson is only used if one of these flags is true,
    ! otherwise its value is irrelevant, set it to 0 to avoid NaN problems
    if ( clubb_config_flags%l_use_C7_Richardson .or. &
         clubb_config_flags%l_use_C11_Richardson ) then

      call compute_Cx_fnc_Richardson( nzm, nzt, ngrdcol, gr,                              & ! In
                                      Lscale_zm, ddzt_umvm_sqd, rho_ds_zm,                & ! In
                                      brunt_vaisala_freq_sqd,                             & ! In
                                      brunt_vaisala_freq_sqd_mixed,                       & ! In
                                      clubb_params,                                       & ! In
                                      clubb_config_flags%l_use_shear_Richardson,          & ! In
                                      clubb_config_flags%l_modify_limiters_for_cnvg_test, & ! In
                                      Cx_fnc_Richardson )                                   ! Out
    else
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm
        do i = 1, ngrdcol
          Cx_fnc_Richardson(i,k) = 0.0
        end do
      end do
      !$acc end parallel loop
    end if

    ! Loop over the 4 main advance subroutines -- advance_xm_wpxp,
    ! advance_wp2_wp3, advance_xp2_xpyp, and advance_windm_edsclrm -- in the
    ! order determined by order_xm_wpxp, order_wp2_wp3, order_xp2_xpyp, and
    ! order_windm.
    do advance_order_loop_iter = 1, 4, 1

      if ( advance_order_loop_iter == order_xm_wpxp ) then

        !----------------------------------------------------------------
        ! Advance rtm/wprtp and thlm/wpthlp one time step
        !----------------------------------------------------------------

        ! Advance the prognostic equations for
        !   the scalar grid means (rtm, thlm, sclrm) and
        !   scalar turbulent fluxes (wprtp, wpthlp, and wpsclrp)
        !   by one time step.
        ! advance_xm_wpxp_bad_wp2 ! Test error comment, DO NOT modify or move
        call advance_xm_wpxp( nzm, nzt, ngrdcol, sclr_dim, sclr_tol, gr, dt_advance, & ! In
                              sigma_sqd_w, wm_zm, wm_zt, wp2, Lscale_zm,             & ! In
                              wp3, Kh_zt, Kh_zm,                                     & ! In
                              stability_correction,                                  & ! In
                              invrs_tau_C6_zm, tau_max_zm, wp2rtp, rtpthvp,          & ! In
                              rtm_forcing, wprtp_forcing, rtm_ref, wp2thlp,          & ! In
                              thlpthvp, thlm_forcing, wpthlp_forcing, thlm_ref,      & ! In
                              rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm,                 & ! In
                              invrs_rho_ds_zt, thv_ds_zm, rtp2, thlp2,               & ! In
                              w_1_zm, w_2_zm, varnce_w_1_zm, varnce_w_2_zm,          & ! In
                              mixt_frac_zm, l_implemented, em, wp2sclrp,             & ! In
                              sclrpthvp, sclrm_forcing, sclrp2, Cx_fnc_Richardson,   & ! In
                              pdf_implicit_coefs_terms,                              & ! In
                              um_forcing, vm_forcing, ug, vg, wpthvp,                & ! In
                              fcor, fcor_y, um_ref, vm_ref, up2, vp2,                & ! In
                              uprcp, vprcp, rc_coef_zm,                              & ! In
                              clubb_params, nu_vert_res_dep, ts_nudge,               & ! In
                              clubb_config_flags%iiPDF_type,                         & ! In
                              clubb_config_flags%penta_solve_method,                 & ! In
                              clubb_config_flags%tridiag_solve_method,               & ! In
                              clubb_config_flags%fill_holes_type,                    & ! In
                              clubb_config_flags%l_predict_upwp_vpwp,                & ! In
                              clubb_config_flags%l_ho_nontrad_coriolis,              & ! In
                              clubb_config_flags%l_ho_trad_coriolis,                 & ! In
                              clubb_config_flags%l_diffuse_rtm_and_thlm,             & ! In
                              clubb_config_flags%l_stability_correct_Kh_N2_zm,       & ! In
                              clubb_config_flags%l_godunov_upwind_wpxp_ta,           & ! In
                              clubb_config_flags%l_upwind_xm_ma,                     & ! In
                              clubb_config_flags%l_uv_nudge,                         & ! In
                              clubb_config_flags%l_tke_aniso,                        & ! In
                              clubb_config_flags%l_diag_Lscale_from_tau,             & ! In
                              clubb_config_flags%l_use_C7_Richardson,                & ! In
                              clubb_config_flags%l_lmm_stepping,                     & ! In
                              clubb_config_flags%l_enable_relaxed_clipping,          & ! In
                              clubb_config_flags%l_linearize_pbl_winds,              & ! In
                              clubb_config_flags%l_mono_flux_lim_thlm,               & ! In
                              clubb_config_flags%l_mono_flux_lim_rtm,                & ! In
                              clubb_config_flags%l_mono_flux_lim_um,                 & ! In
                              clubb_config_flags%l_mono_flux_lim_vm,                 & ! In
                              clubb_config_flags%l_mono_flux_lim_spikefix,           & ! In
                              wprtp_cl_num, wpthlp_cl_num, upwp_cl_num, vpwp_cl_num, & ! InOut
                              stats,                                                 & ! InOut
                              rtm, wprtp, thlm, wpthlp,                              & ! intent(i/o)
                              sclrm, wpsclrp, um, upwp, vm, vpwp,                    & ! intent(i/o)
                              um_pert, vm_pert, upwp_pert, vpwp_pert, err_info )       ! intent(i/o)

        if ( clubb_at_least_debug_level_api( 0 ) ) then
          if ( any(err_info%err_code == clubb_fatal_error) ) then
              write(fstderr,*) err_info%err_header_global
              write(fstderr,*) "Error calling advance_xm_wpxp in advance_clubb_core"
              return
          end if
        end if

        ! Vince Larson clipped rcm in order to prevent rvm < 0.  5 Apr 2008.
        ! This code won't work unless rtm >= 0 !!!
        ! We do not clip rcm_in_layer because rcm_in_layer only influences
        ! radiation, and we do not want to bother recomputing it.  6 Aug 2009
        call clip_rcm( nzt, ngrdcol, rtm,              & ! In
                      'rtm < rcm in advance_xm_wpxp', & ! In
                      rcm )                             ! InOut

      elseif ( advance_order_loop_iter == order_xp2_xpyp ) then

        !----------------------------------------------------------------
        ! Compute some of the variances and covariances.  These include the
        ! variance of total water (rtp2), liquid water potential temperature
        ! (thlp2), their covariance (rtpthlp), and the variance of horizontal
        ! wind (up2 and vp2).  The variance of vertical velocity is computed
        ! in a different section, which will come either earlier or later
        ! depending on the chosen call order.
        !----------------------------------------------------------------

        ! We found that certain cases require a time tendency to run
        ! at shorter timesteps so these are prognosed now.

        ! We found that if we call advance_xp2_xpyp first, we can use a longer timestep.

        ! Advance the prognostic equations
        !   for scalar variances and covariances,
        !   plus the horizontal wind variances by one time step, by one time step.
        call advance_xp2_xpyp( nzm, nzt, ngrdcol, sclr_dim, sclr_tol, gr, sclr_idx, & ! In
                              invrs_tau_xp2_zm, invrs_tau_C4_zm,                   & ! In
                              invrs_tau_C14_zm, wm_zm,                             & ! In
                              rtm, wprtp, thlm, wpthlp, wpthvp, um, vm,            & ! In
                              wp2, wp3, upwp, vpwp,                                & ! In
                              sigma_sqd_w, wprtp2, wpthlp2,                        & ! In
                              wprtpthlp, Kh_zt, rtp2_forcing,                      & ! In
                              thlp2_forcing, rtpthlp_forcing,                      & ! In
                              rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm,               & ! In
                              thv_ds_zm, cloud_frac,                               & ! In
                              pdf_implicit_coefs_terms,                            & ! In
                              dt_advance, fcor_y,                                  & ! In
                              sclrm, wpsclrp,                                      & ! In
                              wpsclrp2, wpsclrprtp, wpsclrpthlp,                   & ! In
                              lhs_splat_wp2,                                       & ! In
                              clubb_params, nu_vert_res_dep,                       & ! In
                              clubb_config_flags%iiPDF_type,                       & ! In
                              clubb_config_flags%tridiag_solve_method,             & ! In
                              clubb_config_flags%fill_holes_type,                  & ! In
                              clubb_config_flags%l_ho_nontrad_coriolis,            & ! In
                              clubb_config_flags%l_min_xp2_from_corr_wx,           & ! In
                              clubb_config_flags%l_C2_cloud_frac,                  & ! In
                              clubb_config_flags%l_upwind_xpyp_ta,                 & ! In
                              clubb_config_flags%l_godunov_upwind_xpyp_ta,         & ! In
                              clubb_config_flags%l_lmm_stepping,                   & ! In
                              l_implemented,                                       & ! In
                              stats,                                               & ! InOut
                              rtp2, thlp2, rtpthlp, up2, vp2,                      & ! InOut
                              sclrp2, sclrprtp, sclrpthlp, err_info )                ! InOut

        if ( clubb_at_least_debug_level_api( 0 ) ) then
          if ( any(err_info%err_code == clubb_fatal_error) ) then
              write(fstderr,*) err_info%err_header_global
              write(fstderr,*) "Error calling advance_xp2_xpyp in advance_clubb_core"
              return
          end if
        end if

        !----------------------------------------------------------------
        ! Covariance clipping for wprtp, wpthlp, wpsclrp, upwp, and vpwp
        ! after subroutine advance_xp2_xpyp updated xp2.
        !----------------------------------------------------------------
        call clip_covars_denom( nzm, ngrdcol, sclr_dim,                        & ! In
                                dt,                                            & ! In
                                rtp2, thlp2, up2, vp2, wp2,                    & ! In
                                sclrp2, clubb_config_flags%l_tke_aniso,        & ! In
                                clubb_config_flags%l_linearize_pbl_winds,      & ! In
                                clubb_config_flags%l_predict_upwp_vpwp,        & ! In
                                stats,                                         & ! InOut
                                wprtp_cl_num, wpthlp_cl_num,                   & ! InOut
                                upwp_cl_num, vpwp_cl_num,                      & ! InOut
                                wprtp, wpthlp, upwp, vpwp, wpsclrp,            & ! InOut
                                upwp_pert, vpwp_pert )                           ! InOut

      elseif ( advance_order_loop_iter == order_wp2_wp3 ) then

        !----------------------------------------------------------------
        ! Advance the 2nd- and 3rd-order moments
        !   of vertical velocity (wp2, wp3) by one timestep.
        !----------------------------------------------------------------

        ! advance_wp2_wp3_bad_wp2 ! Test error comment, DO NOT modify or move
        call advance_wp2_wp3( nzm, nzt, ngrdcol, gr, dt_advance,                    & ! In
                              sfc_elevation, fcor_y, sigma_sqd_w, wm_zm,            & ! In
                              wm_zt,                                                & ! In
                              wpup2, wpvp2, wp2up2, wp2vp2, wp4,                    & ! In
                              wpthvp, wp2thvp, wp2up, um, vm, upwp, vpwp,           & ! In
                              em, Kh_zm, Kh_zt, invrs_tau_C4_zm,                    & ! In
                              invrs_tau_wp3_zt, invrs_tau_C1_zm,                    & ! In
                              rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm,                & ! In
                              invrs_rho_ds_zt, thv_ds_zm,                           & ! In
                              thv_ds_zt, pdf_params%mixt_frac, Cx_fnc_Richardson,   & ! In
                              lhs_splat_wp2, lhs_splat_wp3,                         & ! In
                              pdf_implicit_coefs_terms,                             & ! In
                              wprtp, wpthlp, rtp2, thlp2,                           & ! In
                              clubb_params, nu_vert_res_dep,                        & ! In
                              clubb_config_flags%iiPDF_type,                        & ! In
                              clubb_config_flags%penta_solve_method,                & ! In
                              clubb_config_flags%fill_holes_type,                   & ! In
                              clubb_config_flags%l_min_wp2_from_corr_wx,            & ! In
                              clubb_config_flags%l_upwind_xm_ma,                    & ! In
                              clubb_config_flags%l_tke_aniso,                       & ! In
                              clubb_config_flags%l_standard_term_ta,                & ! In
                              clubb_config_flags%l_partial_upwind_wp3,              & ! In
                              clubb_config_flags%l_damp_wp2_using_em,               & ! In
                              clubb_config_flags%l_use_C11_Richardson,              & ! In
                              clubb_config_flags%l_damp_wp3_Skw_squared,            & ! In
                              clubb_config_flags%l_lmm_stepping,                    & ! In
                              clubb_config_flags%l_use_tke_in_wp3_pr_turb_term,     & ! In
                              clubb_config_flags%l_use_tke_in_wp2_wp3_K_dfsn,       & ! In
                              clubb_config_flags%l_use_wp3_lim_with_smth_Heaviside, & ! In
                              clubb_config_flags%l_wp2_fill_holes_tke,              & ! In
                              clubb_config_flags%l_ho_nontrad_coriolis,             & ! In
                              l_implemented,                                        & ! In
                              stats,                                                & ! InOut
                              up2, vp2, wp2, wp3, err_info )                         ! InOut


        if ( clubb_at_least_debug_level_api( 0 ) ) then
          if ( any(err_info%err_code == clubb_fatal_error) ) then
              write(fstderr,*) err_info%err_header_global
              write(fstderr,*) "Error calling advance_wp2_wp3 in advance_clubb_core"
              return
          end if
        end if

        !----------------------------------------------------------------
        ! Covariance clipping for wprtp, wpthlp, wpsclrp, upwp, and vpwp
        ! after subroutine advance_wp2_wp3 updated wp2.
        !----------------------------------------------------------------
        call clip_covars_denom( nzm, ngrdcol, sclr_dim,                         & ! In
                                dt,                                             & ! In
                                rtp2, thlp2, up2, vp2, wp2,                    & ! In
                                sclrp2, clubb_config_flags%l_tke_aniso,        & ! In
                                clubb_config_flags%l_linearize_pbl_winds,      & ! In
                                clubb_config_flags%l_predict_upwp_vpwp,        & ! In
                                stats,                                          & ! InOut
                                wprtp_cl_num, wpthlp_cl_num,                   & ! InOut
                                upwp_cl_num, vpwp_cl_num,                      & ! InOut
                                wprtp, wpthlp, upwp, vpwp, wpsclrp,            & ! InOut
                                upwp_pert, vpwp_pert )                           ! InOut

      elseif ( advance_order_loop_iter == order_windm ) then

        !----------------------------------------------------------------
        ! Advance the horizontal mean winds (um, vm),
        !   the mean of the eddy-diffusivity scalars (i.e. edsclrm),
        !   and their fluxes (upwp, vpwp, wpedsclrp) by one time step.
        !----------------------------------------------------------------
        call advance_windm_edsclrm( nzm, nzt, ngrdcol, edsclr_dim, gr, dt,      & ! In
                                    wm_zt, Kh_zm, clubb_params,                 & ! In
                                    ug, vg, um_ref, vm_ref,                     & ! In
                                    wp2, up2, vp2, um_forcing, vm_forcing,      & ! In
                                    edsclrm_forcing, p_in_Pa,                   & ! In
                                    rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt,      & ! In
                                    fcor, l_implemented,                        & ! In
                                    nu_vert_res_dep, ts_nudge,                  & ! In
                                    clubb_config_flags%tridiag_solve_method,    & ! In
                                    clubb_config_flags%l_predict_upwp_vpwp,     & ! In
                                    clubb_config_flags%l_upwind_xm_ma,          & ! In
                                    clubb_config_flags%l_uv_nudge,              & ! In
                                    clubb_config_flags%l_tke_aniso,             & ! In
                                    clubb_config_flags%l_lmm_stepping,          & ! In
                                    clubb_config_flags%l_linearize_pbl_winds,   & ! In
                                    clubb_config_flags%l_do_expldiff_rtm_thlm,  & ! In
                                    clubb_config_flags%fill_holes_type,         & ! In
                                    order_xp2_xpyp, order_wp2_wp3, order_windm, & ! In
                                    upwp_cl_num, vpwp_cl_num,                   & ! InOut
                                    stats,                                      & ! InOut
                                    um, vm, thlm, rtm, edsclrm,                & ! InOut
                                    upwp, vpwp, wpedsclrp,                     & ! InOut
                                    um_pert, vm_pert, upwp_pert,               & ! InOut
                                    vpwp_pert, err_info )                        ! InOut

        if ( clubb_at_least_debug_level_api( 0 ) ) then
          if ( any(err_info%err_code == clubb_fatal_error) ) then
              write(fstderr,*) err_info%err_header_global
              write(fstderr,*) "Error calling advance_windm_edsclrm in advance_clubb_core"
              return
          end if
        end if

      endif ! advance_order_loop_iter

    enddo ! advance_order_loop_iter = 1, 4, 1

    !----------------------------------------------------------------
    ! Advance or otherwise calculate <thl'^3>, <rt'^3>, and
    ! <sclr'^3>.
    !----------------------------------------------------------------
    if ( l_advance_xp3 .and. clubb_config_flags%iiPDF_type /= iiPDF_ADG1 ) then

      ! Advance <rt'^3>, <thl'^3>, and <sclr'^3> one model timestep using a
      ! simplified form of the <x'^3> predictive equation.  The simplified
      ! <x'^3> equation can either be advanced from its previous value or
      ! calculated using a steady-state approximation.
      call advance_xp3( nzm, nzt, ngrdcol, sclr_dim, sclr_tol, gr, dt,   & ! In
                        rtm, thlm, rtp2, thlp2, wprtp,                   & ! In
                        wpthlp, wprtp2, wpthlp2, rho_ds_zm,              & ! In
                        invrs_rho_ds_zt, invrs_tau_zt, tau_max_zt,       & ! In
                        sclrm, sclrp2, wpsclrp, wpsclrp2,                & ! In
                        wp2, wp3, upwp, vpwp, up2, vp2,                  & ! In
                        thvm, clubb_params,                              & ! In
                        clubb_config_flags%l_lmm_stepping,               & ! In
                        stats,                                           & ! InOut
                        rtp3, thlp3, sclrp3, up3, vp3 )                    ! InOut

    else

      call diagnose_xp3( nzm, nzt, ngrdcol, sclr_dim, sclr_tol, gr,                 & ! In
                        clubb_config_flags%iiPDF_type, clubb_params,                & ! In
                        wp2, wp3, thvm,                                             & ! In
                        wprtp, wpthlp, rtp2, thlp2, upwp, vpwp, up2, vp2,           & ! In
                        sigma_sqd_w, wpsclrp, sclrp2,                               & ! In
                        stats,                                                      & ! InOut
                        rtp3, thlp3, up3, vp3,                                      & ! InOut
                        sclrp3 )                                                      ! InOut

    end if

    if (      clubb_config_flags%ipdf_call_placement == ipdf_post_advance_fields &
         .or. clubb_config_flags%ipdf_call_placement == ipdf_pre_post_advance_fields ) then

      ! Sample stats in this call to subroutine pdf_closure_driver for
      ! ipdf_post_advance_fields, but not for ipdf_pre_post_advance_fields
      ! because stats were sampled during the first call to subroutine
      ! pdf_closure_driver.
      if ( clubb_config_flags%ipdf_call_placement == ipdf_post_advance_fields ) then
        l_samp_stats_in_pdf_call = .true.
      elseif ( clubb_config_flags%ipdf_call_placement == ipdf_pre_post_advance_fields ) then
        l_samp_stats_in_pdf_call = .false.
      endif

      !########################################################################
      !#######                     CALL CLUBB's PDF                     #######
      !#######   AND OUTPUT PDF PARAMETERS AND INTEGRATED QUANTITITES   #######
      !########################################################################
      ! Given CLUBB's prognosed moments, diagnose CLUBB's PDF parameters
      !   and quantities integrated over that PDF, including
      !   quantities related to clouds, buoyancy, and turbulent advection.
      call compute_sigma_sqd_w( nzm, nzt, ngrdcol, gr,                        & ! In
                                wp3, wp2, thlp2, rtp2,                        & ! In
                                up2, vp2, wpthlp, wprtp, upwp, vpwp,          & ! In
                                clubb_params,                                 & ! In
                                clubb_config_flags%l_predict_upwp_vpwp,       & ! In
                                sigma_sqd_w )                                   ! Out

      call pdf_closure_driver( gr, nzm, nzt, ngrdcol,                       & ! In
                               dt, hydromet_dim, sclr_dim, sclr_tol,        & ! In
                               wprtp, thlm, wpthlp, rtp2, rtp3,             & ! In
                               thlp2, thlp3, rtpthlp, wp2,                  & ! In
                               wp3, wm_zm, wm_zt,                           & ! In
                               um, up2, upwp, up3,                          & ! In
                               vm, vp2, vpwp, vp3,                          & ! In
                               p_in_Pa, exner,                              & ! In
                               thv_ds_zm, thv_ds_zt, rtm_ref,               & ! In
                               wphydrometp,                                 & ! In
                               wp2hmp, rtphmp_zt, thlphmp_zt,               & ! In
                               sclrm, wpsclrp, sclrp2,                      & ! In
                               sclrprtp, sclrpthlp, sclrp3,                 & ! In
                               p_sfc, l_samp_stats_in_pdf_call,             & ! In
                               mixt_frac_max_mag, ts_nudge,                 & ! In
                               rtm_min, rtm_nudge_max_altitude,             & ! In
                               clubb_params,                                & ! In
                               clubb_config_flags%iiPDF_type,               & ! In
                               clubb_config_flags%saturation_formula,       & ! In
                               clubb_config_flags%l_predict_upwp_vpwp,      & ! In
                               clubb_config_flags%l_rtm_nudge,              & ! In
                               clubb_config_flags%l_trapezoidal_rule_zt,    & ! In
                               clubb_config_flags%l_trapezoidal_rule_zm,    & ! In
                               clubb_config_flags%l_call_pdf_closure_twice, & ! In
                               clubb_config_flags%l_use_cloud_cover,        & ! In
                               clubb_config_flags%l_rcm_supersat_adj,       & ! In
                               l_mix_rat_hm,                                & ! In
                               stats,                                       & ! InOut
                               rtm, sigma_sqd_w,                            & ! InOut
                               pdf_implicit_coefs_terms,                    & ! InOut
                               pdf_params, pdf_params_zm, err_info,         & ! InOut
                               rcm, cloud_frac,                             & ! Out
                               ice_supersat_frac, wprcp,                    & ! Out
                               wpthvp, wp2thvp, wp2up,                      & ! Out
                               rtpthvp, thlpthvp,                           & ! Out
                               rcm_in_layer, cloud_cover,                   & ! Out
                               rcp2_zt, thlprcp,                            & ! Out
                               rc_coef_zm, sclrpthvp,                       & ! Out
                               wpup2, wpvp2,                                & ! Out
                               wp2up2, wp2vp2, wp4,                         & ! Out
                               wp2rtp, wprtp2, wp2thlp,                     & ! Out
                               wpthlp2, wprtpthlp,                          & ! Out
                               uprcp, vprcp,                                & ! Out
                               w_up_in_cloud, w_down_in_cloud,              & ! Out
                               cloudy_updraft_frac,                         & ! Out
                               cloudy_downdraft_frac,                       & ! Out
                               wp2sclrp, wpsclrp2,                          & ! Out
                               wpsclrprtp, wpsclrpthlp )                      ! Out

      if ( clubb_at_least_debug_level_api( 0 ) ) then
         if ( any(err_info%err_code == clubb_fatal_error) ) then
            write(fstderr,*) err_info%err_header_global
            write(fstderr,*) "Error calling pdf_closure_driver in advance_clubb_core"
            return
         end if
      end if

    end if ! clubb_config_flags%ipdf_call_placement == ipdf_post_advance_fields
          ! or clubb_config_flags%ipdf_call_placement
          !    == ipdf_pre_post_advance_fields

#ifdef CLUBB_CAM
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzt
      do i = 1, ngrdcol
        qclvar(i,k) = rcp2_zt(i,k)
      end do
    end do
    !$acc end parallel loop
#endif


    !#######################################################################
    !#############            ACCUMULATE STATISTICS            #############
    !#######################################################################

    if ( stats%l_sample ) then

      call stats_accumulate( &
             nzm, nzt, ngrdcol, sclr_dim, edsclr_dim, gr, dt,                             & ! In
             l_implemented, clubb_config_flags%l_host_applies_sfc_fluxes,                 & ! In
             clubb_config_flags%l_stability_correct_tau_zm,                                & ! In
             clubb_params,                                                                 & ! In
             um, vm, upwp, vpwp, up2, vp2,                                                 & ! In
             thlm, rtm, thlm_before, rtm_before, thlm_forcing, rtm_forcing,                & ! In
             wpthlp_sfc, wprtp_sfc, wprtp, wpthlp,                                        & ! In
             wp2, wp3, rtp2, rtp3, thlp2, thlp3,                                           & ! In
             rtpthlp,                                                                       & ! In
             p_in_Pa, exner, rho, rho_zm,                                                  & ! In
             rho_ds_zm, rho_ds_zt, thv_ds_zm, thv_ds_zt,                                   & ! In
             wm_zt, wm_zm, rcm,                                                            & ! In
             cloud_frac,                                                                    & ! In
             thvm, ug, vg,                                                                  & ! In
             ddzt_umvm_sqd, stability_correction,                                          & ! In
             Kh_zt,                                                                         & ! In
             rsat,                                                                          & ! In
             Kh_zm,                                                                         & ! In
             em,                                                                             & ! In
             sclrm, sclrp2,                                                                 & ! In
             sclrprtp, sclrpthlp, sclrm_forcing,                                           & ! In
             wpsclrp, wpedsclrp, edsclrm,                                                  & ! In
             edsclrm_forcing,                                                               & ! In
             clubb_config_flags%saturation_formula,                                         & ! In
             stats )                                                                        ! InOut

      !$acc update host( wp2, vp2, up2, wprtp, wpthlp, rtp2, thlp2, rtpthlp, &
      !$acc              rtm, thlm, um, vm, wp3 )

      call stats_finalize_budget( "wp2_bt", wp2 / dt, stats )
      call stats_finalize_budget( "vp2_bt", vp2 / dt, stats )
      call stats_finalize_budget( "up2_bt", up2 / dt, stats )
      call stats_finalize_budget( "wprtp_bt", wprtp / dt, stats )
      call stats_finalize_budget( "wpthlp_bt", wpthlp / dt, stats )

      call stats_finalize_budget( "rtp2_bt", rtp2 / dt, stats )
      call stats_finalize_budget( "thlp2_bt", thlp2 / dt, stats )
      call stats_finalize_budget( "rtpthlp_bt", rtpthlp / dt, stats )

      call stats_finalize_budget( "rtm_bt", rtm / dt, stats )
      call stats_finalize_budget( "thlm_bt", thlm / dt, stats )
      call stats_finalize_budget( "um_bt", um / dt, stats )
      call stats_finalize_budget( "vm_bt", vm / dt, stats )
      call stats_finalize_budget( "wp3_bt", wp3 / dt, stats )

      if ( clubb_config_flags%l_predict_upwp_vpwp ) then
        !$acc update host( upwp, vpwp )
        call stats_finalize_budget( "upwp_bt", upwp / dt, stats )
        call stats_finalize_budget( "vpwp_bt", vpwp / dt, stats )
      end if

    end if

    if ( clubb_at_least_debug_level_api( 2 ) ) then

      call parameterization_check(                                         &
             nzm, nzt, ngrdcol, sclr_dim, edsclr_dim,                      & ! intent(in)
             thlm_forcing, rtm_forcing, um_forcing,                        & ! intent(in)
             vm_forcing, wm_zm, wm_zt, p_in_Pa,                            & ! intent(in)
             rho_zm, rho, exner, rho_ds_zm,                                & ! intent(in)
             rho_ds_zt, invrs_rho_ds_zm, invrs_rho_ds_zt,                  & ! intent(in)
             thv_ds_zm, thv_ds_zt, wpthlp_sfc, wprtp_sfc, upwp_sfc,        & ! intent(in)
             vpwp_sfc, p_sfc, um, upwp, vm, vpwp, up2, vp2,                & ! intent(in)
             rtm, wprtp, thlm, wpthlp, wp2, wp3,                           & ! intent(in)
             rtp2, thlp2, rtpthlp,                                         & ! intent(in)
             "end of ",                                                    & ! intent(in)
             wpsclrp_sfc, wpedsclrp_sfc, sclrm, wpsclrp,                   & ! intent(in)
             sclrp2,                                                       & ! intent(in)
             sclrprtp, sclrpthlp, sclrm_forcing, edsclrm,                  & ! intent(in)
             edsclrm_forcing,                                              & ! intent(in)
             err_info )                                                      ! intent(inout)

      if ( any(err_info%err_code == clubb_fatal_error) ) then
        write(fstderr,*) "Error occurred during parameterization_check at"// &
                         " end of advance_clubb_core"
        return
      end if

    end if

    !$acc exit data if( sclr_dim > 0 ) &
    !$acc           delete( wp2sclrp, &
    !$acc                   wpsclrp2, wpsclrprtp, wpsclrpthlp )

    !$acc exit data if( edsclr_dim > 0 ) &
    !$acc           delete( wpedsclrp )

    !$acc exit data delete( thvm, &
    !$acc                   wpthlp2, wprtp2, wprtpthlp, Lscale_up, &
    !$acc                   Lscale_zm, Lscale_down, em, tau_zm, &
    !$acc                   ddzt_umvm_sqd, rtm_before, thlm_before, &
    !$acc                   sigma_sqd_w, &
    !$acc                   sigma_sqd_w_zt, em_zt, w_1_zm, w_2_zm, &
    !$acc                   varnce_w_1_zm, &
    !$acc                   varnce_w_2_zm, &
    !$acc                   mixt_frac_zm, rcp2_zt, stability_correction, &
    !$acc                   invrs_tau_C6_zm, invrs_tau_C1_zm, invrs_tau_xp2_zm, &
    !$acc                   invrs_tau_C4_zm, invrs_tau_C14_zm, &
    !$acc                   invrs_tau_zt, invrs_tau_wp3_zt, Cx_fnc_Richardson, &
    !$acc                   brunt_vaisala_freq_sqd, brunt_vaisala_freq_sqd_mixed, &
    !$acc                   brunt_vaisala_freq_sqd_smth, &
    !$acc                   tau_max_zm, tau_max_zt, lhs_splat_wp2, lhs_splat_wp3 )

    return

  end subroutine advance_clubb_core

  !-----------------------------------------------------------------------

  subroutine set_sfc_value_of_flux_profiles( nzm, ngrdcol, sclr_dim, edsclr_dim, gr, &
                                             l_host_applies_sfc_fluxes, l_linearize_pbl_winds, &
                                             wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, &
                                             upwp_sfc_pert, vpwp_sfc_pert, &
                                             wpsclrp_sfc, wpedsclrp_sfc, &
                                             wpthlp, wprtp, upwp, vpwp, &
                                             upwp_pert, vpwp_pert, &
                                             wpsclrp, wpedsclrp )

    ! Description:
    !   Set or clear the surface values of turbulent flux profiles depending on
    !   whether the host model applies surface fluxes outside CLUBB.
    !
    ! References:
    !   None
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd

    use constants_clubb, only: &
        zero

    use grid_class, only: &
        grid

    implicit none

    !--------------------------- Input Variables ---------------------------
    integer, intent(in) :: &
      nzm,        & ! Number of momentum levels
      ngrdcol,    & ! Number of grid columns
      sclr_dim,   & ! Number of passive scalars
      edsclr_dim    ! Number of eddy-diffused scalars

    type(grid), intent(in) :: &
      gr   ! Grid structure

    logical, intent(in) :: &
      l_host_applies_sfc_fluxes,  & ! Whether host model applies surface fluxes outside CLUBB
      l_linearize_pbl_winds         ! Whether to linearize PBL wind equations

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
      wpthlp_sfc,     & ! Surface value of w'thl' [m K/s]
      wprtp_sfc,      & ! Surface value of w'rt' [m kg/(kg s)]
      upwp_sfc,       & ! Surface value of u'w' [m^2/s^2]
      vpwp_sfc,       & ! Surface value of v'w' [m^2/s^2]
      upwp_sfc_pert,  & ! Perturbation surface value of u'w' [m^2/s^2]
      vpwp_sfc_pert     ! Perturbation surface value of v'w' [m^2/s^2]

    real( kind = core_rknd ), dimension(ngrdcol,sclr_dim), intent(in) :: &
      wpsclrp_sfc   ! Surface value of w'sclr' [m sclr/s]

    real( kind = core_rknd ), dimension(ngrdcol,edsclr_dim), intent(in) :: &
      wpedsclrp_sfc   ! Surface value of w'edsclr' [m edsclr/s]

    !----------------------- Input/Output Variables ------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(inout) :: &
      wpthlp,     & ! Turbulent flux of thl [m K/s]
      wprtp,      & ! Turbulent flux of rt [m kg/(kg s)]
      upwp,       & ! Turbulent flux of u [m^2/s^2]
      vpwp,       & ! Turbulent flux of v [m^2/s^2]
      upwp_pert,  & ! Perturbation turbulent flux of u [m^2/s^2]
      vpwp_pert     ! Perturbation turbulent flux of v [m^2/s^2]

    real( kind = core_rknd ), dimension(ngrdcol,nzm,sclr_dim), intent(inout) :: &
      wpsclrp   ! Turbulent flux of passive scalars [m sclr/s]

    real( kind = core_rknd ), dimension(ngrdcol,nzm,edsclr_dim), intent(inout) :: &
      wpedsclrp   ! Turbulent flux of eddy-diffused scalars [m edsclr/s]

    !--------------------------- Local Variables ---------------------------
    integer :: &
      i,     & ! Grid-column loop index
      sclr,  & ! Passive scalar loop index
      edsclr   ! Eddy-diffused scalar loop index

    !----------------------------- Begin Code ------------------------------

    ! SET SURFACE VALUES OF FLUXES (BROUGHT IN)
    ! We only do this for host models that do not apply the flux
    ! elsewhere in the code (e.g. WRF).  In other cases the _sfc variables will
    ! only be used to compute the variance at the surface. -dschanen 8 Sept 2009
    if ( .not. l_host_applies_sfc_fluxes ) then

      !$acc parallel loop gang vector default(present)
      do i = 1, ngrdcol
        wpthlp(i,gr%k_lb_zm) = wpthlp_sfc(i)
        wprtp(i,gr%k_lb_zm)  = wprtp_sfc(i)
        upwp(i,gr%k_lb_zm)   = upwp_sfc(i)
        vpwp(i,gr%k_lb_zm)   = vpwp_sfc(i)
      end do
      !$acc end parallel loop

      if ( l_linearize_pbl_winds ) then
        !$acc parallel loop gang vector default(present)
        do i = 1, ngrdcol
          upwp_pert(i,gr%k_lb_zm) = upwp_sfc_pert(i)
          vpwp_pert(i,gr%k_lb_zm) = vpwp_sfc_pert(i)
        end do
        !$acc end parallel loop
      endif ! l_linearize_pbl_winds

      ! Set fluxes for passive scalars (if enabled)
      if ( sclr_dim > 0 ) then
        !$acc parallel loop gang vector collapse(2) default(present)
        do sclr = 1, sclr_dim
          do i = 1, ngrdcol
            wpsclrp(i,gr%k_lb_zm,sclr) = wpsclrp_sfc(i,sclr)
          end do
        end do
        !$acc end parallel loop
      end if

      if ( edsclr_dim > 0 ) then
        ! wpedsclrp is a local variable that is not saved from
        ! timestep to timestep. Set it to 0 and overwrite the surface.
        wpedsclrp = zero
        !$acc parallel loop gang vector collapse(2) default(present)
        do edsclr = 1, edsclr_dim
          do i = 1, ngrdcol
            wpedsclrp(i,gr%k_lb_zm,edsclr) = wpedsclrp_sfc(i,edsclr)
          end do
        end do
        !$acc end parallel loop
      end if

    else

      !$acc parallel loop gang vector default(present)
      do i = 1, ngrdcol
        wpthlp(i,gr%k_lb_zm) = 0.0_core_rknd
        wprtp(i,gr%k_lb_zm)  = 0.0_core_rknd
        upwp(i,gr%k_lb_zm)   = 0.0_core_rknd
        vpwp(i,gr%k_lb_zm)   = 0.0_core_rknd
      end do
      !$acc end parallel loop

      ! Set fluxes for passive scalars (if enabled)
      if ( sclr_dim > 0 ) then
        !$acc parallel loop gang vector collapse(2) default(present)
        do sclr = 1, sclr_dim
          do i = 1, ngrdcol
            wpsclrp(i,gr%k_lb_zm,sclr) = 0.0_core_rknd
          end do
        end do
        !$acc end parallel loop
      end if

      if ( edsclr_dim > 0 ) then
        ! wpedsclrp is a local variable that is not saved from
        ! timestep to timestep. Set it to 0 and overwrite the surface.
        wpedsclrp = zero
        !$acc parallel loop gang vector collapse(2) default(present)
        do edsclr = 1, edsclr_dim
          do i = 1, ngrdcol
            wpedsclrp(i,gr%k_lb_zm,edsclr) = 0.0_core_rknd
          end do
        end do
        !$acc end parallel loop
      end if

    end if ! ~l_host_applies_sfc_fluxes

  end subroutine set_sfc_value_of_flux_profiles

end module advance_clubb_core_module
