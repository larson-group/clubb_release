
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

  subroutine advance_clubb_core( gr, nzm, nzt, ngrdcol, &           ! intent(in)
               l_implemented, dt, fcor, fcor_y, sfc_elevation, &    ! intent(in)
               hydromet_dim,                                      & ! intent(in)
               sclr_dim, sclr_tol, edsclr_dim, sclr_idx,      &     ! intent(in)
               thlm_forcing, rtm_forcing, um_forcing, vm_forcing, & ! intent(in)
               sclrm_forcing, edsclrm_forcing, wprtp_forcing, &     ! intent(in)
               wpthlp_forcing, rtp2_forcing, thlp2_forcing, &       ! intent(in)
               rtpthlp_forcing, wm_zm, wm_zt, &                     ! intent(in)
               wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, p_sfc, &  ! intent(in)
               wpsclrp_sfc, wpedsclrp_sfc, &                        ! intent(in)
               upwp_sfc_pert, vpwp_sfc_pert, &                      ! intent(in)
               rtm_ref, thlm_ref, um_ref, vm_ref, ug, vg, &         ! Intent(in)
               p_in_Pa, rho_zm, rho, exner, &                       ! intent(in)
               rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &             ! intent(in)
               invrs_rho_ds_zt, thv_ds_zm, thv_ds_zt, &             ! intent(in)
               l_mix_rat_hm, &                                      ! intent(in)
               rfrzm, &                                             ! intent(in)
#ifdef CLUBBND_CAM
               varmu, &                                             ! intent(in)
#endif
               wphydrometp, wp2hmp, rtphmp_zt, thlphmp_zt, &        ! intent(in)
               host_dx, host_dy, &                                  ! intent(in)
               clubb_params, nu_vert_res_dep, lmin, &               ! intent(in)
               mixt_frac_max_mag, T0, ts_nudge, &                   ! Intent(in)
               rtm_min, rtm_nudge_max_altitude, &                   ! Intent(in)
               clubb_config_flags, &                                ! intent(in)
               stats_metadata, &                                    ! intent(in)
               stats_zt, stats_zm, stats_sfc, &                     ! intent(inout)
               um, vm, upwp, vpwp, up2, vp2, up3, vp3, &            ! intent(inout)
               thlm, rtm, wprtp, wpthlp, &                          ! intent(inout)
               wp2, wp3, rtp2, rtp3, thlp2, thlp3, rtpthlp, &       ! intent(inout)
               sclrm,   &                                           ! intent(inout)
#ifdef GFDL
               sclrm_trsport_only,  &  ! h1g, 2010-06-16            ! intent(inout)
#endif
               sclrp2, sclrp3, sclrprtp, sclrpthlp, &               ! intent(inout)
               wpsclrp, edsclrm, &                                  ! intent(inout)
               rcm, cloud_frac, &                                   ! intent(inout)
               wpthvp, wp2thvp, wp2up, rtpthvp, thlpthvp, &         ! intent(inout)
               sclrpthvp, &                                         ! intent(inout)
               wp2rtp, wp2thlp, uprcp, vprcp, rc_coef_zm, wp4, &    ! intent(inout)
               wpup2, wpvp2, wp2up2, wp2vp2, ice_supersat_frac, &   ! intent(inout)
               um_pert, vm_pert, upwp_pert, vpwp_pert, &            ! intent(inout)
               pdf_params, pdf_params_zm, &                         ! intent(inout)
               pdf_implicit_coefs_terms, &                          ! intent(inout)
               err_info, &                                          ! intent(inout)
#ifdef GFDL
               RH_crit, & !h1g, 2010-06-16                          ! intent(inout)
               do_liquid_only_in_clubb, &                           ! intent(in)
#endif
               Kh_zm, Kh_zt, &                                      ! intent(out)
#ifdef CLUBB_CAM
               qclvar, &                                            ! intent(out)
#endif
               thlprcp, wprcp, w_up_in_cloud, w_down_in_cloud, &    ! intent(out)
               cloudy_updraft_frac, cloudy_downdraft_frac, &        ! intent(out)
               rcm_in_layer, cloud_cover, invrs_tau_zm, &           ! intent(out)
               Lscale )                                             ! intent(out)

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
        thl_tol, &
        rt_tol, &
        w_tol, &
        w_tol_sqd, &
        fstderr, &
        zero_threshold, &
        three_halves, &
        one, &
        two, &
        zero, &
        unused_var, &
        grav, &
        eps, &
        min_max_smth_mag

    use parameter_indices, only: &
        nparams,                 & ! Variable(s)
        itaumax,                 &
        ic_K,                    &
        ic_K10,                  &
        ic_K10h,                 &
        imu,                     &
        igamma_coef,             &
        igamma_coefb,            &
        igamma_coefc,            &
        iC_wp2_splat,            &
        ixp3_coef_base,          &
        ixp3_coef_slope,         &
        ilambda0_stability_coef, &
        iup2_sfc_coef,           &
        ia_const,                &
        ia3_coef_min,            &
        ibv_efold

    use parameters_tunable, only: &
        nu_vertical_res_dep    ! Type(s)

    use model_flags, only: &
        clubb_config_flags_type, & ! Type
        l_gamma_Skw, & ! Variable(s)
        l_advance_xp3, &
        iiPDF_ADG1, &
        order_xm_wpxp, &
        order_xp2_xpyp, &
        order_wp2_wp3, &
        order_windm

    use grid_class, only: &
        grid,       & ! Type
        zm2zt_api,  & ! Procedure(s)
        zt2zm_api,  &
        ddzm,       &
        ddzt,       &
        zm2zt2zm

    use numerical_check, only: &
        parameterization_check, & ! Procedure(s)
        calculate_spurious_source

    use pdf_parameter_module, only: &
        pdf_parameter, &
        implicit_coefs_terms

#ifdef GFDL
    use advance_sclrm_Nd_module, only: &  ! h1g, 2010-06-16 begin mod
         advance_sclrm_Nd_diffusion_OG, &
         advance_sclrm_Nd_upwind, &
       advance_sclrm_Nd_semi_implicit     ! h1g, 2010-06-16 end mod
#endif

    use advance_xm_wpxp_module, only: &
        advance_xm_wpxp          ! Compute mean/flux terms

    use advance_xp2_xpyp_module, only: &
        advance_xp2_xpyp     ! Computes variance terms

    use sfc_varnce_module, only:  &
        calc_sfc_varnce ! Procedure

    use mixing_length, only: &
        calc_Lscale_directly,     &  ! for Lscale
        diagnose_Lscale_from_tau, &  ! for Lscale from tau
        set_Lscale_max

    use advance_windm_edsclrm_module, only:  &
        advance_windm_edsclrm  ! Procedure(s)

    use saturation, only:  &
        ! Procedure
        sat_mixrat_liq_api ! Saturation mixing ratio

    use advance_wp2_wp3_module, only:  &
        advance_wp2_wp3 ! Procedure

    use advance_xp3_module, only: &
        advance_xp3    ! Procedure(s)

    use calc_pressure, only: &
        calculate_thvm

    use clubb_precision, only:  &
        core_rknd ! Variable(s)

    use error_code, only: &
        clubb_at_least_debug_level_api,  & ! Procedure
        clubb_fatal_error              ! Constant

    use Skx_module, only: &
        Skx_func,           & ! Procedure(s)
        xp3_LG_2005_ansatz

    use clip_explicit, only: &
        clip_covars_denom, & ! Procedure(s)
        clip_rcm

    use T_in_K_module, only: &
        ! Read values from namelist
        thlm2T_in_K_api ! Procedure

    use sigma_sqd_w_module, only: &
        compute_sigma_sqd_w    ! Procedure(s)

    use stats_clubb_utilities, only: &
        stats_accumulate ! Procedure

    use stats_type_utilities, only:   &
        stat_update_var_pt,   & ! Procedure(s)
        stat_update_var,      &
        stat_begin_update,    &
        stat_begin_update_pt, &
        stat_end_update,      &
        stat_end_update_pt

    use fill_holes, only: &
        fill_holes_vertical_api

    use advance_helper_module, only: &
        calc_stability_correction, & ! Procedure(s)
        compute_Cx_fnc_Richardson, &
        calc_brunt_vaisala_freq_sqd, &
        wp2_term_splat_lhs, &
        wp3_term_splat_lhs, &
        vertical_integral, &
        Lscale_width_vert_avg, &
        smooth_max, &
        calc_Ri_zm, &
        calculate_thlp2_rad, &
        pvertinterp

    use stats_type, only: stats ! Type

    use stats_variables, only: &
      stats_metadata_type

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
      l_implemented    ! True if CLUBB is being run within a large-scale host model,
                       !   rather than a standalone single-column model.

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

#ifdef CLUBBND_CAM
    real( kind = core_rknd ), intent(in), dimension(ngrdcol) :: &
      varmu
#endif

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

    type (stats_metadata_type), intent(in) :: &
      stats_metadata

    !--------------------------- Input/Output Variables ---------------------------
    type (stats), intent(inout), dimension(ngrdcol) :: &
      stats_zt, &
      stats_zm, &
      stats_sfc

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

#ifdef GFDL
    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,nzt,sclr_dim) :: &  ! h1g, 2010-06-16
      sclrm_trsport_only  ! Passive scalar concentration due to pure transport [{units vary}/s]
#endif

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

#ifdef GFDL
    ! hlg, 2010-06-16
    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,nzt, min(1,sclr_dim) , 2) :: &
      RH_crit  ! critical relative humidity for droplet and ice nucleation
! ---> h1g, 2012-06-14
    logical, intent(in)                 ::  do_liquid_only_in_clubb
! <--- h1g, 2012-06-14
#endif

    !--------------------------- Local Variables ---------------------------
    integer :: i, k, edsclr, sclr

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      Skw_zt,       & ! Skewness of w on thermodynamic levels            [-]
      thvm,         & ! Virtual potential temperature                    [K]
      ddzm_thvm_zm    ! d(thvm_zm)/dz, centered over thermodynamic levs. [K/m]

    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
      Skw_zm,       & ! Skewness of w on momentum levels                 [-]
      thvm_zm         ! Virtual potential temperature on momentum levs.  [K]

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      rsat   ! Saturation mixing ratio  ! Brian

    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
      rtprcp, & ! rt'rc'               [kg^2/kg^2]
      rcp2      ! rc'^2                [kg^2/kg^2]

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      wpthlp2,   & ! w'thl'^2    [m K^2/s]
      wprtp2,    & ! w'rt'^2     [m kg^2/kg^2]
      wprtpthlp, & ! w'rt'thl'   [m kg K/kg s]
      wp2rcp       ! w'^2 rc'    [m^2 kg/kg s^2]

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      Lscale_up,   & ! Length scale (upwards component)      [m]
      Lscale_down    ! Length scale (downwards component)    [m]

    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
      wp3_zm,      & ! w'^3                                  [m^3/s^3]
      Lscale_zm      ! Length scale on momentum levels       [m]

    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
      em,     & ! Turbulent Kinetic Energy (TKE)                      [m^2/s^2]
      tau_zm    ! Eddy dissipation time scale on momentum levels      [s]

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      tau_zt    ! Eddy dissipation time scale on thermodynamic levels [s]

    real( kind = core_rknd ), dimension(ngrdcol,nzm,edsclr_dim) :: &
      wpedsclrp   ! w'edsclr'

    real( kind = core_rknd ), dimension(ngrdcol,nzm,sclr_dim) :: &
      sclrprcp       ! sclr'rc'

    real( kind = core_rknd ), dimension(ngrdcol,nzt,sclr_dim) :: &
      wp2sclrp,    & ! w'^2 sclr'
      wpsclrp2,    & ! w'sclr'^2
      wpsclrprtp,  & ! w'sclr'rt'
      wpsclrpthlp    ! w'sclr'thl'

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      wp2_zt,     & ! w'^2 on thermo. grid     [m^2/s^2]
      thlp2_zt,   & ! thl'^2 on thermo. grid   [K^2]
      wpthlp_zt,  & ! w'thl' on thermo. grid   [m K/s]
      wprtp_zt,   & ! w'rt' on thermo. grid    [m kg/(kg s)]
      rtp2_zt,    & ! rt'^2 on therm. grid     [(kg/kg)^2]
      rtpthlp_zt, & ! rt'thl' on thermo. grid  [kg K/kg]
      up2_zt,     & ! u'^2 on thermo. grid     [m^2/s^2]
      vp2_zt,     & ! v'^2 on thermo. grid     [m^2/s^2]
      upwp_zt,    & ! u'w' on thermo. grid     [m^2/s^2]
      vpwp_zt       ! v'w' on thermo. grid     [m^2/s^2]

    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
      wp3_on_wp2,   & ! w'^3 / w'^2 on the zm grid                     [m/s]
      Skw_velocity, & ! Skewness velocity                              [m/s]
      a3_coef         ! The a3 coefficient from CLUBB eqns             [-]

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      wp3_on_wp2_zt, & ! w'^3 / w'^2 on the zt grid [m/s]
      a3_coef_zt       ! The a3 coefficient interpolated to the zt grid [-]

    ! Eric Raut declared this variable solely for output to disk
    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      rc_coef    ! Coefficient of X'r_c' in Eq. (34) on t-levs.  [K/(kg/kg)]

    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
      Km_zm, & ! Eddy diffusivity for momentum on zm grid levels [m^2/s]
      Kmh_zm   ! Eddy diffusivity for thermodynamic variables [m^2/s]

    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
      gamma_Skw_fnc,   & ! Gamma as a function of skewness               [-]
      sigma_sqd_w,     & ! PDF width parameter (momentum levels)         [-]
      sigma_sqd_w_tmp

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      sigma_sqd_w_zt, & ! PDF width parameter (thermodynamic levels)    [-]
      sqrt_em_zt,     & ! sqrt( em ) on zt levels; where em is TKE      [m/s]
      xp3_coef_fnc      ! Coefficient in simple xp3 equation            [-]
!Lscale_weight Uncomment this if you need to use this vairable at some point.

    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
      w_1_zm,        & ! Mean w (1st PDF component)                   [m/s]
      w_2_zm,        & ! Mean w (2nd PDF component)                   [m/s]
      varnce_w_1_zm, & ! Variance of w (1st PDF component)            [m^2/s^2]
      varnce_w_2_zm, & ! Variance of w (2nd PDF component)            [m^2/s^2]
      mixt_frac_zm     ! Weight of 1st PDF component (Sk_w dependent) [-]

    integer :: &
      wprtp_cl_num,   & ! Instance of w'r_t' clipping (1st or 3rd).
      wpthlp_cl_num,  & ! Instance of w'th_l' clipping (1st or 3rd).
      wpsclrp_cl_num, & ! Instance of w'sclr' clipping (1st or 3rd).
      upwp_cl_num,    & ! Instance of u'w' clipping (1st or 2nd).
      vpwp_cl_num       ! Instance of v'w' clipping (1st or 2nd).

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      rcp2_zt,              & ! r_c'^2 (on thermo. grid)             [kg^2/kg^2]
      wpsclrp_zt,           & ! Scalar flux on thermo. levels        [un. vary]
      sclrp2_zt               ! Scalar variance on thermo.levels     [un. vary]

    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
      cloud_frac_zm,        & ! Cloud Fraction on momentum grid      [-]
      ice_supersat_frac_zm, & ! Ice Cloud Fraction on momentum grid  [-]
      rtm_zm,               & ! Total water mixing ratio             [kg/kg]
      thlm_zm,              & ! Liquid potential temperature         [kg/kg]
      rcm_zm                  ! Liquid water mixing ratio on m-levs. [kg/kg]

    real( kind = core_rknd ), dimension(ngrdcol) :: &
      rtm_integral_before, &
      rtm_integral_after, &
      rtm_integral_forcing, &
      rtm_flux_top, &
      rtm_flux_sfc, &
      rtm_spur_src, &
      thlm_integral_before, &
      thlm_integral_after, &
      thlm_integral_forcing, &
      thlm_flux_top, &
      thlm_flux_sfc, &
      thlm_spur_src

    real( kind = core_rknd ), dimension(ngrdcol) :: &
      thlm1000, &
      thlm700                      

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      rcm_supersat_adj, & ! Adjustment to rcm due to spurious supersaturation
      rel_humidity        ! Relative humidity after PDF closure [-]

    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
       stability_correction,         & ! Stability correction factor
       invrs_tau_N2_zm,              & ! Inverse tau with static stability correction applied [1/s]
       invrs_tau_C6_zm,              & ! Inverse tau values used for C6 (pr1) term in wpxp [1/s]
       invrs_tau_C1_zm,              & ! Inverse tau values used for C1 (dp1) term in wp2 [1/s]
       invrs_tau_xp2_zm,             & ! Inverse tau values used for advance_xp2_wpxp [s^-1]
       invrs_tau_N2_iso,             & ! Inverse tau values used for C4 when 
                                       ! l_use_invrs_tau_N2_iso = .true.              [s^-1]
       invrs_tau_C4_zm,              & ! Inverse tau values used for C4 terms         [s^-1]
       invrs_tau_C14_zm,             & ! Inverse tau valuse used for C14 terms        [s^-1]
       invrs_tau_wp2_zm,             & ! Inverse tau values used for advance_wp2_wpxp [s^-1]
       invrs_tau_wpxp_zm,            & ! invrs_tau_C6_zm = invrs_tau_wpxp_zm
       invrs_tau_wp3_zm,             & ! Inverse tau values used for advance_wp3_wp2 [s^-1]
       invrs_tau_no_N2_zm,           & ! One divided by tau (without N2) on zm levels [s^-1]
       invrs_tau_bkgnd,              & ! One divided by tau_wp3 [s^-1]
       invrs_tau_shear,              & ! One divided by tau with stability effects    [s^-1]
       invrs_tau_sfc,                & ! One divided by tau (without N2) on zm levels [s^-1]
       Cx_fnc_Richardson,            & ! Cx_fnc computed from Richardson_num          [-]
       ddzt_um,                      & ! Vertical derivative of um                    [s^-1]
       ddzt_vm,                      & ! Vertical derivative of vm                    [s^-1]
       ddzt_umvm_sqd,                & ! Squared vertical norm of derivative of
                                       ! mean horizontal wind speed                   [s^-2]
       ddzt_umvm_sqd_clipped,        & ! Smooth_maxed version of ddzt_umvm_sqd        [s^-2]
       brunt_vaisala_freq_sqd,       & ! Buoyancy frequency squared, N^2              [s^-2]
       brunt_vaisala_freq_sqd_mixed, & ! A mixture of dry and moist N^2               [s^-2]
       brunt_vaisala_freq_sqd_dry,   & ! dry N^2                                      [s^-2]
       brunt_vaisala_freq_sqd_moist, & ! moist N^2                                    [s^-2]
       brunt_vaisala_freq_sqd_smth,  & ! Mix between dry and moist N^2 that is
                                       ! smoothed in the vertical                     [s^-2]
       brunt_vaisala_freq_sqd_splat, & !                                              [s^-2]
       brunt_vaisala_freq_clipped,   & ! N^2 clipped via smooth_max function          [s^-2]
       Ri_zm                           ! Richardson number                            [-]

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
       invrs_tau_zt,                 & ! Inverse time-scale tau on thermodynamics levels [1/s]
       invrs_tau_wp3_zt,             & ! Inverse tau wp3 at zt levels
       brunt_vaisala_freq_sqd_zt       ! Buoyancy frequency squared on t-levs.        [s^-2]

    real( kind = core_rknd ), parameter :: &
       ufmin = 0.01_core_rknd           ! minimum value of friction velocity     [m/s]

    real( kind = core_rknd ), dimension(ngrdcol) :: &
      Lscale_max    ! Max. allowable mixing length (based on grid box size) [m]

    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
      tau_max_zm    ! Max. allowable eddy dissipation time scale on m-levs  [s]

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      tau_max_zt    ! Max. allowable eddy dissipation time scale on t-levs  [s]

    real( kind = core_rknd ) :: &
      below_grnd_val = 0.01_core_rknd

    real( kind = core_rknd ), dimension(ngrdcol) :: &
      mu

    real( kind = core_rknd ) :: &
      gamma_coef,     & ! CLUBB tunable parameter gamma_coef
      gamma_coefb,    & ! CLUBB tunable parameter gamma_coefb
      gamma_coefc       ! CLUBB tunable parameter gamma_coefc

    ! Flag to sample stats in a particular call to subroutine
    ! pdf_closure_driver.
    logical :: l_samp_stats_in_pdf_call

    ! Flag to determine whether invrs_tau_N2_iso is used in C4 terms.
    ! Important! This flag is only in use when l_diag_Lscale_from_tau = true
    ! Setting l_use_invrs_tau_N2_iso = true will not change anything unless
    ! l_diag_Lscale_from_tau is also true
    logical, parameter :: l_use_invrs_tau_N2_iso = .false.

    logical, parameter :: l_smooth_min_max = .false.  ! whether to apply smooth min 
                                                      ! and max functions

    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
       lhs_splat_wp2    ! LHS coefficient of wp2 splatting term  [1/s]

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
       lhs_splat_wp3    ! LHS coefficient of wp3 splatting term  [1/s]

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
       thlm_ed, &  ! explicit diffusion budget term for thlm
       rtm_ed      ! explicit diffusion budget term for rtm

    ! Variables associated with upgradient momentum contributions due to cumuli
    !real( kind = core_rknd ), dimension(nz) :: &
    !  Km_Skw_factor ! Factor, with value < 1, that reduces eddy diffusivity,
    !                                          Km_zm, in skewed layers
    !real( kind = core_rknd ),parameter :: &
    !  Km_Skw_thresh = zero_threshold, &  ! Value of Skw at which Skw correction kicks in
    !  Km_Skw_factor_efold = 0.5_core_rknd, & ! E-folding rate of exponential Skw correction
    !  Km_Skw_factor_min   = 0.2_core_rknd    ! Minimum value of Km_Skw_factor

    integer :: advance_order_loop_iter

    integer :: smth_type = 2  ! Used for Lscale_width_vert_avg

    !----- Begin Code -----

    !$acc enter data create( Skw_zm, Skw_zt, thvm, thvm_zm, ddzm_thvm_zm, rtprcp, rcp2, &
    !$acc              ddzt_um, ddzt_vm, ddzt_umvm_sqd, ddzt_umvm_sqd_clipped, &
    !$acc              wpthlp2, wprtp2, wprtpthlp, wp2rcp, wp3_zm, Lscale_up, Lscale_zm, &
    !$acc              Lscale_down, em, tau_zm, tau_zt, wp2_zt, thlp2_zt, wpthlp_zt, &
    !$acc              wprtp_zt, rtp2_zt, rtpthlp_zt, up2_zt, vp2_zt, upwp_zt, vpwp_zt, &
    !$acc              Skw_velocity, a3_coef, a3_coef_zt, wp3_on_wp2, wp3_on_wp2_zt, rc_coef, &
    !$acc              Km_zm, Kmh_zm, gamma_Skw_fnc, sigma_sqd_w, sigma_sqd_w_tmp, sigma_sqd_w_zt, &
    !$acc              sqrt_em_zt, xp3_coef_fnc, w_1_zm, w_2_zm, varnce_w_1_zm, varnce_w_2_zm, &
    !$acc              mixt_frac_zm, rcp2_zt, cloud_frac_zm, ice_supersat_frac_zm, rtm_zm, &
    !$acc              thlm_zm, rcm_zm, thlm1000, thlm700, &
    !$acc              rcm_supersat_adj, stability_correction, invrs_tau_N2_zm, &
    !$acc              invrs_tau_C6_zm, invrs_tau_C1_zm, invrs_tau_xp2_zm, invrs_tau_N2_iso, &
    !$acc              invrs_tau_C4_zm, invrs_tau_C14_zm, invrs_tau_wp2_zm, invrs_tau_wpxp_zm, &
    !$acc              invrs_tau_wp3_zm, invrs_tau_no_N2_zm, invrs_tau_bkgnd, invrs_tau_shear, &
    !$acc              invrs_tau_sfc, invrs_tau_zt, invrs_tau_wp3_zt, Cx_fnc_Richardson, &
    !$acc              brunt_vaisala_freq_sqd, brunt_vaisala_freq_sqd_mixed, &
    !$acc              brunt_vaisala_freq_sqd_dry, brunt_vaisala_freq_sqd_moist, &
    !$acc              brunt_vaisala_freq_sqd_splat, brunt_vaisala_freq_sqd_smth, &
    !$acc              brunt_vaisala_freq_sqd_zt, brunt_vaisala_freq_clipped, Ri_zm, Lscale_max, &
    !$acc              tau_max_zm, tau_max_zt, mu, lhs_splat_wp2, lhs_splat_wp3 )

    !$acc enter data if( sclr_dim > 0 ) &
    !$acc            create( sclrprcp, wp2sclrp, &
    !$acc                    wpsclrp2, wpsclrprtp, wpsclrpthlp, wpsclrp_zt, sclrp2_zt )

    !$acc enter data if( edsclr_dim > 0 ) &
    !$acc            create( wpedsclrp )

    if ( clubb_config_flags%l_lmm_stepping ) then
      dt_advance = two * dt
    else
      dt_advance = dt
    end if

    ! Determine the maximum allowable value for Lscale (in meters).
    call set_Lscale_max( ngrdcol, l_implemented, host_dx, host_dy, & ! intent(in)
                         Lscale_max )                                ! intent(out)

    if ( stats_metadata%l_stats .and. stats_metadata%l_stats_samp ) then

      !$acc update host( wm_zt, wm_zm, rho_ds_zt, rtm, gr%dzt, &
      !$acc              rtm, thlm )

      ! Spurious source will only be calculated if rtm_ma and thlm_ma are zero.
      ! Therefore, wm must be zero or l_implemented must be true.
      
      do i = 1, ngrdcol
        if ( l_implemented .or. ( all( abs(wm_zt(i,:)) < eps ) .and. &
             all( abs(wm_zm(i,:)) < eps ) ) ) then
          ! Get the vertical integral of rtm and thlm before this function begins
          ! so that spurious source can be calculated
          rtm_integral_before(i)  &
          = vertical_integral( nzt, rho_ds_zt(i,:), &
                               rtm(i,:), gr%dzt(i,:) )

          thlm_integral_before(i)  &
          = vertical_integral( nzt, rho_ds_zt(i,:), &
                               thlm(i,:), gr%dzt(i,:) )
        end if
      end do
    end if

    !----------------------------------------------------------------
    ! Test input variables
    !----------------------------------------------------------------
    if ( clubb_at_least_debug_level_api( 2 ) ) then

      !$acc update host( thlm_forcing, rtm_forcing, um_forcing, vm_forcing, &
      !$acc              wm_zm, wm_zt, p_in_Pa, rho_zm, rho, exner, rho_ds_zm, &
      !$acc              rho_ds_zt, invrs_rho_ds_zm, invrs_rho_ds_zt, thv_ds_zm, &
      !$acc              thv_ds_zt, wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, p_sfc, &
      !$acc              um, upwp, vm, vpwp, up2, vp2, rtm, wprtp, thlm, wpthlp, &
      !$acc              wp2, wp3, rtp2, thlp2, rtpthlp )

      !$acc update host( wpsclrp_sfc, wpedsclrp_sfc, sclrm, wpsclrp, &
      !$acc              sclrp2, sclrprtp, sclrpthlp, sclrm_forcing ) &
      !$acc if ( sclr_dim > 0 )

      !$acc update host( edsclrm, edsclrm_forcing  ) if ( edsclr_dim > 0 )

      do i = 1, ngrdcol
        call parameterization_check( &
               nzm, nzt, sclr_dim, edsclr_dim, &
               thlm_forcing(i,:), rtm_forcing(i,:), um_forcing(i,:),                     & ! In
               vm_forcing(i,:), wm_zm(i,:), wm_zt(i,:), p_in_Pa(i,:),                    & ! In
               rho_zm(i,:), rho(i,:), exner(i,:), rho_ds_zm(i,:),                        & ! In
               rho_ds_zt(i,:), invrs_rho_ds_zm(i,:), invrs_rho_ds_zt(i,:),               & ! In
               thv_ds_zm(i,:), thv_ds_zt(i,:), wpthlp_sfc(i), wprtp_sfc(i), upwp_sfc(i), & ! In
               vpwp_sfc(i), p_sfc(i), um(i,:), upwp(i,:), vm(i,:), vpwp(i,:), up2(i,:),  & ! In
               vp2(i,:), rtm(i,:), wprtp(i,:), thlm(i,:), wpthlp(i,:), wp2(i,:),         & ! In
               wp3(i,:), rtp2(i,:), thlp2(i,:), rtpthlp(i,:),                            & ! In
               !rcm,                                                                     &
               "beginning of ",                                                          & ! In
               wpsclrp_sfc(i,:), wpedsclrp_sfc(i,:), sclrm(i,:,:), wpsclrp(i,:,:),       & ! In
               sclrp2(i,:,:),                                                            & ! In
               sclrprtp(i,:,:), sclrpthlp(i,:,:), sclrm_forcing(i,:,:), edsclrm(i,:,:),  & ! In
               edsclrm_forcing(i,:,:), &                                                   ! In
               err_info )                                                                  ! Inout

        if ( err_info%err_code(i) == clubb_fatal_error ) then
          write(fstderr,*) err_info%err_header(i)
        endif
      end do

      if ( any(err_info%err_code == clubb_fatal_error) ) then
        write(fstderr,*) "Fatal error detected in parameterization_check when testing input"
        return
      end if

    end if
    !-----------------------------------------------------------------------

    if ( stats_metadata%l_stats_samp ) then

      !$acc update host( rfrzm, wp2, vp2, up2, wprtp, wpthlp, upwp, vpwp, &
      !$acc              rtp2, thlp2, rtpthlp, rtm, thlm, um, vm, wp3 )

      do i = 1, ngrdcol

        call stat_update_var( stats_metadata%irfrzm, rfrzm(i,:), & ! intent(in)
                              stats_zt(i) ) ! intent(inout)

      ! Set up budget stats variables.
        
        !print *, "B stats_zt(i)%accum_field_values", stats_zt(i)%accum_field_values
        !print *, "wp2(i,:) = ", wp2(i,:)

         call stat_begin_update( nzm, stats_metadata%iwp2_bt, wp2(i,:) / dt, & ! intent(in)
                                 stats_zm(i) )           ! intent(inout)
                                 
         !print *, "A stats_zt(i)%accum_field_values", stats_zt(i)%accum_field_values
                                 
                                 
         call stat_begin_update( nzm, stats_metadata%ivp2_bt, vp2(i,:) / dt, & ! intent(in)
                                 stats_zm(i) )           ! intent(inout)
         call stat_begin_update( nzm, stats_metadata%iup2_bt, up2(i,:) / dt, & ! intent(in)
                                 stats_zm(i) )           ! intent(inout)
         call stat_begin_update( nzm, stats_metadata%iwprtp_bt, wprtp(i,:) / dt, & ! intent(in)
                                 stats_zm(i) )               ! intent(inout)
         call stat_begin_update( nzm, stats_metadata%iwpthlp_bt, wpthlp(i,:) / dt, & ! intent(in)
                                 stats_zm(i) )                 ! intent(inout)
         if ( clubb_config_flags%l_predict_upwp_vpwp ) then
            call stat_begin_update( nzm, stats_metadata%iupwp_bt, upwp(i,:) / dt, & ! intent(in)
                                    stats_zm(i) )             ! intent(inout)
            call stat_begin_update( nzm, stats_metadata%ivpwp_bt, vpwp(i,:) / dt, & ! intent(in)
                                    stats_zm(i) )             ! intent(inout)
         endif ! l_predict_upwp_vpwp
         call stat_begin_update( nzm, stats_metadata%irtp2_bt, rtp2(i,:) / dt, & ! intent(in)
                                 stats_zm(i) )             ! intent(inout)
         call stat_begin_update( nzm, stats_metadata%ithlp2_bt, thlp2(i,:) / dt, & ! intent(in)
                                 stats_zm(i) )               ! intent(inout)
         call stat_begin_update( nzm, stats_metadata%irtpthlp_bt, rtpthlp(i,:) / dt, & ! intent(in)
                                 stats_zm(i) )                   ! intent(inout)

         call stat_begin_update( nzt, stats_metadata%irtm_bt, rtm(i,:) / dt, & ! intent(in)
                                 stats_zt(i) )           ! intent(inout)
         call stat_begin_update( nzt, stats_metadata%ithlm_bt, thlm(i,:) / dt, & ! intent(in)
                                 stats_zt(i) )             ! intent(inout)
         call stat_begin_update( nzt, stats_metadata%ium_bt, um(i,:) / dt, & ! intent(in)
                                 stats_zt(i) )         ! intent(inout)
         call stat_begin_update( nzt, stats_metadata%ivm_bt, vm(i,:) / dt, & ! intent(in)
                                 stats_zt(i) )         ! intent(inout)
         call stat_begin_update( nzt, stats_metadata%iwp3_bt, wp3(i,:) / dt, & ! intent(in)
                                 stats_zt(i) )           ! intent(inout)

      end do

    end if

    ! SET SURFACE VALUES OF FLUXES (BROUGHT IN)
    ! We only do this for host models that do not apply the flux
    ! elsewhere in the code (e.g. WRF).  In other cases the _sfc variables will
    ! only be used to compute the variance at the surface. -dschanen 8 Sept 2009
    if ( .not. clubb_config_flags%l_host_applies_sfc_fluxes ) then

      !$acc parallel loop gang vector default(present)
      do i = 1, ngrdcol
        wpthlp(i,gr%k_lb_zm) = wpthlp_sfc(i)
        wprtp(i,gr%k_lb_zm)  = wprtp_sfc(i)
        upwp(i,gr%k_lb_zm)   = upwp_sfc(i)
        vpwp(i,gr%k_lb_zm)   = vpwp_sfc(i)
      end do
      !$acc end parallel loop

      if ( clubb_config_flags%l_linearize_pbl_winds ) then
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
        ! timestep to timestep.  Set it to a value of 0 and then
        ! overwrite the value at the surface.
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
        ! timestep to timestep.  Set it to a value of 0 and then
        ! overwrite the value at the surface.
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

#ifdef CLUBBND_CAM
    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      mu(i) = varmu(i)
    end do
    !$acc end parallel loop
#else
    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      mu(i) = clubb_params(i,imu)
    end do
    !$acc end parallel loop
#endif

    if ( clubb_config_flags%ipdf_call_placement == ipdf_pre_advance_fields &
         .or. clubb_config_flags%ipdf_call_placement &
              == ipdf_pre_post_advance_fields ) then

      ! Sample stats in this call to subroutine pdf_closure_driver for
      ! both of these options (ipdf_pre_advance_fields and
      ! ipdf_pre_post_advance_fields).
      if ( clubb_config_flags%ipdf_call_placement &
           == ipdf_pre_advance_fields ) then
        l_samp_stats_in_pdf_call = .true.
      elseif ( clubb_config_flags%ipdf_call_placement &
               == ipdf_pre_post_advance_fields ) then
        l_samp_stats_in_pdf_call = .true.
      end if

      !########################################################################
      !#######                     CALL CLUBB's PDF                     #######
      !#######   AND OUTPUT PDF PARAMETERS AND INTEGRATED QUANTITITES   #######
      !########################################################################
      call pdf_closure_driver( gr, nzm, nzt, ngrdcol,                       & ! Intent(in)
                               dt, hydromet_dim, sclr_dim, sclr_tol,        & ! Intent(in)
                               wprtp, thlm, wpthlp, rtp2, rtp3,             & ! Intent(in)
                               thlp2, thlp3, rtpthlp, wp2,                  & ! Intent(in)
                               wp3, wm_zm, wm_zt,                           & ! Intent(in)
                               um, up2, upwp, up3,                          & ! Intent(in)
                               vm, vp2, vpwp, vp3,                          & ! Intent(in)
                               p_in_Pa, exner,                              & ! Intent(in)
                               thv_ds_zm, thv_ds_zt, rtm_ref,               & ! Intent(in)
                               wphydrometp,                                 & ! Intent(in)
                               wp2hmp, rtphmp_zt, thlphmp_zt,               & ! Intent(in)
                               sclrm, wpsclrp, sclrp2,                      & ! Intent(in)
                               sclrprtp, sclrpthlp, sclrp3,                 & ! Intent(in)
                               p_sfc, l_samp_stats_in_pdf_call,             & ! Intent(in)
                               mixt_frac_max_mag, ts_nudge,                 & ! Intent(in)
                               rtm_min, rtm_nudge_max_altitude,             & ! Intent(in)
                               clubb_params,                                & ! Intent(in)
                               clubb_config_flags%iiPDF_type,               & ! Intent(in)
                               clubb_config_flags%saturation_formula,       & ! Intent(in)
                               clubb_config_flags%l_predict_upwp_vpwp,      & ! Intent(in)
                               clubb_config_flags%l_rtm_nudge,              & ! Intent(in)
                               clubb_config_flags%l_trapezoidal_rule_zt,    & ! Intent(in)
                               clubb_config_flags%l_trapezoidal_rule_zm,    & ! Intent(in)
                               clubb_config_flags%l_call_pdf_closure_twice, & ! Intent(in)
                               clubb_config_flags%l_use_cloud_cover,        & ! Intent(in)
                               clubb_config_flags%l_rcm_supersat_adj,       & ! Intent(in)
                               l_mix_rat_hm,                                & ! Intent(in)
                               stats_metadata,                              & ! Intent(in)
                               stats_zt, stats_zm,                          & ! Intent(inout)
                               rtm,                                         & ! Intent(inout)
                               pdf_implicit_coefs_terms,                    & ! Intent(inout)
                               pdf_params, pdf_params_zm, err_info,         & ! Intent(inout)
#ifdef GFDL
                               RH_crit(k, : , :),                           & ! Intent(inout)
                               do_liquid_only_in_clubb,                     & ! Intent(in)
#endif
                               rcm, cloud_frac,                             & ! Intent(out)
                               ice_supersat_frac, wprcp,                    & ! Intent(out)
                               sigma_sqd_w, wpthvp, wp2thvp, wp2up,         & ! Intent(out)
                               rtpthvp, thlpthvp, rc_coef,                  & ! Intent(out)
                               rcm_in_layer, cloud_cover,                   & ! Intent(out)
                               rcp2_zt, thlprcp,                            & ! Intent(out)
                               rc_coef_zm, sclrpthvp,                       & ! Intent(out)
                               wpup2, wpvp2,                                & ! Intent(out)
                               wp2up2, wp2vp2, wp4,                         & ! Intent(out)
                               wp2rtp, wprtp2, wp2thlp,                     & ! Intent(out)
                               wpthlp2, wprtpthlp, wp2rcp,                  & ! Intent(out)
                               rtprcp, rcp2,                                & ! Intent(out)
                               uprcp, vprcp,                                & ! Intent(out)
                               w_up_in_cloud, w_down_in_cloud,              & ! Intent(out)
                               cloudy_updraft_frac,                         & ! Intent(out)
                               cloudy_downdraft_frac,                       & ! intent(out)
                               Skw_velocity,                                & ! Intent(out)
                               cloud_frac_zm,                               & ! Intent(out)
                               ice_supersat_frac_zm,                        & ! Intent(out)
                               rtm_zm, thlm_zm, rcm_zm,                     & ! Intent(out)
                               rcm_supersat_adj,                            & ! Intent(out)
                               wp2sclrp, wpsclrp2, sclrprcp,                & ! Intent(out)
                               wpsclrprtp, wpsclrpthlp )                      ! Intent(out)

      if ( clubb_at_least_debug_level_api( 0 ) ) then
        if ( any(err_info%err_code == clubb_fatal_error) ) then
          write(fstderr,*) err_info%err_header_global
          write(fstderr,*) "Error calling pdf_closure_driver in advance_clubb_core"
          return
        endif
      endif
    endif ! clubb_config_flags%ipdf_call_placement == ipdf_pre_advance_fields
          ! or clubb_config_flags%ipdf_call_placement
          !    == ipdf_pre_post_advance_fields

    ! Interpolate wp3 to momentum levels, and wp2 to thermodynamic levels
    ! and then compute Skw for m & t grid.
    ! Positive definite quantity
    wp2_zt(:,:) = zm2zt_api( nzm, nzt, ngrdcol, gr, wp2(:,:), w_tol_sqd )
    wp3_zm(:,:) = zt2zm_api( nzm, nzt, ngrdcol, gr, wp3(:,:) )

    call Skx_func( nzt, ngrdcol, wp2_zt, wp3, &
                   w_tol, clubb_params, &
                   Skw_zt )
                   
    call Skx_func( nzm, ngrdcol, wp2, wp3_zm, &
                   w_tol, clubb_params, &
                   Skw_zm )
   
    if ( clubb_config_flags%ipdf_call_placement == ipdf_post_advance_fields ) then

      ! Calculate sigma_sqd_w here in order to avoid having to pass it in
      ! and out of subroutine advance_clubb_core.
      if ( l_gamma_Skw ) then

        !$acc parallel loop gang vector collapse(2) default(present)
        do k = 1, nzm
          do i = 1, ngrdcol

            gamma_coef = clubb_params(i,igamma_coef)
            gamma_coefb = clubb_params(i,igamma_coefb)
            gamma_coefc = clubb_params(i,igamma_coefc)

            if ( abs(gamma_coef-gamma_coefb) > abs(gamma_coef+gamma_coefb)*eps/2) then
            
              gamma_Skw_fnc(i,k) = gamma_coefb + (gamma_coef-gamma_coefb) &
                    *exp( -(1.0_core_rknd/2.0_core_rknd) * (Skw_zm(i,k)/gamma_coefc)**2 )
            else
              gamma_Skw_fnc(i,k) = gamma_coef
            end if

          end do
        end do
        !$acc end parallel loop

      else

        !$acc parallel loop gang vector collapse(2) default(present)
        do k = 1, nzm
          do i = 1, ngrdcol
            gamma_Skw_fnc(i,k) = clubb_params(i,igamma_coef)
          end do
        end do
        !$acc end parallel loop

      end if

      ! Compute sigma_sqd_w (dimensionless PDF width parameter)
      call compute_sigma_sqd_w( nzm, ngrdcol, &
                                gamma_Skw_fnc, wp2, thlp2, rtp2, &
                                up2, vp2, wpthlp, wprtp, upwp, vpwp, &
                                clubb_config_flags%l_predict_upwp_vpwp, &
                                sigma_sqd_w_tmp )

      ! Smooth in the vertical using interpolation
      sigma_sqd_w(:,:) = zm2zt2zm( nzm, nzt, ngrdcol, gr, sigma_sqd_w_tmp(:,:), &
                                   zero_threshold )

    endif ! clubb_config_flags%ipdf_call_placement == ipdf_post_advance_fields


    ! Compute the a3 coefficient (formula 25 in `Equations for CLUBB')
    ! Note:  a3 has been modified because the wp3 turbulent advection term is
    !        now discretized on its own.  This removes the "- 3" from the end.
!   a3_coef = 3.0_core_rknd * sigma_sqd_w*sigma_sqd_w  &
!      + 6.0_core_rknd*(1.0_core_rknd-sigma_sqd_w)*sigma_sqd_w  &
!      + (1.0_core_rknd-sigma_sqd_w)*(1.0_core_rknd-sigma_sqd_w)

    ! This is a simplified version of the formula above.
    ! Note:  a3 has been modified because the wp3 turbulent advection term is
    !        now discretized on its own.
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzm
      do i = 1, ngrdcol
        a3_coef(i,k) = -2._core_rknd * ( 1._core_rknd - sigma_sqd_w(i,k) )**2 + 3.0_core_rknd
      end do
    end do
    !$acc end parallel loop

    ! We found we obtain fewer spikes in wp3 when we clip a3 to be no greater
    ! than -1.4 -dschanen 4 Jan 2011
    !a3_coef = max( a3_coef, -1.4_core_rknd ) ! Known magic number
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzm
      do i = 1, ngrdcol
        a3_coef(i,k) = max( a3_coef(i,k), clubb_params(i,ia3_coef_min) )
      end do
    end do
    !$acc end parallel loop

    a3_coef_zt(:,:) = zm2zt_api( nzm, nzt, ngrdcol, gr, a3_coef(:,:) )

    ! Interpolate thlp2, rtp2, and rtpthlp to thermodynamic levels.
    thlp2_zt(:,:)   = zm2zt_api( nzm, nzt, ngrdcol, gr, thlp2(:,:), &
                                 thl_tol**2 )  ! Positive def. quantity
    rtp2_zt(:,:)    = zm2zt_api( nzm, nzt, ngrdcol, gr, rtp2(:,:), &
                                 rt_tol**2 )   ! Positive def. quantity
    rtpthlp_zt(:,:) = zm2zt_api( nzm, nzt, ngrdcol, gr, rtpthlp(:,:) )

    ! Compute wp3 / wp2 on zt levels.  Always use the interpolated value in the
    ! denominator since it's less likely to create spikes
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzt
      do i = 1, ngrdcol
        wp3_on_wp2_zt(i,k) = ( wp3(i,k) / max( wp2_zt(i,k), w_tol_sqd ) )
      end do
    end do
    !$acc end parallel loop

    ! Clip wp3_on_wp2_zt if it's too large
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzt
      do i = 1, ngrdcol
        if( wp3_on_wp2_zt(i,k) < 0._core_rknd ) then
          wp3_on_wp2_zt(i,k) = max( -1000._core_rknd, wp3_on_wp2_zt(i,k) )
        else
          wp3_on_wp2_zt(i,k) = min( 1000._core_rknd, wp3_on_wp2_zt(i,k) )
        end if
      end do
    end do
    !$acc end parallel loop

    ! Compute wp3_on_wp2 by interpolating wp3_on_wp2_zt
    wp3_on_wp2(:,:) = zt2zm_api( nzm, nzt, ngrdcol, gr, wp3_on_wp2_zt(:,:) )

    ! Smooth again as above
    wp3_on_wp2_zt(:,:) = zm2zt_api( nzm, nzt, ngrdcol, gr, wp3_on_wp2(:,:) )

    !----------------------------------------------------------------
    ! Compute thvm
    !----------------------------------------------------------------
    call calculate_thvm( nzt, ngrdcol, &
                         thlm, rtm, rcm, exner, thv_ds_zt, &
                         thvm )

    !----------------------------------------------------------------
    ! Compute tke (turbulent kinetic energy)
    !----------------------------------------------------------------
    if ( .not. clubb_config_flags%l_tke_aniso ) then
      ! tke is assumed to be 3/2 of wp2
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

    sqrt_em_zt(:,:) = zm2zt_api( nzm, nzt, ngrdcol, gr, em(:,:), em_min )

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzt
      do i = 1, ngrdcol
        sqrt_em_zt(i,k) = sqrt( sqrt_em_zt(i,k) )
      end do
    end do
    !$acc end parallel loop

    !----------------------------------------------------------------
    ! Compute mixing length and dissipation time
    !----------------------------------------------------------------

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
                                      brunt_vaisala_freq_sqd_dry,                         & ! Out
                                      brunt_vaisala_freq_sqd_moist,                       & ! Out
                                      brunt_vaisala_freq_sqd_smth )                         ! Out

    ! Calculate the norm of the vertical derivative of the mean horizontal wind speed
    ! To feed into the calculation of the Richardson number Ri_zm
    ddzt_um = ddzt( nzm, nzt, ngrdcol, gr, um )
    ddzt_vm = ddzt( nzm, nzt, ngrdcol, gr, vm )

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzm
      do i = 1, ngrdcol
        ddzt_umvm_sqd(i,k) = ddzt_um(i,k)**2 + ddzt_vm(i,k)**2
      end do
    end do
    !$acc end parallel loop

    if ( stats_metadata%l_stats_samp ) then
      !$acc update host(ddzt_umvm_sqd)
      do i = 1, ngrdcol
        call stat_update_var( stats_metadata%iddzt_umvm_sqd, ddzt_umvm_sqd(i,:), & ! intent(in)
                              stats_zm(i) )               ! intent(inout)
      end do
    end if

    ! Calculate Richardson number Ri_zm
    if ( clubb_config_flags%l_modify_limiters_for_cnvg_test ) then

      call calc_Ri_zm(nzm, ngrdcol, brunt_vaisala_freq_sqd_smth, ddzt_umvm_sqd, &
                      0.0_core_rknd, 1.0e-12_core_rknd, Ri_zm )

      Ri_zm = zm2zt2zm( nzm, nzt, ngrdcol, gr, &
                        Ri_zm )

    else ! default method

      if ( l_smooth_min_max ) then

        brunt_vaisala_freq_clipped = smooth_max( nzm, ngrdcol, 1.0e-7_core_rknd, &
                                                 brunt_vaisala_freq_sqd_smth, &
                                                 1.0e-4_core_rknd * min_max_smth_mag )

        ddzt_umvm_sqd_clipped = smooth_max( nzm, ngrdcol, ddzt_umvm_sqd, &
                                        1.0e-7_core_rknd, &
                                        1.0e-6_core_rknd * min_max_smth_mag )

        call calc_Ri_zm(nzm, ngrdcol, brunt_vaisala_freq_clipped, ddzt_umvm_sqd_clipped, &
                        0.0_core_rknd, 0.0_core_rknd, Ri_zm )

      else

        call calc_Ri_zm(nzm, ngrdcol, brunt_vaisala_freq_sqd_smth, ddzt_umvm_sqd, &
                        1.0e-7_core_rknd, 1.0e-7_core_rknd, Ri_zm )

      end if

    end if

    if ( .not. clubb_config_flags%l_diag_Lscale_from_tau ) then ! compute Lscale 1st, using
                                                                ! buoyant parcel calc
      call calc_Lscale_directly ( ngrdcol, nzm, nzt, gr,                             & ! In
                                  l_implemented, p_in_Pa,                            & ! In
                                  exner, rtm, thlm, thvm,                            & ! In
                                  mu, rtp2_zt, thlp2_zt, rtpthlp_zt, pdf_params, em, & ! In
                                  thv_ds_zt, Lscale_max, lmin,                       & ! In
                                  clubb_params,                                      & ! In
                                  clubb_config_flags%saturation_formula,             & ! In
                                  clubb_config_flags%l_Lscale_plume_centered,        & ! In
                                  stats_metadata,                                    & ! In
                                  stats_zt, err_info,                                & ! In/Out
                                  Lscale, Lscale_up, Lscale_down )                     ! Out

      if ( clubb_at_least_debug_level_api( 0 ) ) then
        if ( any(err_info%err_code == clubb_fatal_error) ) then
          write(fstderr,*) err_info%err_header_global
          write(fstderr,*) "Error calling calc_Lscale_directly in advance_clubb_core"
          return
        end if
      end if

      ! Calculate CLUBB's turbulent eddy-turnover time scale as
      !   CLUBB's length scale divided by a velocity scale.
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzt
        do i = 1, ngrdcol
          tau_zt(i,k) = min( Lscale(i,k) / sqrt_em_zt(i,k), clubb_params(i,itaumax) )
        end do
      end do
      !$acc end parallel loop

      Lscale_zm(:,:) = zt2zm_api( nzm, nzt, ngrdcol, gr, Lscale(:,:), zero_threshold )
          
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm
        do i = 1, ngrdcol
          tau_zm(i,k) = min( Lscale_zm(i,k) / sqrt( max( em_min, em(i,k) ) ),  &
                             clubb_params(i,itaumax) )
        end do
      end do
      !$acc end parallel loop

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm
        do i = 1, ngrdcol
          invrs_tau_zm(i,k)      = one / tau_zm(i,k)
          invrs_tau_wp2_zm(i,k)  = invrs_tau_zm(i,k)
          invrs_tau_xp2_zm(i,k)  = invrs_tau_zm(i,k)
          invrs_tau_wpxp_zm(i,k) = invrs_tau_zm(i,k)
          invrs_tau_wp3_zm(i,k)  = invrs_tau_zm(i,k)
          tau_max_zm(i,k) = clubb_params(i,itaumax)
        end do
      end do
      !$acc end parallel loop

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzt
        do i = 1, ngrdcol
          invrs_tau_zt(i,k)     = one / tau_zt(i,k)
          invrs_tau_wp3_zt(i,k) = invrs_tau_zt(i,k)
          tau_max_zt(i,k) = clubb_params(i,itaumax)
        end do
      end do
      !$acc end parallel loop

      ! End Vince Larson's replacement.


    else ! l_diag_Lscale_from_tau = .true., diagnose simple tau and Lscale.

      call diagnose_Lscale_from_tau( nzm, nzt, ngrdcol, gr,                       & ! In
                        upwp_sfc, vpwp_sfc, ddzt_umvm_sqd,                        & ! In
                        ice_supersat_frac,                                        & ! In
                        em, sqrt_em_zt,                                           & ! In
                        ufmin, tau_const,                                         & ! In
                        sfc_elevation, Lscale_max,                                & ! In
                        clubb_params,                                             & ! In
                        stats_metadata,                                           & ! In
                        clubb_config_flags%l_e3sm_config,                         & ! In
                        clubb_config_flags%l_smooth_Heaviside_tau_wpxp,           & ! In
                        brunt_vaisala_freq_sqd_smth, Ri_zm,                       & ! In
                        stats_zm, err_info,                                       & ! Inout
                        invrs_tau_zt, invrs_tau_zm,                               & ! Out
                        invrs_tau_sfc, invrs_tau_no_N2_zm, invrs_tau_bkgnd,       & ! Out
                        invrs_tau_shear, invrs_tau_N2_iso,                        & ! Out
                        invrs_tau_wp2_zm, invrs_tau_xp2_zm,                       & ! Out
                        invrs_tau_wp3_zm, invrs_tau_wp3_zt, invrs_tau_wpxp_zm,    & ! Out
                        tau_max_zm, tau_max_zt, tau_zm, tau_zt,                   & ! Out
                        Lscale, Lscale_up, Lscale_down )                            ! Out

      if ( clubb_at_least_debug_level_api( 0 ) ) then
        if ( any(err_info%err_code == clubb_fatal_error) ) then
          write(fstderr,*) err_info%err_header_global
          write(fstderr, *) "Error calling diagnose_Lscale_from_tau in advance_clubb_core"
          return
        end if
      end if

    end if ! l_diag_Lscale_from_tau



        ! Modification to damp noise in stable region
  ! Vince Larson commented out because it may prevent turbulence from
  !    initiating in unstable regions.  7 Jul 2007
  !       do k = 1, nz
  !         if ( wp2(k) <= 0.005_core_rknd ) then
  !           tau_zt(k) = taumin
  !           tau_zm(k) = taumin
  !         end if
  !       end do
  ! End Vince Larson's commenting.

    !----------------------------------------------------------------
    ! Eddy diffusivity coefficient
    !----------------------------------------------------------------
    ! c_K is 0.548 usually (Duynkerke and Driedonks 1987)
    ! CLUBB uses a smaller value to better fit empirical data.

    ! Calculate CLUBB's eddy diffusivity as
    !   CLUBB's length scale times a velocity scale.
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzt
      do i = 1, ngrdcol
        Kh_zt(i,k) = clubb_params(i,ic_K) * Lscale(i,k) * sqrt_em_zt(i,k)
      end do
    end do
    !$acc end parallel loop

    Lscale_zm(:,:) = zt2zm_api( nzm, nzt, ngrdcol, gr, Lscale(:,:) )

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzm
      do i = 1, ngrdcol
        Kh_zm(i,k) = clubb_params(i,ic_K) * max( Lscale_zm(i,k), zero_threshold )  &
                     * sqrt( max( em(i,k), em_min ) )
      end do
    end do
    !$acc end parallel loop

    ! calculate Brunt-Vaisala frequency used for splatting
    brunt_vaisala_freq_sqd_splat  &
               = Lscale_width_vert_avg( nzm, ngrdcol, gr, smth_type, &
                                        brunt_vaisala_freq_sqd_mixed, Lscale_zm, rho_ds_zm, &
                                        below_grnd_val )

    ! Vertical compression of eddies causes gustiness (increase in up2 and vp2)
    call wp2_term_splat_lhs( nzm, nzt, ngrdcol, gr, clubb_params(:,iC_wp2_splat), & ! Intent(in)
                             brunt_vaisala_freq_sqd_splat,                        & ! Intent(in)
                             lhs_splat_wp2 )                                        ! Intent(out)

    ! Vertical compression of eddies also diminishes w'3
    call wp3_term_splat_lhs( nzm, nzt, ngrdcol, gr, clubb_params(:,iC_wp2_splat), & ! Intent(in)
                             brunt_vaisala_freq_sqd_splat,                        & ! Intent(in)
                             lhs_splat_wp3 )                                        ! Intent(out)

    !----------------------------------------------------------------
    ! Set Surface variances
    !----------------------------------------------------------------
    ! Surface variances should be set here, before the call to either
    ! advance_xp2_xpyp or advance_wp2_wp3.
    ! Surface effects should not be included with any case where the lowest
    ! level is not the ground level.  Brian Griffin.  December 22, 2005.

    ! Diagnose surface variances based on surface fluxes.
    call calc_sfc_varnce( nzm, nzt, ngrdcol, sclr_dim, sclr_idx,          & ! Intent(in)
                          gr, dt, sfc_elevation,                          & ! Intent(in)
                          upwp_sfc, vpwp_sfc, wpthlp, wprtp_sfc,          & ! Intent(in)
                          um, vm, Lscale_up, wpsclrp_sfc,                 & ! Intent(in)
                          lhs_splat_wp2, tau_zm,                          & ! Intent(in)
                          !wp2_splat, tau_zm,                             & ! Intent(in)
                          clubb_config_flags%l_vary_convect_depth,        & ! Intent(in)
                          T0,                                             & ! Intent(in)
                          clubb_params(:,iup2_sfc_coef),                  & ! Intent(in)
                          clubb_params(:,ia_const),                       & ! Intent(in)
                          stats_metadata,                                 & ! Intent(in)
                          stats_zm,                                       & ! Intent(inout)
                          wp2, up2, vp2,                                  & ! Intent(inout)
                          thlp2, rtp2, rtpthlp,                           & ! Intent(inout)
                          sclrp2, sclrprtp, sclrpthlp,                    & ! Intent(inout)
                          err_info )                                        ! Intent(inout)

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

    if ( stats_metadata%l_stats_samp ) then

      !$acc update host( rtm, rcm, thlm, exner, p_in_Pa )

      do i = 1, ngrdcol
        call stat_update_var( stats_metadata%irvm, rtm(i,:) - rcm(i,:), & !intent(in)
                              stats_zt(i) )               !intent(inout)
      end do

      ! Output relative humidity (q/q∗ where q∗ is the saturation mixing ratio over liquid)
      ! Added an extra check for stats_metadata%irel_humidity > 0; otherwise, 
      ! if both stats_metadata%irsat = 0 and
      ! stats_metadata%irel_humidity = 0, rsat is not computed, leading to a 
      ! floating-point exception
      ! when stat_update_var is called for rel_humidity.  ldgrant
      if ( stats_metadata%irel_humidity > 0 ) then
        
        rsat = sat_mixrat_liq_api( nzt, ngrdcol, gr, p_in_Pa, &
                                   thlm2T_in_K_api( nzt, ngrdcol, thlm, exner, rcm ), &
                                   clubb_config_flags%saturation_formula )

        ! Recompute rsat and rel_humidity. They might have changed.
        do i = 1, ngrdcol
          rel_humidity(i,:) = (rtm(i,:) - rcm(i,:)) / rsat(i,:)

          call stat_update_var( stats_metadata%irel_humidity, rel_humidity(i,:), & ! intent(in)
                                stats_zt(i))                                       ! intent(inout)
        end do
      end if ! stats_metadata%irel_humidity > 0
    end if ! stats_metadata%l_stats_samp

    !----------------------------------------------------------------
    ! Advance rtm/wprtp and thlm/wpthlp one time step
    !----------------------------------------------------------------
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

    ! Here we determine if we're using tau_zm or tau_N2_zm, which is tau
    ! that has been stability corrected for stably stratified regions.
    ! -dschanen 7 Nov 2014
    if ( clubb_config_flags%l_stability_correct_tau_zm ) then

      ! Determine stability correction factor
      call calc_stability_correction( nzm, nzt, ngrdcol, gr,                              & ! In
                                      thlm, Lscale_zm, em,                                & ! In
                                      exner, rtm, rcm,                                    & ! In
                                      p_in_Pa, thvm, ice_supersat_frac,                   & ! In
                                      clubb_params(:,ilambda0_stability_coef),            & ! In
                                      clubb_params(:,ibv_efold), T0,                      & ! In
                                      clubb_config_flags%saturation_formula,              & ! In
                                      clubb_config_flags%l_brunt_vaisala_freq_moist,      & ! In
                                      clubb_config_flags%l_use_thvm_in_bv_freq,           & ! In
                                      clubb_config_flags%l_modify_limiters_for_cnvg_test, & ! In
                                      stability_correction )                                ! Out

      if ( stats_metadata%l_stats_samp ) then
        !$acc update host( stability_correction )
        do i = 1, ngrdcol
          call stat_update_var( stats_metadata%istability_correction, & ! In
                                stability_correction(i,:), & ! In
                                stats_zm(i) ) ! In/Out
        end do
      end if

      ! Determine the static stability corrected version of tau_zm
      ! Create a damping time scale that is more strongly damped at the
      ! altitudes where the Brunt-Vaisala frequency (N^2) is large.
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm
        do i = 1, ngrdcol
          invrs_tau_N2_zm(i,k) = invrs_tau_zm(i,k) * stability_correction(i,k)
          invrs_tau_C6_zm(i,k) = invrs_tau_N2_zm(i,k)
          invrs_tau_C1_zm(i,k) = invrs_tau_N2_zm(i,k)
        end do
      end do
      !$acc end parallel loop
    else
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm
        do i = 1, ngrdcol
          invrs_tau_N2_zm(i,k) = unused_var
          invrs_tau_C6_zm(i,k) = invrs_tau_wpxp_zm(i,k)
          invrs_tau_C1_zm(i,k) = invrs_tau_wp2_zm(i,k)
        end do
      end do
      !$acc end parallel loop
    end if ! l_stability_correction

    ! Set invrs_tau variables for C4 and C14
      !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzm
      do i = 1, ngrdcol
        invrs_tau_C14_zm(i,k) = invrs_tau_wp2_zm(i,k)
      end do
    end do
    !$acc end parallel loop

    if ( .not. clubb_config_flags%l_diag_Lscale_from_tau .and. l_use_invrs_tau_N2_iso) then
      write(fstderr,*) "Error! l_use_invrs_tau_N2_iso is not used when "// &
                       "l_diag_Lscale_from_tau=false."// &
                       "If you want to use Lscale code, go to file "// &
                       "src/CLUBB_core/advance_clubb_core_module.F90 and "// &
                       "change l_use_invrs_tau_N2_iso to false"
      error stop
    end if

    if ( .not. l_use_invrs_tau_N2_iso ) then
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm
        do i = 1, ngrdcol
          invrs_tau_C4_zm(i,k) = invrs_tau_wp2_zm(i,k)
        end do
      end do
      !$acc end parallel loop
    else
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm
        do i = 1, ngrdcol
          invrs_tau_C4_zm(i,k) = invrs_tau_N2_iso(i,k)
        end do
      end do
      !$acc end parallel loop
    end if

    if ( stats_metadata%l_stats_samp ) then

      !$acc update host( invrs_tau_zm, invrs_tau_xp2_zm, invrs_tau_wp2_zm, invrs_tau_wpxp_zm, &
      !$acc              Ri_zm, invrs_tau_wp3_zm, invrs_tau_no_N2_zm, invrs_tau_bkgnd, &
      !$acc              invrs_tau_sfc, invrs_tau_shear, brunt_vaisala_freq_sqd, &
      !$acc              brunt_vaisala_freq_sqd_splat, brunt_vaisala_freq_sqd_mixed, &
      !$acc              brunt_vaisala_freq_sqd_moist, brunt_vaisala_freq_sqd_dry, &
      !$acc              brunt_vaisala_freq_sqd_smth )

      do i = 1, ngrdcol

        call stat_update_var(stats_metadata%iinvrs_tau_zm, & ! intent(in)
                             invrs_tau_zm(i,:), & ! intent(in)
                             stats_zm(i))                      ! intent(inout)
        call stat_update_var(stats_metadata%iinvrs_tau_xp2_zm, & ! intent(in)
                             invrs_tau_xp2_zm(i,:), & ! intent(in)
                             stats_zm(i))                              ! intent(inout)
        call stat_update_var(stats_metadata%iinvrs_tau_wp2_zm, & ! intent(in)
                             invrs_tau_wp2_zm(i,:), & ! intent(in)
                             stats_zm(i))                              ! intent(inout)
        call stat_update_var(stats_metadata%iinvrs_tau_wpxp_zm, & ! intent(in)
                             invrs_tau_wpxp_zm(i,:), & ! intent(in)
                             stats_zm(i))                                ! intent(inout)
        call stat_update_var(stats_metadata%iRi_zm, & ! intent(in)
                             Ri_zm(i,:), & ! intent(in)
                             stats_zm(i))                  ! intent(inout)
        call stat_update_var(stats_metadata%iinvrs_tau_wp3_zm, & ! intent(in)
                             invrs_tau_wp3_zm(i,:), & ! intent(in)
                             stats_zm(i))                                ! intent(inout)

        if ( clubb_config_flags%l_diag_Lscale_from_tau ) then
          call stat_update_var(stats_metadata%iinvrs_tau_no_N2_zm, & ! intent(in)
                               invrs_tau_no_N2_zm(i,:), & ! intent(in)
                               stats_zm(i))                                  ! intent(inout)
          call stat_update_var(stats_metadata%iinvrs_tau_bkgnd, & ! intent(in)
                               invrs_tau_bkgnd(i,:), & ! intent(in)
                               stats_zm(i))                            ! intent(inout)
          call stat_update_var(stats_metadata%iinvrs_tau_sfc, & ! intent(in)
                               invrs_tau_sfc(i,:), & ! intent(in)
                               stats_zm(i))                        ! intent(inout)
          call stat_update_var(stats_metadata%iinvrs_tau_shear, & ! intent(in)
                               invrs_tau_shear(i,:), & ! intent(in)
                               stats_zm(i))                            ! intent(inout)
        end if
        call stat_update_var(stats_metadata%ibrunt_vaisala_freq_sqd, & ! intent(in)
                             brunt_vaisala_freq_sqd(i,:), & ! intent(in)
                             stats_zm(i))
        call stat_update_var(stats_metadata%ibrunt_vaisala_freq_sqd_splat, & ! intent(in)
                             brunt_vaisala_freq_sqd_splat(i,:), & ! intent(in)
                             stats_zm(i))                                       ! intent(inout)
        call stat_update_var(stats_metadata%ibrunt_vaisala_freq_sqd_mixed, & ! intent(in)
                             brunt_vaisala_freq_sqd_mixed(i,:), & ! intent(in)
                             stats_zm(i))                                          ! intent(inout)
        call stat_update_var(stats_metadata%ibrunt_vaisala_freq_sqd_moist, & ! intent(in)
                             brunt_vaisala_freq_sqd_moist(i,:), & ! intent(in)
                             stats_zm(i))                                          ! intent(inout)
        call stat_update_var(stats_metadata%ibrunt_vaisala_freq_sqd_dry, & ! intent(in)
                             brunt_vaisala_freq_sqd_dry(i,:), & ! intent(in)
                             stats_zm(i))                                          ! intent(inout)
        call stat_update_var(stats_metadata%ibrunt_vaisala_freq_sqd_smth, & ! intent(in)
                             brunt_vaisala_freq_sqd_smth(i,:), & ! intent(in)
                             stats_zm(i))                                          ! intent(inout)

      end do
    end if

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

      ! Advance the prognostic equations for
      !   the scalar grid means (rtm, thlm, sclrm) and
      !   scalar turbulent fluxes (wprtp, wpthlp, and wpsclrp)
      !   by one time step.
      ! advance_xm_wpxp_bad_wp2 ! Test error comment, DO NOT modify or move
      call advance_xm_wpxp( nzm, nzt, ngrdcol, sclr_dim, sclr_tol, gr, dt_advance, & ! intent(in)  
                            sigma_sqd_w, wm_zm, wm_zt, wp2, Lscale_zm,             & ! intent(in)
                            wp3_on_wp2, wp3_on_wp2_zt, Kh_zt, Kh_zm,               & ! intent(in)
                            stability_correction,                                  & ! intent(in)
                            invrs_tau_C6_zm, tau_max_zm, Skw_zm, wp2rtp, rtpthvp,  & ! intent(in)
                            rtm_forcing, wprtp_forcing, rtm_ref, wp2thlp,          & ! intent(in)
                            thlpthvp, thlm_forcing, wpthlp_forcing, thlm_ref,      & ! intent(in)
                            rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm,                 & ! intent(in)
                            invrs_rho_ds_zt, thv_ds_zm, rtp2, thlp2,               & ! intent(in)
                            w_1_zm, w_2_zm, varnce_w_1_zm, varnce_w_2_zm,          & ! intent(in)
                            mixt_frac_zm, l_implemented, em, wp2sclrp,             & ! intent(in)
                            sclrpthvp, sclrm_forcing, sclrp2, Cx_fnc_Richardson,   & ! intent(in)
                            pdf_implicit_coefs_terms,                              & ! intent(in)
                            um_forcing, vm_forcing, ug, vg, wpthvp,                & ! intent(in)
                            fcor, fcor_y, um_ref, vm_ref, up2, vp2,                & ! intent(in)
                            uprcp, vprcp, rc_coef_zm,                              & ! intent(in)
                            clubb_params, nu_vert_res_dep, ts_nudge,               & ! intent(in)
                            clubb_config_flags%iiPDF_type,                         & ! intent(in)
                            clubb_config_flags%penta_solve_method,                 & ! intent(in)
                            clubb_config_flags%tridiag_solve_method,               & ! intent(in)
                            clubb_config_flags%fill_holes_type,                    & ! intent(in)
                            clubb_config_flags%l_predict_upwp_vpwp,                & ! intent(in)
                            clubb_config_flags%l_nontraditional_Coriolis,          & ! intent(in)
                            clubb_config_flags%l_traditional_Coriolis,             & ! intent(in)
                            clubb_config_flags%l_diffuse_rtm_and_thlm,             & ! intent(in)
                            clubb_config_flags%l_stability_correct_Kh_N2_zm,       & ! intent(in)
                            clubb_config_flags%l_godunov_upwind_wpxp_ta,           & ! intent(in)
                            clubb_config_flags%l_upwind_xm_ma,                     & ! intent(in)
                            clubb_config_flags%l_uv_nudge,                         & ! intent(in)
                            clubb_config_flags%l_tke_aniso,                        & ! intent(in)
                            clubb_config_flags%l_diag_Lscale_from_tau,             & ! intent(in)
                            clubb_config_flags%l_use_C7_Richardson,                & ! intent(in)
                            clubb_config_flags%l_lmm_stepping,                     & ! intent(in)
                            clubb_config_flags%l_enable_relaxed_clipping,          & ! intent(in)
                            clubb_config_flags%l_linearize_pbl_winds,              & ! intent(in)
                            clubb_config_flags%l_mono_flux_lim_thlm,               & ! intent(in)
                            clubb_config_flags%l_mono_flux_lim_rtm,                & ! intent(in)
                            clubb_config_flags%l_mono_flux_lim_um,                 & ! intent(in)
                            clubb_config_flags%l_mono_flux_lim_vm,                 & ! intent(in)
                            clubb_config_flags%l_mono_flux_lim_spikefix,           & ! intent(in)
                            order_xm_wpxp, order_xp2_xpyp, order_wp2_wp3,          & ! intent(in)
                            stats_metadata,                                        & ! intent(in)
                            stats_zt, stats_zm, stats_sfc,                         & ! intent(i/o)
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
                     rcm )                             ! In/Out

#ifdef GFDL
      do i = 1, ngrdcol
        call advance_sclrm_Nd_diffusion_OG( dt, &  ! h1g, 2012-06-16     ! In
                                            sclrm(i,:,:), sclrm_trsport_only(i,:,:), & ! In/Out
                                            Kh_zm(i,:),  cloud_frac(i,:) )         ! In
      end do
#endif

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
      call advance_xp2_xpyp( nzm, nzt, ngrdcol, sclr_dim, sclr_tol, gr, sclr_idx, & ! intent(in)
                             invrs_tau_xp2_zm, invrs_tau_C4_zm,                   & ! intent(in)
                             invrs_tau_C14_zm, wm_zm,                             & ! intent(in)
                             rtm, wprtp, thlm, wpthlp, wpthvp, um, vm,            & ! intent(in)
                             wp2, wp2_zt, wp3, upwp, vpwp,                        & ! intent(in)
                             sigma_sqd_w, wprtp2, wpthlp2,                        & ! intent(in)
                             wprtpthlp, Kh_zt, rtp2_forcing,                      & ! intent(in)
                             thlp2_forcing, rtpthlp_forcing,                      & ! intent(in)
                             rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm,               & ! intent(in)
                             thv_ds_zm, cloud_frac,                               & ! intent(in)
                             wp3_on_wp2, wp3_on_wp2_zt,                           & ! intent(in)
                             pdf_implicit_coefs_terms,                            & ! intent(in)
                             dt_advance, fcor_y,                                   & ! intent(in)
                             sclrm, wpsclrp,                                      & ! intent(in)
                             wpsclrp2, wpsclrprtp, wpsclrpthlp,                   & ! intent(in)
                             lhs_splat_wp2,                                       & ! intent(in)
                             clubb_params, nu_vert_res_dep,                       & ! intent(in)
                             clubb_config_flags%iiPDF_type,                       & ! intent(in)
                             clubb_config_flags%tridiag_solve_method,             & ! intent(in)
                             clubb_config_flags%fill_holes_type,                  & ! intent(in)
                             clubb_config_flags%l_predict_upwp_vpwp,              & ! intent(in)
                             clubb_config_flags%l_nontraditional_Coriolis,        & ! intent(in)
                             clubb_config_flags%l_min_xp2_from_corr_wx,           & ! intent(in)
                             clubb_config_flags%l_C2_cloud_frac,                  & ! intent(in)
                             clubb_config_flags%l_upwind_xpyp_ta,                 & ! intent(in)
                             clubb_config_flags%l_godunov_upwind_xpyp_ta,         & ! intent(in)
                             clubb_config_flags%l_lmm_stepping,                   & ! intent(in)
                             stats_metadata,                                      & ! intent(in)
                             stats_zt, stats_zm, stats_sfc,                       & ! intent(inout)
                             rtp2, thlp2, rtpthlp, up2, vp2,                      & ! intent(inout)
                             sclrp2, sclrprtp, sclrpthlp, err_info )                ! intent(inout)
      
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
      if ( order_xp2_xpyp < order_xm_wpxp &
           .and. order_xp2_xpyp < order_wp2_wp3 ) then
         wprtp_cl_num   = 1 ! First instance of w'r_t' clipping.
         wpthlp_cl_num  = 1 ! First instance of w'th_l' clipping.
         wpsclrp_cl_num = 1 ! First instance of w'sclr' clipping.
         if ( clubb_config_flags%l_predict_upwp_vpwp ) then
            upwp_cl_num = 1 ! First instance of u'w' clipping.
            vpwp_cl_num = 1 ! First instance of v'w' clipping.
         endif
      elseif ( order_xp2_xpyp > order_xm_wpxp &
               .and. order_xp2_xpyp > order_wp2_wp3 ) then
         wprtp_cl_num   = 3 ! Third instance of w'r_t' clipping.
         wpthlp_cl_num  = 3 ! Third instance of w'th_l' clipping.
         wpsclrp_cl_num = 3 ! Third instance of w'sclr' clipping.
         if ( clubb_config_flags%l_predict_upwp_vpwp ) then
            upwp_cl_num = 3 ! Third instance of u'w' clipping.
            vpwp_cl_num = 3 ! Third instance of v'w' clipping.
         endif
      else
         wprtp_cl_num   = 2 ! Second instance of w'r_t' clipping.
         wpthlp_cl_num  = 2 ! Second instance of w'th_l' clipping.
         wpsclrp_cl_num = 2 ! Second instance of w'sclr' clipping.
         if ( clubb_config_flags%l_predict_upwp_vpwp ) then
            upwp_cl_num = 2 ! Second instance of u'w' clipping.
            vpwp_cl_num = 2 ! Second instance of v'w' clipping.
         endif
      endif

      if ( .not. clubb_config_flags%l_predict_upwp_vpwp ) then
         if ( order_xp2_xpyp < order_wp2_wp3 &
              .and. order_xp2_xpyp < order_windm ) then
            upwp_cl_num = 1 ! First instance of u'w' clipping.
            vpwp_cl_num = 1 ! First instance of v'w' clipping.
         elseif ( order_xp2_xpyp > order_wp2_wp3 &
                  .and. order_xp2_xpyp > order_windm ) then
            upwp_cl_num = 3 ! Third instance of u'w' clipping.
            vpwp_cl_num = 3 ! Third instance of v'w' clipping.
         else
            upwp_cl_num = 2 ! Second instance of u'w' clipping.
            vpwp_cl_num = 2 ! Second instance of v'w' clipping.
         endif ! l_predict_upwp_vpwp
      endif

      call clip_covars_denom( nzm, ngrdcol, sclr_dim, dt,                   & ! intent(in)
                              rtp2, thlp2, up2, vp2, wp2,                   & ! intent(in)
                              sclrp2, wprtp_cl_num, wpthlp_cl_num,          & ! intent(in)
                              wpsclrp_cl_num, upwp_cl_num, vpwp_cl_num,     & ! intent(in)
                              clubb_config_flags%l_predict_upwp_vpwp,       & ! intent(in)
                              clubb_config_flags%l_tke_aniso,               & ! intent(in)
                              clubb_config_flags%l_linearize_pbl_winds,     & ! intent(in)
                              stats_metadata,                               & ! intent(in)
                              stats_zm,                                     & ! intent(inout)
                              wprtp, wpthlp, upwp, vpwp, wpsclrp,           & ! intent(inout)
                              upwp_pert, vpwp_pert )                          ! intent(inout)
      
     elseif ( advance_order_loop_iter == order_wp2_wp3 ) then

      !----------------------------------------------------------------
      ! Advance the 2nd- and 3rd-order moments
      !   of vertical velocity (wp2, wp3) by one timestep.
      !----------------------------------------------------------------

      ! advance_wp2_wp3_bad_wp2 ! Test error comment, DO NOT modify or move
      call advance_wp2_wp3( nzm, nzt, ngrdcol, gr, dt_advance,                    & ! intent(in)
                            sfc_elevation, fcor_y, sigma_sqd_w, wm_zm,            & ! intent(in)
                            wm_zt, a3_coef, a3_coef_zt, wp3_on_wp2,               & ! intent(in)
                            wpup2, wpvp2, wp2up2, wp2vp2, wp4,                    & ! intent(in)
                            wpthvp, wp2thvp, wp2up, um, vm, upwp, vpwp,           & ! intent(in)
                            em, Kh_zm, Kh_zt, invrs_tau_C4_zm,                    & ! intent(in)
                            invrs_tau_wp3_zt, invrs_tau_C1_zm, Skw_zm,            & ! intent(in)
                            Skw_zt, rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm,        & ! intent(in)
                            invrs_rho_ds_zt, thv_ds_zm,                           & ! intent(in)
                            thv_ds_zt, pdf_params%mixt_frac, Cx_fnc_Richardson,   & ! intent(in)
                            lhs_splat_wp2, lhs_splat_wp3,                         & ! intent(in)
                            pdf_implicit_coefs_terms,                             & ! intent(in)
                            wprtp, wpthlp, rtp2, thlp2,                           & ! intent(in)
                            clubb_params, nu_vert_res_dep,                        & ! intent(in)
                            clubb_config_flags%iiPDF_type,                        & ! intent(in)
                            clubb_config_flags%penta_solve_method,                & ! intent(in)
                            clubb_config_flags%fill_holes_type,                   & ! intent(in)
                            clubb_config_flags%l_min_wp2_from_corr_wx,            & ! intent(in)
                            clubb_config_flags%l_upwind_xm_ma,                    & ! intent(in)
                            clubb_config_flags%l_tke_aniso,                       & ! intent(in)
                            clubb_config_flags%l_standard_term_ta,                & ! intent(in)
                            clubb_config_flags%l_partial_upwind_wp3,              & ! intent(in)
                            clubb_config_flags%l_damp_wp2_using_em,               & ! intent(in)
                            clubb_config_flags%l_use_C11_Richardson,              & ! intent(in)
                            clubb_config_flags%l_damp_wp3_Skw_squared,            & ! intent(in)
                            clubb_config_flags%l_lmm_stepping,                    & ! intent(in)
                            clubb_config_flags%l_use_tke_in_wp3_pr_turb_term,     & ! intent(in)
                            clubb_config_flags%l_use_tke_in_wp2_wp3_K_dfsn,       & ! intent(in)
                            clubb_config_flags%l_use_wp3_lim_with_smth_Heaviside, & ! intent(in)
                            clubb_config_flags%l_wp2_fill_holes_tke,              & ! intent(in)
                            clubb_config_flags%l_nontraditional_Coriolis,         & ! intent(in)
                            stats_metadata,                                       & ! intent(in)
                            stats_zt, stats_zm, stats_sfc,                        & ! intent(inout)
                            up2, vp2, wp2, wp3, wp3_zm, wp2_zt, err_info )          ! intent(inout)

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

      if ( order_wp2_wp3 < order_xm_wpxp &
           .and. order_wp2_wp3 < order_xp2_xpyp ) then
         wprtp_cl_num   = 1 ! First instance of w'r_t' clipping.
         wpthlp_cl_num  = 1 ! First instance of w'th_l' clipping.
         wpsclrp_cl_num = 1 ! First instance of w'sclr' clipping.
         if ( clubb_config_flags%l_predict_upwp_vpwp ) then
            upwp_cl_num = 1 ! First instance of u'w' clipping.
            vpwp_cl_num = 1 ! First instance of v'w' clipping.
         endif
      elseif ( order_wp2_wp3 > order_xm_wpxp &
               .and. order_wp2_wp3 > order_xp2_xpyp ) then
         wprtp_cl_num   = 3 ! Third instance of w'r_t' clipping.
         wpthlp_cl_num  = 3 ! Third instance of w'th_l' clipping.
         wpsclrp_cl_num = 3 ! Third instance of w'sclr' clipping.
         if ( clubb_config_flags%l_predict_upwp_vpwp ) then
            upwp_cl_num = 3 ! Third instance of u'w' clipping.
            vpwp_cl_num = 3 ! Third instance of v'w' clipping.
         endif
      else
         wprtp_cl_num   = 2 ! Second instance of w'r_t' clipping.
         wpthlp_cl_num  = 2 ! Second instance of w'th_l' clipping.
         wpsclrp_cl_num = 2 ! Second instance of w'sclr' clipping.
         if ( clubb_config_flags%l_predict_upwp_vpwp ) then
            upwp_cl_num = 2 ! Second instance of u'w' clipping.
            vpwp_cl_num = 2 ! Second instance of v'w' clipping.
         endif
      endif
    
      if ( .not. clubb_config_flags%l_predict_upwp_vpwp ) then
         if ( order_wp2_wp3 < order_xp2_xpyp &
              .and. order_wp2_wp3 < order_windm ) then
            upwp_cl_num = 1 ! First instance of u'w' clipping.
            vpwp_cl_num = 1 ! First instance of v'w' clipping.
         elseif ( order_wp2_wp3 > order_xp2_xpyp &
                  .and. order_wp2_wp3 > order_windm ) then
            upwp_cl_num = 3 ! Third instance of u'w' clipping.
            vpwp_cl_num = 3 ! Third instance of v'w' clipping.
         else
            upwp_cl_num = 2 ! Second instance of u'w' clipping.
            vpwp_cl_num = 2 ! Second instance of v'w' clipping.
         endif ! l_predict_upwp_vpwp
      endif

      call clip_covars_denom( nzm, ngrdcol, sclr_dim, dt,                   & ! intent(in)
                              rtp2, thlp2, up2, vp2, wp2,                   & ! intent(in)
                              sclrp2, wprtp_cl_num, wpthlp_cl_num,          & ! intent(in)
                              wpsclrp_cl_num, upwp_cl_num, vpwp_cl_num,     & ! intent(in)
                              clubb_config_flags%l_predict_upwp_vpwp,       & ! intent(in)
                              clubb_config_flags%l_tke_aniso,               & ! intent(in)
                              clubb_config_flags%l_linearize_pbl_winds,     & ! intent(in)
                              stats_metadata,                               & ! intent(in)
                              stats_zm,                                     & ! intent(inout)
                              wprtp, wpthlp, upwp, vpwp, wpsclrp,           & ! intent(inout)
                              upwp_pert, vpwp_pert )                          ! intent(inout)

     elseif ( advance_order_loop_iter == order_windm ) then

      !----------------------------------------------------------------
      ! Advance the horizontal mean winds (um, vm),
      !   the mean of the eddy-diffusivity scalars (i.e. edsclrm),
      !   and their fluxes (upwp, vpwp, wpedsclrp) by one time step.
      !----------------------------------------------------------------
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm
        do i = 1, ngrdcol
          Km_zm(i,k) = Kh_zm(i,k) * clubb_params(i,ic_K10)   ! Coefficient for momentum

          Kmh_zm(i,k) = Kh_zm(i,k) * clubb_params(i,ic_K10h) ! Coefficient for thermo
        end do
      end do
      !$acc end parallel loop

      if ( edsclr_dim > 1 .and. clubb_config_flags%l_do_expldiff_rtm_thlm ) then
        !$acc parallel loop gang vector collapse(2) default(present)
        do k = 1, nzt
          do i = 1, ngrdcol
            edsclrm(i,k,edsclr_dim-1) = thlm(i,k)
            edsclrm(i,k,edsclr_dim) = rtm(i,k)
          end do
        end do
        !$acc end parallel loop
      end if

      call advance_windm_edsclrm( nzm, nzt, ngrdcol, edsclr_dim, gr, dt,      & ! intent(in)
                                  wm_zt, Km_zm, Kmh_zm,                       & ! intent(in)
                                  ug, vg, um_ref, vm_ref,                     & ! intent(in)
                                  wp2, up2, vp2, um_forcing, vm_forcing,      & ! intent(in)
                                  edsclrm_forcing,                            & ! intent(in)
                                  rho_ds_zm, invrs_rho_ds_zt,                 & ! intent(in)
                                  fcor, l_implemented,                        & ! intent(in)
                                  nu_vert_res_dep, ts_nudge,                  & ! intent(in)
                                  clubb_config_flags%tridiag_solve_method,    & ! intent(in)
                                  clubb_config_flags%l_predict_upwp_vpwp,     & ! intent(in)
                                  clubb_config_flags%l_upwind_xm_ma,          & ! intent(in)
                                  clubb_config_flags%l_uv_nudge,              & ! intent(in)
                                  clubb_config_flags%l_tke_aniso,             & ! intent(in)
                                  clubb_config_flags%l_lmm_stepping,          & ! intent(in)
                                  clubb_config_flags%l_linearize_pbl_winds,   & ! intent(in)
                                  order_xp2_xpyp, order_wp2_wp3, order_windm, & ! intent(in)
                                  stats_metadata,                             & ! intent(in)
                                  stats_zt, stats_zm, stats_sfc,              & ! intent(inout)
                                  um, vm, edsclrm,                            & ! intent(inout)
                                  upwp, vpwp, wpedsclrp,                      & ! intent(inout)
                                  um_pert, vm_pert, upwp_pert,                & ! intent(inout)
                                  vpwp_pert, err_info )                         ! intent(inout)

      if ( clubb_at_least_debug_level_api( 0 ) ) then
         if ( any(err_info%err_code == clubb_fatal_error) ) then
            write(fstderr,*) err_info%err_header_global
            write(fstderr,*) "Error calling advance_windm_edsclrm in advance_clubb_core"
            return
         end if
      end if

      if ( edsclr_dim > 1 .and. clubb_config_flags%l_do_expldiff_rtm_thlm ) then

        call pvertinterp( nzt, ngrdcol, gr,                 & ! intent(in)
                          p_in_Pa, 70000.0_core_rknd, thlm, & ! intent(in)
                          thlm700 )                           ! intent(out)

        call pvertinterp( nzt, ngrdcol, gr,                   & ! intent(in)
                          p_in_Pa, 100000.0_core_rknd, thlm,  & ! intent(in)
                          thlm1000 )                            ! intent(out)

        if ( stats_metadata%l_stats_samp ) then

          !$acc update host( edsclrm, thlm, rtm, thlm700, thlm1000 )

          ! thlm_ed and rtm_ed are budget terms intended to track the effect of explicit diffusion
          do k = 1, nzt
            do i = 1, ngrdcol         
              if ( thlm700(i) - thlm1000(i) < 20.0_core_rknd ) then
                thlm_ed(i,k) = ( edsclrm(i,k,edsclr_dim-1) - thlm(i,k) ) / dt
                rtm_ed(i,k)  = ( edsclrm(i,k,edsclr_dim)   - rtm(i,k) )  / dt
              else
                thlm_ed(i,k) = zero
                rtm_ed(i,k)  = zero
              end if
            end do
          end do

        end if

        !$acc parallel loop gang vector collapse(2) default(present)
        do k = 1, nzt
          do i = 1, ngrdcol         
            if ( thlm700(i) - thlm1000(i) < 20.0_core_rknd ) then
              thlm(i,k) = edsclrm(i,k,edsclr_dim-1)
              rtm(i,k)  = edsclrm(i,k,edsclr_dim)
            end if
          end do
        end do
        !$acc end parallel loop
        
      end if

      if ( stats_metadata%l_stats_samp ) then

        if ( edsclr_dim <= 1 .or. .not. clubb_config_flags%l_do_expldiff_rtm_thlm ) then
          ! thlm_ed and rtm_ed are budget terms intended to track the effect of explicit diffusion
          thlm_ed(:,:) = zero
          rtm_ed(:,:)  = zero
        end if

        do i = 1, ngrdcol
          call stat_update_var( stats_metadata%ithlm_ed, thlm_ed(i,:),  & ! intent(in)
                                stats_zt(i) )                             ! intent(inout)
          call stat_update_var( stats_metadata%irtm_ed, rtm_ed(i,:),    & ! intent(in)
                                stats_zt(i) )                             ! intent(inout)
        end do
      end if

      ! Eric Raut: this seems dangerous to call without any attached flag.
      ! Hence the preprocessor.
#ifdef CLUBB_CAM
      do edsclr=1,edsclr_dim
        ! upper_hf_level = nzt since we are filling the zt levels
        call fill_holes_vertical_api( nzt, ngrdcol, zero_threshold, 1, nzt, & ! In
                                      gr%dzt, rho_ds_zt, gr%grid_dir_indx,  & ! In
                                      clubb_config_flags%fill_holes_type,   & ! In
                                      edsclrm(:,:,edsclr) )                   ! InOut
      enddo
#endif

     endif ! advance_order_loop_iter

    enddo ! advance_order_loop_iter = 1, 4, 1

    !----------------------------------------------------------------
    ! Advance or otherwise calculate <thl'^3>, <rt'^3>, and
    ! <sclr'^3>.
    !----------------------------------------------------------------
    if ( l_advance_xp3 &
         .and. clubb_config_flags%iiPDF_type /= iiPDF_ADG1 ) then

      ! Advance <rt'^3>, <thl'^3>, and <sclr'^3> one model timestep using a
      ! simplified form of the <x'^3> predictive equation.  The simplified
      ! <x'^3> equation can either be advanced from its previous value or
      ! calculated using a steady-state approximation.
      call advance_xp3( nzm, nzt, ngrdcol, sclr_dim, sclr_tol, gr, dt, & ! Intent(in)
                        rtm, thlm, rtp2, thlp2, wprtp,                 & ! Intent(in)
                        wpthlp, wprtp2, wpthlp2, rho_ds_zm,            & ! Intent(in)
                        invrs_rho_ds_zt, invrs_tau_zt, tau_max_zt,     & ! Intent(in)
                        sclrm, sclrp2, wpsclrp, wpsclrp2,              & ! Intent(in)
                        clubb_config_flags%l_lmm_stepping,             & ! intent(in)
                        stats_metadata,                                & ! intent(in)
                        stats_zt,                                      & ! intent(inout)
                        rtp3, thlp3, sclrp3 )                            ! Intent(inout)

      ! Use a modified form of the Larson and Golaz (2005) ansatz for the
      ! ADG1 PDF to calculate <u'^3> and <v'^3> for another type of PDF.
      call Skx_func( nzt, ngrdcol, wp2_zt, wp3, &
                     w_tol, clubb_params, &
                     Skw_zt )

      upwp_zt(:,:) = zm2zt_api( nzm, nzt, ngrdcol, gr, upwp(:,:) )
      vpwp_zt(:,:) = zm2zt_api( nzm, nzt, ngrdcol, gr, vpwp(:,:) )
      ! Positive def. quantity
      up2_zt(:,:)  = zm2zt_api( nzm, nzt, ngrdcol, gr, up2(:,:), w_tol_sqd )
      ! Positive def. quantity
      vp2_zt(:,:)  = zm2zt_api( nzm, nzt, ngrdcol, gr, vp2(:,:), w_tol_sqd )

      thvm_zm(:,:)                  = zt2zm_api( nzm, nzt, ngrdcol, gr, thvm(:,:), zero_threshold )
      ddzm_thvm_zm(:,:)             = ddzm( nzm, nzt, ngrdcol, gr, thvm_zm(:,:) )
      brunt_vaisala_freq_sqd_zt(:,:) = max( ( grav / thvm(:,:) ) * ddzm_thvm_zm(:,:), zero )

      ! The xp3_coef_fnc is used in place of sigma_sqd_w_zt when the ADG1 PDF
      ! is not being used.  The xp3_coef_fnc provides some extra tunability to
      ! the simple xp3 equation.
      ! When xp3_coef_fnc goes to 0, the value of Skx goes to the smallest
      ! magnitude permitted by the function.  When xp3_coef_fnc goes to 1, the
      ! magnitude of Skx becomes huge.
      do k = 1, nzt
        do i = 1, ngrdcol
          xp3_coef_fnc(i,k) = clubb_params(i,ixp3_coef_base) &
                              + ( one - clubb_params(i,ixp3_coef_slope) ) &
                                * ( one - exp( brunt_vaisala_freq_sqd_zt(i,k) &
                                               / clubb_params(i,ixp3_coef_slope) ) )
        end do
      end do

      call xp3_LG_2005_ansatz( nzt, ngrdcol, Skw_zt, upwp_zt, wp2_zt, &
                               up2_zt, xp3_coef_fnc, &
                               clubb_params, w_tol, &
                               up3 )

      call xp3_LG_2005_ansatz( nzt, ngrdcol, Skw_zt, vpwp_zt, wp2_zt, &
                               vp2_zt, xp3_coef_fnc, &
                               clubb_params, w_tol, &
                               vp3 )

    else ! .not. l_advance_xp3 .or. clubb_config_flags%iiPDF_type = iiPDF_ADG1

      ! The ADG1 PDF must use this option.
      call Skx_func( nzt, ngrdcol, wp2_zt, wp3, &
                     w_tol, clubb_params, &
                     Skw_zt )

      wpthlp_zt(:,:) = zm2zt_api( nzm, nzt, ngrdcol, gr, wpthlp(:,:) )
      wprtp_zt(:,:)  = zm2zt_api( nzm, nzt, ngrdcol, gr, wprtp(:,:) )
      ! Positive def. quantity
      thlp2_zt(:,:)  = zm2zt_api( nzm, nzt, ngrdcol, gr, thlp2(:,:), thl_tol**2 )
      ! Positive def. quantity
      rtp2_zt(:,:)   = zm2zt_api( nzm, nzt, ngrdcol, gr, rtp2(:,:), rt_tol**2 )

      upwp_zt(:,:) = zm2zt_api( nzm, nzt, ngrdcol, gr, upwp(:,:) )
      vpwp_zt(:,:) = zm2zt_api( nzm, nzt, ngrdcol, gr, vpwp(:,:) )
      ! Positive def. quantity
      up2_zt(:,:)  = zm2zt_api( nzm, nzt, ngrdcol, gr, up2(:,:), w_tol_sqd )
      ! Positive def. quantity
      vp2_zt(:,:)  = zm2zt_api( nzm, nzt, ngrdcol, gr, vp2(:,:), w_tol_sqd )

      if ( clubb_config_flags%iiPDF_type == iiPDF_ADG1 ) then

        ! Use the Larson and Golaz (2005) ansatz for the ADG1 PDF to
        ! calculate <rt'^3>, <thl'^3>, <u'^3>, <v'^3>, and <sclr'^3>.
        sigma_sqd_w_zt(:,:) = zm2zt_api( nzm, nzt, ngrdcol, gr, sigma_sqd_w(:,:), zero_threshold )

        call xp3_LG_2005_ansatz( nzt, ngrdcol, Skw_zt, wpthlp_zt, wp2_zt, &
                                 thlp2_zt, sigma_sqd_w_zt, &
                                 clubb_params, thl_tol, &
                                 thlp3 )

        call xp3_LG_2005_ansatz( nzt, ngrdcol, Skw_zt, wprtp_zt, wp2_zt, &
                                 rtp2_zt, sigma_sqd_w_zt, &
                                 clubb_params, rt_tol, &
                                 rtp3 )

        call xp3_LG_2005_ansatz( nzt, ngrdcol, Skw_zt, upwp_zt, wp2_zt, &
                                 up2_zt, sigma_sqd_w_zt, &
                                 clubb_params, w_tol, &
                                 up3 )

        call xp3_LG_2005_ansatz( nzt, ngrdcol, Skw_zt, vpwp_zt, wp2_zt, &
                                 vp2_zt, sigma_sqd_w_zt, &
                                 clubb_params, w_tol, &
                                 vp3 )

        do sclr = 1, sclr_dim, 1
          
          wpsclrp_zt(:,:) = zm2zt_api( nzm, nzt, ngrdcol, gr, wpsclrp(:,:,sclr), &
                                       sclr_tol(sclr)**2 )
          sclrp2_zt(:,:)  = zm2zt_api( nzm, nzt, ngrdcol, gr, sclrp2(:,:,sclr), &
                                       sclr_tol(sclr)**2 )

          call xp3_LG_2005_ansatz( nzt, ngrdcol, Skw_zt, wpsclrp_zt, wp2_zt, &
                                   sclrp2_zt, sigma_sqd_w_zt, &
                                   clubb_params, sclr_tol(sclr), &
                                   sclrp3 )

        enddo ! sclr = 1, sclr_dim

      else ! clubb_config_flags%iiPDF_type /= iiPDF_ADG1

        ! Use a modified form of the Larson and Golaz (2005) ansatz for the
        ! ADG1 PDF to calculate <u'^3> and <v'^3> for another type of PDF.
        thvm_zm(:,:)                   = zt2zm_api( nzm, nzt, ngrdcol, gr, thvm(:,:), &
                                                    zero_threshold )
        ddzm_thvm_zm(:,:)              = ddzm( nzm, nzt, ngrdcol, gr, thvm_zm(:,:) )
        brunt_vaisala_freq_sqd_zt(:,:) = max( ( grav / thvm(:,:) ) * ddzm_thvm_zm(:,:), zero )
        
        
        ! Initialize sigma_sqd_w_zt to zero so we don't break output
        do k = 1, nzt
          do i = 1, ngrdcol
            sigma_sqd_w_zt(i,k) = zero
          end do
        end do

        ! The xp3_coef_fnc is used in place of sigma_sqd_w_zt when the
        ! ADG1 PDF is not being used.  The xp3_coef_fnc provides some extra
        ! tunability to the simple xp3 equation.
        ! When xp3_coef_fnc goes to 0, the value of Skx goes to the smallest
        ! magnitude permitted by the function.  When xp3_coef_fnc goes to 1,
        ! the magnitude of Skx becomes huge.
        ! The value of Skx becomes large near cloud top, where there is a
        ! higher degree of static stability.  The exp{ } portion of the
        ! xp3_coef_fnc allows the xp3_coef_fnc to become larger in regions
        ! of high static stability, producing larger magnitude values of Skx.
        do k = 1, nzt
          do i = 1, ngrdcol
            xp3_coef_fnc(i,k) = clubb_params(i,ixp3_coef_base) &
              + ( one - clubb_params(i,ixp3_coef_slope) ) &
                * ( one - exp( brunt_vaisala_freq_sqd_zt(i,k) / clubb_params(i,ixp3_coef_slope) ) )
          end do
        end do
        
        call xp3_LG_2005_ansatz( nzt, ngrdcol, Skw_zt, wpthlp_zt, wp2_zt, &
                                 thlp2_zt, xp3_coef_fnc, &
                                 clubb_params, thl_tol, &
                                 thlp3 )

        call xp3_LG_2005_ansatz( nzt, ngrdcol, Skw_zt, wprtp_zt, wp2_zt, &
                                 rtp2_zt, xp3_coef_fnc, &
                                 clubb_params, rt_tol, &
                                 rtp3 )

        call xp3_LG_2005_ansatz( nzt, ngrdcol, Skw_zt, upwp_zt, wp2_zt, &
                                 up2_zt, xp3_coef_fnc, &
                                 clubb_params, w_tol, &
                                 up3 )

        call xp3_LG_2005_ansatz( nzt, ngrdcol, Skw_zt, vpwp_zt, wp2_zt, &
                                 vp2_zt, xp3_coef_fnc, &
                                 clubb_params, w_tol, &
                                 vp3 )

        do sclr = 1, sclr_dim, 1
          
          wpsclrp_zt(:,:) = zm2zt_api( nzm, nzt, ngrdcol, gr, wpsclrp(:,:,sclr) )
          sclrp2_zt(:,:)  = zm2zt_api( nzm, nzt, ngrdcol, gr, sclrp2(:,:,sclr), sclr_tol(sclr)**2 )

          call xp3_LG_2005_ansatz( nzt, ngrdcol, Skw_zt(:,:), wpsclrp_zt(:,:), wp2_zt(:,:), &
                                   sclrp2_zt(:,:), xp3_coef_fnc(:,:), &
                                   clubb_params, sclr_tol(sclr), &
                                   sclrp3(:,:,sclr) )
        end do ! sclr = 1, sclr_dim

      end if ! clubb_config_flags%iiPDF_type == iiPDF_ADG1

    end if ! l_advance_xp3 .and. clubb_config_flags%iiPDF_type /= iiPDF_ADG1

    if ( clubb_config_flags%ipdf_call_placement == ipdf_post_advance_fields &
         .or. clubb_config_flags%ipdf_call_placement &
              == ipdf_pre_post_advance_fields ) then

      ! Sample stats in this call to subroutine pdf_closure_driver for
      ! ipdf_post_advance_fields, but not for ipdf_pre_post_advance_fields
      ! because stats were sampled during the first call to subroutine
      ! pdf_closure_driver.
      if ( clubb_config_flags%ipdf_call_placement &
          == ipdf_post_advance_fields ) then
        l_samp_stats_in_pdf_call = .true.
      elseif ( clubb_config_flags%ipdf_call_placement &
              == ipdf_pre_post_advance_fields ) then
        l_samp_stats_in_pdf_call = .false.
      endif

      !########################################################################
      !#######                     CALL CLUBB's PDF                     #######
      !#######   AND OUTPUT PDF PARAMETERS AND INTEGRATED QUANTITITES   #######
      !########################################################################
      ! Given CLUBB's prognosed moments, diagnose CLUBB's PDF parameters
      !   and quantities integrated over that PDF, including
      !   quantities related to clouds, buoyancy, and turbulent advection.
      call pdf_closure_driver( gr, nzm, nzt, ngrdcol,                       & ! Intent(in)
                               dt, hydromet_dim, sclr_dim, sclr_tol,        & ! Intent(in)
                               wprtp, thlm, wpthlp, rtp2, rtp3,             & ! Intent(in)
                               thlp2, thlp3, rtpthlp, wp2,                  & ! Intent(in)
                               wp3, wm_zm, wm_zt,                           & ! Intent(in)
                               um, up2, upwp, up3,                          & ! Intent(in)
                               vm, vp2, vpwp, vp3,                          & ! Intent(in)
                               p_in_Pa, exner,                              & ! Intent(in)
                               thv_ds_zm, thv_ds_zt, rtm_ref,               & ! Intent(in)
                               wphydrometp,                                 & ! Intent(in)
                               wp2hmp, rtphmp_zt, thlphmp_zt,               & ! Intent(in)
                               sclrm, wpsclrp, sclrp2,                      & ! Intent(in)
                               sclrprtp, sclrpthlp, sclrp3,                 & ! Intent(in)
                               p_sfc, l_samp_stats_in_pdf_call,             & ! Intent(in)
                               mixt_frac_max_mag, ts_nudge,                 & ! Intent(in)
                               rtm_min, rtm_nudge_max_altitude,             & ! Intent(in)
                               clubb_params,                                & ! Intent(in)
                               clubb_config_flags%iiPDF_type,               & ! Intent(in)
                               clubb_config_flags%saturation_formula,       & ! Intent(in)
                               clubb_config_flags%l_predict_upwp_vpwp,      & ! Intent(in)
                               clubb_config_flags%l_rtm_nudge,              & ! Intent(in)
                               clubb_config_flags%l_trapezoidal_rule_zt,    & ! Intent(in)
                               clubb_config_flags%l_trapezoidal_rule_zm,    & ! Intent(in)
                               clubb_config_flags%l_call_pdf_closure_twice, & ! Intent(in)
                               clubb_config_flags%l_use_cloud_cover,        & ! Intent(in)
                               clubb_config_flags%l_rcm_supersat_adj,       & ! Intent(in)
                               l_mix_rat_hm,                                & ! Intent(in)
                               stats_metadata,                              & ! Intent(in)
                               stats_zt, stats_zm,                          & ! Intent(inout)
                               rtm,                                         & ! Intent(inout)
                               pdf_implicit_coefs_terms,                    & ! Intent(inout)
                               pdf_params, pdf_params_zm, err_info,         & ! Intent(inout)
#ifdef GFDL
                               RH_crit(k, : , :),                           & ! Intent(inout)
                               do_liquid_only_in_clubb,                     & ! Intent(in)
#endif
                               rcm, cloud_frac,                             & ! Intent(out)
                               ice_supersat_frac, wprcp,                    & ! Intent(out)
                               sigma_sqd_w, wpthvp, wp2thvp, wp2up,         & ! Intent(out)
                               rtpthvp, thlpthvp, rc_coef,                  & ! Intent(out)
                               rcm_in_layer, cloud_cover,                   & ! Intent(out)
                               rcp2_zt, thlprcp,                            & ! Intent(out)
                               rc_coef_zm, sclrpthvp,                       & ! Intent(out)
                               wpup2, wpvp2,                                & ! Intent(out)
                               wp2up2, wp2vp2, wp4,                         & ! Intent(out)
                               wp2rtp, wprtp2, wp2thlp,                     & ! Intent(out)
                               wpthlp2, wprtpthlp, wp2rcp,                  & ! Intent(out)
                               rtprcp, rcp2,                                & ! Intent(out)
                               uprcp, vprcp,                                & ! Intent(out)
                               w_up_in_cloud, w_down_in_cloud,              & ! Intent(out)
                               cloudy_updraft_frac,                         & ! Intent(out)
                               cloudy_downdraft_frac,                       & ! intent(out)
                               Skw_velocity,                                & ! Intent(out)
                               cloud_frac_zm,                               & ! Intent(out)
                               ice_supersat_frac_zm,                        & ! Intent(out)
                               rtm_zm, thlm_zm, rcm_zm,                     & ! Intent(out)
                               rcm_supersat_adj,                            & ! Intent(out)
                               wp2sclrp, wpsclrp2, sclrprcp,                & ! Intent(out)
                               wpsclrprtp, wpsclrpthlp )                      ! Intent(out)

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

    if ( stats_metadata%l_stats_samp ) then

      !$acc update host( wp2, vp2, up2, wprtp, wpthlp, upwp, vpwp, rtp2, thlp2, &
      !$acc              rtpthlp, rtm, thlm, um, vm, wp3, &
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
      !$acc              um, vm, upwp, vpwp, up2, vp2, &
      !$acc              thlm, rtm, wprtp, wpthlp, &
      !$acc              wp2, wp3, rtp2, rtp3, thlp2, thlp3, rtpthlp, &
      !$acc              wpthvp, wp2thvp, rtpthvp, thlpthvp, &
      !$acc              p_in_Pa, exner, rho, rho_zm, &
      !$acc              rho_ds_zm, rho_ds_zt, thv_ds_zm, thv_ds_zt, &
      !$acc              wm_zt, wm_zm, rcm, wprcp, rc_coef_zm, rc_coef, &
      !$acc              rcm_zm, rtm_zm, thlm_zm, cloud_frac, ice_supersat_frac, &
      !$acc              cloud_frac_zm, ice_supersat_frac_zm, rcm_in_layer, &
      !$acc              cloud_cover, rcm_supersat_adj, sigma_sqd_w, &
      !$acc              thvm, ug, vg, Lscale, wpthlp2, wp2thlp, wprtp2, wp2rtp, &
      !$acc              Lscale_up, Lscale_down, tau_zt, Kh_zt, wp2rcp, &
      !$acc              wprtpthlp, sigma_sqd_w_zt, wp2_zt, thlp2_zt, &
      !$acc              wpthlp_zt, wprtp_zt, rtp2_zt, rtpthlp_zt, up2_zt, &
      !$acc              vp2_zt, upwp_zt, vpwp_zt, wpup2, wpvp2, &
      !$acc              wp2up2, wp2vp2, wp4, &
      !$acc              tau_zm, Kh_zm, thlprcp, &
      !$acc              rtprcp, rcp2, em, a3_coef, a3_coef_zt, &
      !$acc              wp3_zm, wp3_on_wp2, wp3_on_wp2_zt, Skw_velocity, &
      !$acc              w_up_in_cloud, w_down_in_cloud, &
      !$acc              cloudy_updraft_frac, cloudy_downdraft_frac )

      !$acc if( clubb_config_flags%l_call_pdf_closure_twice ) &
      !$acc update host( pdf_params_zm%w_1, pdf_params_zm%w_2, &
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
      !$acc              pdf_params_zm%ice_supersat_frac_2 )

      !$acc update host( sclrm, sclrp2, sclrprtp, sclrpthlp, sclrm_forcing, &
      !$acc              sclrpthvp, wpsclrp, sclrprcp, wp2sclrp, wpsclrp2, &
      !$acc              wpsclrprtp, wpsclrpthlp, wpedsclrp ) &
      !$acc if ( sclr_dim > 0 )

      !$acc update host( edsclrm, edsclrm_forcing  ) if ( edsclr_dim > 0 )

      do i = 1, ngrdcol

        call stat_end_update( nzm, stats_metadata%iwp2_bt, wp2(i,:) / dt, & ! intent(in)
                              stats_zm(i) )           ! intent(inout)
        call stat_end_update( nzm, stats_metadata%ivp2_bt, vp2(i,:) / dt, & ! intent(in)
                              stats_zm(i) )           ! intent(inout)
        call stat_end_update( nzm, stats_metadata%iup2_bt, up2(i,:) / dt, & ! intent(in)
                              stats_zm(i) )           ! intent(inout)
        call stat_end_update( nzm, stats_metadata%iwprtp_bt, wprtp(i,:) / dt, & ! intent(in)
                              stats_zm(i) )               ! intent(inout)
        call stat_end_update( nzm, stats_metadata%iwpthlp_bt, wpthlp(i,:) / dt, & ! intent(in)
                              stats_zm(i) )                 ! intent(inout)
        if ( clubb_config_flags%l_predict_upwp_vpwp ) then
           call stat_end_update( nzm, stats_metadata%iupwp_bt, upwp(i,:) / dt, & ! intent(in)
                                 stats_zm(i) )             ! intent(inout)
           call stat_end_update( nzm, stats_metadata%ivpwp_bt, vpwp(i,:) / dt, & ! intent(in)
                                 stats_zm(i) )             ! intent(inout)
        endif ! l_predict_upwp_vpwp
        call stat_end_update( nzm, stats_metadata%irtp2_bt, rtp2(i,:) / dt, & ! intent(in)
                              stats_zm(i) )             ! intent(inout)
        call stat_end_update( nzm, stats_metadata%ithlp2_bt, thlp2(i,:) / dt, & ! intent(in)
                              stats_zm(i) )               ! intent(inout)
        call stat_end_update( nzm, stats_metadata%irtpthlp_bt, rtpthlp(i,:) / dt, & ! intent(in)
                              stats_zm(i) )                   ! intent(inout)
 
        call stat_end_update( nzt, stats_metadata%irtm_bt, rtm(i,:) / dt, & ! intent(in)
                              stats_zt(i) )           ! intent(inout)
        call stat_end_update( nzt, stats_metadata%ithlm_bt, thlm(i,:) / dt, & ! intent(in)
                              stats_zt(i) )             ! intent(inout)
        call stat_end_update( nzt, stats_metadata%ium_bt, um(i,:) / dt, & ! intent(in)
                              stats_zt(i) )         ! intent(inout)
        call stat_end_update( nzt, stats_metadata%ivm_bt, vm(i,:) / dt, & ! intent(in)
                              stats_zt(i) )         ! intent(inout)
        call stat_end_update( nzt, stats_metadata%iwp3_bt, wp3(i,:) / dt, & ! intent(in)
                              stats_zt(i) )           ! intent(inout)
      end do

      if ( stats_metadata%iwpthlp_zt > 0 ) then
        wpthlp_zt(:,:)  = zm2zt_api( nzm, nzt, ngrdcol, gr, wpthlp(:,:) )
      end if

      if ( stats_metadata%iwprtp_zt > 0 ) then
        wprtp_zt(:,:)   = zm2zt_api( nzm, nzt, ngrdcol, gr, wprtp(:,:) )
      end if

      if ( stats_metadata%iup2_zt > 0 ) then
        up2_zt(:,:) = zm2zt_api( nzm, nzt, ngrdcol, gr, up2(:,:), w_tol_sqd )
      end if

      if (stats_metadata%ivp2_zt > 0 ) then
        vp2_zt(:,:) = zm2zt_api( nzm, nzt, ngrdcol, gr, vp2(:,:), w_tol_sqd )
      end if

      if ( stats_metadata%iupwp_zt > 0 ) then
        upwp_zt(:,:) = zm2zt_api( nzm, nzt, ngrdcol, gr, upwp(:,:) )
      end if

      if ( stats_metadata%ivpwp_zt > 0 ) then
        vpwp_zt(:,:) = zm2zt_api( nzm, nzt, ngrdcol, gr, vpwp(:,:) )
      end if

      do i = 1, ngrdcol

        call stats_accumulate( &
               nzm, nzt, i, sclr_dim, edsclr_dim, gr%invrs_dzm(i,:), gr%zt(i,:),          & ! In
               gr%grid_dir * gr%dzm(i,:), gr%grid_dir * gr%dzt(i,:), dt,                  & ! In
               um(i,:), vm(i,:), upwp(i,:), vpwp(i,:), up2(i,:), vp2(i,:),                & ! In
               thlm(i,:), rtm(i,:), wprtp(i,:), wpthlp(i,:),                              & ! In
               wp2(i,:), wp3(i,:), rtp2(i,:), rtp3(i,:), thlp2(i,:), thlp3(i,:),          & ! In
               rtpthlp(i,:),                                                              & ! In
               wpthvp(i,:), wp2thvp(i,:), wp2up(i,:), rtpthvp(i,:), thlpthvp(i,:),        & ! In
               p_in_Pa(i,:), exner(i,:), rho(i,:), rho_zm(i,:),                           & ! In
               rho_ds_zm(i,:), rho_ds_zt(i,:), thv_ds_zm(i,:), thv_ds_zt(i,:),            & ! In
               wm_zt(i,:), wm_zm(i,:), rcm(i,:), wprcp(i,:), rc_coef(i,:),                & ! In
               rc_coef_zm(i,:),                                                           & ! In
               rcm_zm(i,:), rtm_zm(i,:), thlm_zm(i,:), cloud_frac(i,:),                   & ! In
               ice_supersat_frac(i,:),                                                    & ! In
               cloud_frac_zm(i,:), ice_supersat_frac_zm(i,:), rcm_in_layer(i,:),          & ! In
               cloud_cover(i,:), rcm_supersat_adj(i,:), sigma_sqd_w(i,:),                 & ! In
               thvm(i,:), ug(i,:), vg(i,:), Lscale(i,:), wpthlp2(i,:), wp2thlp(i,:),      & ! In
               wprtp2(i,:), wp2rtp(i,:),                                                  & ! In
               Lscale_up(i,:), Lscale_down(i,:), tau_zt(i,:), Kh_zt(i,:), wp2rcp(i,:),    & ! In
               wprtpthlp(i,:), sigma_sqd_w_zt(i,:), rsat(i,:), wp2_zt(i,:),               & ! In
               thlp2_zt(i,:),                                                             & ! In
               wpthlp_zt(i,:), wprtp_zt(i,:), rtp2_zt(i,:), rtpthlp_zt(i,:),              & ! In
               up2_zt(i,:),                                                               & ! In
               vp2_zt(i,:), upwp_zt(i,:), vpwp_zt(i,:), wpup2(i,:), wpvp2(i,:),           & ! In
               wp2up2(i,:), wp2vp2(i,:), wp4(i,:),                                        & ! In
               tau_zm(i,:), Kh_zm(i,:), thlprcp(i,:),                                     & ! In
               rtprcp(i,:), rcp2(i,:), em(i,:), a3_coef(i,:), a3_coef_zt(i,:),            & ! In
               wp3_zm(i,:), wp3_on_wp2(i,:), wp3_on_wp2_zt(i,:), Skw_velocity(i,:),       & ! In
               w_up_in_cloud(i,:), w_down_in_cloud(i,:),                                  & ! In
               cloudy_updraft_frac(i,:), cloudy_downdraft_frac(i,:),                      & ! In
               pdf_params, pdf_params_zm,                                                 & ! In
               sclrm(i,:,:), sclrp2(i,:,:),                                               & ! In
               sclrprtp(i,:,:), sclrpthlp(i,:,:), sclrm_forcing(i,:,:), sclrpthvp(i,:,:), & ! In
               wpsclrp(i,:,:), sclrprcp(i,:,:), wp2sclrp(i,:,:), wpsclrp2(i,:,:),         & ! In
               wpsclrprtp(i,:,:),                                                         & ! In
               wpsclrpthlp(i,:,:), wpedsclrp(i,:,:), edsclrm(i,:,:),                      & ! In
               edsclrm_forcing(i,:,:),                                                    & ! In
               clubb_config_flags%saturation_formula,                                     & ! In
               clubb_config_flags%l_call_pdf_closure_twice,                               & ! In
               stats_metadata,                                                            & ! In
               stats_zt(i), stats_zm(i), stats_sfc(i) )                                     ! In/Out
      end do
    endif ! stats_metadata%l_stats_samp

    if ( clubb_at_least_debug_level_api( 2 ) ) then

      !$acc update host( thlm_forcing, rtm_forcing, um_forcing, vm_forcing, &
      !$acc              wm_zm, wm_zt, p_in_Pa, rho_zm, rho, exner, rho_ds_zm, &
      !$acc              rho_ds_zt, invrs_rho_ds_zm, invrs_rho_ds_zt, thv_ds_zm, &
      !$acc              thv_ds_zt, wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, p_sfc, &
      !$acc              um, upwp, vm, vpwp, up2, vp2, rtm, wprtp, thlm, wpthlp, &
      !$acc              wp2, wp3, rtp2, thlp2, rtpthlp )

      !$acc update host( wpsclrp_sfc, wpedsclrp_sfc, sclrm, wpsclrp, &
      !$acc              sclrp2, sclrprtp, sclrpthlp, sclrm_forcing ) &
      !$acc if ( sclr_dim > 0 )

      !$acc update host( edsclrm, edsclrm_forcing  ) if ( edsclr_dim > 0 )

      do i = 1, ngrdcol
        call parameterization_check( &
             nzm, nzt, sclr_dim, edsclr_dim, &
             thlm_forcing(i,:), rtm_forcing(i,:), um_forcing(i,:),                     & ! In
             vm_forcing(i,:), wm_zm(i,:), wm_zt(i,:), p_in_Pa(i,:),                    & ! In
             rho_zm(i,:), rho(i,:), exner(i,:), rho_ds_zm(i,:),                        & ! In
             rho_ds_zt(i,:), invrs_rho_ds_zm(i,:), invrs_rho_ds_zt(i,:),               & ! In
             thv_ds_zm(i,:), thv_ds_zt(i,:), wpthlp_sfc(i), wprtp_sfc(i), upwp_sfc(i), & ! In
             vpwp_sfc(i), p_sfc(i), um(i,:), upwp(i,:), vm(i,:), vpwp(i,:), up2(i,:),  & ! In
             vp2(i,:), rtm(i,:), wprtp(i,:), thlm(i,:), wpthlp(i,:), wp2(i,:),         & ! In
             wp3(i,:), rtp2(i,:), thlp2(i,:), rtpthlp(i,:),                            & ! In
             !rcm,                                                                     &
             "end of ",                                                                & ! In
             wpsclrp_sfc(i,:), wpedsclrp_sfc(i,:), sclrm(i,:,:), wpsclrp(i,:,:),       & ! In
             sclrp2(i,:,:),                                                            & ! In
             sclrprtp(i,:,:), sclrpthlp(i,:,:), sclrm_forcing(i,:,:), edsclrm(i,:,:),  & ! In
             edsclrm_forcing(i,:,:), &                                                   ! In
             err_info )                                                                  ! Inout
        if ( err_info%err_code(i) == clubb_fatal_error ) then
          write(fstderr, *) err_info%err_header(i)
        endif

      end do

      if ( any(err_info%err_code == clubb_fatal_error) ) then
        write(fstderr,*) "Error occurred during parameterization_check at"// &
                         " end of advance_clubb_core"
        return
      end if

    end if

    if ( stats_metadata%l_stats .and. stats_metadata%l_stats_samp ) then

      !$acc update host( wm_zt, wm_zm, rho_ds_zm, wprtp, wprtp_sfc, rho_ds_zt, &
      !$acc              rtm, rtm_forcing, thlm, thlm_forcing )

      ! Spurious source will only be calculated if rtm_ma and thlm_ma are zero.
      ! Therefore, wm must be zero or l_implemented must be true.
      do i = 1, ngrdcol
        if ( l_implemented .or. &
             (all( abs(wm_zt(i,:)) < eps ) .and. all( abs(wm_zm(i,:)) < eps ))) then
          ! Calculate the spurious source for rtm
          rtm_flux_top(i) = rho_ds_zm(i,gr%k_ub_zm) * wprtp(i,gr%k_ub_zm)

          if ( .not. clubb_config_flags%l_host_applies_sfc_fluxes ) then
            rtm_flux_sfc(i) = rho_ds_zm(i,gr%k_lb_zm) * wprtp_sfc(i)
          else
            rtm_flux_sfc(i) = 0.0_core_rknd
          end if

          rtm_integral_after(i)  &
          = vertical_integral( nzt, rho_ds_zt(i,:), &
                               rtm(i,:), gr%dzt(i,:) )

          rtm_integral_forcing(i)  &
          = vertical_integral( nzt, rho_ds_zt(i,:), &
                               rtm_forcing(i,:), gr%dzt(i,:) )

          rtm_spur_src(i)  &
          = calculate_spurious_source( rtm_integral_after(i), &
                                       rtm_integral_before(i), &
                                       rtm_flux_top(i), rtm_flux_sfc(i), &
                                       rtm_integral_forcing(i), &
                                       dt )

          ! Calculate the spurious source for thlm
          thlm_flux_top(i) = rho_ds_zm(i,gr%k_ub_zm) * wpthlp(i,gr%k_ub_zm)

          if ( .not. clubb_config_flags%l_host_applies_sfc_fluxes ) then
            thlm_flux_sfc(i) = rho_ds_zm(i,gr%k_lb_zm) * wpthlp_sfc(i)
          else
            thlm_flux_sfc(i) = 0.0_core_rknd
          end if

          thlm_integral_after(i)  &
          = vertical_integral( nzt, rho_ds_zt(i,:), &
                               thlm(i,:), gr%dzt(i,:) )

          thlm_integral_forcing(i)  &
          = vertical_integral( nzt, rho_ds_zt(i,:), &
                               thlm_forcing(i,:), gr%dzt(i,:) )

          thlm_spur_src(i)  &
          = calculate_spurious_source( thlm_integral_after(i), &
                                       thlm_integral_before(i), &
                                       thlm_flux_top(i), thlm_flux_sfc(i), &
                                       thlm_integral_forcing(i), &
                                       dt )
        else ! If l_implemented is false, we don't want spurious source output
          rtm_spur_src(i) = -9999.0_core_rknd
          thlm_spur_src(i) = -9999.0_core_rknd
        end if
      end do

      ! Write the var to stats
      do i = 1, ngrdcol
        call stat_update_var_pt( stats_metadata%irtm_spur_src, 1, & ! intent(in)
                                 rtm_spur_src(i),                 & ! intent(in)
                                 stats_sfc(i) )                     ! intent(inout)
        call stat_update_var_pt( stats_metadata%ithlm_spur_src, 1, & ! intent(in)
                                 thlm_spur_src(i),                 & ! intent(in)
                                 stats_sfc(i) )                      ! intent(inout)
      end do
    end if

    !$acc exit data if( sclr_dim > 0 ) &
    !$acc           delete( sclrprcp, wp2sclrp, &
    !$acc                   wpsclrp2, wpsclrprtp, wpsclrpthlp, wpsclrp_zt, sclrp2_zt )

    !$acc exit data if( edsclr_dim > 0 ) &
    !$acc           delete( wpedsclrp )

    !$acc exit data delete( Skw_zm, Skw_zt, thvm, thvm_zm, ddzm_thvm_zm, rtprcp, rcp2, &
    !$acc                   wpthlp2, wprtp2, wprtpthlp, wp2rcp, wp3_zm, Lscale_up, &
    !$acc                   Lscale_zm, Lscale_down, em, tau_zm, tau_zt, sigma_sqd_w_zt, &
    !$acc                   ddzt_um, ddzt_vm, ddzt_umvm_sqd, ddzt_umvm_sqd_clipped, &
    !$acc                   wp2_zt, thlp2_zt, wpthlp_zt, wprtp_zt, &
    !$acc                   rtp2_zt, rtpthlp_zt, up2_zt, vp2_zt, upwp_zt, vpwp_zt, &
    !$acc                   Skw_velocity, a3_coef, a3_coef_zt, wp3_on_wp2, wp3_on_wp2_zt, &
    !$acc                   rc_coef, Km_zm, Kmh_zm, gamma_Skw_fnc, sigma_sqd_w, sigma_sqd_w_tmp, &
    !$acc                   sqrt_em_zt, xp3_coef_fnc, w_1_zm, w_2_zm, varnce_w_1_zm, &
    !$acc                   varnce_w_2_zm, &
    !$acc                   mixt_frac_zm, rcp2_zt, cloud_frac_zm, ice_supersat_frac_zm, rtm_zm, &
    !$acc                   thlm_zm, rcm_zm, thlm1000, thlm700, &
    !$acc                   rcm_supersat_adj, stability_correction, invrs_tau_N2_zm, &
    !$acc                   invrs_tau_C6_zm, invrs_tau_C1_zm, invrs_tau_xp2_zm, invrs_tau_N2_iso, &
    !$acc                   invrs_tau_C4_zm, invrs_tau_C14_zm, invrs_tau_wp2_zm, &
    !$acc                   invrs_tau_wpxp_zm, &
    !$acc                   invrs_tau_wp3_zm, invrs_tau_no_N2_zm, invrs_tau_bkgnd, &
    !$acc                   invrs_tau_shear, &
    !$acc                   invrs_tau_sfc, invrs_tau_zt, invrs_tau_wp3_zt, Cx_fnc_Richardson, &
    !$acc                   brunt_vaisala_freq_sqd, brunt_vaisala_freq_sqd_mixed, &
    !$acc                   brunt_vaisala_freq_sqd_dry, brunt_vaisala_freq_sqd_moist, &
    !$acc                   brunt_vaisala_freq_sqd_splat, brunt_vaisala_freq_sqd_smth, &
    !$acc                   brunt_vaisala_freq_sqd_zt, brunt_vaisala_freq_clipped, Ri_zm, &
    !$acc                   Lscale_max, tau_max_zm, tau_max_zt, mu, lhs_splat_wp2, lhs_splat_wp3 )
    
    return

  end subroutine advance_clubb_core

end module advance_clubb_core_module
