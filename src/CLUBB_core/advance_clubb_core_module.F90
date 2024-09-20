
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
      ipdf_pre_post_advance_fields, & ! Call both before and after advancing
                                      ! predictive fields
      saturation_bolton, &
      saturation_gfdl, &
      saturation_flatau, &
      saturation_lookup
                             

  implicit none

  public ::  &
    check_clubb_settings, &
    advance_clubb_core, &
    cleanup_clubb_core, &
    set_Lscale_max, &
    calculate_thlp2_rad

  private ! Default Scope

  ! Advance subroutine ordering variables
  integer, parameter, private :: &
    order_xm_wpxp = 1, &
    order_xp2_xpyp = 2, &
    order_wp2_wp3 = 3, &
    order_windm = 4

  contains

  !-----------------------------------------------------------------------

  !#######################################################################
  !#######################################################################
  ! If you change the argument list of advance_clubb_core you also have to
  ! change the calls to this function in the host models CAM, WRF, SAM
  ! and GFDL.
  !#######################################################################
  !#######################################################################
  subroutine advance_clubb_core( gr, nzm, nzt, ngrdcol, &           ! intent(in)
               l_implemented, dt, fcor, sfc_elevation,            & ! intent(in)
               hydromet_dim,                                      & ! intent(in)
               sclr_dim, sclr_tol, edsclr_dim, sclr_idx,      & ! intent(in)
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
               hydromet, l_mix_rat_hm, &                            ! intent(in)
               rfrzm, radf, &                                       ! intent(in)
#ifdef CLUBBND_CAM
               varmu, &                                             ! intent(in)
#endif
               wphydrometp, wp2hmp, rtphmp_zt, thlphmp_zt, &        ! intent(in)
               host_dx, host_dy, &                                  ! intent(in)
               clubb_params, nu_vert_res_dep, lmin, &               ! intent(in)
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
               wpthvp, wp2thvp, rtpthvp, thlpthvp, &                ! intent(inout)
               sclrpthvp, &                                         ! intent(inout)
               wp2rtp, wp2thlp, uprcp, vprcp, rc_coef_zm, wp4, &    ! intent(inout)
               wpup2, wpvp2, wp2up2, wp2vp2, ice_supersat_frac, &   ! intent(inout)
               um_pert, vm_pert, upwp_pert, vpwp_pert, &            ! intent(inout)
               pdf_params, pdf_params_zm, &                         ! intent(inout)
               pdf_implicit_coefs_terms, &                          ! intent(inout)
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
               err_code_out )                                       ! intent(out)

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
        eps

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
        iiPDF_ADG1

    use grid_class, only: &
        grid, & ! Type
        zm2zt,  & ! Procedure(s)
        zt2zm, &
        ddzm, &
        ddzt, &
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
        calc_Lscale_directly,  &  ! for Lscale
        diagnose_Lscale_from_tau  ! for Lscale from tau

    use advance_windm_edsclrm_module, only:  &
        advance_windm_edsclrm  ! Procedure(s)

    use saturation, only:  &
        ! Procedure
        sat_mixrat_liq ! Saturation mixing ratio

    use advance_wp2_wp3_module, only:  &
        advance_wp2_wp3 ! Procedure

    use advance_xp3_module, only: &
        advance_xp3    ! Procedure(s)

    use calc_pressure, only: &
        calculate_thvm

    use clubb_precision, only:  &
        core_rknd ! Variable(s)

    use error_code, only: &
        clubb_at_least_debug_level,  & ! Procedure
        err_code,                    & ! Error Indicator
        clubb_no_error, &              ! Constant
        clubb_fatal_error              ! Constant

    use Skx_module, only: &
        Skx_func,           & ! Procedure(s)
        xp3_LG_2005_ansatz

    use clip_explicit, only: &
        clip_covars_denom ! Procedure(s)

    use T_in_K_module, only: &
        ! Read values from namelist
        thlm2T_in_K ! Procedure

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
        fill_holes_vertical

    use advance_helper_module, only: &
        calc_stability_correction, & ! Procedure(s)
        compute_Cx_fnc_Richardson, &
        calc_brunt_vaisala_freq_sqd, &
        wp2_term_splat_lhs, &
        wp3_term_splat_lhs, &
        vertical_integral, &
        Lscale_width_vert_avg, &
        smooth_max, &
        calc_Ri_zm

    use interpolation, only: &
        pvertinterp

    use stats_type, only: stats ! Type
    
    use pdf_parameter_module, only: &
      copy_single_pdf_params_to_multi, &
      copy_multi_pdf_params_to_single, &
      init_pdf_params

    use stats_variables, only: &
      stats_metadata_type

    use array_index, only: &
      sclr_idx_type

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

    type (grid), target, intent(in) :: gr

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

    real( kind = core_rknd ), dimension(ngrdcol,nzt,hydromet_dim), intent(in) :: &
      hydromet           ! Array of hydrometeors                [units vary]

    logical, dimension(hydromet_dim), intent(in) :: &
      l_mix_rat_hm   ! if true, then the quantity is a hydrometeor mixing ratio

    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: &
      radf          ! Buoyancy production at cloud top due to longwave radiative cooling [m^2/s^3]

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
      lmin    ! Min. value for the length scale    [m]

    type( clubb_config_flags_type ), intent(in) :: &
      clubb_config_flags ! Derived type holding all configurable CLUBB flags

    type (stats_metadata_type), intent(in) :: &
      stats_metadata

    !--------------------------- Input/Output Variables ---------------------------
    type (stats), target, intent(inout), dimension(ngrdcol) :: &
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
      pdf_params,    & ! Fortran structure of PDF parameters on thermodynamic levels    [units vary]
      pdf_params_zm    ! Fortran structure of PDF parameters on momentum levels        [units vary]

    type(implicit_coefs_terms), intent(inout) :: &
      pdf_implicit_coefs_terms    ! Implicit coefs / explicit terms [units vary]

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
      Lscale,      & ! Length scale                          [m]
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

    ! "base" smoothing magnitude before scaling for the respective data structure.
    ! See https://github.com/larson-group/clubb/issues/965#issuecomment-1119816722
    ! for a plot on how output behaves with varying min_max_smth_mag
    real( kind = core_rknd ), parameter :: &
      min_max_smth_mag = 1.0e-9_core_rknd

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

    integer, intent(out) :: &
      err_code_out  ! Error code indicator

    type(pdf_parameter) :: pdf_params_single_col(ngrdcol), &
                           pdf_params_zm_single_col(ngrdcol)

    integer :: advance_order_loop_iter

    integer :: smth_type = 2  ! Used for Lscale_width_vert_avg

    !----- Begin Code -----

    !$acc enter data create( Skw_zm, Skw_zt, thvm, thvm_zm, ddzm_thvm_zm, rtprcp, rcp2, &
    !$acc              ddzt_um, ddzt_vm, ddzt_umvm_sqd, ddzt_umvm_sqd_clipped, &
    !$acc              wpthlp2, wprtp2, wprtpthlp, wp2rcp, wp3_zm, Lscale, Lscale_up, Lscale_zm, &
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

    err_code_out = clubb_no_error  ! Initialize to no error value

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
    if ( clubb_at_least_debug_level( 2 ) ) then

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
               edsclrm_forcing(i,:,:) )                                                    ! In

      end do

      if ( err_code == clubb_fatal_error ) then
        write(fstderr,*) "Fatal error when testing input"
        err_code_out = err_code
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
        wpthlp(i,1) = wpthlp_sfc(i)
        wprtp(i,1)  = wprtp_sfc(i)
        upwp(i,1)   = upwp_sfc(i)
        vpwp(i,1)   = vpwp_sfc(i)
      end do
      !$acc end parallel loop

      if ( clubb_config_flags%l_linearize_pbl_winds ) then
        !$acc parallel loop gang vector default(present)
        do i = 1, ngrdcol
          upwp_pert(i,1) = upwp_sfc_pert(i)
          vpwp_pert(i,1) = vpwp_sfc_pert(i)
        end do
        !$acc end parallel loop
      endif ! l_linearize_pbl_winds

      ! Set fluxes for passive scalars (if enabled)
      if ( sclr_dim > 0 ) then
        !$acc parallel loop gang vector collapse(2) default(present)
        do sclr = 1, sclr_dim
          do i = 1, ngrdcol
            wpsclrp(i,1,sclr)   = wpsclrp_sfc(i,sclr)
          end do
        end do
        !$acc end parallel loop
      end if

      if ( edsclr_dim > 0 ) then
        !$acc parallel loop gang vector collapse(2) default(present)
        do edsclr = 1, edsclr_dim
          do i = 1, ngrdcol
            wpedsclrp(i,1,edsclr) = wpedsclrp_sfc(i,edsclr)
          end do
        end do
        !$acc end parallel loop
      end if

    else

      !$acc parallel loop gang vector default(present)
      do i = 1, ngrdcol
        wpthlp(i,1) = 0.0_core_rknd
        wprtp(i,1)  = 0.0_core_rknd
        upwp(i,1)   = 0.0_core_rknd
        vpwp(i,1)   = 0.0_core_rknd
      end do
      !$acc end parallel loop

      ! Set fluxes for passive scalars (if enabled)
      if ( sclr_dim > 0 ) then
        !$acc parallel loop gang vector collapse(2) default(present)
        do sclr = 1, sclr_dim
          do i = 1, ngrdcol
            wpsclrp(i,1,sclr) = 0.0_core_rknd
          end do
        end do
        !$acc end parallel loop
      end if

      if ( edsclr_dim > 0 ) then
        !$acc parallel loop gang vector collapse(2) default(present)
        do edsclr = 1, edsclr_dim
          do i = 1, ngrdcol
            wpedsclrp(i,1,edsclr) = 0.0_core_rknd
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
                               ! rfrzm, hydromet,                             &
                               wphydrometp,                                 & ! Intent(in)
                               wp2hmp, rtphmp_zt, thlphmp_zt,               & ! Intent(in)
                               sclrm, wpsclrp, sclrp2,                      & ! Intent(in)
                               sclrprtp, sclrpthlp, sclrp3,                 & ! Intent(in)
                               p_sfc, l_samp_stats_in_pdf_call,             & ! Intent(in)
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
                               pdf_params, pdf_params_zm,                   & ! Intent(inout)
#ifdef GFDL
                               RH_crit(k, : , :),                           & ! Intent(inout)
                               do_liquid_only_in_clubb,                     & ! Intent(in)
#endif
                               rcm, cloud_frac,                             & ! Intent(out)
                               ice_supersat_frac, wprcp,                    & ! Intent(out)
                               sigma_sqd_w, wpthvp, wp2thvp,                & ! Intent(out)
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
      
    endif ! clubb_config_flags%ipdf_call_placement == ipdf_pre_advance_fields
          ! or clubb_config_flags%ipdf_call_placement
          !    == ipdf_pre_post_advance_fields

    ! Interpolate wp3 to momentum levels, and wp2 to thermodynamic levels
    ! and then compute Skw for m & t grid.
    wp2_zt(:,:) = zm2zt( nzm, nzt, ngrdcol, gr, wp2(:,:), &
                         w_tol_sqd )  ! Positive definite quantity
    wp3_zm(:,:) = zt2zm( nzm, nzt, ngrdcol, gr, wp3(:,:) )

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

    a3_coef_zt(:,:) = zm2zt( nzm, nzt, ngrdcol, gr, &
                             a3_coef(:,:) )

    ! Interpolate thlp2, rtp2, and rtpthlp to thermodynamic levels.
    thlp2_zt(:,:)   = zm2zt( nzm, nzt, ngrdcol, gr, &
                             thlp2(:,:), &
                             thl_tol**2 )  ! Positive def. quantity
    rtp2_zt(:,:)    = zm2zt( nzm, nzt, ngrdcol, gr, &
                             rtp2(:,:), &
                             rt_tol**2 )   ! Positive def. quantity
    rtpthlp_zt(:,:) = zm2zt( nzm, nzt, ngrdcol, gr, &
                             rtpthlp(:,:) )

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
    wp3_on_wp2(:,:) = zt2zm( nzm, nzt, ngrdcol, gr, wp3_on_wp2_zt(:,:) )

    ! Smooth again as above
    wp3_on_wp2_zt(:,:) = zm2zt( nzm, nzt, ngrdcol, gr, &
                                wp3_on_wp2(:,:) )

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

    sqrt_em_zt(:,:) = zm2zt( nzm, nzt, ngrdcol, gr, em(:,:), &
                             em_min )

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
                                      clubb_params(:,ibv_efold),                          & ! In
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
                                  stats_zt,                                          & ! In/Out
                                  Lscale, Lscale_up, Lscale_down )                     ! Out

      if ( clubb_at_least_debug_level( 0 ) ) then
        if ( err_code == clubb_fatal_error ) then
          err_code_out = err_code
          write(fstderr,*) "Error calling calc_Lscale_directly"
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

      tau_zm(:,:) = zt2zm( nzm, nzt, ngrdcol, gr, Lscale(:,:), zero_threshold )
          
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm
        do i = 1, ngrdcol
          tau_zm(i,k) = min( tau_zm(i,k) / sqrt( max( em_min, em(i,k) ) ),  &
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
                        stats_zm,                                                 & ! Inout
                        invrs_tau_zt, invrs_tau_zm,                               & ! Out
                        invrs_tau_sfc, invrs_tau_no_N2_zm, invrs_tau_bkgnd,       & ! Out
                        invrs_tau_shear, invrs_tau_N2_iso,                        & ! Out
                        invrs_tau_wp2_zm, invrs_tau_xp2_zm,                       & ! Out
                        invrs_tau_wp3_zm, invrs_tau_wp3_zt, invrs_tau_wpxp_zm,    & ! Out
                        tau_max_zm, tau_max_zt, tau_zm, tau_zt,                   & ! Out
                        Lscale, Lscale_up, Lscale_down )                            ! Out

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

    Lscale_zm(:,:) = zt2zm( nzm, nzt, ngrdcol, gr, Lscale(:,:) )

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
                          clubb_params(:,iup2_sfc_coef),                  & ! Intent(in)
                          clubb_params(:,ia_const),                       & ! Intent(in)
                          stats_metadata,                                 & ! Intent(in)
                          stats_zm,                                       & ! Intent(inout)
                          wp2, up2, vp2,                                  & ! Intent(inout)
                          thlp2, rtp2, rtpthlp,                           & ! Intent(inout)
                          sclrp2, sclrprtp, sclrpthlp )                     ! Intent(inout)

    if ( clubb_at_least_debug_level( 0 ) ) then
      if ( err_code == clubb_fatal_error ) then
        err_code_out = err_code
        write(fstderr, *) "Error calling calc_sfc_varnce"
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
        
        rsat = sat_mixrat_liq( nzt, ngrdcol, p_in_Pa, &
                               thlm2T_in_K( nzt, ngrdcol, thlm, exner, rcm ), &
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
      w_1_zm(:,:)        = zt2zm( nzm, nzt, ngrdcol, gr, pdf_params%w_1(:,:) )
      w_2_zm(:,:)        = zt2zm( nzm, nzt, ngrdcol, gr, pdf_params%w_2(:,:) )
      varnce_w_1_zm(:,:) = zt2zm( nzm, nzt, ngrdcol, gr, pdf_params%varnce_w_1(:,:) )
      varnce_w_2_zm(:,:) = zt2zm( nzm, nzt, ngrdcol, gr, pdf_params%varnce_w_2(:,:) )
      mixt_frac_zm(:,:)  = zt2zm( nzm, nzt, ngrdcol, gr, pdf_params%mixt_frac(:,:) )
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
                                      clubb_params(:,ibv_efold),                          & ! In
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
                                      stats_metadata,                                     & ! In
                                      stats_zm,                                           & ! In/Out
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
                            fcor, um_ref, vm_ref, up2, vp2,                        & ! intent(in)
                            uprcp, vprcp, rc_coef_zm,                              & ! intent(in)
                            clubb_params, nu_vert_res_dep,                         & ! intent(in)
                            clubb_config_flags%iiPDF_type,                         & ! intent(in)
                            clubb_config_flags%penta_solve_method,                 & ! intent(in)
                            clubb_config_flags%tridiag_solve_method,               & ! intent(in)
                            clubb_config_flags%l_predict_upwp_vpwp,                & ! intent(in)
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
                            um_pert, vm_pert, upwp_pert, vpwp_pert )                 ! intent(i/o)

      if ( clubb_at_least_debug_level( 0 ) ) then
         if ( err_code == clubb_fatal_error ) then
            err_code_out = err_code
            write(fstderr,*) "Error calling advance_xm_wpxp"
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
                             dt_advance,                                          & ! intent(in)
                             sclrm, wpsclrp,                                      & ! intent(in)
                             wpsclrp2, wpsclrprtp, wpsclrpthlp,                   & ! intent(in)
                             lhs_splat_wp2,                                       & ! intent(in)
                             clubb_params, nu_vert_res_dep,                       & ! intent(in)
                             clubb_config_flags%iiPDF_type,                       & ! intent(in)
                             clubb_config_flags%tridiag_solve_method,             & ! intent(in)
                             clubb_config_flags%l_predict_upwp_vpwp,              & ! intent(in)
                             clubb_config_flags%l_min_xp2_from_corr_wx,           & ! intent(in)
                             clubb_config_flags%l_C2_cloud_frac,                  & ! intent(in)
                             clubb_config_flags%l_upwind_xpyp_ta,                 & ! intent(in)
                             clubb_config_flags%l_godunov_upwind_xpyp_ta,         & ! intent(in)
                             clubb_config_flags%l_lmm_stepping,                   & ! intent(in)
                             stats_metadata,                                      & ! In
                             stats_zt, stats_zm, stats_sfc,                       & ! intent(inout)
                             rtp2, thlp2, rtpthlp, up2, vp2,                      & ! intent(inout)
                             sclrp2, sclrprtp, sclrpthlp)                           ! intent(inout)
      
      if ( clubb_at_least_debug_level( 0 ) ) then
         if ( err_code == clubb_fatal_error ) then
            err_code_out = err_code
            write(fstderr,*) "Error calling advance_xp2_xpyp"
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
                            sfc_elevation, sigma_sqd_w, wm_zm,                    & ! intent(in)
                            wm_zt, a3_coef, a3_coef_zt, wp3_on_wp2,               & ! intent(in)
                            wpup2, wpvp2, wp2up2, wp2vp2, wp4,                    & ! intent(in)
                            wpthvp, wp2thvp, um, vm, upwp, vpwp,                  & ! intent(in)
                            up2, vp2, em, Kh_zm, Kh_zt, invrs_tau_C4_zm,          & ! intent(in)
                            invrs_tau_wp3_zt, invrs_tau_C1_zm, Skw_zm,            & ! intent(in)
                            Skw_zt, rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm,        & ! intent(in)
                            invrs_rho_ds_zt, radf, thv_ds_zm,                     & ! intent(in)
                            thv_ds_zt, pdf_params%mixt_frac, Cx_fnc_Richardson,   & ! intent(in)
                            lhs_splat_wp2, lhs_splat_wp3,                         & ! intent(in)
                            pdf_implicit_coefs_terms,                             & ! intent(in)
                            wprtp, wpthlp, rtp2, thlp2,                           & ! intent(in)
                            clubb_params, nu_vert_res_dep,                        & ! intent(in)
                            clubb_config_flags%iiPDF_type,                        & ! intent(in)
                            clubb_config_flags%penta_solve_method,                & ! intent(in)
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
                            stats_metadata,                                       & ! intent(in)
                            stats_zt, stats_zm, stats_sfc,                        & ! intent(inout)
                            wp2, wp3, wp3_zm, wp2_zt )                              ! intent(inout)

      if ( clubb_at_least_debug_level( 0 ) ) then
         if ( err_code == clubb_fatal_error ) then
            err_code_out = err_code
            write(fstderr,*) "Error calling advance_wp2_wp3"
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
                                  nu_vert_res_dep,                            & ! intent(in)
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
                                  um_pert, vm_pert, upwp_pert, vpwp_pert )      ! intent(inout)

      if ( edsclr_dim > 1 .and. clubb_config_flags%l_do_expldiff_rtm_thlm ) then

        call pvertinterp( nzt, ngrdcol,                     & ! intent(in)
                          p_in_Pa, 70000.0_core_rknd, thlm, & ! intent(in)
                          thlm700 )                           ! intent(out)

        call pvertinterp( nzt, ngrdcol,                       & ! intent(in)
                          p_in_Pa, 100000.0_core_rknd, thlm,  & ! intent(in)
                          thlm1000 )                            ! intent(out)

        !$acc parallel loop gang vector collapse(2) default(present)
        do k = 1, nzt
          do i = 1, ngrdcol         
            if ( thlm700(i) - thlm1000(i) < 20.0_core_rknd ) then
              thlm(i,k) = edsclrm(i,k,edsclr_dim-1)
              rtm(i,k) = edsclrm(i,k,edsclr_dim)
            end if
          end do
        end do
        !$acc end parallel loop

      end if

      ! Eric Raut: this seems dangerous to call without any attached flag.
      ! Hence the preprocessor.
#ifdef CLUBB_CAM
      do edsclr=1,edsclr_dim
        ! upper_hf_level = nzt since we are filling the zt levels
        call fill_holes_vertical( nzt, ngrdcol, zero_threshold, 1, nzt, & ! In
                                  gr%dzt, rho_ds_zt,                    & ! In
                                  edsclrm(:,:,edsclr) )                    ! InOut
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

      upwp_zt(:,:) = zm2zt( nzm, nzt, ngrdcol, gr, &
                            upwp(:,:) )
      vpwp_zt(:,:) = zm2zt( nzm, nzt, ngrdcol, gr, &
                            vpwp(:,:) )
      up2_zt(:,:)  = zm2zt( nzm, nzt, ngrdcol, gr, up2(:,:), &
                            w_tol_sqd ) ! Positive def. quantity
      vp2_zt(:,:)  = zm2zt( nzm, nzt, ngrdcol, gr, vp2(:,:), &
                            w_tol_sqd ) ! Positive def. quantity

      thvm_zm(:,:)                   = zt2zm( nzm, nzt, ngrdcol, gr, thvm(:,:), zero_threshold )
      ddzm_thvm_zm(:,:)              = ddzm( nzm, nzt, ngrdcol, gr, thvm_zm(:,:) )
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

      wpthlp_zt(:,:) = zm2zt( nzm, nzt, ngrdcol, gr, &
                              wpthlp(:,:) )
      wprtp_zt(:,:)  = zm2zt( nzm, nzt, ngrdcol, gr, &
                              wprtp(:,:) )
      thlp2_zt(:,:)  = zm2zt( nzm, nzt, ngrdcol, gr, thlp2(:,:), &
                              thl_tol**2 ) ! Positive def. quantity
      rtp2_zt(:,:)   = zm2zt( nzm, nzt, ngrdcol, gr, rtp2(:,:), &
                              rt_tol**2 )   ! Positive def. quantity

      upwp_zt(:,:) = zm2zt( nzm, nzt, ngrdcol, gr, &
                            upwp(:,:) )
      vpwp_zt(:,:) = zm2zt( nzm, nzt, ngrdcol, gr, &
                            vpwp(:,:) )
      up2_zt(:,:)  = zm2zt( nzm, nzt, ngrdcol, gr, up2(:,:), &
                            w_tol_sqd ) ! Positive def. quantity
      vp2_zt(:,:)  = zm2zt( nzm, nzt, ngrdcol, gr, vp2(:,:), &
                            w_tol_sqd ) ! Positive def. quantity

      if ( clubb_config_flags%iiPDF_type == iiPDF_ADG1 ) then

        ! Use the Larson and Golaz (2005) ansatz for the ADG1 PDF to
        ! calculate <rt'^3>, <thl'^3>, <u'^3>, <v'^3>, and <sclr'^3>.
        sigma_sqd_w_zt(:,:) = zm2zt( nzm, nzt, ngrdcol, gr, sigma_sqd_w(:,:), &
                                     zero_threshold )

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
          
          wpsclrp_zt(:,:) = zm2zt( nzm, nzt, ngrdcol, gr, wpsclrp(:,:,sclr), &
                                   sclr_tol(sclr)**2 )
          sclrp2_zt(:,:)  = zm2zt( nzm, nzt, ngrdcol, gr, sclrp2(:,:,sclr), &
                                   sclr_tol(sclr)**2 )

          call xp3_LG_2005_ansatz( nzt, ngrdcol, Skw_zt, wpsclrp_zt, wp2_zt, &
                                   sclrp2_zt, sigma_sqd_w_zt, &
                                   clubb_params, sclr_tol(sclr), &
                                   sclrp3 )

        enddo ! sclr = 1, sclr_dim

      else ! clubb_config_flags%iiPDF_type /= iiPDF_ADG1

        ! Use a modified form of the Larson and Golaz (2005) ansatz for the
        ! ADG1 PDF to calculate <u'^3> and <v'^3> for another type of PDF.
        thvm_zm(:,:)                   = zt2zm( nzm, nzt, ngrdcol, gr, thvm(:,:), zero_threshold )
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
          
          wpsclrp_zt(:,:) = zm2zt( nzm, nzt, ngrdcol, gr, &
                                   wpsclrp(:,:,sclr) )
          sclrp2_zt(:,:)  = zm2zt( nzm, nzt, ngrdcol, gr, sclrp2(:,:,sclr), &
                                   sclr_tol(sclr)**2 )

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
                               ! rfrzm, hydromet,                             &
                               wphydrometp,                                 & ! Intent(in)
                               wp2hmp, rtphmp_zt, thlphmp_zt,               & ! Intent(in)
                               sclrm, wpsclrp, sclrp2,                      & ! Intent(in)
                               sclrprtp, sclrpthlp, sclrp3,                 & ! Intent(in)
                               p_sfc, l_samp_stats_in_pdf_call,             & ! Intent(in)
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
                               pdf_params, pdf_params_zm,                   & ! Intent(inout)
#ifdef GFDL
                               RH_crit(k, : , :),                           & ! Intent(inout)
                               do_liquid_only_in_clubb,                     & ! Intent(in)
#endif
                               rcm, cloud_frac,                             & ! Intent(out)
                               ice_supersat_frac, wprcp,                    & ! Intent(out)
                               sigma_sqd_w, wpthvp, wp2thvp,                & ! Intent(out)
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
      !$acc              pdf_params_zm%ice_supersat_frac_2, &
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
        wpthlp_zt(:,:)  = zm2zt( nzm, nzt, ngrdcol, gr, &
                                 wpthlp(:,:) )
      end if

      if ( stats_metadata%iwprtp_zt > 0 ) then
        wprtp_zt(:,:)   = zm2zt( nzm, nzt, ngrdcol, gr, &
                                 wprtp(:,:) )
      end if

      if ( stats_metadata%iup2_zt > 0 ) then
        up2_zt(:,:) = zm2zt( nzm, nzt, ngrdcol, gr, up2(:,:), &
                             w_tol_sqd )
      end if

      if (stats_metadata%ivp2_zt > 0 ) then
        vp2_zt(:,:) = zm2zt( nzm, nzt, ngrdcol, gr, vp2(:,:), &
                             w_tol_sqd )
      end if

      if ( stats_metadata%iupwp_zt > 0 ) then
        upwp_zt(:,:) = zm2zt( nzm, nzt, ngrdcol, gr, &
                              upwp(:,:) )
      end if

      if ( stats_metadata%ivpwp_zt > 0 ) then
        vpwp_zt(:,:) = zm2zt( nzm, nzt, ngrdcol, gr, &
                              vpwp(:,:) )
      end if
      
      do i = 1, ngrdcol
        
        ! Allocate arrays in single column versions of pdf_params
        call init_pdf_params( nzt, 1, pdf_params_single_col(i) )
        call init_pdf_params( nzm, 1, pdf_params_zm_single_col(i) )
        
        ! Copy multicolumn pdf_params to single column version  
        call copy_multi_pdf_params_to_single( pdf_params, i, &
                                              pdf_params_single_col(i) )
                                              
        call copy_multi_pdf_params_to_single( pdf_params_zm, i, &
                                              pdf_params_zm_single_col(i) )
        
        call stats_accumulate( &
               nzm, nzt, sclr_dim, edsclr_dim,                                            & ! In
               gr%invrs_dzm(i,:), gr%zt(i,:), gr%dzm(i,:), gr%dzt(i,:), dt,               & ! In
               um(i,:), vm(i,:), upwp(i,:), vpwp(i,:), up2(i,:), vp2(i,:),                & ! In
               thlm(i,:), rtm(i,:), wprtp(i,:), wpthlp(i,:),                              & ! In
               wp2(i,:), wp3(i,:), rtp2(i,:), rtp3(i,:), thlp2(i,:), thlp3(i,:),          & ! In
               rtpthlp(i,:),                                                              & ! In
               wpthvp(i,:), wp2thvp(i,:), rtpthvp(i,:), thlpthvp(i,:),                    & ! In
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
               pdf_params_single_col(i), pdf_params_zm_single_col(i),                     & ! In
               sclrm(i,:,:), sclrp2(i,:,:),                                               & ! In
               sclrprtp(i,:,:), sclrpthlp(i,:,:), sclrm_forcing(i,:,:), sclrpthvp(i,:,:), & ! In
               wpsclrp(i,:,:), sclrprcp(i,:,:), wp2sclrp(i,:,:), wpsclrp2(i,:,:),         & ! In
               wpsclrprtp(i,:,:),                                                         & ! In
               wpsclrpthlp(i,:,:), wpedsclrp(i,:,:), edsclrm(i,:,:),                      & ! In
               edsclrm_forcing(i,:,:),                                                    & ! In
               clubb_config_flags%saturation_formula,                                     & ! In
               stats_metadata,                                                            & ! In
               stats_zt(i), stats_zm(i), stats_sfc(i) )                                     ! In/Out
      end do
    endif ! stats_metadata%l_stats_samp

    if ( clubb_at_least_debug_level( 2 ) ) then

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
             edsclrm_forcing(i,:,:) )                                                    ! In
      end do

      if ( err_code == clubb_fatal_error ) then
        write(fstderr,*) "Error occurred during parameterization_check at"// &
                         " end of advance_clubb_core"
        err_code_out = err_code
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
          rtm_flux_top(i) = rho_ds_zm(i,nzm) * wprtp(i,nzm)

          if ( .not. clubb_config_flags%l_host_applies_sfc_fluxes ) then
            rtm_flux_sfc(i) = rho_ds_zm(i,1) * wprtp_sfc(i)
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
          thlm_flux_top(i) = rho_ds_zm(i,nzm) * wpthlp(i,nzm)

          if ( .not. clubb_config_flags%l_host_applies_sfc_fluxes ) then
            thlm_flux_sfc(i) = rho_ds_zm(i,1) * wpthlp_sfc(i)
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
    !$acc                   wpthlp2, wprtp2, wprtpthlp, wp2rcp, wp3_zm, Lscale, Lscale_up, &
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

  !=============================================================================
  subroutine pdf_closure_driver( gr, nzm, nzt, ngrdcol,                   & ! Intent(in)
                                 dt, hydromet_dim, sclr_dim, sclr_tol,    & ! Intent(in)
                                 wprtp, thlm, wpthlp, rtp2, rtp3,         & ! Intent(in)
                                 thlp2, thlp3, rtpthlp, wp2,              & ! Intent(in)
                                 wp3, wm_zm, wm_zt,                       & ! Intent(in)
                                 um, up2, upwp, up3,                      & ! Intent(in)
                                 vm, vp2, vpwp, vp3,                      & ! Intent(in)
                                 p_in_Pa, exner,                          & ! Intent(in)
                                 thv_ds_zm, thv_ds_zt, rtm_ref,           & ! Intent(in)
!                                rfrzm, hydromet,                         &
                                 wphydrometp,                             & ! Intent(in)
                                 wp2hmp, rtphmp_zt, thlphmp_zt,           & ! Intent(in)
                                 sclrm, wpsclrp, sclrp2,                  & ! Intent(in)
                                 sclrprtp, sclrpthlp, sclrp3,             & ! Intent(in)
                                 p_sfc, l_samp_stats_in_pdf_call,         & ! Intent(in)
                                 clubb_params,                            & ! Intent(in)
                                 iiPDF_type,                              & ! Intent(in)
                                 saturation_formula,                      & ! Intent(in)
                                 l_predict_upwp_vpwp,                     & ! Intent(in)
                                 l_rtm_nudge,                             & ! Intent(in)
                                 l_trapezoidal_rule_zt,                   & ! Intent(in)
                                 l_trapezoidal_rule_zm,                   & ! Intent(in)
                                 l_call_pdf_closure_twice,                & ! Intent(in)
                                 l_use_cloud_cover,                       & ! Intent(in)
                                 l_rcm_supersat_adj,                      & ! Intent(in
                                 l_mix_rat_hm,                            & ! Intent(in)
                                 stats_metadata,                          & ! Intent(in)
                                 stats_zt, stats_zm,                      & ! Intent(inout)
                                 rtm,                                     & ! Intent(inout)
                                 pdf_implicit_coefs_terms,                & ! Intent(inout)
                                 pdf_params, pdf_params_zm,               & ! Intent(inout)
#ifdef GFDL
                                 RH_crit(k, : , :),                       & ! Intent(inout)
                                 do_liquid_only_in_clubb,                 & ! Intent(in)
#endif
                                 rcm, cloud_frac,                         & ! Intent(out)
                                 ice_supersat_frac, wprcp,                & ! Intent(out)
                                 sigma_sqd_w, wpthvp, wp2thvp,            & ! Intent(out)
                                 rtpthvp, thlpthvp, rc_coef,              & ! Intent(out)
                                 rcm_in_layer, cloud_cover,               & ! Intent(out)
                                 rcp2_zt, thlprcp,                        & ! Intent(out)
                                 rc_coef_zm, sclrpthvp,                   & ! Intent(out)
                                 wpup2, wpvp2,                            & ! Intent(out)
                                 wp2up2, wp2vp2, wp4,                     & ! Intent(out)
                                 wp2rtp, wprtp2, wp2thlp,                 & ! Intent(out)
                                 wpthlp2, wprtpthlp, wp2rcp,              & ! Intent(out)
                                 rtprcp, rcp2,                            & ! Intent(out)
                                 uprcp, vprcp,                            & ! Intent(out)
                                 w_up_in_cloud, w_down_in_cloud,          & ! Intent(out)
                                 cloudy_updraft_frac,                     & ! Intent(out)
                                 cloudy_downdraft_frac,                   & ! intent(out)
                                 Skw_velocity,                            & ! Intent(out)
                                 cloud_frac_zm,                           & ! Intent(out)
                                 ice_supersat_frac_zm,                    & ! Intent(out)
                                 rtm_zm, thlm_zm, rcm_zm,                 & ! Intent(out)
                                 rcm_supersat_adj,                        & ! Intent(out)
                                 wp2sclrp, wpsclrp2, sclrprcp,            & ! Intent(out)
                                 wpsclrprtp, wpsclrpthlp )                  ! Intent(out)

    use grid_class, only: &
        grid, & ! Type
        zt2zm, & ! Procedure(s)
        zm2zt, &
        zm2zt2zm

    use constants_clubb, only: &
        one_half,       & ! Variable(s)
        w_tol,          & 
        w_tol_sqd,      &
        rt_tol,         &
        thl_tol,        &
        p0,             &
        kappa,          &
        fstderr,        &
        zero,           &
        zero_threshold, &
        eps

    use pdf_parameter_module, only: &
        pdf_parameter,        & ! Variable Type
        implicit_coefs_terms, &  ! Variable Type
        init_pdf_implicit_coefs_terms ! Procedure

    use parameters_model, only: &
        ts_nudge,               & ! Variable(s)
        rtm_min,                &
        rtm_nudge_max_altitude

    use parameter_indices, only: &
        nparams,         & ! Variable(s)
        igamma_coef,     &
        igamma_coefb,    &
        igamma_coefc

    use pdf_closure_module, only: &
        pdf_closure ! Procedure(s)

    use Skx_module, only: &
        Skx_func    ! Procedure(s)

    use sigma_sqd_w_module, only: &
        compute_sigma_sqd_w    ! Procedure(s)

    use pdf_utilities, only: &
        compute_mean_binormal    ! Procedure(s)

    use T_in_K_module, only: &
        thlm2T_in_K    ! Procedure(s)

    use saturation, only:  &
        sat_mixrat_liq    ! Procedure(s)

    use model_flags, only: &
        l_gamma_Skw,      & ! Variable(s)
        iiPDF_new,        & ! new PDF
        iiPDF_new_hybrid    ! new hybrid PDF

    use error_code, only: &
        clubb_at_least_debug_level,  & ! Procedure
        err_code,                    & ! Error Indicator
        clubb_fatal_error              ! Constant

    use stats_type_utilities, only: &
        stat_update_var,    & ! Procedure(s)
        stat_update_var_pt

    use stats_variables, only: &
        stats_metadata_type

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    use stats_type, only: stats ! Type

    implicit none

    !------------------------------- Input Variables -------------------------------
    type (grid), target, intent(in) :: &
      gr

    integer, intent(in) :: &
      nzm, &
      nzt, &
      ngrdcol

    real( kind = core_rknd ), intent(in) ::  &
      dt  ! Current timestep duration    [s]

    integer, intent(in) :: &
      hydromet_dim,   & ! Total number of hydrometeor species       [#]
      sclr_dim          ! Number of passive scalars                 [#]

    real( kind = core_rknd ), intent(in), dimension(sclr_dim) :: & 
      sclr_tol          ! Threshold(s) on the passive scalars  [units vary]

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) ::  &
      !rtm,       & ! total water mixing ratio, r_t (thermo. levels) [kg/kg]
      thlm,      & ! liq. water pot. temp., th_l (thermo. levels)   [K]
      rtp3,      & ! r_t'^3 (thermodynamic levels)                  [(kg/kg)^3]
      thlp3,     & ! th_l'^3 (thermodynamic levels)                 [K^3]
      wp3,       & ! w'^3 (thermodynamic levels)                    [m^3/s^3]
      wm_zt,     & ! w mean wind component on thermo. levels        [m/s]
      p_in_Pa,   & ! Air pressure (thermodynamic levels)            [Pa]
      exner,     & ! Exner function (thermodynamic levels)          [-]
      thv_ds_zt, & ! Dry, base-state theta_v on thermo. levs.       [K]
      rtm_ref!,  & ! Initial total water mixing ratio               [kg/kg]
      !rfrzm        ! Total ice-phase water mixing ratio             [kg/kg]

    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) ::  &
      wprtp,     & ! w' r_t' (momentum levels)                      [(kg/kg)m/s]
      wpthlp,    & ! w' th_l' (momentum levels)                     [(m/s) K]
      rtp2,      & ! r_t'^2 (momentum levels)                       [(kg/kg)^2]
      thlp2,     & ! th_l'^2 (momentum levels)                      [K^2]
      rtpthlp,   & ! r_t' th_l' (momentum levels)                   [(kg/kg) K]
      wp2,       & ! w'^2 (momentum levels)                         [m^2/s^2]
      wm_zm,     & ! w mean wind component on momentum levels       [m/s]
      thv_ds_zm    ! Dry, base-state theta_v on momentum levs.      [K]

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) ::  &
      um,          & ! Grid-mean eastward wind     [m/s]
      up3,         & ! u'^3                        [(m/s)^3]
      vm,          & ! Grid-mean northward wind    [m/s]
      vp3            ! v'^3                        [(m/s)^3]

    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) ::  &
      up2,         & ! u'^2                        [(m/s)^2]
      upwp,        & ! u'w'                        [(m/s)^2]
      vp2,         & ! v'^2                        [(m/s)^2]
      vpwp           ! v'w'                        [(m/s)^2]

    ! Hydrometeor variables
    !real( kind = core_rknd ), dimension(gr%nz,hydromet_dim), intent(in) :: &
    !hydromet       ! Mean of hydrometeor fields               [units vary]

    real( kind = core_rknd ), dimension(ngrdcol,nzt,hydromet_dim), intent(in) :: &
      wp2hmp,      & ! Third-order moment:  < w'^2 hm' >    [(m/s)^2 <hm units>]
      rtphmp_zt,   & ! Covariance of rt and hm (on t-levs.) [(kg/kg) <hm units>]
      thlphmp_zt     ! Covariance of thl and hm (on t-levs.)      [K <hm units>]

    real( kind = core_rknd ), dimension(ngrdcol,nzm,hydromet_dim), intent(in) :: &
      wphydrometp    ! Covariance of w and a hydrometeor      [(m/s) <hm units>]

    ! Passive scalar variables
    real( kind = core_rknd ), dimension(ngrdcol,nzt,sclr_dim), intent(in) :: &
      sclrm,     & ! Passive scalar mean (thermo. levels) [units vary]
      sclrp3       ! sclr'^3 (thermodynamic levels)       [{units vary}^3]

    real( kind = core_rknd ), dimension(ngrdcol,nzm,sclr_dim), intent(in) :: &
      wpsclrp,   & ! w'sclr' (momentum levels)            [{units vary} m/s]
      sclrp2,    & ! sclr'^2 (momentum levels)            [{units vary}^2]
      sclrprtp,  & ! sclr'rt' (momentum levels)           [{units vary} (kg/kg)]
      sclrpthlp    ! sclr'thl' (momentum levels)          [{units vary} K]

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
      p_sfc        ! Pressure at surface                  [Pa]

    logical, intent(in) :: &
      l_samp_stats_in_pdf_call    ! Sample stats in this call to this subroutine

    real( kind = core_rknd ), dimension(ngrdcol,nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

    integer, intent(in) :: &
      iiPDF_type,      & ! Selected option for the two-component normal (double
                         ! Gaussian) PDF type to use for the w, rt, and theta-l (or
                         ! w, chi, and eta) portion of CLUBB's multivariate,
                         ! two-component PDF.
      saturation_formula ! Integer that stores the saturation formula to be used

    logical, dimension(hydromet_dim), intent(in) :: &
      l_mix_rat_hm   ! if true, then the quantity is a hydrometeor mixing ratio

    logical, intent(in) :: &
      l_predict_upwp_vpwp,      & ! Flag to predict <u'w'> and <v'w'> along with <u> and <v>
                                  ! alongside the advancement of <rt>, <w'rt'>, <thl>, <wpthlp>,
                                  ! <sclr>, and <w'sclr'> in subroutine advance_xm_wpxp.
                                  ! Otherwise, <u'w'> and <v'w'> are still approximated by eddy
                                  ! diffusivity when <u> and <v> are advanced in subroutine
                                  ! advance_windm_edsclrm.
      l_rtm_nudge,              & ! For rtm nudging
      l_trapezoidal_rule_zt,    & ! If true, the trapezoidal rule is called for the
                                  ! thermodynamic-level variables output from pdf_closure.
      l_trapezoidal_rule_zm,    & ! If true, the trapezoidal rule is called for three
                                  ! momentum-level variables – wpthvp, thlpthvp, and rtpthvp -
                                  ! output from pdf_closure.
      l_call_pdf_closure_twice, & ! This logical flag determines whether or not to call subroutine
                                  ! pdf_closure twice.  If true, pdf_closure is called first on
                                  ! thermodynamic levels and then on momentum levels so that each
                                  ! variable is computed on its native level.  If false,
                                  ! pdf_closure is only called on thermodynamic levels, and
                                  ! variables which belong on momentum levels are interpolated.
      l_use_cloud_cover,        & ! Use cloud_cover and rcm_in_layer to help boost cloud_frac and
                                  ! rcm to help increase cloudiness at coarser grid resolutions.
      l_rcm_supersat_adj          ! Add excess supersaturated vapor to cloud water

    type (stats_metadata_type), intent(in) :: &
      stats_metadata

    !------------------------------- InOut Variables -------------------------------
    type (stats), target, dimension(ngrdcol), intent(inout) :: &
      stats_zt, &
      stats_zm

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(inout) ::  &
      rtm    ! total water mixing ratio, r_t (thermo. levels) [kg/kg]

    type(implicit_coefs_terms), intent(inout) :: &
      pdf_implicit_coefs_terms    ! Implicit coefs / explicit terms [units vary]

    ! Variable being passed back to and out of advance_clubb_core.
    type(pdf_parameter), intent(inout) :: &
      pdf_params,    & ! PDF parameters                           [units vary]
      pdf_params_zm    ! PDF parameters                           [units vary]

#ifdef GFDL
    ! hlg, 2010-06-16
    real( kind = core_rknd ), dimension(ngrdcol,nz, min(1,sclr_dim) , 2), intent(inout) :: &
      RH_crit  ! critical relative humidity for droplet and ice nucleation
! ---> h1g, 2012-06-14
    logical, intent(in)                 ::  do_liquid_only_in_clubb
! <--- h1g, 2012-06-14
#endif

    !------------------------------- Output Variables -------------------------------
    ! Variables being passed back to and out of advance_clubb_core.
    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) ::  &
      rcm,               & ! mean r_c (thermodynamic levels)        [kg/kg]
      cloud_frac,        & ! cloud fraction (thermodynamic levels)  [-]
      ice_supersat_frac, & ! ice supersat. frac. (thermo. levels)   [-]
      wp2thvp,           & ! < w'^2 th_v' > (thermodynamic levels)  [m^2/s^2 K]
      rc_coef,           & ! Coefficient of X'r_c' (thermo. levs.)  [K/(kg/kg)]
      rcm_in_layer,      & ! rcm in cloud layer                     [kg/kg]
      cloud_cover,       & ! cloud cover                            [-]
      rcp2_zt              ! r_c'^2 (on thermo. grid)               [kg^2/kg^2]

    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(out) ::  &
      wprcp,             & ! < w'r_c' > (momentum levels)           [m/s kg/kg]
      sigma_sqd_w,       & ! PDF width parameter (momentum levels)  [-]
      wpthvp,            & ! < w' th_v' > (momentum levels)         [kg/kg K]
      rtpthvp,           & ! < r_t' th_v' > (momentum levels)       [kg/kg K]
      thlpthvp,          & ! < th_l' th_v' > (momentum levels)      [K^2]
      thlprcp,           & ! < th_l' r_c' > (momentum levels)       [K kg/kg]
      rc_coef_zm           ! Coefficient of X'r_c' on m-levs.       [K/(kg/kg)]

    ! Variable being passed back to and out of advance_clubb_core.
    real( kind = core_rknd ), dimension(ngrdcol,nzm,sclr_dim), intent(out) :: &
      sclrpthvp    ! < sclr' th_v' > (momentum levels)   [units vary]

    ! Variables being passed back to only advance_clubb_core (for statistics).
    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) ::  &
      wpup2,     & ! < w'u'^2 > (thermodynamic levels)        [m^3/s^3]
      wpvp2,     & ! < w'v'^2 > (thermodynamic levels)        [m^3/s^3]
      wp2rtp,    & ! < w'^2 r_t' > (thermodynamic levels)     [m^2/s^2 kg/kg]
      wprtp2,    & ! < w' r_t'^2 > (thermodynamic levels)     [m/s kg^2/kg^2]
      wp2thlp,   & ! < w'^2 th_l' > (thermodynamic levels)    [m^2/s^2 K]
      wpthlp2,   & ! < w' th_l'^2 > (thermodynamic levels)    [m/s K^2]
      wprtpthlp, & ! < w' r_t' th_l' > (thermodynamic levels) [m/s kg/kg K]
      wp2rcp       ! < w'^2 r_c' > (thermodynamic levels)     [m^2/s^2 kg/kg]

    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(out) ::  &
      wp2up2,    & ! < w'^2u'^2 > (momentum levels)           [m^4/s^4]
      wp2vp2,    & ! < w'^2v'^2 > (momentum levels)           [m^4/s^4]
      wp4,       & ! < w'^4 > (momentum levels)               [m^4/s^4]
      rtprcp,    & ! < r_t' r_c' > (momentum levels)          [kg^2/kg^2]
      rcp2         ! Variance of r_c (momentum levels)        [kg^2/kg^2]

    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(out) ::  &
      uprcp,                 & ! < u' r_c' >                [(m kg)/(s kg)]
      vprcp                    ! < v' r_c' >                [(m kg)/(s kg)]

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) ::  &
      w_up_in_cloud,         & ! mean cloudy updraft vel    [m/s]
      w_down_in_cloud,       & ! mean cloudy downdraft vel  [m/s]
      cloudy_updraft_frac,   & ! cloudy updraft fraction    [-]
      cloudy_downdraft_frac    ! cloudy downdraft fraction  [-]

    ! Variables being passed back to only advance_clubb_core (for statistics).
    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(out) ::  &
      Skw_velocity,         & ! Skewness velocity                        [m/s]
      cloud_frac_zm,        & ! Cloud Fraction on momentum levels        [-]
      ice_supersat_frac_zm, & ! Ice supersat. frac. on momentum levels   [-]
      rtm_zm,               & ! Total water mixing ratio at mom. levs.   [kg/kg]
      thlm_zm,              & ! Liquid water pot. temp. at mom. levs.    [K]
      rcm_zm                  ! rcm at momentum levels                   [kg/kg]

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) ::  &
      rcm_supersat_adj        ! Adjust. to rcm due to spurious supersat. [kg/kg]

    real( kind = core_rknd ), dimension(ngrdcol,nzt,sclr_dim), intent(out) :: &
      wp2sclrp,    & ! < w'^2 sclr' > (thermodynamic levels)      [units vary]
      wpsclrp2,    & ! < w' sclr'^2 > (thermodynamic levels)      [units vary]
      wpsclrprtp,  & ! < w' sclr' r_t' > (thermodynamic levels)   [units vary]
      wpsclrpthlp    ! < w' sclr' th_l' > (thermodynamic levels)  [units vary]

    real( kind = core_rknd ), dimension(ngrdcol,nzm,sclr_dim), intent(out) :: &
      sclrprcp       ! < sclr' r_c' > (momentum levels)           [units vary]

    !------------------------------- Local Variables -------------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      wp2_zt,           & ! wp2 interpolated to thermodynamic levels   [m^2/s^2]
      rtp2_zt,          & ! rtp2 interpolated to thermodynamic levels  [kg^2/kg^2]
      thlp2_zt,         & ! thlp2 interpolated to thermodynamic levels [K^2]
      wprtp_zt,         & ! wprtp interpolated to thermodynamic levels [m/s kg/kg]
      wpthlp_zt,        & ! wpthlp interpolated to thermodynamic levs. [m/s K]
      rtpthlp_zt,       & ! rtpthlp interp. to thermodynamic levels    [kg/kg K]
      up2_zt,           & ! up2 interpolated to thermodynamic levels   [m^2/s^2]
      vp2_zt,           & ! vp2 interpolated to thermodynamic levels   [m^2/s^2]
      upwp_zt,          & ! upwp interpolated to thermodynamic levels  [m^2/s^2]
      vpwp_zt,          & ! vpwp interpolated to thermodynamic levels  [m^2/s^2]
      gamma_Skw_fnc_zt, & ! Gamma as a function of skewness (t-levs.)  [-]
      sigma_sqd_w_zt,   & ! PDF width parameter (thermodynamic levels) [-]
      Skw_zt,           & ! Skewness of w on thermodynamic levels      [-]
      Skrt_zt,          & ! Skewness of rt on thermodynamic levels     [-]
      Skthl_zt,         & ! Skewness of thl on thermodynamic levels    [-]
      Sku_zt,           & ! Skewness of u on thermodynamic levels      [-]
      Skv_zt              ! Skewness of v on thermodynamic levels      [-]

    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
      wp3_zm,           & ! wp3 interpolated to momentum levels        [m^3/s^3]
      rtp3_zm,          & ! rtp3 interpolated to momentum levels       [kg^3/kg^3]
      thlp3_zm,         & ! thlp3 interpolated to momentum levels      [K^3]
      up3_zm,           & ! up3 interpolated to momentum levels        [m^3/s^3]
      vp3_zm,           & ! vp3 interpolated to momentum levels        [m^3/s^3]
      gamma_Skw_fnc,    & ! Gamma as a function of skewness            [-]
      Skw_zm,           & ! Skewness of w on momentum levels           [-]
      Skrt_zm,          & ! Skewness of rt on momentum levels          [-]
      Skthl_zm,         & ! Skewness of thl on momentum levels         [-]
      Sku_zm,           & ! Skewness of u on momentum levels           [-]
      Skv_zm              ! Skewness of v on momentum levels           [-]

    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
      w_up_in_cloud_zm,         & ! Avg. cloudy updraft velocity; m-levs   [m/s]
      w_down_in_cloud_zm,       & ! Avg. cloudy downdraft velocity; m-levs [m/s]
      cloudy_updraft_frac_zm,   & ! cloudy updraft fraction; m-levs        [-]
      cloudy_downdraft_frac_zm    ! cloudy downdraft fraction; m-levs      [-]

    ! Interpolated values for optional second call to PDF closure.
    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
      p_in_Pa_zm, & ! Pressure interpolated to momentum levels  [Pa]
      exner_zm      ! Exner interpolated to momentum levels     [-]

    real( kind = core_rknd ), dimension(ngrdcol,nzt,hydromet_dim) :: &
      wphydrometp_zt    ! Covariance of w and hm (on t-levs.) [(m/s) <hm units>]

    real( kind = core_rknd ), dimension(ngrdcol,nzm,hydromet_dim) :: &
      wp2hmp_zm,      & ! Moment <w'^2 hm'> (on m-levs.)    [(m/s)^2 <hm units>]
      rtphmp,         & ! Covariance of rt and hm           [(kg/kg) <hm units>]
      thlphmp           ! Covariance of thl and hm                [K <hm units>]

    real( kind = core_rknd ), dimension(ngrdcol,nzt,sclr_dim) :: &
      wpsclrp_zt,   & ! w' sclr' interpolated to thermo. levels
      sclrp2_zt,    & ! sclr'^2 interpolated to thermo. levels
      sclrprtp_zt,  & ! sclr' r_t' interpolated to thermo. levels
      sclrpthlp_zt, & ! sclr' th_l' interpolated thermo. levels
      Sksclr_zt       ! Skewness of sclr on thermodynamic levels      [-]

    real( kind = core_rknd ), dimension(ngrdcol,nzm,sclr_dim) :: &
      sclrp3_zm,    & ! sclr'^3 interpolated to momentum levels
      Sksclr_zm       ! Skewness of sclr on momentum levels           [-]

    ! These local variables are declared because they originally belong on the
    ! momentum grid levels, but pdf_closure outputs them on the thermodynamic
    ! grid levels.
    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
      wpup2_zm,    & ! w'u'^2 (on momentum grid)        [m^3/s^3]
      wpvp2_zm       ! w'v'^2 (on momentum grid)        [m^3/s^3]

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      wp2up2_zt,   & ! w'^2u'^2 (on thermo. grid)       [m^4/s^4]
      wp2vp2_zt,   & ! w'^2v'^2 (on thermo. grid)       [m^4/s^4]
      wp4_zt,      & ! w'^4 (on thermo. grid)           [m^4/s^4]
      wpthvp_zt,   & ! Buoyancy flux (on thermo. grid)  [(K m)/s]
      rtpthvp_zt,  & ! r_t' th_v' (on thermo. grid)     [(kg K)/kg]
      thlpthvp_zt, & ! th_l' th_v' (on thermo. grid)    [K^2]
      wprcp_zt,    & ! w' r_c' (on thermo. grid)        [(m kg)/(s kg)]
      rtprcp_zt,   & ! r_t' r_c' (on thermo. grid)      [(kg^2)/(kg^2)]
      thlprcp_zt,  & ! th_l' r_c' (on thermo. grid)     [(K kg)/kg]
      uprcp_zt,    & ! u' r_c' (on thermo. grid)        [(m kg)/(s kg)]
      vprcp_zt       ! v' r_c' (on thermo. grid)        [(m kg)/(s kg)]

    real( kind = core_rknd ), dimension(ngrdcol,nzt,sclr_dim) :: &
      sclrpthvp_zt, & ! sclr'th_v' (on thermo. grid)
      sclrprcp_zt     ! sclr'rc' (on thermo. grid)

    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
      wprtp2_zm,    & ! < w' r_t'^2 > on momentum levels      [m/s kg^2/kg^2]
      wp2rtp_zm,    & ! < w'^2 r_t' > on momentum levels      [m^2/s^2 kg/kg]
      wpthlp2_zm,   & ! < w' th_l'^2 > on momentum levels     [m/s K^2]
      wp2thlp_zm,   & ! < w'^2 th_l' > on momentum levels     [m^2/s^2 K]
      wprtpthlp_zm, & ! < w' r_t' th_l' > on momentum levels  [m/s kg/kg K]
      wp2thvp_zm,   & ! < w'^2 th_v' > on momentum levels     [m^2/s^2 K]
      wp2rcp_zm       ! < w'^2 r_c' > on momentum levles      [m^2/s^2 kg/kg]

    real( kind = core_rknd ), dimension(ngrdcol,nzm,sclr_dim) :: &
      wpsclrprtp_zm,  & ! w'sclr'rt' on momentum grid
      wpsclrp2_zm,    & ! w'sclr'^2 on momentum grid
      wpsclrpthlp_zm, & ! w'sclr'thl' on momentum grid
      wp2sclrp_zm,    & ! w'^2 sclr' on momentum grid
      sclrm_zm          ! Passive scalar mean on momentum grid

    ! Output from new PDF for recording statistics.
    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      F_w,   & ! Parameter for the spread of the PDF component means of w    [-]
      F_rt,  & ! Parameter for the spread of the PDF component means of rt   [-]
      F_thl    ! Parameter for the spread of the PDF component means of thl  [-]

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      min_F_w,   & ! Minimum allowable value of parameter F_w      [-]
      max_F_w,   & ! Maximum allowable value of parameter F_w      [-]
      min_F_rt,  & ! Minimum allowable value of parameter F_rt     [-]
      max_F_rt,  & ! Maximum allowable value of parameter F_rt     [-]
      min_F_thl, & ! Minimum allowable value of parameter F_thl    [-]
      max_F_thl    ! Maximum allowable value of parameter F_thl    [-]

    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
      F_w_zm,       &
      F_rt_zm,      &
      F_thl_zm,     &
      min_F_w_zm,   &
      max_F_w_zm,   &
      min_F_rt_zm,  &
      max_F_rt_zm,  &
      min_F_thl_zm, &
      max_F_thl_zm

    type(implicit_coefs_terms) :: &
      pdf_implicit_coefs_terms_zm

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      rsat,             & ! Saturation mixing ratio from mean rt and thl.
      rel_humidity        ! Relative humidity after PDF closure [-]

    real( kind = core_rknd ) :: &
      gamma_coef,     & ! CLUBB tunable parameter gamma_coef
      gamma_coefb,    & ! CLUBB tunable parameter gamma_coefb
      gamma_coefc       ! CLUBB tunable parameter gamma_coefc
      
    real( kind = core_rknd ), dimension(ngrdcol, nzm) :: &
      um_zm, &
      vm_zm, &
      sigma_sqd_w_tmp

    real( kind = core_rknd ), dimension(ngrdcol, nzt) :: &
      T_in_K

    logical :: l_spur_supersat   ! Spurious supersaturation?

    integer :: i, k, j, sclr

    !-------------------------------------- Begin Code --------------------------------------

    !$acc enter data create( wp2_zt,wp3_zm, rtp2_zt,rtp3_zm, thlp2_zt,  thlp3_zm, &
    !$acc                    wprtp_zt, wpthlp_zt, rtpthlp_zt, up2_zt, up3_zm, &
    !$acc                    vp2_zt, vp3_zm, upwp_zt, vpwp_zt, gamma_Skw_fnc, &
    !$acc                    gamma_Skw_fnc_zt,sigma_sqd_w_zt,  Skw_zt, Skw_zm, &
    !$acc                    Skrt_zt, Skrt_zm, Skthl_zt, Skthl_zm, Sku_zt, &
    !$acc                    Sku_zm, Skv_zt, Skv_zm, wp2up2_zt, &
    !$acc                    wp2vp2_zt, wp4_zt, wpthvp_zt, rtpthvp_zt, thlpthvp_zt, &
    !$acc                    wprcp_zt, rtprcp_zt, thlprcp_zt, uprcp_zt, vprcp_zt, &
    !$acc                    rsat, rel_humidity, um_zm, vm_zm, T_in_K, sigma_sqd_w_tmp )

    !$acc enter data if( l_call_pdf_closure_twice ) &
    !$acc            create( w_up_in_cloud_zm, wpup2_zm, wpvp2_zm, &
    !$acc                    w_down_in_cloud_zm, cloudy_updraft_frac_zm,  &
    !$acc                    cloudy_downdraft_frac_zm, p_in_Pa_zm, exner_zm, &
    !$acc                    wprtp2_zm, wp2rtp_zm, wpthlp2_zm, &
    !$acc                    wp2thlp_zm, wprtpthlp_zm, wp2thvp_zm, wp2rcp_zm )

    !$acc enter data if( sclr_dim > 0 ) &
    !$acc            create( wpsclrp_zt, sclrp2_zt, sclrp3_zm, sclrprtp_zt, sclrpthlp_zt, &
    !$acc                    Sksclr_zt, Sksclr_zm, sclrpthvp_zt, sclrprcp_zt, wpsclrprtp_zm, &
    !$acc                    wpsclrp2_zm, wpsclrpthlp_zm, wp2sclrp_zm, sclrm_zm )

    !$acc enter data if( hydromet_dim > 0 ) create( wphydrometp_zt, wp2hmp_zm, rtphmp, thlphmp )

    !---------------------------------------------------------------------------
    ! Interpolate wp3, rtp3, thlp3, up3, vp3, and sclrp3 to momentum levels, and
    ! wp2, rtp2, thlp2, up2, vp2, and sclrp2 to thermodynamic levels, and then
    ! compute Skw, Skrt, Skthl, Sku, Skv, and Sksclr for both the momentum and
    ! thermodynamic grid levels.
    !---------------------------------------------------------------------------
    wp2_zt(:,:)   = zm2zt( nzm, nzt, ngrdcol, gr, wp2(:,:), &
                           w_tol_sqd ) ! Positive definite quantity
    wp3_zm(:,:)   = zt2zm( nzm, nzt, ngrdcol, gr, wp3(:,:) )
    thlp2_zt(:,:) = zm2zt( nzm, nzt, ngrdcol, gr, thlp2(:,:), &
                           thl_tol**2 ) ! Positive definite quantity
    thlp3_zm(:,:) = zt2zm( nzm, nzt, ngrdcol, gr, thlp3(:,:) )
    rtp2_zt(:,:)  = zm2zt( nzm, nzt, ngrdcol, gr, rtp2(:,:), &
                           rt_tol**2 ) ! Positive definite quantity
    rtp3_zm(:,:)  = zt2zm( nzm, nzt, ngrdcol, gr, rtp3(:,:) )
    up2_zt(:,:)   = zm2zt( nzm, nzt, ngrdcol, gr, up2(:,:), &
                           w_tol_sqd ) ! Positive definite quantity
    up3_zm(:,:)   = zt2zm( nzm, nzt, ngrdcol, gr, up3(:,:) )
    vp2_zt(:,:)   = zm2zt( nzm, nzt, ngrdcol, gr, vp2(:,:), &
                           w_tol_sqd ) ! Positive definite quantity
    vp3_zm(:,:)   = zt2zm( nzm, nzt, ngrdcol, gr, vp3(:,:) )

    do sclr = 1, sclr_dim, 1
      sclrp2_zt(:,:,sclr) = zm2zt( nzm, nzt, ngrdcol, gr, sclrp2(:,:,sclr), &
                                   sclr_tol(sclr)**2 ) ! Pos. def. quantity
      sclrp3_zm(:,:,sclr) = zt2zm( nzm, nzt, ngrdcol, gr, sclrp3(:,:,sclr) )
    end do ! sclr = 1, sclr_dim, 1

    call Skx_func( nzt, ngrdcol, wp2_zt, wp3, &
                   w_tol, clubb_params, &
                   Skw_zt )
                   
    call Skx_func( nzm, ngrdcol, wp2, wp3_zm, &
                   w_tol, clubb_params, &
                   Skw_zm )    
                   
    call Skx_func( nzt, ngrdcol, thlp2_zt, thlp3, &
                   thl_tol, clubb_params, &
                   Skthl_zt )  
                   
    call Skx_func( nzm, ngrdcol, thlp2, thlp3_zm, &
                   thl_tol, clubb_params, &
                   Skthl_zm )  
                   
    call Skx_func( nzt, ngrdcol, rtp2_zt, rtp3, &
                   rt_tol, clubb_params, &
                   Skrt_zt )   
                   
    call Skx_func( nzm, ngrdcol, rtp2, rtp3_zm, &
                   rt_tol, clubb_params, &
                   Skrt_zm )   
                   
    call Skx_func( nzt, ngrdcol, up2_zt, up3, &
                   w_tol, clubb_params, &
                   Sku_zt )   
                                      
    call Skx_func( nzm, ngrdcol, up2, up3_zm, &
                   w_tol, clubb_params, &
                   Sku_zm )   
                   
    call Skx_func( nzt, ngrdcol, vp2_zt, vp3, &
                   w_tol, clubb_params, &
                   Skv_zt )   
                   
    call Skx_func( nzm, ngrdcol, vp2, vp3_zm, &
                   w_tol, clubb_params, &
                   Skv_zm )      

    do sclr = 1, sclr_dim
      
      call Skx_func( nzt, ngrdcol, sclrp2_zt(:,:,sclr), sclrp3(:,:,sclr), &
                     sclr_tol(sclr), clubb_params, &
                     Sksclr_zt(:,:,sclr) )   
                     
      call Skx_func( nzm, ngrdcol, sclrp2(:,:,sclr), sclrp3_zm(:,:,sclr), &
                     sclr_tol(sclr), clubb_params, &
                     Sksclr_zm(:,:,sclr) )  
                      
    end do ! sclr = 1, sclr_dim, 1

    if ( stats_metadata%l_stats_samp .and. l_samp_stats_in_pdf_call ) then

      !$acc update host( Skw_zt, Skw_zm, Skthl_zt, Skrt_zt, Skrt_zm, Skthl_zm )

      do i = 1, ngrdcol
        call stat_update_var( stats_metadata%iSkw_zt, Skw_zt(i,:), & ! In
                              stats_zt(i) ) ! In/Out
        call stat_update_var( stats_metadata%iSkw_zm, Skw_zm(i,:), &
                              stats_zm(i) ) ! In/Out
        call stat_update_var( stats_metadata%iSkthl_zt, Skthl_zt(i,:), &
                              stats_zt(i) ) ! In/Out
        call stat_update_var( stats_metadata%iSkthl_zm, Skthl_zm(i,:), &
                              stats_zm(i) ) ! In/Out
        call stat_update_var( stats_metadata%iSkrt_zt, Skrt_zt(i,:), &
                              stats_zt(i) ) ! In/Out
        call stat_update_var( stats_metadata%iSkrt_zm, Skrt_zm(i,:), &
                              stats_zm(i) ) ! In/Out
      end do
    end if

    ! The right hand side of this conjunction is only for reducing cpu time,
    ! since the more complicated formula is mathematically equivalent
    if ( l_gamma_Skw ) then

      !----------------------------------------------------------------
      ! Compute gamma as a function of Skw  - 14 April 06 dschanen
      !----------------------------------------------------------------
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm
        do i = 1, ngrdcol

          gamma_coef  = clubb_params(i,igamma_coef)
          gamma_coefb = clubb_params(i,igamma_coefb)
          gamma_coefc = clubb_params(i,igamma_coefc)

          if ( abs( gamma_coef - gamma_coefb ) > abs( gamma_coef + gamma_coefb ) * eps/2 ) then

            gamma_Skw_fnc(i,k) = gamma_coefb &
                                 + ( gamma_coef - gamma_coefb ) &
                                   * exp( -one_half * ( Skw_zm(i,k) / gamma_coefc )**2 )

          else

            gamma_Skw_fnc(i,k) = gamma_coef

          end if

        end do
      end do
      !$acc end parallel loop

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzt
        do i = 1, ngrdcol

          gamma_coef  = clubb_params(i,igamma_coef)
          gamma_coefb = clubb_params(i,igamma_coefb)
          gamma_coefc = clubb_params(i,igamma_coefc)

          if ( abs( gamma_coef - gamma_coefb ) > abs( gamma_coef + gamma_coefb ) * eps/2 ) then

            gamma_Skw_fnc_zt(i,k) = gamma_coefb &
                                    + ( gamma_coef - gamma_coefb ) &
                                      * exp( -one_half * ( Skw_zt(i,k) / gamma_coefc )**2 )

          else

            gamma_Skw_fnc_zt(i,k) = gamma_coef

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

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzt
        do i = 1, ngrdcol
          gamma_Skw_fnc_zt(i,k) = clubb_params(i,igamma_coef)
        end do
      end do
      !$acc end parallel loop

    end if

    if ( stats_metadata%l_stats_samp .and. l_samp_stats_in_pdf_call ) then
      !$acc update host(gamma_Skw_fnc)      
      do i = 1, ngrdcol
        call stat_update_var( stats_metadata%igamma_Skw_fnc, gamma_Skw_fnc(i,:), & ! intent(in)
                              stats_zm(i) )                       ! intent(inout)
      end do
    endif

    ! Compute sigma_sqd_w (dimensionless PDF width parameter)
    call compute_sigma_sqd_w( nzm, ngrdcol, &
                              gamma_Skw_fnc, wp2, thlp2, rtp2, &
                              up2, vp2, wpthlp, wprtp, upwp, vpwp, &
                              l_predict_upwp_vpwp, &
                              sigma_sqd_w_tmp )

    ! Smooth in the vertical using interpolation
    sigma_sqd_w(:,:) = zm2zt2zm( nzm, nzt, ngrdcol, gr, sigma_sqd_w_tmp(:,:), &
                                 zero_threshold ) ! Pos. def. quantity


    ! Interpolate the the stats_zt grid
    sigma_sqd_w_zt(:,:) = zm2zt( nzm, nzt, ngrdcol, gr, sigma_sqd_w(:,:), &
                                 zero_threshold )  ! Pos. def. quantity

    !---------------------------------------------------------------------------
    ! Interpolate thlp2, rtp2, and rtpthlp to thermodynamic levels,
    !---------------------------------------------------------------------------

    ! Interpolate variances to the stats_zt grid (statistics and closure)
    rtp2_zt(:,:)    = zm2zt( nzm, nzt, ngrdcol, gr, rtp2(:,:), &
                             rt_tol**2 )   ! Positive def. quantity
    thlp2_zt(:,:)   = zm2zt( nzm, nzt, ngrdcol, gr, thlp2(:,:), &
                             thl_tol**2 ) ! Positive def. quantity
    up2_zt(:,:)     = zm2zt( nzm, nzt, ngrdcol, gr, up2(:,:), &
                             w_tol_sqd )    ! Positive def. quantity
    vp2_zt(:,:)     = zm2zt( nzm, nzt, ngrdcol, gr, vp2(:,:), &
                             w_tol_sqd )    ! Positive def. quantity
    wprtp_zt(:,:)   = zm2zt( nzm, nzt, ngrdcol, gr, &
                             wprtp(:,:) )
    wpthlp_zt(:,:)  = zm2zt( nzm, nzt, ngrdcol, gr, &
                             wpthlp(:,:) )
    rtpthlp_zt(:,:) = zm2zt( nzm, nzt, ngrdcol, gr, &
                             rtpthlp(:,:) )
    upwp_zt(:,:)    = zm2zt( nzm, nzt, ngrdcol, gr, &
                             upwp(:,:) )
    vpwp_zt(:,:)    = zm2zt( nzm, nzt, ngrdcol, gr, &
                             vpwp(:,:) )

    ! Compute skewness velocity for stats output purposes
    if ( stats_metadata%iSkw_velocity > 0 ) then
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm
        do i = 1, ngrdcol
          Skw_velocity(i,k) = ( 1.0_core_rknd / ( 1.0_core_rknd - sigma_sqd_w(i,k) ) ) &
                       * ( wp3_zm(i,k) / max( wp2(i,k), w_tol_sqd ) )
        end do
      end do
      !$acc end parallel loop
    end if

    !----------------------------------------------------------------
    ! Call closure scheme
    !----------------------------------------------------------------

    ! Put passive scalar input on the t grid for the PDF
    do sclr = 1, sclr_dim
      wpsclrp_zt(:,:,sclr)   = zm2zt( nzm, nzt, ngrdcol, gr, &
                                      wpsclrp(:,:,sclr) )
      sclrp2_zt(:,:,sclr)    = zm2zt( nzm, nzt, ngrdcol, gr, sclrp2(:,:,sclr), &
                                      sclr_tol(sclr)**2 ) ! Pos. def. quantity
      sclrprtp_zt(:,:,sclr)  = zm2zt( nzm, nzt, ngrdcol, gr, &
                                      sclrprtp(:,:,sclr) )
      sclrpthlp_zt(:,:,sclr) = zm2zt( nzm, nzt, ngrdcol, gr, &
                                      sclrpthlp(:,:,sclr) )
    end do ! sclr = 1, sclr_dim, 1

    ! Interpolate hydrometeor mixed moments to momentum levels.
    do j = 1, hydromet_dim
      wphydrometp_zt(:,:,j) = zm2zt( nzm, nzt, ngrdcol, gr, &
                                     wphydrometp(:,:,j) )
    end do ! i = 1, hydromet_dim, 1

    call pdf_closure( nzt, ngrdcol, sclr_dim, sclr_tol,    & ! intent(in)
           hydromet_dim, p_in_Pa, exner, thv_ds_zt,        & ! intent(in)
           wm_zt, wp2_zt, wp3,                             & ! intent(in)
           Skw_zt, Skthl_zt, Skrt_zt, Sku_zt, Skv_zt,      & ! intent(in)
           rtm, rtp2_zt, wprtp_zt,                         & ! intent(in)
           thlm, thlp2_zt, wpthlp_zt,                      & ! intent(in)
           um, up2_zt, upwp_zt,                            & ! intent(in)
           vm, vp2_zt, vpwp_zt,                            & ! intent(in)
           rtpthlp_zt,                                     & ! intent(in)
           sclrm, wpsclrp_zt, sclrp2_zt,                   & ! intent(in)
           sclrprtp_zt, sclrpthlp_zt, Sksclr_zt,           & ! intent(in)
           gamma_Skw_fnc_zt,                               & ! intent(in)
#ifdef GFDL
           RH_crit,                                        & ! intent(inout)
           do_liquid_only_in_clubb,                        & ! intent(in)
#endif
           wphydrometp_zt, wp2hmp,                         & ! intent(in)
           rtphmp_zt, thlphmp_zt,                          & ! intent(in)
           clubb_params,                                   & ! intent(in)
           saturation_formula,                             & ! intent(in)
           stats_metadata,                                 & ! intent(in)
           iiPDF_type,                                     & ! intent(in)
           l_mix_rat_hm,                                   & ! intent(in)
           sigma_sqd_w_zt,                                 & ! intent(inout)
           pdf_params, pdf_implicit_coefs_terms,           & ! intent(inout)
           wpup2, wpvp2,                                   & ! intent(out)
           wp2up2_zt, wp2vp2_zt, wp4_zt,                   & ! intent(out)
           wprtp2, wp2rtp,                                 & ! intent(out)
           wpthlp2, wp2thlp, wprtpthlp,                    & ! intent(out)
           cloud_frac, ice_supersat_frac,                  & ! intent(out)
           rcm, wpthvp_zt, wp2thvp, rtpthvp_zt,            & ! intent(out)
           thlpthvp_zt, wprcp_zt, wp2rcp, rtprcp_zt,       & ! intent(out)
           thlprcp_zt, rcp2_zt,                            & ! intent(out)
           uprcp_zt, vprcp_zt,                             & ! intent(out)
           w_up_in_cloud, w_down_in_cloud,                 & ! intent(out)
           cloudy_updraft_frac, cloudy_downdraft_frac,     & ! intent(out)
           F_w, F_rt, F_thl,                               & ! intent(out)
           min_F_w, max_F_w,                               & ! intent(out)
           min_F_rt, max_F_rt,                             & ! intent(out)
           min_F_thl, max_F_thl,                           & ! intent(out)
           wpsclrprtp, wpsclrp2, sclrpthvp_zt,             & ! intent(out)
           wpsclrpthlp, sclrprcp_zt, wp2sclrp,             & ! intent(out)
           rc_coef                                         ) ! intent(out)

    ! Subroutine may produce NaN values, and if so, return
    if ( clubb_at_least_debug_level( 0 ) ) then
       if ( err_code == clubb_fatal_error ) then
          write(fstderr,*) "After pdf_closure"
          return
       endif
    endif

    ! Stats output
    if ( stats_metadata%l_stats_samp .and. l_samp_stats_in_pdf_call ) then

      do i = 1, ngrdcol
        call stat_update_var( stats_metadata%iF_w, F_w(i,:), & ! intent(in)
                              stats_zt(i) )   ! intent(inout)
        call stat_update_var( stats_metadata%iF_rt, F_rt(i,:), & ! intent(in)
                              stats_zt(i) )     ! intent(inout)
        call stat_update_var( stats_metadata%iF_thl, F_thl(i,:), & ! intent(in)
                              stats_zt(i) )       ! intent(inout)
        call stat_update_var( stats_metadata%imin_F_w, min_F_w(i,:), & ! intent(in)
                              stats_zt(i) )           ! intent(inout)
        call stat_update_var( stats_metadata%imax_F_w, max_F_w(i,:), & ! intent(in)
                              stats_zt(i) )           ! intent(inout)
        call stat_update_var( stats_metadata%imin_F_rt, min_F_rt(i,:), & ! intent(in)
                              stats_zt(i) )             ! intent(inout)
        call stat_update_var( stats_metadata%imax_F_rt, max_F_rt(i,:), & ! intent(in)
                              stats_zt(i) )             ! intent(inout)
        call stat_update_var( stats_metadata%imin_F_thl, min_F_thl(i,:), & ! intent(in)
                              stats_zt(i) )               ! intent(inout)
        call stat_update_var( stats_metadata%imax_F_thl, max_F_thl(i,:), & ! intent(in)
                              stats_zt(i) )               ! intent(inout)
      end do
    end if

    if( l_rtm_nudge ) then
      ! Nudge rtm to prevent excessive drying
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzt
        do i = 1, ngrdcol
          if ( rtm(i,k) < rtm_min .and. gr%zt(i,k) < rtm_nudge_max_altitude ) then
            rtm(i,k) = rtm(i,k) + (rtm_ref(i,k) - rtm(i,k)) * ( dt / ts_nudge )
          end if
        end do
      end do
      !$acc end parallel loop
    end if

    if ( l_call_pdf_closure_twice ) then

      ! Call pdf_closure a second time on momentum levels, to
      ! output (rather than interpolate) the variables which
      ! belong on the momentum levels.

      ! Interpolate sclrm to the momentum level for use in
      ! the second call to pdf_closure
      do sclr = 1, sclr_dim
        ! Clip if extrap. causes sclrm_zm to be less than sclr_tol
        sclrm_zm(:,:,sclr) = zt2zm( nzm, nzt, ngrdcol, gr, sclrm(:,:,sclr), sclr_tol(sclr) )
      end do ! sclr = 1, sclr_dim

      ! Interpolate pressure, p_in_Pa, to momentum levels.
      ! Since the surface (or model lower boundary) is located at momentum level
      ! k = 1, the pressure there is p_sfc.
      p_in_Pa_zm(:,:) = zt2zm( nzm, nzt, ngrdcol, gr, p_in_Pa(:,:) )

      !$acc parallel loop gang vector default(present)
      do i = 1, ngrdcol
        p_in_Pa_zm(i,1) = p_sfc(i)

        ! Clip pressure if the extrapolation leads to a negative value of pressure
        p_in_Pa_zm(i,nzm) = max( p_in_Pa_zm(i,nzm), 0.5_core_rknd*p_in_Pa(i,nzt) )
      end do
      !$acc end parallel loop

      ! Set exner at momentum levels, exner_zm, based on p_in_Pa_zm.
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm
        do i = 1, ngrdcol
          exner_zm(i,k) = (p_in_Pa_zm(i,k)/p0)**kappa
        end do
      end do
      !$acc end parallel loop

      ! Clip if extrapolation at the top level causes rtm_zm to be < rt_tol
      rtm_zm(:,:) = zt2zm( nzm, nzt, ngrdcol, gr, rtm(:,:), rt_tol )

        ! Clip if extrapolation at the top level causes thlm_zm to be < thl_tol
      thlm_zm(:,:) = zt2zm( nzm, nzt, ngrdcol, gr, thlm(:,:), thl_tol )

      ! Interpolate hydrometeor mixed moments to momentum levels.
      do j = 1, hydromet_dim
        rtphmp(:,:,j)    = zt2zm( nzm, nzt, ngrdcol, gr, rtphmp_zt(:,:,j) )
        thlphmp(:,:,j)   = zt2zm( nzm, nzt, ngrdcol, gr, thlphmp_zt(:,:,j) )
        wp2hmp_zm(:,:,j) = zt2zm( nzm, nzt, ngrdcol, gr, wp2hmp(:,:,j) )
      end do ! i = 1, hydromet_dim, 1
      
      um_zm(:,:) = zt2zm( nzm, nzt, ngrdcol, gr, um(:,:) )
      vm_zm(:,:) = zt2zm( nzm, nzt, ngrdcol, gr, vm(:,:) )
      
      ! pdf_implicit_coefs_terms is only used in the iiPDF_new and iiPDF_new_hybrid closures.
      ! So we only need to initialize our local _zm version if we're working with one of those.
      if ( iiPDF_type == iiPDF_new .or. iiPDF_type == iiPDF_new_hybrid ) then
        call init_pdf_implicit_coefs_terms( nzm, ngrdcol, sclr_dim, &     ! Intent(in)
                                            pdf_implicit_coefs_terms_zm ) ! Intent(out)
      end if 

      ! Call pdf_closure to output the variables which belong on the momentum grid.
      call pdf_closure( nzm, ngrdcol, sclr_dim, sclr_tol,          & ! intent(in)
             hydromet_dim, p_in_Pa_zm, exner_zm, thv_ds_zm,        & ! intent(in)
             wm_zm, wp2, wp3_zm,                                   & ! intent(in)
             Skw_zm, Skthl_zm, Skrt_zm, Sku_zm, Skv_zm,            & ! intent(in)
             rtm_zm, rtp2, wprtp,                                  & ! intent(in)
             thlm_zm, thlp2, wpthlp,                               & ! intent(in)
             um_zm, up2, upwp,                                     & ! intent(in)
             vm_zm, vp2, vpwp,                                     & ! intent(in)
             rtpthlp,                                              & ! intent(in)
             sclrm_zm, wpsclrp, sclrp2,                            & ! intent(in)
             sclrprtp, sclrpthlp, Sksclr_zm,                       & ! intent(in)
             gamma_Skw_fnc,                                        & ! intent(in)
#ifdef GFDL
             RH_crit,                                              & ! intent(inout)
             do_liquid_only_in_clubb,                              & ! intent(in)
#endif
             wphydrometp, wp2hmp_zm,                               & ! intent(in)
             rtphmp, thlphmp,                                      & ! intent(in)
             clubb_params,                                         & ! intent(in)
             saturation_formula,                                   & ! intent(in)
             stats_metadata,                                       & ! intent(in)
             iiPDF_type,                                           & ! intent(in)
             l_mix_rat_hm,                                         & ! intent(in)
             sigma_sqd_w,                                          & ! intent(inout)
             pdf_params_zm, pdf_implicit_coefs_terms_zm,           & ! intent(inout)
             wpup2_zm, wpvp2_zm,                                   & ! intent(out)
             wp2up2, wp2vp2, wp4,                                  & ! intent(out)
             wprtp2_zm, wp2rtp_zm,                                 & ! intent(out)
             wpthlp2_zm, wp2thlp_zm, wprtpthlp_zm,                 & ! intent(out)
             cloud_frac_zm, ice_supersat_frac_zm,                  & ! intent(out)
             rcm_zm, wpthvp, wp2thvp_zm, rtpthvp,                  & ! intent(out)
             thlpthvp, wprcp, wp2rcp_zm, rtprcp,                   & ! intent(out)
             thlprcp, rcp2,                                        & ! intent(out)
             uprcp, vprcp,                                         & ! intent(out)
             w_up_in_cloud_zm, w_down_in_cloud_zm,                 & ! intent(out)
             cloudy_updraft_frac_zm, cloudy_downdraft_frac_zm,     & ! intent(out)
             F_w_zm, F_rt_zm, F_thl_zm,                            & ! intent(out)
             min_F_w_zm, max_F_w_zm,                               & ! intent(out)
             min_F_rt_zm, max_F_rt_zm,                             & ! intent(out)
             min_F_thl_zm, max_F_thl_zm,                           & ! intent(out)
             wpsclrprtp_zm, wpsclrp2_zm, sclrpthvp,                & ! intent(out)
             wpsclrpthlp_zm, sclrprcp, wp2sclrp_zm,                & ! intent(out)
             rc_coef_zm                                            ) ! intent(out)

      ! Subroutine may produce NaN values, and if so, return
      if ( clubb_at_least_debug_level( 0 ) ) then
         if ( err_code == clubb_fatal_error ) then
            write(fstderr,*) "After pdf_closure"
            return
         endif
      endif

    else ! l_call_pdf_closure_twice is false
      
      ! Interpolate momentum variables output from the first call to
      ! pdf_closure back to momentum grid.
      wp4(:,:) = zt2zm( nzm, nzt, ngrdcol, gr, wp4_zt(:,:), zero_threshold )  ! Pos. def. quantity

      !$acc parallel loop gang vector default(present)
      do i = 1, ngrdcol
        ! Since top momentum level is higher than top thermo level,
        ! set variables at top momentum level to 0.
        wp4(i,nzm) = zero
        ! Set wp4 to 0 at the lowest momentum level (momentum level 1).
        ! The value of wp4 at momentum level 1 is found by interpolation of
        ! the values produced by the PDF for wp4_zt at thermodynamic levels
        ! 1 and 2.  This value is unreliable at thermodynamic level 1.
        wp4(i,1) = zero
      end do
      !$acc end parallel loop

#ifndef CLUBB_CAM
      ! CAM-CLUBB needs cloud water variance thus always compute this
      if ( stats_metadata%ircp2 > 0 ) then
#endif
        rcp2(:,:) = zt2zm( nzm, nzt, ngrdcol, gr, rcp2_zt(:,:), &
                           zero_threshold )  ! Pos. def. quantity
#ifndef CLUBB_CAM
        !$acc parallel loop gang vector default(present) 
        do i = 1, ngrdcol
          rcp2(i,nzm) = zero
        end do
        !$acc end parallel loop
      endif
#endif

      wpthvp(:,:)      = zt2zm( nzm, nzt, ngrdcol, gr, wpthvp_zt(:,:) )
      thlpthvp(:,:)    = zt2zm( nzm, nzt, ngrdcol, gr, thlpthvp_zt(:,:) )
      rtpthvp(:,:)     = zt2zm( nzm, nzt, ngrdcol, gr, rtpthvp_zt(:,:) )
      wprcp(:,:)       = zt2zm( nzm, nzt, ngrdcol, gr, wprcp_zt(:,:) )
      rc_coef_zm(:,:)  = zt2zm( nzm, nzt, ngrdcol, gr, rc_coef(:,:) )
      rtprcp(:,:)      = zt2zm( nzm, nzt, ngrdcol, gr, rtprcp_zt(:,:) )
      thlprcp(:,:)     = zt2zm( nzm, nzt, ngrdcol, gr, thlprcp_zt(:,:) )
      uprcp(:,:)       = zt2zm( nzm, nzt, ngrdcol, gr, uprcp_zt(:,:) )
      vprcp(:,:)       = zt2zm( nzm, nzt, ngrdcol, gr, vprcp_zt(:,:) )
      wp2up2(:,:)      = zt2zm( nzm, nzt, ngrdcol, gr, wp2up2_zt(:,:) )
      wp2vp2(:,:)      = zt2zm( nzm, nzt, ngrdcol, gr, wp2vp2_zt(:,:) )

      !$acc parallel loop gang vector default(present) 
      do i = 1, ngrdcol 
        wpthvp(i,nzm)     = 0.0_core_rknd
        thlpthvp(i,nzm)   = 0.0_core_rknd
        rtpthvp(i,nzm)    = 0.0_core_rknd
        wprcp(i,nzm)      = 0.0_core_rknd
        rc_coef_zm(i,nzm) = 0.0_core_rknd
        rtprcp(i,nzm)     = 0.0_core_rknd
        thlprcp(i,nzm)    = 0.0_core_rknd
        uprcp(i,nzm)      = 0.0_core_rknd
        vprcp(i,nzm)      = 0.0_core_rknd
        wp2up2(i,nzm)     = 0.0_core_rknd
        wp2vp2(i,nzm)     = 0.0_core_rknd
      end do
      !$acc end parallel loop

      ! Initialize variables to avoid uninitialized variables.
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm
        do i = 1, ngrdcol
          cloud_frac_zm(i,k)        = 0.0_core_rknd
          ice_supersat_frac_zm(i,k) = 0.0_core_rknd
          rcm_zm(i,k)               = 0.0_core_rknd
          rtm_zm(i,k)               = 0.0_core_rknd
          thlm_zm(i,k)              = 0.0_core_rknd
        end do
      end do
      !$acc end parallel loop

      ! Interpolate passive scalars back onto the m grid
      do sclr = 1, sclr_dim
        sclrpthvp(:,:,sclr)       = zt2zm( nzm, nzt, ngrdcol, gr, sclrpthvp_zt(:,:,sclr) )
        sclrprcp(:,:,sclr)        = zt2zm( nzm, nzt, ngrdcol, gr, sclrprcp_zt(:,:,sclr) )

        !$acc parallel loop gang vector default(present)
        do i = 1, ngrdcol
          sclrpthvp(i,nzm,sclr) = 0.0_core_rknd
          sclrprcp(i,nzm,sclr)  = 0.0_core_rknd
        end do
        !$acc end parallel loop

      end do ! sclr=1, sclr_dim

    end if ! l_call_pdf_closure_twice
    
    if ( stats_metadata%l_stats_samp .and. l_samp_stats_in_pdf_call ) then
      !$acc update host( uprcp, vprcp )
      do i = 1, ngrdcol
        call stat_update_var( stats_metadata%iuprcp,  uprcp(i,:),  & ! intent(in)
                              stats_zm(i) )           ! intent(inout)
        call stat_update_var( stats_metadata%ivprcp,  vprcp(i,:),  & ! intent(in)
                              stats_zm(i) )           ! intent(inout)
      end do
    end if
    
    ! If l_trapezoidal_rule_zt is true, call trapezoidal_rule_zt for
    ! thermodynamic-level variables output from pdf_closure.
    ! ldgrant June 2009
    if ( l_trapezoidal_rule_zt ) then
      call trapezoidal_rule_zt( nzm, nzt, ngrdcol, sclr_dim, gr,             & ! intent(in)
                                l_call_pdf_closure_twice,                    & ! intent(in)
                                stats_metadata,                              & ! intent(in)
                                wprtp2, wpthlp2,                             & ! intent(inout)
                                wprtpthlp, cloud_frac, ice_supersat_frac,    & ! intent(inout)
                                rcm, wp2thvp, wpsclrprtp, wpsclrp2,          & ! intent(inout)
                                wpsclrpthlp,                                 & ! intent(inout)
                                wprtp2_zm, wpthlp2_zm,                       & ! intent(inout)
                                wprtpthlp_zm, cloud_frac_zm,                 & ! intent(inout)
                                ice_supersat_frac_zm, rcm_zm, wp2thvp_zm,    & ! intent(inout)
                                wpsclrprtp_zm, wpsclrp2_zm, wpsclrpthlp_zm )   ! intent(inout)
    else ! l_trapezoidal_rule_zt
      cloud_frac_zm = zt2zm( nzm, nzt, ngrdcol, gr, cloud_frac )
      ! Since top momentum level is higher than top thermo. level,
      ! set variables at top momentum level to 0.
      !$acc parallel loop gang vector default(present)
      do i = 1, ngrdcol
        cloud_frac_zm(i,nzm)  = 0.0_core_rknd
      end do
      !$acc end parallel loop
    end if ! l_trapezoidal_rule_zt

    ! If l_trapezoidal_rule_zm is true, call trapezoidal_rule_zm for
    ! the important momentum-level variabes output from pdf_closure.
    ! ldgrant Feb. 2010
    if ( l_trapezoidal_rule_zm ) then
      call trapezoidal_rule_zm( nzm, nzt, ngrdcol, gr,              & ! intent(in)
                                wpthvp_zt, thlpthvp_zt, rtpthvp_zt, & ! intent(in)
                                wpthvp, thlpthvp, rtpthvp )           ! intent(inout)
    end if ! l_trapezoidal_rule_zm


    ! Vince Larson clipped rcm in order to prevent rvm < 0.  5 Apr 2008.
    ! This code won't work unless rtm >= 0 !!!
    ! We do not clip rcm_in_layer because rcm_in_layer only influences
    ! radiation, and we do not want to bother recomputing it.
    ! Code is duplicated from below to ensure that relative humidity
    ! is calculated properly.  3 Sep 2009
    call clip_rcm( nzt, ngrdcol, rtm,             & ! intent(in)
                   'rtm < rcm after pdf_closure', & ! intent(in)
                   rcm )                            ! intent(inout)

    ! Compute variables cloud_cover and rcm_in_layer.
    ! Added July 2009
    call compute_cloud_cover( gr, nzt, ngrdcol,             & ! intent(in)
                              pdf_params, cloud_frac, rcm, & ! intent(in)
                              cloud_cover, rcm_in_layer )    ! intent(out)

    if ( l_use_cloud_cover ) then
      ! Use cloud_cover and rcm_in_layer to help boost cloud_frac and rcm to help
      ! increase cloudiness at coarser grid resolutions.
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzt
        do i = 1, ngrdcol
          cloud_frac(i,k) = cloud_cover(i,k)
          rcm(i,k) = rcm_in_layer(i,k)
        end do
      end do
    !$acc end parallel loop
    end if

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzt
      do i = 1, ngrdcol
        ! Clip cloud fraction here if it still exceeds 1.0 due to round off
        cloud_frac(i,k) = min( 1.0_core_rknd, cloud_frac(i,k) )
        ! Ditto with ice cloud fraction
        ice_supersat_frac(i,k) = min( 1.0_core_rknd, ice_supersat_frac(i,k) )
      end do
    end do
    !$acc end parallel loop

    T_in_K = thlm2T_in_K( nzt, ngrdcol, thlm, exner, rcm )
    rsat = sat_mixrat_liq( nzt, ngrdcol, p_in_Pa, T_in_K, saturation_formula )

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzt
      do i = 1, ngrdcol
        rel_humidity(i,k) = (rtm(i,k) - rcm(i,k)) / rsat(i,k)
        rcm_supersat_adj(i,k) = zero
      end do
    end do
    !$acc end parallel loop
      
    if ( l_rcm_supersat_adj ) then
      ! +PAB mods, take remaining supersaturation that may exist
      !   after CLUBB PDF call and add it to rcm.  Supersaturation
      !   may exist after PDF call due to issues with calling PDF on the
      !   thermo grid and momentum grid and the interpolation between the two
      l_spur_supersat = .false.
      

      !$acc parallel loop gang vector collapse(2) default(present) reduction(.or.:l_spur_supersat)
      do k = 1, nzt
        do i = 1, ngrdcol
          if (rel_humidity(i,k) > 1.0_core_rknd) then
            rcm_supersat_adj(i,k) = (rtm(i,k) - rcm(i,k)) - rsat(i,k)
            rcm(i,k) = rcm(i,k) + rcm_supersat_adj(i,k)
            l_spur_supersat = .true.
          end if
        end do
      end do
      !$acc end parallel loop

      if ( clubb_at_least_debug_level( 1 ) .and. l_spur_supersat ) then
        write(fstderr,*) 'Warning: spurious supersaturation was removed after pdf_closure!'
      end if

    end if ! l_rcm_supersat_adj

    !$acc exit data delete( wp2_zt,wp3_zm, rtp2_zt,rtp3_zm, thlp2_zt,  thlp3_zm, &
    !$acc                   wprtp_zt, wpthlp_zt, rtpthlp_zt, up2_zt, up3_zm, &
    !$acc                   vp2_zt, vp3_zm, upwp_zt, vpwp_zt, gamma_Skw_fnc, &
    !$acc                   gamma_Skw_fnc_zt,sigma_sqd_w_zt,  Skw_zt, Skw_zm, &
    !$acc                   Skrt_zt, Skrt_zm, Skthl_zt, Skthl_zm, Sku_zt, &
    !$acc                   Sku_zm, Skv_zt, Skv_zm, wp2up2_zt, &
    !$acc                   wp2vp2_zt, wp4_zt, wpthvp_zt, rtpthvp_zt, thlpthvp_zt, &
    !$acc                   wprcp_zt, rtprcp_zt, thlprcp_zt, uprcp_zt, vprcp_zt, &
    !$acc                   rsat, rel_humidity, um_zm, vm_zm, T_in_K, sigma_sqd_w_tmp )

    !$acc exit data if( l_call_pdf_closure_twice ) &
    !$acc           delete( w_up_in_cloud_zm, wpup2_zm, wpvp2_zm, &
    !$acc                   w_down_in_cloud_zm, cloudy_updraft_frac_zm,  &
    !$acc                   cloudy_downdraft_frac_zm, p_in_Pa_zm, exner_zm, &
    !$acc                   wprtp2_zm, wp2rtp_zm, wpthlp2_zm, &
    !$acc                   wp2thlp_zm, wprtpthlp_zm, wp2thvp_zm, wp2rcp_zm )

    !$acc exit data if( sclr_dim > 0 ) &
    !$acc           delete( wpsclrp_zt, sclrp2_zt, sclrp3_zm, sclrprtp_zt, sclrpthlp_zt, &
    !$acc                   Sksclr_zt, Sksclr_zm, sclrpthvp_zt, sclrprcp_zt, wpsclrprtp_zm, &
    !$acc                   wpsclrp2_zm, wpsclrpthlp_zm, wp2sclrp_zm, sclrm_zm )

    !$acc exit data if( hydromet_dim > 0 ) delete( wphydrometp_zt, wp2hmp_zm, rtphmp, thlphmp )

    return

  end subroutine pdf_closure_driver

  !=============================================================================
  subroutine check_clubb_settings( nzmax, params,       & ! intent(in)
                                   l_implemented,       & ! intent(in)
                                   l_input_fields,      & ! intent(in)
                                   clubb_config_flags,  & ! intent(in)
                                   err_code_out )         ! intent(out)

      ! Description:
      !   Subroutine to set up the model for execution.
      !
      ! References:
      !   None
      !---------------------------------------------------------------------

      use parameter_indices, only:  &
          nparams,      & ! Variable(s)
          iC1,          & ! Constant(s)
          iC1b,         &
          iC2rt,        &
          iC2thl,       &
          iC2rtthl,     &
          iC6rt,        &
          iC6rtb,       &
          iC6thl,       &
          iC6thlb,      &
          iC14

      use parameters_tunable, only: &
          setup_parameters,    & ! Procedure
          nu_vertical_res_dep    ! Type(s)

      use constants_clubb, only:  &
          fstderr, &  ! Variable(s)
          one, &
          eps

      use error_code, only: &
          clubb_at_least_debug_level,  & ! Procedures
          initialize_error_headers,    &
          err_code,                    & ! Error Indicator
          clubb_no_error, &              ! Constant
          clubb_fatal_error              ! Constant

      use model_flags, only: &
          clubb_config_flags_type, & ! Type
          iiPDF_ADG1,       & ! Variable(s)
          iiPDF_ADG2,       &
          iiPDF_3D_Luhar,   &
          iiPDF_new,        &
          iiPDF_TSDADG,     &
          iiPDF_LY93,       &
          iiPDF_new_hybrid, &
          l_explicit_turbulent_adv_wpxp
          
#ifdef CLUBB_GPU
      use model_flags, only: &
          lapack
#endif

      use clubb_precision, only: &
          core_rknd ! Variable(s)

      use sponge_layer_damping, only: &
          thlm_sponge_damp_settings,    & ! Variable(s)
          rtm_sponge_damp_settings,     &
          uv_sponge_damp_settings,      &
          wp2_sponge_damp_settings,     &
          wp3_sponge_damp_settings,     &
          up2_vp2_sponge_damp_settings


      implicit none

      !---------------------- Input Variables ----------------------

      ! Grid definition
      integer, intent(in) :: nzmax  ! Vertical grid levels            [#]
      !                      Only true when used in a host model
      !                      CLUBB determines what nzmax should be
      !                      given zm_init and zm_top when
      !                      running in standalone mode.

      real( kind = core_rknd ), intent(in), dimension(nparams) :: &
        params  ! Including C1, nu1, nu2, etc.

      ! Flags
      logical, intent(in) ::  &
        l_implemented, & ! Flag for whether CLUBB is being run in a host model
        l_input_fields   ! Flag for whether LES input fields are being used

      type(clubb_config_flags_type), intent(in) :: &
        clubb_config_flags

      !---------------------- Output Variables ----------------------
      integer, intent(out) :: &
        err_code_out  ! Error code indicator

      !---------------------- Begin Code ----------------------

      err_code_out = clubb_no_error ! Initialize to no error value
      call initialize_error_headers

#ifdef CLUBB_GPU
      if ( clubb_config_flags%penta_solve_method == lapack ) then
        write(fstderr,*) "WARNING: The penta-diagonal lapack solver is not GPU accelerated"
        write(fstderr,*) " Set penta_solve_method = 2, to use an accelerated penta-diagonal solver"
      end if

      if ( clubb_config_flags%tridiag_solve_method == lapack ) then
        write(fstderr,*) "WARNING: The tri-diagonal lapack solver is not GPU accelerated"
        write(fstderr,*) " Set tridiag_solve_method = 2, to use an accelerated tri-diagonal solver"
      end if
#endif

      ! Sanity check
      if ( clubb_at_least_debug_level( 0 ) ) then

        if ( clubb_config_flags%l_damp_wp2_using_em .and. &
           (abs(params(iC1) - params(iC14)) > abs(params(iC1) + params(iC14)) / 2 * eps .or. &
             clubb_config_flags%l_stability_correct_tau_zm) ) then
          write(fstderr,*) "l_damp_wp2_using_em = T requires C1=C14 and" &
                            // " l_stability_correct_tau_zm = F"
          write(fstderr,*) "C1 = ", params(iC1)
          write(fstderr,*) "C14 = ", params(iC14)
          write(fstderr,*) "l_stability_correct_tau_zm = ", &
                           clubb_config_flags%l_stability_correct_tau_zm
          write(fstderr,*) "Fatal error in check_clubb_settings"
          err_code_out = clubb_fatal_error
          return
        end if

      end if

      ! Sanity check for the saturation formula
      select case ( clubb_config_flags%saturation_formula )
      case ( saturation_bolton )
        ! Using the Bolton 1980 approximations for SVP over vapor/ice

      case ( saturation_gfdl )
        ! Using the Flatau, et al. polynomial approximation for SVP over vapor/ice

      case ( saturation_flatau )   ! h1g, 2010-06-16
        ! Using the GFDL SVP formula (Goff-Gratch)

        ! Add new saturation formulas after this

      case ( saturation_lookup )
        ! Using the lookup table

      case default
        write(fstderr,*) "Unknown approx. of saturation vapor pressure: ", &
           clubb_config_flags%saturation_formula
        write(fstderr,*) "Fatal error in check_clubb_settings"
        err_code_out = clubb_fatal_error
        return
      end select

      ! Check for the type of two component normal (double Gaussian) PDF being
      ! used for w, rt, and theta-l (or w, chi, and eta).
      if ( clubb_config_flags%iiPDF_type < iiPDF_ADG1 &
           .or. clubb_config_flags%iiPDF_type > iiPDF_new_hybrid ) then
         write(fstderr,*) "Unknown type of double Gaussian PDF selected: ", &
                          clubb_config_flags%iiPDF_type
         write(fstderr,*) "iiPDF_type = ", clubb_config_flags%iiPDF_type
         write(fstderr,*) "Fatal error in check_clubb_settings"
         err_code_out = clubb_fatal_error
         return
      endif ! iiPDF_type < iiPDF_ADG1 or iiPDF_type > iiPDF_lY93

      ! The ADG2 and 3D Luhar PDFs can only be used as part of input fields.
      if ( clubb_config_flags%iiPDF_type == iiPDF_ADG2 ) then
         if ( .not. l_input_fields ) then
            write(fstderr,*) "The ADG2 PDF can only be used with" &
                             // " input fields (l_input_fields = .true.)."
            write(fstderr,*) "iiPDF_type = ", clubb_config_flags%iiPDF_type
            write(fstderr,*) "l_input_fields = ", l_input_fields
            write(fstderr,*) "Fatal error in check_clubb_settings"
            err_code_out = clubb_fatal_error
            return
         endif ! .not. l_input_fields
      endif ! iiPDF_type == iiPDF_ADG2

      if ( clubb_config_flags%iiPDF_type == iiPDF_3D_Luhar ) then
         if ( .not. l_input_fields ) then
            write(fstderr,*) "The 3D Luhar PDF can only be used with" &
                             // " input fields (l_input_fields = .true.)."
            write(fstderr,*) "iiPDF_type = ", clubb_config_flags%iiPDF_type
            write(fstderr,*) "l_input_fields = ", l_input_fields
            write(fstderr,*) "Fatal error in check_clubb_settings"
            err_code_out = clubb_fatal_error
            return
         endif ! .not. l_input_fields
      endif ! iiPDF_type == iiPDF_3D_Luhar

      ! This also currently applies to the new PDF until it has been fully
      ! implemented.
      if ( clubb_config_flags%iiPDF_type == iiPDF_new ) then
         if ( .not. l_input_fields ) then
            write(fstderr,*) "The new PDF can only be used with" &
                             // " input fields (l_input_fields = .true.)."
            write(fstderr,*) "iiPDF_type = ", clubb_config_flags%iiPDF_type
            write(fstderr,*) "l_input_fields = ", l_input_fields
            write(fstderr,*) "Fatal error in check_clubb_settings"
            err_code_out = clubb_fatal_error
            return
         endif ! .not. l_input_fields
      endif ! iiPDF_type == iiPDF_new

      ! This also currently applies to the TSDADG PDF until it has been fully
      ! implemented.
      if ( clubb_config_flags%iiPDF_type == iiPDF_TSDADG ) then
         if ( .not. l_input_fields ) then
            write(fstderr,*) "The new TSDADG PDF can only be used with" &
                             // " input fields (l_input_fields = .true.)."
            write(fstderr,*) "iiPDF_type = ", clubb_config_flags%iiPDF_type
            write(fstderr,*) "l_input_fields = ", l_input_fields
            write(fstderr,*) "Fatal error in check_clubb_settings"
            err_code_out = clubb_fatal_error
            return
         endif ! .not. l_input_fields
      endif ! iiPDF_type == iiPDF_TSDADG

      ! This also applies to Lewellen and Yoh (1993).
      if ( clubb_config_flags%iiPDF_type == iiPDF_LY93 ) then
         if ( .not. l_input_fields ) then
            write(fstderr,*) "The Lewellen and Yoh PDF can only be used with" &
                             // " input fields (l_input_fields = .true.)."
            write(fstderr,*) "iiPDF_type = ", clubb_config_flags%iiPDF_type
            write(fstderr,*) "l_input_fields = ", l_input_fields
            write(fstderr,*) "Fatal error in check_clubb_settings"
            err_code_out = clubb_fatal_error
            return
         endif ! .not. l_input_fields
      endif ! iiPDF_type == iiPDF_LY93

      ! Check the option for the placement of the call to CLUBB's PDF.
      if ( clubb_config_flags%ipdf_call_placement < ipdf_pre_advance_fields &
           .or. clubb_config_flags%ipdf_call_placement > ipdf_pre_post_advance_fields ) then
         write(fstderr,*) "Invalid option selected for ipdf_call_placement: ", &
                          clubb_config_flags%ipdf_call_placement
         write(fstderr,*) "Fatal error in check_clubb_settings"
         err_code_out = clubb_fatal_error
         return
      endif

      ! The l_predict_upwp_vpwp flag requires that the ADG1 PDF is used
      ! implicitly in subroutine advance_xm_wpxp.
      if ( clubb_config_flags%l_predict_upwp_vpwp ) then

         ! When l_predict_upwp_vpwp is enabled, the
         ! l_explicit_turbulent_adv_wpxp flag must be turned off.
         ! Otherwise, explicit turbulent advection would require PDF parameters
         ! for u and v to be calculated in PDF closure.  These would be needed
         ! to calculate integrated fields such as wp2up, etc.
         if ( l_explicit_turbulent_adv_wpxp ) then
            write(fstderr,*) "The l_explicit_turbulent_adv_wpxp option" &
                             // " is not currently set up for use with the" &
                             // " l_predict_upwp_vpwp code."
            write(fstderr,*) "Fatal error in check_clubb_settings"
            err_code_out = clubb_fatal_error
            return
         endif ! l_explicit_turbulent_adv_wpxp

         ! When l_predict_upwp_vpwp is enabled, the PDF type must be set to
         ! the ADG1 PDF or the new hybrid PDF.  The other PDFs are not currently
         ! set up to calculate variables needed for implicit or semi-implicit
         ! turbulent advection, such as coef_wp2up_implicit, etc.
         if ( ( clubb_config_flags%iiPDF_type /= iiPDF_ADG1 ) &
              .and. ( clubb_config_flags%iiPDF_type /= iiPDF_new_hybrid ) ) then
            write(fstderr,*) "Currently, only the ADG1 PDF and the new hybrid" &
                             // " PDF are set up for use with the" &
                             // " l_predict_upwp_vpwp code."
            write(fstderr,*) "Fatal error in check_clubb_settings"
            err_code_out = clubb_fatal_error
            return
         endif ! iiPDF_type /= iiPDF_ADG1

      endif ! l_predict_upwp_vpwp

      ! The flags l_min_xp2_from_corr_wx and l_enable_relaxed_clipping must
      ! have opposite values.
      if ( ( clubb_config_flags%l_min_xp2_from_corr_wx ) &
         .and. ( clubb_config_flags%l_enable_relaxed_clipping ) ) then
         write(fstderr,*) "Invalid configuration: l_min_xp2_from_corr_wx = T " &
                          // "and l_enable_relaxed_clipping = T"
         write(fstderr,*) "They must have opposite values"
         write(fstderr,*) "Fatal error in check_clubb_settings"
         err_code_out = clubb_fatal_error
         return
      elseif ( ( .not. clubb_config_flags%l_min_xp2_from_corr_wx ) &
               .and. ( .not. clubb_config_flags%l_enable_relaxed_clipping ) ) then
         write(fstderr,*) "Invalid configuration: l_min_xp2_from_corr_wx = F " &
                          // "and l_enable_relaxed_clipping = F"
         write(fstderr,*) "They must have opposite values"
         write(fstderr,*) "Fatal error in check_clubb_settings"
         !err_code_out = clubb_fatal_error
         !return
      endif

      ! Checking for the code that orders CLUBB's advance_ subroutines
      if ( order_xm_wpxp < 1 .or. order_xm_wpxp > 4 ) then
         write(fstderr,*) "The variable order_xm_wpxp must have a value " &
                          // "between 1 and 4"
         write(fstderr,*) "order_xm_wpxp = ", order_xm_wpxp
         write(fstderr,*) "Fatal error in check_clubb_settings"
         err_code_out = clubb_fatal_error
         return
      elseif ( order_xm_wpxp == order_wp2_wp3 &
               .or. order_xm_wpxp == order_xp2_xpyp &
               .or. order_xm_wpxp == order_windm ) then
         write(fstderr,*) "The variable order_xm_wpxp has the same value " &
                          // "as another order_ variable.  Please give each " &
                          // "order index a unique value."
         write(fstderr,*) "order_xm_wpxp = ", order_xm_wpxp
         write(fstderr,*) "order_wp2_wp3 = ", order_wp2_wp3
         write(fstderr,*) "order_xp2_xpyp = ", order_xp2_xpyp
         write(fstderr,*) "order_windm = ", order_windm
         write(fstderr,*) "Fatal error in check_clubb_settings"
         err_code_out = clubb_fatal_error
         return
      endif

      if ( order_wp2_wp3 < 1 .or. order_wp2_wp3 > 4 ) then
         write(fstderr,*) "The variable order_wp2_wp3 must have a value " &
                          // "between 1 and 4"
         write(fstderr,*) "order_wp2_wp3 = ", order_wp2_wp3
         write(fstderr,*) "Fatal error in check_clubb_settings"
         err_code_out = clubb_fatal_error
         return
      elseif ( order_wp2_wp3 == order_xm_wpxp &
               .or. order_wp2_wp3 == order_xp2_xpyp &
               .or. order_wp2_wp3 == order_windm ) then
         write(fstderr,*) "The variable order_wp2_wp3 has the same value " &
                          // "as another order_ variable.  Please give each " &
                          // "order index a unique value."
         write(fstderr,*) "order_wp2_wp3 = ", order_wp2_wp3
         write(fstderr,*) "order_xm_wpxp = ", order_xm_wpxp
         write(fstderr,*) "order_xp2_xpyp = ", order_xp2_xpyp
         write(fstderr,*) "order_windm = ", order_windm
         write(fstderr,*) "Fatal error in check_clubb_settings"
         err_code_out = clubb_fatal_error
         return
      endif

      if ( order_xp2_xpyp < 1 .or. order_xp2_xpyp > 4 ) then
         write(fstderr,*) "The variable order_xp2_xpyp must have a value " &
                          // "between 1 and 4"
         write(fstderr,*) "order_xp2_xpyp = ", order_xp2_xpyp
         write(fstderr,*) "Fatal error in check_clubb_settings"
         err_code_out = clubb_fatal_error
         return
      elseif ( order_xp2_xpyp == order_wp2_wp3 &
               .or. order_xp2_xpyp == order_xm_wpxp &
               .or. order_xp2_xpyp == order_windm ) then
         write(fstderr,*) "The variable order_xp2_xpyp has the same value " &
                          // "as another order_ variable.  Please give each " &
                          // "order index a unique value."
         write(fstderr,*) "order_xp2_xpyp = ", order_xp2_xpyp
         write(fstderr,*) "order_wp2_wp3 = ", order_wp2_wp3
         write(fstderr,*) "order_xm_wpxp = ", order_xm_wpxp
         write(fstderr,*) "order_windm = ", order_windm
         write(fstderr,*) "Fatal error in check_clubb_settings"
         err_code_out = clubb_fatal_error
         return
      endif

      if ( order_windm < 1 .or. order_windm > 4 ) then
         write(fstderr,*) "The variable order_windm must have a value " &
                          // "between 1 and 4"
         write(fstderr,*) "order_windm = ", order_windm
         write(fstderr,*) "Fatal error in check_clubb_settings"
         err_code_out = clubb_fatal_error
         return
      elseif ( order_windm == order_wp2_wp3 &
               .or. order_windm == order_xp2_xpyp &
               .or. order_windm == order_xm_wpxp ) then
         write(fstderr,*) "The variable order_windm has the same value " &
                          // "as another order_ variable.  Please give each " &
                          // "order index a unique value."
         write(fstderr,*) "order_windm = ", order_windm
         write(fstderr,*) "order_wp2_wp3 = ", order_wp2_wp3
         write(fstderr,*) "order_xp2_xpyp = ", order_xp2_xpyp
         write(fstderr,*) "order_xm_wpxp = ", order_xm_wpxp
         write(fstderr,*) "Fatal error in check_clubb_settings"
         err_code_out = clubb_fatal_error
         return
      endif

      ! Checking that when the l_diag_Lscale_from_tau is enabled, the
      ! relevant Cx tunable parameters are all set to a value of 1 (as
      ! you're supposed to tune the C_invrs_tau_ parameters instead).
      if ( clubb_config_flags%l_diag_Lscale_from_tau ) then

         ! Note: someday when we can successfully run with all these parameters
         ! having a value of 1, the "Warning" messages should be removed and the
         ! "Fatal error" messages should be uncommented.

         ! C1 must have a value of 1
         if ( params(iC1) > one .or. params(iC1) < one ) then
            write(fstderr,*) "When the l_diag_Lscale_from_tau flag is " &
                             // "enabled, C1 must have a value of 1."
            write(fstderr,*) "C1 = ", params(iC1)
            write(fstderr,*) "Warning in check_clubb_settings"
            !write(fstderr,*) "Fatal error in check_clubb_settings"
            !err_code_out = clubb_fatal_error
         endif ! C1 check

         ! C1b must have a value of 1
         if ( params(iC1b) > one .or. params(iC1b) < one ) then
            write(fstderr,*) "When the l_diag_Lscale_from_tau flag is " &
                             // "enabled, C1b must have a value of 1."
            write(fstderr,*) "C1b = ", params(iC1b)
            write(fstderr,*) "Warning in check_clubb_settings"
            !write(fstderr,*) "Fatal error in check_clubb_settings"
            !err_code_out = clubb_fatal_error
         endif ! C1b check

         ! C2rt must have a value of 1
         if ( params(iC2rt) > one .or. params(iC2rt) < one ) then
            write(fstderr,*) "When the l_diag_Lscale_from_tau flag is " &
                             // "enabled, C2rt must have a value of 1."
            write(fstderr,*) "C2rt = ", params(iC2rt)
            write(fstderr,*) "Warning in check_clubb_settings"
            !write(fstderr,*) "Fatal error in check_clubb_settings"
            !err_code_out = clubb_fatal_error
         endif ! C2rt check

         ! C2thl must have a value of 1
         if ( params(iC2thl) > one .or. params(iC2thl) < one ) then
            write(fstderr,*) "When the l_diag_Lscale_from_tau flag is " &
                             // "enabled, C2thl must have a value of 1."
            write(fstderr,*) "C2thl = ", params(iC2thl)
            write(fstderr,*) "Warning in check_clubb_settings"
            !write(fstderr,*) "Fatal error in check_clubb_settings"
            !err_code_out = clubb_fatal_error
         endif ! C2thl check

         ! C2rtthl must have a value of 1
         if ( params(iC2rtthl) > one .or. params(iC2rtthl) < one ) then
            write(fstderr,*) "When the l_diag_Lscale_from_tau flag is " &
                             // "enabled, C2rtthl must have a value of 1."
            write(fstderr,*) "C2rtthl = ", params(iC2rtthl)
            write(fstderr,*) "Warning in check_clubb_settings"
            !write(fstderr,*) "Fatal error in check_clubb_settings"
            !err_code_out = clubb_fatal_error
         endif ! C2rtthl check

         ! C6rt must have a value of 1
         if ( params(iC6rt) > one .or. params(iC6rt) < one ) then
            write(fstderr,*) "When the l_diag_Lscale_from_tau flag is " &
                             // "enabled, C6rt must have a value of 1."
            write(fstderr,*) "C6rt = ", params(iC6rt)
            write(fstderr,*) "Warning in check_clubb_settings"
            !write(fstderr,*) "Fatal error in check_clubb_settings"
            !err_code_out = clubb_fatal_error
         endif ! C6rt check

         ! C6rtb must have a value of 1
         if ( params(iC6rtb) > one .or. params(iC6rtb) < one ) then
            write(fstderr,*) "When the l_diag_Lscale_from_tau flag is " &
                             // "enabled, C6rtb must have a value of 1."
            write(fstderr,*) "C6rtb = ", params(iC6rtb)
            write(fstderr,*) "Warning in check_clubb_settings"
            !write(fstderr,*) "Fatal error in check_clubb_settings"
            !err_code_out = clubb_fatal_error
         endif ! C6rtb check

         ! C6thl must have a value of 1
         if ( params(iC6thl) > one .or. params(iC6thl) < one ) then
            write(fstderr,*) "When the l_diag_Lscale_from_tau flag is " &
                             // "enabled, C6thl must have a value of 1."
            write(fstderr,*) "C6thl = ", params(iC6thl)
            write(fstderr,*) "Warning in check_clubb_settings"
            !write(fstderr,*) "Fatal error in check_clubb_settings"
            !err_code_out = clubb_fatal_error
         endif ! C6thl check

         ! C6thlb must have a value of 1
         if ( params(iC6thlb) > one .or. params(iC6thlb) < one ) then
            write(fstderr,*) "When the l_diag_Lscale_from_tau flag is " &
                             // "enabled, C6thlb must have a value of 1."
            write(fstderr,*) "C6thlb = ", params(iC6thlb)
            write(fstderr,*) "Warning in check_clubb_settings"
            !write(fstderr,*) "Fatal error in check_clubb_settings"
            !err_code_out = clubb_fatal_error
         endif ! C6thlb check

         ! C14 must have a value of 1
         if ( params(iC14) > one .or. params(iC14) < one ) then
            write(fstderr,*) "When the l_diag_Lscale_from_tau flag is " &
                             // "enabled, C14 must have a value of 1."
            write(fstderr,*) "C14 = ", params(iC14)
            write(fstderr,*) "Warning in check_clubb_settings"
            !write(fstderr,*) "Fatal error in check_clubb_settings"
            !err_code_out = clubb_fatal_error
         endif ! C14 check

      endif ! l_diag_Lscale_from_tau

      if ( l_implemented ) then

        if ( clubb_config_flags%l_rtm_nudge ) then
          write(fstderr,*) "l_rtm_nudge must be set to .false. when " &
                           // "l_implemented = .true."
          write(fstderr,*) "Fatal error in check_clubb_settings"
          err_code_out = clubb_fatal_error
        end if

        if ( clubb_config_flags%l_uv_nudge ) then
          write(fstderr,*) "l_rtm_nudge must be set to .false. when " &
                           // "l_implemented = .true."
          write(fstderr,*) "Fatal error in check_clubb_settings"
          err_code_out = clubb_fatal_error
        end if

        if ( thlm_sponge_damp_settings%l_sponge_damping ) then
          write(fstderr,*) "thlm_sponge_damp_settings%l_sponge_damping " &
                           // "must be set to .false. when  l_implemented = .true."
          write(fstderr,*) "Fatal error in check_clubb_settings"
          err_code_out = clubb_fatal_error
        end if

        if ( rtm_sponge_damp_settings%l_sponge_damping ) then
          write(fstderr,*) "rtm_sponge_damp_settings%l_sponge_damping " &
                           // "must be set to .false. when  l_implemented = .true."
          write(fstderr,*) "Fatal error in check_clubb_settings"
          err_code_out = clubb_fatal_error
        end if

        if ( uv_sponge_damp_settings%l_sponge_damping ) then
          write(fstderr,*) "uv_sponge_damp_settings%l_sponge_damping " &
                           // "must be set to .false. when  l_implemented = .true."
          write(fstderr,*) "Fatal error in check_clubb_settings"
          err_code_out = clubb_fatal_error
        end if

        if ( wp2_sponge_damp_settings%l_sponge_damping ) then
          write(fstderr,*) "wp2_sponge_damp_settings%l_sponge_damping " &
                           // "must be set to .false. when  l_implemented = .true."
          write(fstderr,*) "Fatal error in check_clubb_settings"
          err_code_out = clubb_fatal_error
        end if

        if ( wp3_sponge_damp_settings%l_sponge_damping ) then
          write(fstderr,*) "wp3_sponge_damp_settings%l_sponge_damping " &
                           // "must be set to .false. when  l_implemented = .true."
          write(fstderr,*) "Fatal error in check_clubb_settings"
          err_code_out = clubb_fatal_error
        end if

        if ( up2_vp2_sponge_damp_settings%l_sponge_damping ) then
          write(fstderr,*) "up2_vp2_sponge_damp_settings%l_sponge_damping " &
                           // "must be set to .false. when  l_implemented = .true."
          write(fstderr,*) "Fatal error in check_clubb_settings"
          err_code_out = clubb_fatal_error
        end if

      end if

      err_code = err_code_out

      return

    end subroutine check_clubb_settings

    !----------------------------------------------------------------------------
    subroutine cleanup_clubb_core( gr )

      ! Description:
      !   Frees memory used by the model itself.
      !
      ! References:
      !   None
      !---------------------------------------------------------------------------

      use grid_class, only: &
        cleanup_grid, & ! Procedure
        grid            ! Type

      implicit none

      type(grid), target, intent(inout) :: gr

      !----- Begin Code -----

      ! De-allocate the arrays for the grid
      call cleanup_grid( gr ) ! intent(in)

      return
    end subroutine cleanup_clubb_core

    !-----------------------------------------------------------------------
    subroutine trapezoidal_rule_zt( nzm, nzt, ngrdcol, sclr_dim, gr,             & ! intent(in)
                                    l_call_pdf_closure_twice,                    & ! intent(in)
                                    stats_metadata,                              & ! intent(in)
                                    wprtp2, wpthlp2,                             & ! intent(inout)
                                    wprtpthlp, cloud_frac, ice_supersat_frac,    & ! intent(inout)
                                    rcm, wp2thvp, wpsclrprtp, wpsclrp2,          & ! intent(inout)
                                    wpsclrpthlp,                                 & ! intent(inout)
                                    wprtp2_zm, wpthlp2_zm,                       & ! intent(inout)
                                    wprtpthlp_zm, cloud_frac_zm,                 & ! intent(inout)
                                    ice_supersat_frac_zm, rcm_zm, wp2thvp_zm,    & ! intent(inout)
                                    wpsclrprtp_zm, wpsclrp2_zm, wpsclrpthlp_zm )   ! intent(inout)
                 
      !
      ! Description:
      !   This subroutine takes the output variables on the thermo.
      !   grid and either: interpolates them to the momentum grid, or uses the
      !   values output from the second call to pdf_closure on momentum levels if
      !   l_call_pdf_closure_twice is true.  It then calls the function
      !   trapezoid_zt to recompute the variables on the thermo. grid.
      !
      !   ldgrant June 2009
      !
      ! Note:
      !   The argument variables in the last 5 lines of the subroutine
      !   (wprtp2_zm through pdf_params_zm) are declared intent(inout) because
      !   if l_call_pdf_closure_twice is true, these variables will already have
      !   values from pdf_closure on momentum levels and will not be altered in
      !   this subroutine.  However, if l_call_pdf_closure_twice is false, these
      !   variables will not have values yet and will be interpolated to
      !   momentum levels in this subroutine.
      ! References:
      !   None
      !-----------------------------------------------------------------------

      use grid_class, only: &
          grid, & ! Type
          zt2zm ! Procedure

      use pdf_parameter_module, only: &
          pdf_parameter ! Derived data type

      use clubb_precision, only: &
          core_rknd ! Variable(s)

      use stats_variables, only: &
          stats_metadata_type

      implicit none

      !------------------------ Input variables ------------------------
      integer, intent(in) :: &
        nzm, &
        nzt, &
        ngrdcol, &
        sclr_dim

      type (grid), target, intent(in) :: &
        gr

      logical, intent(in) :: &
        l_call_pdf_closure_twice

      type (stats_metadata_type), intent(in) :: &
        stats_metadata

      !------------------------ Input/Output variables ------------------------
      ! Thermodynamic level variables output from the first call to pdf_closure
      real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(inout) :: &
        wprtp2,             & ! w'rt'^2                   [m kg^2/kg^2]
        wpthlp2,            & ! w'thl'^2                  [m K^2/s]
        wprtpthlp,          & ! w'rt'thl'                 [m kg K/kg s]
        cloud_frac,         & ! Cloud Fraction            [-]
        ice_supersat_frac,  & ! Ice Cloud Fraction        [-]
        rcm,                & ! Liquid water mixing ratio [kg/kg]
        wp2thvp               ! w'^2 th_v'                [m^2 K/s^2]

      real( kind = core_rknd ), dimension(ngrdcol,nzt,sclr_dim), intent(inout) :: &
        wpsclrprtp,  & ! w'sclr'rt'
        wpsclrp2,    & ! w'sclr'^2
        wpsclrpthlp    ! w'sclr'thl'

      ! Thermo. level variables brought to momentum levels either by
      ! interpolation (in subroutine trapezoidal_rule_zt) or by
      ! the second call to pdf_closure (in subroutine advance_clubb_core)
      real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(inout) :: &
        wprtp2_zm,            & ! w'rt'^2 on momentum grid                   [m kg^2/kg^2]
        wpthlp2_zm,           & ! w'thl'^2 on momentum grid                  [m K^2/s]
        wprtpthlp_zm,         & ! w'rt'thl' on momentum grid                 [m kg K/kg s]
        cloud_frac_zm,        & ! Cloud Fraction on momentum grid            [-]
        ice_supersat_frac_zm, & ! Ice Cloud Fraction on momentum grid        [-]
        rcm_zm,               & ! Liquid water mixing ratio on momentum grid [kg/kg]
        wp2thvp_zm              ! w'^2 th_v' on momentum grid                [m^2 K/s^2]

      real( kind = core_rknd ), dimension(ngrdcol,nzm,sclr_dim), intent(inout) :: &
        wpsclrprtp_zm,  & ! w'sclr'rt' on momentum grid
        wpsclrp2_zm,    & ! w'sclr'^2 on momentum grid
        wpsclrpthlp_zm    ! w'sclr'thl' on momentum grid

      !------------------------ Local variables ------------------------

      integer :: i, sclr

      !----------------------- Begin Code -----------------------------

      ! Store components of pdf_params in the locally declared variables
      ! We only apply the trapezoidal rule to these when
      ! l_apply_rule_to_pdf_params is true.  This is because when we apply the
      ! rule to the final result of pdf_closure rather than the intermediate
      ! results it can lead to an inconsistency in how we determine which
      ! PDF component a point is in and whether the point is in or out of cloud,
      ! which is turn will break the latin hypercube code that samples
      ! preferentially in cloud. -dschanen 13 Feb 2012


      ! If l_call_pdf_closure_twice is true, the _zm variables already have
      ! values from the second call to pdf_closure in advance_clubb_core.
      ! If it is false, the variables are interpolated to the _zm levels.
      if ( .not. l_call_pdf_closure_twice ) then

        ! Interpolate thermodynamic variables to the momentum grid.
        wprtp2_zm                   = zt2zm( nzm, nzt, ngrdcol, gr, wprtp2 )
        wpthlp2_zm                  = zt2zm( nzm, nzt, ngrdcol, gr, wpthlp2 )
        wprtpthlp_zm                = zt2zm( nzm, nzt, ngrdcol, gr, wprtpthlp )
        cloud_frac_zm               = zt2zm( nzm, nzt, ngrdcol, gr, cloud_frac )
        ice_supersat_frac_zm        = zt2zm( nzm, nzt, ngrdcol, gr, ice_supersat_frac )
        rcm_zm                      = zt2zm( nzm, nzt, ngrdcol, gr, rcm )
        wp2thvp_zm                  = zt2zm( nzm, nzt, ngrdcol, gr, wp2thvp )

        ! Since top momentum level is higher than top thermo. level,
        ! set variables at top momentum level to 0.
        !$acc parallel loop gang vector default(present)
        do i = 1, ngrdcol
          wprtp2_zm(i,nzm)             = 0.0_core_rknd
          wpthlp2_zm(i,nzm)            = 0.0_core_rknd
          wprtpthlp_zm(i,nzm)          = 0.0_core_rknd
          cloud_frac_zm(i,nzm)         = 0.0_core_rknd
          ice_supersat_frac_zm(i,nzm)  = 0.0_core_rknd
          rcm_zm(i,nzm)                = 0.0_core_rknd
          wp2thvp_zm(i,nzm)            = 0.0_core_rknd
        end do
        !$acc end parallel loop

        do sclr = 1, sclr_dim
          wpsclrprtp_zm(:,:,sclr)   = zt2zm( nzm, nzt, ngrdcol, gr, wpsclrprtp(:,:,sclr) )
          wpsclrp2_zm(:,:,sclr)     = zt2zm( nzm, nzt, ngrdcol, gr, wpsclrp2(:,:,sclr) )
          wpsclrpthlp_zm(:,:,sclr)  = zt2zm( nzm, nzt, ngrdcol, gr, wpsclrpthlp(:,:,sclr) )

          !$acc parallel loop gang vector default(present)
          do i = 1, ngrdcol
            wpsclrprtp_zm(i,nzm,sclr)  = 0.0_core_rknd
            wpsclrp2_zm(i,nzm,sclr)    = 0.0_core_rknd
            wpsclrpthlp_zm(i,nzm,sclr) = 0.0_core_rknd
          end do
          !$acc end parallel loop
        end do ! sclr = 1, sclr_dim

      end if ! .not. l_call_pdf_closure_twice

      if ( stats_metadata%l_stats ) then
        ! Use the trapezoidal rule to recompute the variables on the stats_zt level
        if ( stats_metadata%iwprtp2 > 0 ) then
          call calc_trapezoid_zt( nzm, nzt, ngrdcol, gr, &
                                  wprtp2, wprtp2_zm, &
                                  wprtp2 )
        end if
        if ( stats_metadata%iwpthlp2 > 0 ) then
          call calc_trapezoid_zt( nzm, nzt, ngrdcol, gr, &
                                  wpthlp2, wpthlp2_zm, &
                                  wpthlp2 )
        end if
        if ( stats_metadata%iwprtpthlp > 0 ) then
          call calc_trapezoid_zt( nzm, nzt, ngrdcol, gr, &
                                  wprtpthlp, wprtpthlp_zm, &
                                  wprtpthlp )
        end if

        do sclr = 1, sclr_dim
          if ( stats_metadata%iwpsclrprtp(sclr) > 0 ) then
            call calc_trapezoid_zt( nzm, nzt, ngrdcol, gr, &
                                    wpsclrprtp(:,:,sclr), wpsclrprtp_zm(:,:,sclr), &
                                    wpsclrprtp(:,:,sclr) )
          end if
          if ( stats_metadata%iwpsclrpthlp(sclr) > 0 ) then
            call calc_trapezoid_zt( nzm, nzt, ngrdcol, gr, &
                                    wpsclrpthlp(:,:,sclr), wpsclrpthlp_zm(:,:,sclr), &
                                    wpsclrpthlp(:,:,sclr) )
          end if
          if ( stats_metadata%iwpsclrp2(sclr) > 0 ) then
            call calc_trapezoid_zt( nzm, nzt, ngrdcol,  gr, &
                                    wpsclrp2(:,:,sclr), wpsclrp2_zm(:,:,sclr), &
                                    wpsclrp2(:,:,sclr) )
          end if
          
        end do ! sclr = 1, sclr_dim
      end if ! l_stats

      call calc_trapezoid_zt( nzm, nzt, ngrdcol, gr, &
                              cloud_frac, cloud_frac_zm, &
                              cloud_frac )
                              
      call calc_trapezoid_zt( nzm, nzt, ngrdcol, gr, &
                              ice_supersat_frac, ice_supersat_frac_zm, &
                              ice_supersat_frac )

      call calc_trapezoid_zt( nzm, nzt, ngrdcol, gr, &
                              rcm, rcm_zm, &
                              rcm )

      call calc_trapezoid_zt( nzm, nzt, ngrdcol, gr, &
                              wp2thvp, wp2thvp_zm, &
                              wp2thvp )

      ! End of trapezoidal rule

      return
    end subroutine trapezoidal_rule_zt
    
    !-----------------------------------------------------------------------
    subroutine trapezoidal_rule_zm( nzm, nzt, ngrdcol, gr,              & ! intent(in)
                                    wpthvp_zt, thlpthvp_zt, rtpthvp_zt, & ! intent(in)
                                    wpthvp, thlpthvp, rtpthvp )           ! intent(inout)
      !
      ! Description:
      !   This subroutine recomputes three variables on the
      !   momentum grid from pdf_closure -- wpthvp, thlpthvp, and
      !   rtpthvp -- by calling the function trapezoid_zm.  Only these three
      !   variables are used in this subroutine because they are the only
      !   pdf_closure momentum variables used elsewhere in CLUBB.
      !
      !   The _zt variables are output from the first call to pdf_closure.
      !   The _zm variables are output from the second call to pdf_closure
      !   on the momentum levels.
      !   This is done before the call to this subroutine.
      !
      !   ldgrant Feb. 2010
      !
      !  References:
      !    None
      !-----------------------------------------------------------------------

      use grid_class, only: grid

      use clubb_precision, only: &
        core_rknd ! variable(s)

      implicit none

      ! ----------------------- Input variables -----------------------
      integer, intent(in) :: &
        nzm, &
        nzt, &
        ngrdcol

      type (grid), target, intent(in) :: gr
    
      real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
        wpthvp_zt,   & ! Buoyancy flux (on thermo. grid)  [(K m)/s]
        thlpthvp_zt, & ! th_l' th_v' (on thermo. grid)    [K^2]
        rtpthvp_zt     ! r_t' th_v' (on thermo. grid)     [(kg K)/kg]

      ! ----------------------- Input/Output variables -----------------------
      real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(inout) :: &
        wpthvp,   & ! Buoyancy flux   [(K m)/s]
        thlpthvp, & ! th_l' th_v'     [K^2]
        rtpthvp     ! r_t' th_v'      [(kg K)/kg]

      ! ----------------------- Begin Code -----------------------

      ! Use the trapezoidal rule to recompute the variables on the zm level
      call calc_trapezoid_zm( nzm, nzt, ngrdcol, gr, wpthvp, wpthvp_zt, & ! Intent(in) 
                              wpthvp )                                    ! Intent(out)
                         
      call calc_trapezoid_zm( nzm, nzt, ngrdcol, gr, thlpthvp, thlpthvp_zt, & ! Intent(in)
                              thlpthvp )                                      ! Intent(out)
                         
      call calc_trapezoid_zm( nzm, nzt, ngrdcol, gr, rtpthvp, rtpthvp_zt, & ! Intent(in)
                              rtpthvp )                                     ! Intent(out)

      return
    end subroutine trapezoidal_rule_zm

    !-----------------------------------------------------------------------
    subroutine calc_trapezoid_zt( nzm, nzt, ngrdcol, gr, &
                                  variable_zt, variable_zm, &
                                  trapezoid_zt )
      !
      ! Description:
      !   Function which uses the trapezoidal rule from calculus
      !   to recompute the values for the variables on the thermo. grid which
      !   are output from the first call to pdf_closure in module clubb_core.
      !
      !   ldgrant June 2009
      !--------------------------------------------------------------------

      use grid_class, only: grid

      use clubb_precision, only: &
        core_rknd ! Variable(s)

      implicit none

      ! ---------------- Input Variables ----------------
      integer, intent(in) :: &
        nzm, &
        nzt, &
        ngrdcol

      type (grid), target, intent(in) :: gr
      
      real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
        variable_zt    ! Variable on the zt grid

      real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: &
        variable_zm    ! Variable on the zm grid

      ! ---------------- Output Variable ----------------
      real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) :: &
        trapezoid_zt

      ! ---------------- Local Variables ----------------
      integer :: i, k ! Loop index

      ! ---------------- Begin Code ----------------

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzt
        do i = 1, ngrdcol
          ! Trapezoidal rule from calculus
          trapezoid_zt(i,k) =  0.5_core_rknd * ( variable_zm(i,k+1) + variable_zt(i,k) ) &
                               * ( gr%zm(i,k+1) - gr%zt(i,k) ) * gr%invrs_dzt(i,k) &
                               + 0.5_core_rknd * ( variable_zt(i,k) + variable_zm(i,k) ) &
                                 * ( gr%zt(i,k) - gr%zm(i,k) ) * gr%invrs_dzt(i,k)
        end do
      end do ! k = 2, gr%nzt
      !$acc end parallel loop

      return
    end subroutine calc_trapezoid_zt

    !-----------------------------------------------------------------------
    subroutine calc_trapezoid_zm( nzm, nzt, ngrdcol, gr, variable_zm, variable_zt, &
                                  trapezoid_zm )
      !
      ! Description:
      !   Function which uses the trapezoidal rule from calculus
      !   to recompute the values for the important variables on the momentum
      !   grid which are output from pdf_closure in module clubb_core.
      !   These momentum variables only include wpthvp, thlpthvp, and rtpthvp.
      !
      !   ldgrant Feb. 2010
      !--------------------------------------------------------------------

      use grid_class, only: grid

      use clubb_precision, only: &
          core_rknd ! Variable(s)

      implicit none

      ! -------------------- Input Variables --------------------
      integer, intent(in) :: &
        nzm, &
        nzt, &
        ngrdcol

      type (grid), target, intent(in) :: gr
      
      real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: &
        variable_zm    ! Variable on the zm grid

      real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
        variable_zt    ! Variable on the zt grid

      ! -------------------- Output Variable --------------------
      real( kind = core_rknd ), intent(out), dimension(ngrdcol,nzm) :: &
        trapezoid_zm
 
      ! -------------------- Local Variables --------------------
      integer :: i, k ! Loop index

      ! -------------------- Begin Code --------------------

      ! Boundary conditions: trapezoidal rule not valid at top zm level, nzmax.
      ! Trapezoidal rule also not used at zm level 1.
      !$acc parallel loop gang vector default(present) 
      do i = 1, ngrdcol
        trapezoid_zm(i,1)   = variable_zm(i,1)
        trapezoid_zm(i,nzm) = variable_zm(i,nzm)
      end do
      !$acc end parallel loop

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 2, nzm-1
        do i = 1, ngrdcol
          ! Trapezoidal rule from calculus
          trapezoid_zm(i,k) =  0.5_core_rknd * ( variable_zt(i,k) + variable_zm(i,k) ) &
                               * ( gr%zt(i,k) - gr%zm(i,k) ) * gr%invrs_dzm(i,k) &
                               + 0.5_core_rknd * ( variable_zm(i,k) + variable_zt(i,k-1) ) &
                                 * ( gr%zm(i,k) - gr%zt(i,k-1) ) * gr%invrs_dzm(i,k)
        end do
      end do 
      !$acc end parallel loop

      return
    end subroutine calc_trapezoid_zm

    !-----------------------------------------------------------------------
    subroutine compute_cloud_cover( gr, nzt, ngrdcol, &
                                    pdf_params, cloud_frac, rcm, & ! intent(in)
                                    cloud_cover, rcm_in_layer )    ! intent(out)
      !
      ! Description:
      !   Subroutine to compute cloud cover (the amount of sky
      !   covered by cloud) and rcm in layer (liquid water mixing ratio in
      !   the portion of the grid box filled by cloud).
      !
      ! References:
      !   Definition of 's' comes from:
      !   ``The Gaussian Cloud Model Relations'' G. L. Mellor (1977)
      !   JAS, Vol. 34, pp. 356--358.
      !
      ! Notes:
      !   Added July 2009
      !---------------------------------------------------------------------

      use constants_clubb, only: &
          rc_tol, & ! Variable(s)
          fstderr, &
          unused_var

      use grid_class, only: grid

      use pdf_parameter_module, only: &
          pdf_parameter ! Derived data type

      use clubb_precision, only: &
          core_rknd ! Variable(s)

      use error_code, only: &
        clubb_at_least_debug_level,  & ! Procedure
        err_code,                    & ! Error Indicator
        clubb_fatal_error              ! Constant

      implicit none

      !------------------------ Input variables ------------------------
      integer, intent(in) :: &
        ngrdcol,  & ! Number of grid columns
        nzt         ! Number of thermodynamic vertical levels

      type (grid), target, intent(in) :: gr

      real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
        cloud_frac, & ! Cloud fraction             [-]
        rcm           ! Liquid water mixing ratio  [kg/kg]

      type (pdf_parameter), intent(in) :: &
        pdf_params    ! PDF Parameters  [units vary]

      !------------------------ Output variables ------------------------
      real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) :: &
        cloud_cover,  & ! Cloud cover                               [-]
        rcm_in_layer    ! Liquid water mixing ratio in cloud layer  [kg/kg]

      !------------------------ Local variables ------------------------
      real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
        chi_mean,              & ! Mean extended cloud water mixing ratio of the
                                 ! two Gaussian distributions
        vert_cloud_frac_upper, & ! Fraction of cloud in top half of grid box
        vert_cloud_frac_lower, & ! Fraction of cloud in bottom half of grid box
        vert_cloud_frac          ! Fraction of cloud filling the grid box in the vertical

      integer :: i, k, kp1, km1

      !------------------------ Begin code ------------------------

      !$acc enter data create( chi_mean, vert_cloud_frac_upper, &
      !$acc                    vert_cloud_frac_lower, vert_cloud_frac )

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzt
        do i = 1, ngrdcol

          chi_mean(i,k) =      pdf_params%mixt_frac(i,k)  * pdf_params%chi_1(i,k) + &
                      (1.0_core_rknd-pdf_params%mixt_frac(i,k)) * pdf_params%chi_2(i,k)

        end do
      end do
      !$acc end parallel loop

      !$acc parallel loop gang vector collapse(2) default(present) reduction(.or.:err_code)
      do k = 1, nzt-1
        do i = 1, ngrdcol

          km1 = max( k-1, 1 )
          kp1 = min( k+1, nzt )

          if ( rcm(i,k) < rc_tol ) then ! No cloud at this level

            cloud_cover(i,k)  = cloud_frac(i,k)
            rcm_in_layer(i,k) = rcm(i,k)

          else if ( ( rcm(i,kp1) >= rc_tol ) .and. ( rcm(i,km1) >= rc_tol ) ) then
            ! There is cloud above and below,
            !   so assume cloud fills grid box from top to bottom

            cloud_cover(i,k) = cloud_frac(i,k)
            rcm_in_layer(i,k) = rcm(i,k)

          else if ( ( rcm(i,kp1) < rc_tol ) .or. ( rcm(i,km1) < rc_tol) ) then
            ! Cloud may fail to reach gridbox top or base or both

            ! First let the cloud fill the entire grid box, then overwrite
            ! vert_cloud_frac_upper(k) and/or vert_cloud_frac_lower(k)
            ! for a cloud top, cloud base, or one-point cloud.
            vert_cloud_frac_upper(i,k) = 0.5_core_rknd
            vert_cloud_frac_lower(i,k) = 0.5_core_rknd

            if ( rcm(i,kp1) < rc_tol ) then ! Cloud top

              vert_cloud_frac_upper(i,k) = &
                       ( ( 0.5_core_rknd / gr%invrs_dzm(i,kp1) ) / ( gr%zm(i,kp1) - gr%zt(i,k) ) ) &
                       * ( rcm(i,k) / ( rcm(i,k) + abs( chi_mean(i,kp1) ) ) )

              vert_cloud_frac_upper(i,k) = min( 0.5_core_rknd, vert_cloud_frac_upper(i,k) )

              ! Make the transition in cloudiness more gradual than using
              ! the above min statement alone.
              vert_cloud_frac_upper(i,k) = vert_cloud_frac_upper(i,k) + &
                ( ( rcm(i,kp1)/rc_tol )*( 0.5_core_rknd -vert_cloud_frac_upper(i,k) ) )

            else

              vert_cloud_frac_upper(i,k) = 0.5_core_rknd

            end if

            if ( rcm(i,km1) < rc_tol ) then ! Cloud base

              vert_cloud_frac_lower(i,k) = &
                       ( ( 0.5_core_rknd / gr%invrs_dzm(i,k) ) / ( gr%zt(i,k) - gr%zm(i,k) ) ) &
                       * ( rcm(i,k) / ( rcm(i,k) + abs( chi_mean(i,km1) ) ) )

              vert_cloud_frac_lower(i,k) = min( 0.5_core_rknd, vert_cloud_frac_lower(i,k) )

              ! Make the transition in cloudiness more gradual than using
              ! the above min statement alone.
              vert_cloud_frac_lower(i,k) = vert_cloud_frac_lower(i,k) + &
                ( ( rcm(i,km1)/rc_tol )*( 0.5_core_rknd -vert_cloud_frac_lower(i,k) ) )

            else

              vert_cloud_frac_lower(i,k) = 0.5_core_rknd

            end if

            vert_cloud_frac(i,k) = &
              vert_cloud_frac_upper(i,k) + vert_cloud_frac_lower(i,k)

            vert_cloud_frac(i,k) = &
              max( cloud_frac(i,k), min( 1.0_core_rknd, vert_cloud_frac(i,k) ) )

            cloud_cover(i,k)  = cloud_frac(i,k) / vert_cloud_frac(i,k)
            rcm_in_layer(i,k) = rcm(i,k) / vert_cloud_frac(i,k)

          else

            ! This case should not be entered
            cloud_cover(i,k) = unused_var
            rcm_in_layer(i,k) = unused_var
            err_code = clubb_fatal_error

          end if ! rcm(k) < rc_tol
          
        end do
      end do ! k = 1, gr%nzt-1, 1
      !$acc end parallel loop

      !$acc parallel loop gang vector default(present)
      do i = 1, ngrdcol
        cloud_cover(i,nzt) = cloud_frac(i,nzt)
        rcm_in_layer(i,nzt) = rcm(i,nzt)
      end do
      !$acc end parallel loop

      if ( clubb_at_least_debug_level( 0 ) ) then
        if ( err_code == clubb_fatal_error ) then

          !$acc update host( pdf_params%mixt_frac, pdf_params%chi_1, pdf_params%chi_2, &
          !$acc              cloud_frac, rcm )

          write(fstderr,*)  &
             "ERROR: compute_cloud_cover entered a conditional case it should not have "

          write(fstderr,*) "pdf_params%mixt_frac = ", pdf_params%mixt_frac
          write(fstderr,*) "pdf_params%chi_1 = ", pdf_params%chi_1
          write(fstderr,*) "pdf_params%chi_2 = ", pdf_params%chi_2
          write(fstderr,*) "cloud_frac = ", cloud_frac
          write(fstderr,*) "rcm = ", rcm
        end if
      end if 

      !$acc exit data delete( chi_mean, vert_cloud_frac_upper, &
      !$acc                   vert_cloud_frac_lower, vert_cloud_frac )

      return

    end subroutine compute_cloud_cover

    !-----------------------------------------------------------------------
    subroutine clip_rcm ( nzt, ngrdcol, rtm, & ! intent(in)
                          message,           & ! intent(in)
                          rcm )                ! intent(inout)
      !
      ! Description:
      !   Subroutine that reduces cloud water (rcm) whenever
      !   it exceeds total water (rtm = vapor + liquid).
      !   This avoids negative values of rvm = water vapor mixing ratio.
      !   However, it will not ensure that rcm <= rtm if rtm <= 0.
      !
      ! References:
      !   None
      !---------------------------------------------------------------------

      use error_code, only: &
        clubb_at_least_debug_level  ! Procedure

      use constants_clubb, only: &
        fstderr, & ! Variable(s)
        zero_threshold

      use clubb_precision, only: &
        core_rknd ! Variable(s)

      implicit none

      ! Input variables
      integer, intent(in) :: &
        nzt, &
        ngrdcol
    
      real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
        rtm           ! Total water mixing ratio             [kg/kg]

      character(len= * ), intent(in) :: message

      real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(inout) :: &
        rcm           ! Cloud water mixing ratio  [kg/kg]

      integer :: i, k

      ! ------------ Begin code ---------------

      !$acc data copyin( rtm ) &
      !$acc        copy( rcm )

      if ( clubb_at_least_debug_level( 3 ) ) then

        !$acc update host( rcm, rtm )

        do k = 1, nzt
          do i = 1, ngrdcol 
            
            if ( rtm(i,k) < rcm(i,k) ) then

              write(fstderr,*) message, ' at k=', k, ' at i=', i, 'rcm(k) = ', rcm(i,k), &
                'rtm(k) = ', rtm(i,k), '.',  ' Clipping rcm.'

            end if ! rtm(k) < rcm(k)
            
          end do
        end do
      end if ! clubb_at_least_debug_level( 3 )

      ! Vince Larson clipped rcm in order to prevent rvm < 0.  5 Apr 2008.
      ! This code won't work unless rtm >= 0 !!!
      ! We do not clip rcm_in_layer because rcm_in_layer only influences
      ! radiation, and we do not want to bother recomputing it.  6 Aug 2009
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzt
        do i = 1, ngrdcol 
          if ( rtm(i,k) < rcm(i,k) ) then
            rcm(i,k) = max( zero_threshold, rtm(i,k) - epsilon( rtm(i,k) ) )
          end if ! rtm(k) < rcm(k)
        end do
      end do
      !$acc end parallel loop

      !$acc end data

      return
    end subroutine clip_rcm

    !-----------------------------------------------------------------------------
    subroutine set_Lscale_max( ngrdcol, l_implemented, host_dx, host_dy, &
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

      use clubb_precision, only: &
        core_rknd ! Variable(s)

      implicit none

      !----------------------- Input Variables -----------------------
      integer, intent(in) :: &
        ngrdcol
      
      logical, intent(in) :: &
        l_implemented     ! Flag to see if CLUBB is running on it's own,
                          ! or if it's implemented as part of a host model.

      real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
        host_dx, & ! Host model's east-west horizontal grid spacing     [m]
        host_dy    ! Host model's north-south horizontal grid spacing   [m]

      !----------------------- Output Variable -----------------------
      real( kind = core_rknd ), dimension(ngrdcol), intent(out) :: &
        Lscale_max    ! Maximum allowable value for Lscale   [m]

      !----------------------- Local Variable -----------------------
      integer :: i

      !----------------------- Begin Code-----------------------

      ! Determine the maximum allowable value for Lscale (in meters).
      if ( l_implemented ) then
        !$acc parallel loop gang vector default(present)
        do i = 1, ngrdcol
          Lscale_max(i) = 0.25_core_rknd * min( host_dx(i), host_dy(i) )
        end do
        !$acc end parallel loop
      else
        !$acc parallel loop gang vector default(present)
        do i = 1, ngrdcol
          Lscale_max(i) = 1.0e5_core_rknd
        end do
        !$acc end parallel loop
      end if

      return
    end subroutine set_Lscale_max

!===============================================================================
  subroutine calculate_thlp2_rad &
                  ( nzm, rcm_zm, thlprcp, radht_zm, & ! Intent(in)
                    clubb_params,                  & ! Intent(in)
                    thlp2_forcing )                  ! Intent(inout)

  ! Description:
  !   Computes the contribution of radiative cooling to thlp2

  ! References:
  !   See clubb:ticket:632
  !----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd                     ! Constant(s)

    use constants_clubb, only: &
        two, &
        rc_tol

    use parameter_indices, only: &
        nparams, & ! Variable(s)
        ithlp2_rad_coef

    implicit none

  ! Input Variables
    integer, intent(in) :: &
      nzm                   ! Number of momentum vertical levels             [-]

    real( kind = core_rknd ), dimension(nzm), intent(in) :: &
      rcm_zm, &             ! Cloud water mixing ratio on momentum grid      [kg/kg]
      thlprcp, &            ! thl'rc'                                        [K kg/kg]
      radht_zm              ! SW + LW heating rate (on momentum grid)        [K/s]

    real( kind = core_rknd ), dimension(nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

  ! Input/Output Variables
    real( kind = core_rknd ), dimension(nzm), intent(inout) :: &
      thlp2_forcing         ! <th_l'^2> forcing (momentum levels)            [K^2/s]

  ! Local Variables
    integer :: &
      k                     ! Loop iterator                                  [-]

  !----------------------------------------------------------------------


    do k = 1, nzm

       if ( rcm_zm(k) > rc_tol ) then

          thlp2_forcing(k) &
          = thlp2_forcing(k) + &
            clubb_params(ithlp2_rad_coef) * ( two ) * radht_zm(k) / rcm_zm(k) * thlprcp(k)

       end if

    end do


    return
  end subroutine calculate_thlp2_rad


    !-----------------------------------------------------------------------
end module advance_clubb_core_module
