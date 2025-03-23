module generalized_grid_test

  implicit none

  public :: clubb_generalized_grid_testing
  private :: check_flipped_results 

  contains

  !-----------------------------------------------------------------------
  subroutine clubb_generalized_grid_testing &
             ( gr, gr_desc, nzm, nzt, ngrdcol, &                          ! Intent(in)
               l_implemented, dt, fcor, sfc_elevation, &                  ! Intent(in)
               hydromet_dim, &                                            ! intent(in)
               sclr_dim, sclr_tol, edsclr_dim, sclr_idx, &                ! intent(in)
               thlm_forcing, rtm_forcing, um_forcing, vm_forcing, &       ! Intent(in)
               sclrm_forcing, edsclrm_forcing, wprtp_forcing, &           ! Intent(in)
               wpthlp_forcing, rtp2_forcing, thlp2_forcing, &             ! Intent(in)
               rtpthlp_forcing, wm_zm, wm_zt, &                           ! Intent(in)
               wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, p_sfc, &        ! Intent(in)
               wpsclrp_sfc, wpedsclrp_sfc,  &                             ! Intent(in)
               upwp_sfc_pert, vpwp_sfc_pert, &                            ! intent(in)
               rtm_ref, thlm_ref, um_ref, vm_ref, ug, vg, &               ! Intent(in)
               p_in_Pa, rho_zm, rho, exner, &                             ! Intent(in)
               rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &                   ! Intent(in)
               invrs_rho_ds_zt, thv_ds_zm, thv_ds_zt, &                   ! Intent(in) 
               l_mix_rat_hm, &                                            ! Intent(in)
               rfrzm, wphydrometp, &                                      ! Intent(in)
               wp2hmp, rtphmp_zt, thlphmp_zt, &                           ! Intent(in)
               host_dx, host_dy, &                                        ! Intent(in)
               clubb_params, nu_vert_res_dep, lmin, &                     ! Intent(in)
               clubb_config_flags, &                                      ! Intent(in)
               stats_metadata, &                                          ! Intent(in)
               stats_zt, stats_zm, stats_sfc, &                           ! intent(inout)
               um, vm, upwp, vpwp, up2, vp2, up3, vp3, &                  ! Intent(inout)
               thlm, rtm, wprtp, wpthlp, &                                ! Intent(inout)
               wp2, wp3, rtp2, rtp3, thlp2, thlp3, rtpthlp, &             ! Intent(inout)
               sclrm, sclrp2, sclrp3, sclrprtp, sclrpthlp, &              ! Intent(inout)
               wpsclrp, edsclrm, err_code, &                              ! Intent(inout)
               rcm, cloud_frac, &                                         ! Intent(inout)
               wpthvp, wp2thvp, rtpthvp, thlpthvp, &                      ! Intent(inout)
               sclrpthvp, &                                               ! Intent(inout)
               wp2rtp, wp2thlp, uprcp, vprcp, rc_coef_zm, wp4, &          ! intent(inout)
               wpup2, wpvp2, wp2up2, wp2vp2, ice_supersat_frac, &         ! intent(inout)
               um_pert, vm_pert, upwp_pert, vpwp_pert, &                  ! intent(inout)
               pdf_params, pdf_params_zm, &                               ! Intent(inout)
               pdf_implicit_coefs_terms, &                                ! intent(inout)
               Kh_zm, Kh_zt, &                                            ! intent(out)
               thlprcp, wprcp, w_up_in_cloud, w_down_in_cloud, &          ! Intent(out)
               cloudy_updraft_frac, cloudy_downdraft_frac, &              ! Intent(out)
               rcm_in_layer, cloud_cover, invrs_tau_zm, &                 ! Intent(out)
               Lscale )                                                   ! Intent(out)

    use grid_class, only: &
        grid, & ! Type(s)
        flip    ! Procedure(s)

    use clubb_api_module, only: &
        advance_clubb_core_api

    use clubb_precision, only: &
        core_rknd

    use array_index, only: &
        sclr_idx_type

    use parameter_indices, only: &
        nparams

    use parameters_tunable, only: &
        nu_vertical_res_dep

    use pdf_parameter_module, only: &
        pdf_parameter,                 &
        implicit_coefs_terms,          &
        init_pdf_params,               &
        init_pdf_implicit_coefs_terms

    use model_flags, only: &
        clubb_config_flags_type

    use stats_type, only: &
        stats ! Type

    use stats_variables, only: &
        stats_metadata_type

    use error_code, only: &
        clubb_generalized_grd_test_err, &
        clubb_no_error, &
        clubb_fatal_error

    implicit none

    !------------------------ Input Variables ----------------------------
    type(grid), intent(in) :: &
      gr,      & ! Grid type variable for ascending grid
      gr_desc    ! Grid type variable for descending grid

    integer, intent(in) :: &
      nzm,     &
      nzt,     &
      ngrdcol

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

    logical, dimension(hydromet_dim), intent(in) :: &
      l_mix_rat_hm   ! if true, then the quantity is a hydrometeor mixing ratio

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

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) :: &
      Lscale     ! Length scale         [m]

    integer, intent(inout) :: err_code  ! Diagnostic, for if some calculation goes amiss.

#ifdef GFDL
    ! hlg, 2010-06-16
    real( kind = core_rknd ), intent(inOUT), dimension(ngrdcol,nzt, min(1,sclr_dim) , 2) :: &
      RH_crit  ! critical relative humidity for droplet and ice nucleation
    logical, intent(in)                 ::  do_liquid_only_in_clubb
#endif

    !------------------------ Local Variables ------------------------

    real( kind = core_rknd ), dimension(ngrdcol,nzt) ::  &
      thlm_forcing_flip,    & ! theta_l forcing (thermodynamic levels)    [K/s]
      rtm_forcing_flip,     & ! r_t forcing (thermodynamic levels)        [(kg/kg)/s]
      um_forcing_flip,      & ! u wind forcing (thermodynamic levels)     [m/s/s]
      vm_forcing_flip,      & ! v wind forcing (thermodynamic levels)     [m/s/s]
      wm_zt_flip,           & ! w mean wind component on thermo. levels   [m/s]
      rho_flip,             & ! Air density on thermodynamic levels       [kg/m^3]
      rho_ds_zt_flip,       & ! Dry, static density on thermo. levels     [kg/m^3]
      invrs_rho_ds_zt_flip, & ! Inv. dry, static density @ thermo. levs.  [m^3/kg]
      thv_ds_zt_flip,       & ! Dry, base-state theta_v on thermo. levs.  [K]
      rfrzm_flip              ! Total ice-phase water mixing ratio        [kg/kg]

    real( kind = core_rknd ), dimension(ngrdcol,nzm) ::  &
      wprtp_forcing_flip,   & ! <w'r_t'> forcing (momentum levels)    [m*K/s^2]
      wpthlp_forcing_flip,  & ! <w'th_l'> forcing (momentum levels)   [m*(kg/kg)/s^2]
      rtp2_forcing_flip,    & ! <r_t'^2> forcing (momentum levels)    [(kg/kg)^2/s]
      thlp2_forcing_flip,   & ! <th_l'^2> forcing (momentum levels)   [K^2/s]
      rtpthlp_forcing_flip, & ! <r_t'th_l'> forcing (momentum levels) [K*(kg/kg)/s]
      wm_zm_flip,           & ! w mean wind component on momentum levels  [m/s]
      rho_zm_flip,          & ! Air density on momentum levels            [kg/m^3]
      rho_ds_zm_flip,       & ! Dry, static density on momentum levels    [kg/m^3]
      invrs_rho_ds_zm_flip, & ! Inv. dry, static density @ momentum levs. [m^3/kg]
      thv_ds_zm_flip          ! Dry, base-state theta_v on momentum levs. [K]

    real( kind = core_rknd ), dimension(ngrdcol,nzm, hydromet_dim) :: &
      wphydrometp_flip    ! Covariance of w and a hydrometeor   [(m/s) <hm units>]

    real( kind = core_rknd ), dimension(ngrdcol,nzt, hydromet_dim) :: &
      wp2hmp_flip,      & ! Third moment: <w'^2> * <hydro.'>    [(m/s)^2 <hm units>]
      rtphmp_zt_flip,   & ! Covariance of rt and a hydrometeor  [(kg/kg) <hm units>]
      thlphmp_zt_flip     ! Covariance of thl and a hydrometeor [K <hm units>]

    real( kind = core_rknd ), dimension(ngrdcol,nzt,sclr_dim) :: &
      sclrm_forcing_flip    ! Passive scalar forcing         [{units vary}/s]

    real( kind = core_rknd ), dimension(ngrdcol,nzt,edsclr_dim) :: &
      edsclrm_forcing_flip  ! Eddy passive scalar forcing    [{units vary}/s]

    real( kind = core_rknd ), dimension(ngrdcol,nzt) ::  &
      rtm_ref_flip,  & ! Initial total water mixing ratio             [kg/kg]
      thlm_ref_flip, & ! Initial liquid water potential temperature   [K]
      um_ref_flip,   & ! Initial u wind; Michael Falk                 [m/s]
      vm_ref_flip,   & ! Initial v wind; Michael Falk                 [m/s]
      ug_flip,       & ! u geostrophic wind                           [m/s]
      vg_flip          ! v geostrophic wind                           [m/s]

    type (stats_metadata_type) :: &
      stats_metadata_flip

    real( kind = core_rknd ), dimension(ngrdcol,nzt) ::  &
      um_flip,      & ! u mean wind component (thermodynamic levels)   [m/s]
      vm_flip,      & ! v mean wind component (thermodynamic levels)   [m/s]
      up3_flip,     & ! u'^3 (thermodynamic levels)                    [m^3/s^3]
      vp3_flip,     & ! v'^3 (thermodynamic levels)                    [m^3/s^3]
      rtm_flip,     & ! total water mixing ratio, r_t (thermo. levels) [kg/kg]
      thlm_flip,    & ! liq. water pot. temp., th_l (thermo. levels)   [K]
      rtp3_flip,    & ! r_t'^3 (thermodynamic levels)                  [(kg/kg)^3]
      thlp3_flip,   & ! th_l'^3 (thermodynamic levels)                 [K^3]
      wp3_flip        ! w'^3 (thermodynamic levels)                    [m^3/s^3]

    real( kind = core_rknd ), dimension(ngrdcol,nzm) ::  &
      upwp_flip,    & ! u'w' (momentum levels)                         [m^2/s^2]
      vpwp_flip,    & ! v'w' (momentum levels)                         [m^2/s^2]
      up2_flip,     & ! u'^2 (momentum levels)                         [m^2/s^2]
      vp2_flip,     & ! v'^2 (momentum levels)                         [m^2/s^2]
      wprtp_flip,   & ! w' r_t' (momentum levels)                      [(kg/kg) m/s]
      wpthlp_flip,  & ! w' th_l' (momentum levels)                     [(m/s) K]
      rtp2_flip,    & ! r_t'^2 (momentum levels)                       [(kg/kg)^2]
      thlp2_flip,   & ! th_l'^2 (momentum levels)                      [K^2]
      rtpthlp_flip, & ! r_t' th_l' (momentum levels)                   [(kg/kg) K]
      wp2_flip        ! w'^2 (momentum levels)                         [m^2/s^2]

    real( kind = core_rknd ), dimension(ngrdcol,nzt,sclr_dim) :: &
      sclrm_flip,     & ! Passive scalar mean (thermo. levels) [units vary]
      sclrp3_flip       ! sclr'^3 (thermodynamic levels)       [{units vary}^3]

    real( kind = core_rknd ), dimension(ngrdcol,nzm,sclr_dim) :: &
      wpsclrp_flip,   & ! w'sclr' (momentum levels)            [{units vary} m/s]
      sclrp2_flip,    & ! sclr'^2 (momentum levels)            [{units vary}^2]
      sclrprtp_flip,  & ! sclr'rt' (momentum levels)           [{units vary} (kg/kg)]
      sclrpthlp_flip    ! sclr'thl' (momentum levels)          [{units vary} K]

    real( kind = core_rknd ), dimension(ngrdcol,nzt) ::  &
      p_in_Pa_flip, & ! Air pressure (thermodynamic levels)       [Pa]
      exner_flip      ! Exner function (thermodynamic levels)     [-]

    real( kind = core_rknd ), dimension(ngrdcol,nzt) ::  &
      rcm_flip,        & ! cloud water mixing ratio, r_c (thermo. levels) [kg/kg]
      cloud_frac_flip, & ! cloud fraction (thermodynamic levels)          [-]
      wp2thvp_flip       ! < w'^2 th_v' > (thermodynamic levels)          [m^2/s^2 K]

    real( kind = core_rknd ), dimension(ngrdcol,nzm) ::  &
      wpthvp_flip,     & ! < w' th_v' > (momentum levels)                 [kg/kg K]
      rtpthvp_flip,    & ! < r_t' th_v' > (momentum levels)               [kg/kg K]
      thlpthvp_flip      ! < th_l' th_v' > (momentum levels)              [K^2]

    real( kind = core_rknd ), dimension(ngrdcol,nzm,sclr_dim) :: &
      sclrpthvp_flip     ! < sclr' th_v' > (momentum levels)   [units vary]

    real( kind = core_rknd ), dimension(ngrdcol,nzt) ::  &
      wp2rtp_flip,            & ! w'^2 rt' (thermodynamic levels)      [m^2/s^2 kg/kg]
      wp2thlp_flip,           & ! w'^2 thl' (thermodynamic levels)     [m^2/s^2 K]
      wpup2_flip,             & ! w'u'^2 (thermodynamic levels)        [m^3/s^3]
      wpvp2_flip,             & ! w'v'^2 (thermodynamic levels)        [m^3/s^3]
      ice_supersat_frac_flip    ! ice cloud fraction (thermo. levels)  [-]

    real( kind = core_rknd ), dimension(ngrdcol,nzm) ::  &
      uprcp_flip,             & ! < u' r_c' > (momentum levels)        [(m/s)(kg/kg)]
      vprcp_flip,             & ! < v' r_c' > (momentum levels)        [(m/s)(kg/kg)]
      rc_coef_zm_flip,        & ! Coef of X'r_c' in Eq. (34) (m-levs.) [K/(kg/kg)]
      wp4_flip,               & ! w'^4 (momentum levels)               [m^4/s^4]
      wp2up2_flip,            & ! w'^2 u'^2 (momentum levels)          [m^4/s^4]
      wp2vp2_flip               ! w'^2 v'^2 (momentum levels)          [m^4/s^4]

    type(pdf_parameter) :: &
      pdf_params_flip,    & ! PDF parameters (thermodynamic levels)    [units vary]
      pdf_params_zm_flip    ! PDF parameters on momentum levels        [units vary]

    type(implicit_coefs_terms) :: &
      pdf_implicit_coefs_terms_flip    ! Implicit coefs / explicit terms [units vary]

#ifdef GFDL
    real( kind = core_rknd ), dimension(ngrdcol,nzt,sclr_dim) :: &  ! h1g, 2010-06-16
      sclrm_trsport_only  ! Passive scalar concentration due to pure transport [{units vary}/s]
#endif

    real( kind = core_rknd ), dimension(ngrdcol,nzt,edsclr_dim) :: &
      edsclrm_flip   ! Eddy passive scalar mean (thermo. levels)   [units vary]

    real( kind = core_rknd ), dimension(ngrdcol,nzt) ::  &
      rcm_in_layer_flip, & ! rcm in cloud layer                              [kg/kg]
      cloud_cover_flip     ! cloud cover                                     [-]

    ! Variables that need to be output for use in host models
    real( kind = core_rknd ), dimension(ngrdcol,nzm) ::  &
      wprcp_flip,                 & ! w'r_c' (momentum levels)              [(kg/kg) m/s]
      invrs_tau_zm_flip             ! One divided by tau on zm levels       [1/s]

    real( kind = core_rknd ), dimension(ngrdcol,nzt) ::  &
      w_up_in_cloud_flip,         & ! Average cloudy updraft velocity       [m/s]
      w_down_in_cloud_flip,       & ! Average cloudy downdraft velocity     [m/s]
      cloudy_updraft_frac_flip,   & ! cloudy updraft fraction               [-]
      cloudy_downdraft_frac_flip    ! cloudy downdraft fraction             [-]

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      Kh_zt_flip    ! Eddy diffusivity coefficient on thermodynamic levels   [m^2/s]

    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
      Kh_zm_flip    ! Eddy diffusivity coefficient on momentum levels        [m^2/s]

#ifdef CLUBB_CAM
    real( kind = core_rknd), dimension(ngrdcol,nzt) :: &
      qclvar_flip        ! cloud water variance
#endif

    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
      thlprcp_flip    ! thl'rc'              [K kg/kg]

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      Lscale_flip     ! Length scale         [m]

    integer :: &
      i, sclr, edsclr, hm_idx

    logical :: l_differences = .false.


      ! Set up "flipped" variables for call to descending grid.
      do i = 1, ngrdcol

         thlm_forcing_flip(i,:) = flip( thlm_forcing(i,:), nzt )
         rtm_forcing_flip(i,:) = flip( rtm_forcing(i,:), nzt )
         um_forcing_flip(i,:) = flip( um_forcing(i,:), nzt )
         vm_forcing_flip(i,:) = flip( vm_forcing(i,:), nzt )
         wprtp_forcing_flip(i,:) = flip( wprtp_forcing(i,:), nzm )
         wpthlp_forcing_flip(i,:) = flip( wpthlp_forcing(i,:), nzm )
         rtp2_forcing_flip(i,:) = flip( rtp2_forcing(i,:), nzm )
         thlp2_forcing_flip(i,:) = flip( thlp2_forcing(i,:), nzm )
         rtpthlp_forcing_flip(i,:) = flip( rtpthlp_forcing(i,:), nzm )
         wm_zm_flip(i,:) = flip( wm_zm(i,:), nzm )
         wm_zt_flip(i,:) = flip( wm_zt(i,:), nzt )
         rtm_ref_flip(i,:) = flip( rtm_ref(i,:), nzt )
         thlm_ref_flip(i,:) = flip( thlm_ref(i,:), nzt )
         um_ref_flip(i,:) = flip( um_ref(i,:), nzt )
         vm_ref_flip(i,:) = flip( vm_ref(i,:), nzt )
         ug_flip(i,:) = flip( ug(i,:), nzt )
         vg_flip(i,:) = flip( vg(i,:), nzt )
         p_in_Pa_flip(i,:) = flip( p_in_Pa(i,:), nzt )
         rho_zm_flip(i,:) = flip( rho_zm(i,:), nzm )
         rho_flip(i,:) = flip( rho(i,:), nzt )
         exner_flip(i,:) = flip( exner(i,:), nzt )
         rho_ds_zm_flip(i,:) = flip( rho_ds_zm(i,:), nzm )
         rho_ds_zt_flip(i,:) = flip( rho_ds_zt(i,:), nzt )
         invrs_rho_ds_zm_flip(i,:) = flip( invrs_rho_ds_zm(i,:), nzm )
         invrs_rho_ds_zt_flip(i,:) = flip( invrs_rho_ds_zt(i,:), nzt )
         thv_ds_zt_flip(i,:) = flip( thv_ds_zt(i,:), nzt )
         thv_ds_zm_flip(i,:) = flip( thv_ds_zm(i,:), nzm )
         rfrzm_flip(i,:) = flip( rfrzm(i,:), nzt )

         um_flip(i,:) = flip( um(i,:), nzt )
         vm_flip(i,:) = flip( vm(i,:), nzt )
         upwp_flip(i,:) = flip( upwp(i,:), nzm )
         vpwp_flip(i,:) = flip( vpwp(i,:), nzm )
         up2_flip(i,:) = flip( up2(i,:), nzm )
         vp2_flip(i,:) = flip( vp2(i,:), nzm )
         up3_flip(i,:) = flip( up3(i,:), nzt )
         vp3_flip(i,:) = flip( vp3(i,:), nzt )
         thlm_flip(i,:) = flip( thlm(i,:), nzt )
         rtm_flip(i,:) = flip( rtm(i,:), nzt )
         wpthlp_flip(i,:) = flip( wpthlp(i,:), nzm )
         wprtp_flip(i,:) = flip( wprtp(i,:), nzm )
         wp2_flip(i,:) = flip( wp2(i,:), nzm )
         wp3_flip(i,:) = flip( wp3(i,:), nzt )
         rtp2_flip(i,:) = flip( rtp2(i,:), nzm )
         rtp3_flip(i,:) = flip( rtp3(i,:), nzt )
         thlp2_flip(i,:) = flip( thlp2(i,:), nzm )
         thlp3_flip(i,:) = flip( thlp3(i,:), nzt )
         rtpthlp_flip(i,:) = flip( rtpthlp(i,:), nzm )
         rcm_flip(i,:) = flip( rcm(i,:), nzt )
         cloud_frac_flip(i,:) = flip( cloud_frac(i,:), nzt )
         wpthvp_flip(i,:) = flip( wpthvp(i,:), nzm )
         wp2thvp_flip(i,:) = flip( wp2thvp(i,:), nzt )
         rtpthvp_flip(i,:) = flip( rtpthvp(i,:), nzm )
         thlpthvp_flip(i,:) = flip( thlpthvp(i,:), nzm )
         wp2rtp_flip(i,:) = flip( wp2rtp(i,:), nzt )
         wp2thlp_flip(i,:) = flip( wp2thlp(i,:), nzt )
         uprcp_flip(i,:) = flip( uprcp(i,:), nzm )
         vprcp_flip(i,:) = flip( vprcp(i,:), nzm )
         rc_coef_zm_flip(i,:) = flip( rc_coef_zm(i,:), nzm )
         wp4_flip(i,:) = flip( wp4(i,:), nzm )
         wpup2_flip(i,:) = flip( wpup2(i,:), nzt )
         wpvp2_flip(i,:) = flip( wpvp2(i,:), nzt )
         wp2up2_flip(i,:) = flip( wp2up2(i,:), nzm )
         wp2vp2_flip(i,:) = flip( wp2vp2(i,:), nzm )
         ice_supersat_frac_flip(i,:) = flip( ice_supersat_frac(i,:), nzt )

         ! sclr variables
         do sclr = 1, sclr_dim
            sclrm_forcing_flip(i,:,sclr) = flip( sclrm_forcing(i,:,sclr), nzt )
            sclrm_flip(i,:,sclr) = flip( sclrm(i,:,sclr), nzt )
            wpsclrp_flip(i,:,sclr) = flip( wpsclrp(i,:,sclr), nzm )
            sclrprtp_flip(i,:,sclr) = flip( sclrprtp(i,:,sclr), nzm )
            sclrpthlp_flip(i,:,sclr) = flip( sclrpthlp(i,:,sclr), nzm )
            sclrpthvp_flip(i,:,sclr) = flip( sclrpthvp(i,:,sclr), nzm )
            sclrp2_flip(i,:,sclr) = flip( sclrp2(i,:,sclr), nzm )
            sclrp3_flip(i,:,sclr) = flip( sclrp3(i,:,sclr), nzt )
         enddo

         ! edsclr variables
         do edsclr = 1, edsclr_dim
            edsclrm_forcing_flip(i,:,edsclr) = flip( edsclrm_forcing(i,:,edsclr), nzt )
            edsclrm_flip(i,:,edsclr) = flip( edsclrm(i,:,edsclr), nzt )
         enddo

         ! hydrometeor variables
         do hm_idx = 1, hydromet_dim
            wphydrometp_flip(i,:,hm_idx) = flip( wphydrometp(i,:,hm_idx), nzm )
            wp2hmp_flip(i,:,hm_idx) = flip( wp2hmp(i,:,hm_idx), nzt )
            rtphmp_zt_flip(i,:,hm_idx) = flip( rtphmp_zt(i,:,hm_idx), nzt )
            thlphmp_zt_flip(i,:,hm_idx) = flip( thlphmp_zt(i,:,hm_idx), nzt )
         enddo ! hm_idx = 1, hydromet_dim

         ! Allocate space and initialize flipped pdf parameter terms
         call init_pdf_params( gr%nzt, ngrdcol, pdf_params_flip )
         call init_pdf_params( gr%nzm, ngrdcol, pdf_params_zm_flip )

         call init_pdf_implicit_coefs_terms( gr%nzt, ngrdcol, sclr_dim, &   ! Intent(in)
                                             pdf_implicit_coefs_terms_flip ) ! Intent(out)

         ! pdf_params
         pdf_params_flip%w_1(i,:) = flip( pdf_params%w_1(i,:), nzt )
         pdf_params_flip%w_2(i,:) = flip( pdf_params%w_2(i,:), nzt )
         pdf_params_flip%varnce_w_1(i,:) = flip( pdf_params%varnce_w_1(i,:), nzt )
         pdf_params_flip%varnce_w_2(i,:) = flip( pdf_params%varnce_w_2(i,:), nzt )
         pdf_params_flip%rt_1(i,:) = flip( pdf_params%rt_1(i,:), nzt )
         pdf_params_flip%rt_2(i,:) = flip( pdf_params%rt_2(i,:), nzt )
         pdf_params_flip%varnce_rt_1(i,:) = flip( pdf_params%varnce_rt_1(i,:), nzt )
         pdf_params_flip%varnce_rt_2(i,:) = flip( pdf_params%varnce_rt_2(i,:), nzt )
         pdf_params_flip%thl_1(i,:) = flip( pdf_params%thl_1(i,:), nzt )
         pdf_params_flip%thl_2(i,:) = flip( pdf_params%thl_2(i,:), nzt )
         pdf_params_flip%varnce_thl_1(i,:) = flip( pdf_params%varnce_thl_1(i,:), nzt )
         pdf_params_flip%varnce_thl_2(i,:) = flip( pdf_params%varnce_thl_2(i,:), nzt )
         pdf_params_flip%corr_w_rt_1(i,:) = flip( pdf_params%corr_w_rt_1(i,:), nzt )
         pdf_params_flip%corr_w_rt_2(i,:) = flip( pdf_params%corr_w_rt_2(i,:), nzt )
         pdf_params_flip%corr_w_thl_1(i,:) = flip( pdf_params%corr_w_thl_1(i,:), nzt )
         pdf_params_flip%corr_w_thl_2(i,:) = flip( pdf_params%corr_w_thl_2(i,:), nzt )
         pdf_params_flip%corr_rt_thl_1(i,:) = flip( pdf_params%corr_rt_thl_1(i,:), nzt )
         pdf_params_flip%corr_rt_thl_2(i,:) = flip( pdf_params%corr_rt_thl_2(i,:), nzt )
         pdf_params_flip%alpha_thl(i,:) = flip( pdf_params%alpha_thl(i,:), nzt )
         pdf_params_flip%alpha_rt(i,:) = flip( pdf_params%alpha_rt(i,:), nzt )
         pdf_params_flip%crt_1(i,:) = flip( pdf_params%crt_1(i,:), nzt )
         pdf_params_flip%crt_2(i,:) = flip( pdf_params%crt_2(i,:), nzt )
         pdf_params_flip%cthl_1(i,:) = flip( pdf_params%cthl_1(i,:), nzt )
         pdf_params_flip%cthl_2(i,:) = flip( pdf_params%cthl_2(i,:), nzt )
         pdf_params_flip%chi_1(i,:) = flip( pdf_params%chi_1(i,:), nzt )
         pdf_params_flip%chi_2(i,:) = flip( pdf_params%chi_2(i,:), nzt )
         pdf_params_flip%stdev_chi_1(i,:) = flip( pdf_params%stdev_chi_1(i,:), nzt )
         pdf_params_flip%stdev_chi_2(i,:) = flip( pdf_params%stdev_chi_2(i,:), nzt )
         pdf_params_flip%stdev_eta_1(i,:) = flip( pdf_params%stdev_eta_1(i,:), nzt )
         pdf_params_flip%stdev_eta_2(i,:) = flip( pdf_params%stdev_eta_2(i,:), nzt )
         pdf_params_flip%covar_chi_eta_1(i,:) = flip( pdf_params%covar_chi_eta_1(i,:), nzt )
         pdf_params_flip%covar_chi_eta_2(i,:) = flip( pdf_params%covar_chi_eta_2(i,:), nzt )
         pdf_params_flip%corr_w_chi_1(i,:) = flip( pdf_params%corr_w_chi_1(i,:), nzt )
         pdf_params_flip%corr_w_chi_2(i,:) = flip( pdf_params%corr_w_chi_2(i,:), nzt )
         pdf_params_flip%corr_w_eta_1(i,:) = flip( pdf_params%corr_w_eta_1(i,:), nzt )
         pdf_params_flip%corr_w_eta_2(i,:) = flip( pdf_params%corr_w_eta_2(i,:), nzt )
         pdf_params_flip%corr_chi_eta_1(i,:) = flip( pdf_params%corr_chi_eta_1(i,:), nzt )
         pdf_params_flip%corr_chi_eta_2(i,:) = flip( pdf_params%corr_chi_eta_2(i,:), nzt )
         pdf_params_flip%rsatl_1(i,:) = flip( pdf_params%rsatl_1(i,:), nzt )
         pdf_params_flip%rsatl_2(i,:) = flip( pdf_params%rsatl_2(i,:), nzt )
         pdf_params_flip%rc_1(i,:) = flip( pdf_params%rc_1(i,:), nzt )
         pdf_params_flip%rc_2(i,:) = flip( pdf_params%rc_2(i,:), nzt )
         pdf_params_flip%cloud_frac_1(i,:) = flip( pdf_params%cloud_frac_1(i,:), nzt )
         pdf_params_flip%cloud_frac_2(i,:) = flip( pdf_params%cloud_frac_2(i,:), nzt )
         pdf_params_flip%mixt_frac(i,:) = flip( pdf_params%mixt_frac(i,:), nzt )
         pdf_params_flip%ice_supersat_frac_1(i,:) = flip( pdf_params%ice_supersat_frac_1(i,:), nzt )
         pdf_params_flip%ice_supersat_frac_2(i,:) = flip( pdf_params%ice_supersat_frac_2(i,:), nzt )

         ! pdf_params_zm
         pdf_params_zm_flip%w_1(i,:) = flip( pdf_params_zm%w_1(i,:), nzm )
         pdf_params_zm_flip%w_2(i,:) = flip( pdf_params_zm%w_2(i,:), nzm )
         pdf_params_zm_flip%varnce_w_1(i,:) = flip( pdf_params_zm%varnce_w_1(i,:), nzm )
         pdf_params_zm_flip%varnce_w_2(i,:) = flip( pdf_params_zm%varnce_w_2(i,:), nzm )
         pdf_params_zm_flip%rt_1(i,:) = flip( pdf_params_zm%rt_1(i,:), nzm )
         pdf_params_zm_flip%rt_2(i,:) = flip( pdf_params_zm%rt_2(i,:), nzm )
         pdf_params_zm_flip%varnce_rt_1(i,:) = flip( pdf_params_zm%varnce_rt_1(i,:), nzm )
         pdf_params_zm_flip%varnce_rt_2(i,:) = flip( pdf_params_zm%varnce_rt_2(i,:), nzm )
         pdf_params_zm_flip%thl_1(i,:) = flip( pdf_params_zm%thl_1(i,:), nzm )
         pdf_params_zm_flip%thl_2(i,:) = flip( pdf_params_zm%thl_2(i,:), nzm )
         pdf_params_zm_flip%varnce_thl_1(i,:) = flip( pdf_params_zm%varnce_thl_1(i,:), nzm )
         pdf_params_zm_flip%varnce_thl_2(i,:) = flip( pdf_params_zm%varnce_thl_2(i,:), nzm )
         pdf_params_zm_flip%corr_w_rt_1(i,:) = flip( pdf_params_zm%corr_w_rt_1(i,:), nzm )
         pdf_params_zm_flip%corr_w_rt_2(i,:) = flip( pdf_params_zm%corr_w_rt_2(i,:), nzm )
         pdf_params_zm_flip%corr_w_thl_1(i,:) = flip( pdf_params_zm%corr_w_thl_1(i,:), nzm )
         pdf_params_zm_flip%corr_w_thl_2(i,:) = flip( pdf_params_zm%corr_w_thl_2(i,:), nzm )
         pdf_params_zm_flip%corr_rt_thl_1(i,:) = flip( pdf_params_zm%corr_rt_thl_1(i,:), nzm )
         pdf_params_zm_flip%corr_rt_thl_2(i,:) = flip( pdf_params_zm%corr_rt_thl_2(i,:), nzm )
         pdf_params_zm_flip%alpha_thl(i,:) = flip( pdf_params_zm%alpha_thl(i,:), nzm )
         pdf_params_zm_flip%alpha_rt(i,:) = flip( pdf_params_zm%alpha_rt(i,:), nzm )
         pdf_params_zm_flip%crt_1(i,:) = flip( pdf_params_zm%crt_1(i,:), nzm )
         pdf_params_zm_flip%crt_2(i,:) = flip( pdf_params_zm%crt_2(i,:), nzm )
         pdf_params_zm_flip%cthl_1(i,:) = flip( pdf_params_zm%cthl_1(i,:), nzm )
         pdf_params_zm_flip%cthl_2(i,:) = flip( pdf_params_zm%cthl_2(i,:), nzm )
         pdf_params_zm_flip%chi_1(i,:) = flip( pdf_params_zm%chi_1(i,:), nzm )
         pdf_params_zm_flip%chi_2(i,:) = flip( pdf_params_zm%chi_2(i,:), nzm )
         pdf_params_zm_flip%stdev_chi_1(i,:) = flip( pdf_params_zm%stdev_chi_1(i,:), nzm )
         pdf_params_zm_flip%stdev_chi_2(i,:) = flip( pdf_params_zm%stdev_chi_2(i,:), nzm )
         pdf_params_zm_flip%stdev_eta_1(i,:) = flip( pdf_params_zm%stdev_eta_1(i,:), nzm )
         pdf_params_zm_flip%stdev_eta_2(i,:) = flip( pdf_params_zm%stdev_eta_2(i,:), nzm )
         pdf_params_zm_flip%covar_chi_eta_1(i,:) = flip( pdf_params_zm%covar_chi_eta_1(i,:), nzm )
         pdf_params_zm_flip%covar_chi_eta_2(i,:) = flip( pdf_params_zm%covar_chi_eta_2(i,:), nzm )
         pdf_params_zm_flip%corr_w_chi_1(i,:) = flip( pdf_params_zm%corr_w_chi_1(i,:), nzm )
         pdf_params_zm_flip%corr_w_chi_2(i,:) = flip( pdf_params_zm%corr_w_chi_2(i,:), nzm )
         pdf_params_zm_flip%corr_w_eta_1(i,:) = flip( pdf_params_zm%corr_w_eta_1(i,:), nzm )
         pdf_params_zm_flip%corr_w_eta_2(i,:) = flip( pdf_params_zm%corr_w_eta_2(i,:), nzm )
         pdf_params_zm_flip%corr_chi_eta_1(i,:) = flip( pdf_params_zm%corr_chi_eta_1(i,:), nzm )
         pdf_params_zm_flip%corr_chi_eta_2(i,:) = flip( pdf_params_zm%corr_chi_eta_2(i,:), nzm )
         pdf_params_zm_flip%rsatl_1(i,:) = flip( pdf_params_zm%rsatl_1(i,:), nzm )
         pdf_params_zm_flip%rsatl_2(i,:) = flip( pdf_params_zm%rsatl_2(i,:), nzm )
         pdf_params_zm_flip%rc_1(i,:) = flip( pdf_params_zm%rc_1(i,:), nzm )
         pdf_params_zm_flip%rc_2(i,:) = flip( pdf_params_zm%rc_2(i,:), nzm )
         pdf_params_zm_flip%cloud_frac_1(i,:) = flip( pdf_params_zm%cloud_frac_1(i,:), nzm )
         pdf_params_zm_flip%cloud_frac_2(i,:) = flip( pdf_params_zm%cloud_frac_2(i,:), nzm )
         pdf_params_zm_flip%mixt_frac(i,:) = flip( pdf_params_zm%mixt_frac(i,:), nzm )
         pdf_params_zm_flip%ice_supersat_frac_1(i,:) = flip( pdf_params_zm%ice_supersat_frac_1(i,:), nzm )
         pdf_params_zm_flip%ice_supersat_frac_2(i,:) = flip( pdf_params_zm%ice_supersat_frac_2(i,:), nzm )

         ! pdf_implicit_coefs_terms
         pdf_implicit_coefs_terms_flip%coef_wp4_implicit(i,:) &
            = flip( pdf_implicit_coefs_terms%coef_wp4_implicit(i,:), nzt )
         pdf_implicit_coefs_terms_flip%coef_wp2rtp_implicit(i,:) & 
            = flip( pdf_implicit_coefs_terms%coef_wp2rtp_implicit(i,:), nzt )
         pdf_implicit_coefs_terms_flip%term_wp2rtp_explicit(i,:) &
            = flip( pdf_implicit_coefs_terms%term_wp2rtp_explicit(i,:), nzt )
         pdf_implicit_coefs_terms_flip%coef_wp2thlp_implicit(i,:) &
            = flip( pdf_implicit_coefs_terms%coef_wp2thlp_implicit(i,:), nzt )
         pdf_implicit_coefs_terms_flip%term_wp2thlp_explicit(i,:) &
            = flip( pdf_implicit_coefs_terms%term_wp2thlp_explicit(i,:), nzt )
         pdf_implicit_coefs_terms_flip%coef_wp2up_implicit(i,:) &
            = flip( pdf_implicit_coefs_terms%coef_wp2up_implicit(i,:), nzt )
         pdf_implicit_coefs_terms_flip%term_wp2up_explicit(i,:) &
            = flip( pdf_implicit_coefs_terms%term_wp2up_explicit(i,:), nzt )
         pdf_implicit_coefs_terms_flip%coef_wp2vp_implicit(i,:) &
            = flip( pdf_implicit_coefs_terms%coef_wp2vp_implicit(i,:), nzt )
         pdf_implicit_coefs_terms_flip%term_wp2vp_explicit(i,:) &
            = flip( pdf_implicit_coefs_terms%term_wp2vp_explicit(i,:), nzt )
         pdf_implicit_coefs_terms_flip%coef_wprtp2_implicit(i,:) &
            = flip( pdf_implicit_coefs_terms%coef_wprtp2_implicit(i,:), nzt )
         pdf_implicit_coefs_terms_flip%term_wprtp2_explicit(i,:) &
            = flip( pdf_implicit_coefs_terms%term_wprtp2_explicit(i,:), nzt )
         pdf_implicit_coefs_terms_flip%coef_wpthlp2_implicit(i,:) &
            = flip( pdf_implicit_coefs_terms%coef_wpthlp2_implicit(i,:), nzt )
         pdf_implicit_coefs_terms_flip%term_wpthlp2_explicit(i,:) &
            = flip( pdf_implicit_coefs_terms%term_wpthlp2_explicit(i,:), nzt )
         pdf_implicit_coefs_terms_flip%coef_wprtpthlp_implicit(i,:) &
            = flip( pdf_implicit_coefs_terms%coef_wprtpthlp_implicit(i,:), nzt )
         pdf_implicit_coefs_terms_flip%term_wprtpthlp_explicit(i,:) &
            = flip( pdf_implicit_coefs_terms%term_wprtpthlp_explicit(i,:), nzt )
         pdf_implicit_coefs_terms_flip%coef_wpup2_implicit(i,:) &
            = flip( pdf_implicit_coefs_terms%coef_wpup2_implicit(i,:), nzt )
         pdf_implicit_coefs_terms_flip%term_wpup2_explicit(i,:) &
            = flip( pdf_implicit_coefs_terms%term_wpup2_explicit(i,:), nzt )
         pdf_implicit_coefs_terms_flip%coef_wpvp2_implicit(i,:) &
            = flip( pdf_implicit_coefs_terms%coef_wpvp2_implicit(i,:), nzt )
         pdf_implicit_coefs_terms_flip%term_wpvp2_explicit(i,:) &
            = flip( pdf_implicit_coefs_terms%term_wpvp2_explicit(i,:), nzt )
         if ( sclr_dim > 0 ) then
            do sclr = 1, sclr_dim
               pdf_implicit_coefs_terms_flip%coef_wp2sclrp_implicit(i,:,sclr) &
                  = flip( pdf_implicit_coefs_terms%coef_wp2sclrp_implicit(i,:,sclr), nzt )
               pdf_implicit_coefs_terms_flip%term_wp2sclrp_explicit(i,:,sclr) &
                  = flip( pdf_implicit_coefs_terms%term_wp2sclrp_explicit(i,:,sclr), nzt )
               pdf_implicit_coefs_terms_flip%coef_wpsclrp2_implicit(i,:,sclr) &
                  = flip( pdf_implicit_coefs_terms%coef_wpsclrp2_implicit(i,:,sclr), nzt )
               pdf_implicit_coefs_terms_flip%term_wpsclrp2_explicit(i,:,sclr) &
                  = flip( pdf_implicit_coefs_terms%term_wpsclrp2_explicit(i,:,sclr), nzt )
               pdf_implicit_coefs_terms_flip%coef_wprtpsclrp_implicit(i,:,sclr) &
                  = flip( pdf_implicit_coefs_terms%coef_wprtpsclrp_implicit(i,:,sclr), nzt )
               pdf_implicit_coefs_terms_flip%term_wprtpsclrp_explicit(i,:,sclr) &
                  = flip( pdf_implicit_coefs_terms%term_wprtpsclrp_explicit(i,:,sclr), nzt )
               pdf_implicit_coefs_terms_flip%coef_wpthlpsclrp_implicit(i,:,sclr) &
                  = flip( pdf_implicit_coefs_terms%coef_wpthlpsclrp_implicit(i,:,sclr), nzt )
               pdf_implicit_coefs_terms_flip%term_wpthlpsclrp_explicit(i,:,sclr) &
                  = flip( pdf_implicit_coefs_terms%term_wpthlpsclrp_explicit(i,:,sclr), nzt )
            enddo
         endif ! sclr_dim > 0

         Kh_zm_flip(i,:) = flip( Kh_zm(i,:), nzm )
         Kh_zt_flip(i,:) = flip( Kh_zt(i,:), nzt )
         thlprcp_flip(i,:) = flip( thlprcp(i,:), nzm )
         wprcp_flip(i,:) = flip( wprcp(i,:), nzm )
         w_up_in_cloud_flip(i,:) = flip( w_up_in_cloud(i,:), nzt )
         w_down_in_cloud_flip(i,:) = flip( w_down_in_cloud(i,:), nzt )
         cloudy_updraft_frac_flip(i,:) = flip( cloudy_updraft_frac(i,:), nzt )
         cloudy_downdraft_frac_flip(i,:) = flip( cloudy_downdraft_frac(i,:), nzt )
         rcm_in_layer_flip(i,:) = flip( rcm_in_layer(i,:), nzt )
         cloud_cover_flip(i,:) = flip( cloud_cover(i,:), nzt )
         invrs_tau_zm_flip(i,:) = flip( invrs_tau_zm(i,:), nzm )
         Lscale_flip(i,:) = flip( Lscale(i,:), nzt )

      enddo ! i = 1, ngrdcol

      ! Set statistical sampling to false for descending grid to avoid
      ! the error from having too many statistical samples.
      stats_metadata_flip%l_stats_samp = .false.


      ! Call advance_clubb_core_api for the ascending grid direction
      call advance_clubb_core_api( &
              gr, nzm, nzt, ngrdcol, &                             ! Intent(in)
              l_implemented, dt, fcor, sfc_elevation, &            ! Intent(in)
              hydromet_dim, &                                      ! intent(in)
              sclr_dim, sclr_tol, edsclr_dim, sclr_idx, &          ! intent(in)
              thlm_forcing, rtm_forcing, um_forcing, vm_forcing, & ! Intent(in)
              sclrm_forcing, edsclrm_forcing, wprtp_forcing, &     ! Intent(in)
              wpthlp_forcing, rtp2_forcing, thlp2_forcing, &       ! Intent(in)
              rtpthlp_forcing, wm_zm, wm_zt, &                     ! Intent(in)
              wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, p_sfc, &  ! Intent(in)
              wpsclrp_sfc, wpedsclrp_sfc,  &                       ! Intent(in)
              upwp_sfc_pert, vpwp_sfc_pert, &                      ! intent(in)
              rtm_ref, thlm_ref, um_ref, vm_ref, ug, vg, &         ! Intent(in)
              p_in_Pa, rho_zm, rho, exner, &                       ! Intent(in)
              rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &             ! Intent(in)
              invrs_rho_ds_zt, thv_ds_zm, thv_ds_zt, &             ! Intent(in) 
              l_mix_rat_hm, &                                      ! Intent(in)
              rfrzm, wphydrometp, &                                ! Intent(in)
              wp2hmp, rtphmp_zt, thlphmp_zt, &                     ! Intent(in)
              host_dx, host_dy, &                                  ! Intent(in)
              clubb_params, nu_vert_res_dep, lmin, &               ! Intent(in)
              clubb_config_flags, &                                ! Intent(in)
              stats_metadata, &                                    ! Intent(in)
              stats_zt, stats_zm, stats_sfc, &                     ! intent(inout)
              um, vm, upwp, vpwp, up2, vp2, up3, vp3, &            ! Intent(inout)
              thlm, rtm, wprtp, wpthlp, &                          ! Intent(inout)
              wp2, wp3, rtp2, rtp3, thlp2, thlp3, rtpthlp, &       ! Intent(inout)
              sclrm, sclrp2, sclrp3, sclrprtp, sclrpthlp, &        ! Intent(inout)
              wpsclrp, edsclrm, err_code, &                        ! Intent(inout)
              rcm, cloud_frac, &                                   ! Intent(inout)
              wpthvp, wp2thvp, rtpthvp, thlpthvp, &                ! Intent(inout)
              sclrpthvp, &                                         ! Intent(inout)
              wp2rtp, wp2thlp, uprcp, vprcp, rc_coef_zm, wp4, &    ! intent(inout)
              wpup2, wpvp2, wp2up2, wp2vp2, ice_supersat_frac, &   ! intent(inout)
              um_pert, vm_pert, upwp_pert, vpwp_pert, &            ! intent(inout)
              pdf_params, pdf_params_zm, &                         ! Intent(inout)
              pdf_implicit_coefs_terms, &                          ! intent(inout)
              Kh_zm, Kh_zt, &                                      ! intent(out)
              thlprcp, wprcp, w_up_in_cloud, w_down_in_cloud, &    ! Intent(out)
              cloudy_updraft_frac, cloudy_downdraft_frac, &        ! Intent(out)
              rcm_in_layer, cloud_cover, invrs_tau_zm, &           ! Intent(out)
              Lscale )                                             ! Intent(out)


      ! In the case of a fatal error during the 1st call (ascending grid) to 
      ! advance_clubb_core_api, reset the error code.
      !
      ! Otherwise, a mismatch between the ascending grid and the descending grid
      ! will be caused by the fact that the 2nd call to advance_clubb_core_api
      ! will be returned from right away because err_code is already set to 
      ! clubb_fatal_error upon entering advance_clubb_core_api. Thus, there
      ! likely won't be as many variables calculated during the descending call
      ! as there were during the ascending call. This causes a mismatch and a
      ! false failure of the generalized grid test.
      !
      ! In the case of a fatal error, as long as the generalized grid code is
      ! working properly, the fatal error will occur in the exact same spot
      ! in the call to the descending grid. The fatal error code will be output
      ! (and not overwritten this time) and that will cause the code to exit
      ! the run and stop. However, the results between the ascending grid and
      ! the descending grid will still match.
      if ( err_code == clubb_fatal_error ) then
         ! Reset error code
         err_code = clubb_no_error
      endif


      ! Call advance_clubb_core_api for the descending grid direction
      ! All variables with a vertical dimension should be "flip" variables
      ! in this call.
      call advance_clubb_core_api( &
              gr_desc, nzm, nzt, ngrdcol, &                                         ! Intent(in)
              l_implemented, dt, fcor, sfc_elevation, &                             ! Intent(in)
              hydromet_dim, &                                                       ! intent(in)
              sclr_dim, sclr_tol, edsclr_dim, sclr_idx, &                           ! intent(in)
              thlm_forcing_flip, rtm_forcing_flip, um_forcing_flip, vm_forcing_flip, & ! Intent(in)
              sclrm_forcing_flip, edsclrm_forcing_flip, wprtp_forcing_flip, &       ! Intent(in)
              wpthlp_forcing_flip, rtp2_forcing_flip, thlp2_forcing_flip, &         ! Intent(in)
              rtpthlp_forcing_flip, wm_zm_flip, wm_zt_flip, &                       ! Intent(in)
              wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, p_sfc, &                   ! Intent(in)
              wpsclrp_sfc, wpedsclrp_sfc,  &                                        ! Intent(in)
              upwp_sfc_pert, vpwp_sfc_pert, &                                       ! intent(in)
              rtm_ref_flip, thlm_ref_flip, um_ref_flip, vm_ref_flip, ug_flip, vg_flip, & ! Intent(in)
              p_in_Pa_flip, rho_zm_flip, rho_flip, exner_flip, &                    ! Intent(in)
              rho_ds_zm_flip, rho_ds_zt_flip, invrs_rho_ds_zm_flip, &               ! Intent(in)
              invrs_rho_ds_zt_flip, thv_ds_zm_flip, thv_ds_zt_flip, &               ! Intent(in) 
              l_mix_rat_hm, &                                                       ! Intent(in)
              rfrzm_flip, wphydrometp_flip, &                                       ! Intent(in)
              wp2hmp_flip, rtphmp_zt_flip, thlphmp_zt_flip, &                       ! Intent(in)
              host_dx, host_dy, &                                                   ! Intent(in)
              clubb_params, nu_vert_res_dep, lmin, &                                ! Intent(in)
              clubb_config_flags, &                                                 ! Intent(in)
              stats_metadata_flip, &                                                ! Intent(in)
              stats_zt, stats_zm, stats_sfc, &                                      ! intent(inout)
              um_flip, vm_flip, upwp_flip, vpwp_flip, up2_flip, vp2_flip, up3_flip, vp3_flip, & ! Intent(inout)
              thlm_flip, rtm_flip, wprtp_flip, wpthlp_flip, &                       ! Intent(inout)
              wp2_flip, wp3_flip, rtp2_flip, rtp3_flip, thlp2_flip, thlp3_flip, rtpthlp_flip, & ! Intent(inout)
              sclrm_flip, sclrp2_flip, sclrp3_flip, sclrprtp_flip, sclrpthlp_flip, & ! Intent(inout)
              wpsclrp_flip, edsclrm_flip, err_code, &                               ! Intent(inout)
              rcm_flip, cloud_frac_flip, &                                          ! Intent(inout)
              wpthvp_flip, wp2thvp_flip, rtpthvp_flip, thlpthvp_flip, &             ! Intent(inout)
              sclrpthvp_flip, &                                                     ! Intent(inout)
              wp2rtp_flip, wp2thlp_flip, uprcp_flip, vprcp_flip, rc_coef_zm_flip, wp4_flip, & ! intent(inout)
              wpup2_flip, wpvp2_flip, wp2up2_flip, wp2vp2_flip, ice_supersat_frac_flip, & ! intent(inout)
              um_pert, vm_pert, upwp_pert, vpwp_pert, &                             ! intent(inout)
              pdf_params_flip, pdf_params_zm_flip, &                                ! Intent(inout)
              pdf_implicit_coefs_terms_flip, &                                      ! intent(inout)
              Kh_zm_flip, Kh_zt_flip, &                                             ! intent(out)
              thlprcp_flip, wprcp_flip, w_up_in_cloud_flip, w_down_in_cloud_flip, & ! Intent(out)
              cloudy_updraft_frac_flip, cloudy_downdraft_frac_flip, &               ! Intent(out)
              rcm_in_layer_flip, cloud_cover_flip, invrs_tau_zm_flip, &             ! Intent(out)
              Lscale_flip )                                                         ! Intent(out)


      ! Compare the ascending grid variables to the descending grid variables
      ! rtm
      call check_flipped_results( "rtm", rtm, rtm_flip, nzt, ngrdcol, &
                                  l_differences )
      ! wprtp
      call check_flipped_results( "wprtp", wprtp, wprtp_flip, nzm, ngrdcol, &
                                  l_differences )
      ! thlm
      call check_flipped_results( "thlm", thlm, thlm_flip, nzt, ngrdcol, &
                                  l_differences )
      ! wpthlp
      call check_flipped_results( "wpthlp", wpthlp, wpthlp_flip, nzm, ngrdcol, &
                                  l_differences )
      ! wp2
      call check_flipped_results( "wp2", wp2, wp2_flip, nzm, ngrdcol, &
                                  l_differences )
      ! wp3
      call check_flipped_results( "wp3", wp3, wp3_flip, nzt, ngrdcol, &
                                  l_differences )
      ! rtp2
      call check_flipped_results( "rtp2", rtp2, rtp2_flip, nzm, ngrdcol, &
                                  l_differences )
      ! thlp2
      call check_flipped_results( "thlp2", thlp2, thlp2_flip, nzm, ngrdcol, &
                                  l_differences )
      ! rtpthlp
      call check_flipped_results( "rtpthlp", rtpthlp, rtpthlp_flip, nzm, ngrdcol, &
                                  l_differences )
      ! rtp3
      call check_flipped_results( "rtp3", rtp3, rtp3_flip, nzt, ngrdcol, &
                                  l_differences )
      ! thlp3
      call check_flipped_results( "thlp3", thlp3, thlp3_flip, nzt, ngrdcol, &
                                  l_differences )
      ! um
      call check_flipped_results( "um", um, um_flip, nzt, ngrdcol, &
                                  l_differences )
      ! vm
      call check_flipped_results( "vm", vm, vm_flip, nzt, ngrdcol, &
                                  l_differences )
      ! upwp
      call check_flipped_results( "upwp", upwp, upwp_flip, nzm, ngrdcol, &
                                  l_differences )
      ! vpwp
      call check_flipped_results( "vpwp", vpwp, vpwp_flip, nzm, ngrdcol, &
                                  l_differences )
      ! up2
      call check_flipped_results( "up2", up2, up2_flip, nzm, ngrdcol, &
                                  l_differences )
      ! vp2
      call check_flipped_results( "vp2", vp2, vp2_flip, nzm, ngrdcol, &
                                  l_differences )
      ! up3
      call check_flipped_results( "up3", up3, up3_flip, nzt, ngrdcol, &
                                  l_differences )
      ! vp3
      call check_flipped_results( "vp3", vp3, vp3_flip, nzt, ngrdcol, &
                                  l_differences )
      ! rcm
      call check_flipped_results( "rcm", rcm, rcm_flip, nzt, ngrdcol, &
                                  l_differences )
      ! cloud_frac
      call check_flipped_results( "cloud_frac", cloud_frac, cloud_frac_flip, nzt, ngrdcol, &
                                  l_differences )
      ! wpthvp
      call check_flipped_results( "wpthvp", wpthvp, wpthvp_flip, nzm, ngrdcol, &
                                  l_differences )
      ! wp2thvp
      call check_flipped_results( "wp2thvp", wp2thvp, wp2thvp_flip, nzt, ngrdcol, &
                                  l_differences )
      ! rtpthvp
      call check_flipped_results( "rtpthvp", rtpthvp, rtpthvp_flip, nzm, ngrdcol, &
                                  l_differences )
      ! thlpthvp
      call check_flipped_results( "thlpthvp", thlpthvp, thlpthvp_flip, nzm, ngrdcol, &
                                  l_differences )
      ! wp2rtp
      call check_flipped_results( "wp2rtp", wp2rtp, wp2rtp_flip, nzt, ngrdcol, &
                                  l_differences )
      ! wp2thlp
      call check_flipped_results( "wp2thlp", wp2thlp, wp2thlp_flip, nzt, ngrdcol, &
                                  l_differences )
      ! uprcp
      call check_flipped_results( "uprcp", uprcp, uprcp_flip, nzm, ngrdcol, &
                                  l_differences )
      ! vprcp
      call check_flipped_results( "vprcp", vprcp, vprcp_flip, nzm, ngrdcol, &
                                  l_differences )
      ! rc_coef_zm
      call check_flipped_results( "rc_coef_zm", rc_coef_zm, rc_coef_zm_flip, nzm, ngrdcol, &
                                  l_differences )
      ! wp4
      call check_flipped_results( "wp4", wp4, wp4_flip, nzm, ngrdcol, &
                                  l_differences )
      ! wpup2
      call check_flipped_results( "wpup2", wpup2, wpup2_flip, nzt, ngrdcol, &
                                  l_differences )
      ! wpvp2
      call check_flipped_results( "wpvp2", wpvp2, wpvp2_flip, nzt, ngrdcol, &
                                  l_differences )
      ! wp2up2
      call check_flipped_results( "wp2up2", wp2up2, wp2up2_flip, nzm, ngrdcol, &
                                  l_differences )
      ! wp2vp2
      call check_flipped_results( "wp2vp2", wp2vp2, wp2vp2_flip, nzm, ngrdcol, &
                                  l_differences )
      ! ice_supersat_frac
      call check_flipped_results( "ice_supersat_frac", ice_supersat_frac, &
                                  ice_supersat_frac_flip, nzt, ngrdcol, &
                                  l_differences )
      ! Kh_zm
      call check_flipped_results( "Kh_zm", Kh_zm, Kh_zm_flip, nzm, ngrdcol, &
                                  l_differences )
      ! Kh_zt
      call check_flipped_results( "Kh_zt", Kh_zt, Kh_zt_flip, nzt, ngrdcol, &
                                  l_differences )
      ! thlprcp
      call check_flipped_results( "thlprcp", thlprcp, thlprcp_flip, nzm, ngrdcol, &
                                  l_differences )
      ! wprcp
      call check_flipped_results( "wprcp", wprcp, wprcp_flip, nzm, ngrdcol, &
                                  l_differences )
      ! w_up_in_cloud
      call check_flipped_results( "w_up_in_cloud", w_up_in_cloud, &
                                  w_up_in_cloud_flip, nzt, ngrdcol, &
                                  l_differences )
      ! w_down_in_cloud
      call check_flipped_results( "w_down_in_cloud", w_down_in_cloud, &
                                  w_down_in_cloud_flip, nzt, ngrdcol, &
                                  l_differences )
      ! cloudy_updraft_frac
      call check_flipped_results( "cloudy_updraft_frac", cloudy_updraft_frac, &
                                  cloudy_updraft_frac_flip, nzt, ngrdcol, &
                                  l_differences )
      ! cloudy_downdraft_frac
      call check_flipped_results( "cloudy_downdraft_frac", cloudy_downdraft_frac, &
                                  cloudy_downdraft_frac_flip, nzt, ngrdcol, &
                                  l_differences )
      ! rcm_in_layer
      call check_flipped_results( "rcm_in_layer", rcm_in_layer, rcm_in_layer_flip, nzt, ngrdcol, &
                                  l_differences )
      ! cloud_cover
      call check_flipped_results( "cloud_cover", cloud_cover, cloud_cover_flip, nzt, ngrdcol, &
                                  l_differences )
      ! invrs_tau_zm
      call check_flipped_results( "invrs_tau_zm", invrs_tau_zm, invrs_tau_zm_flip, nzm, ngrdcol, &
                                  l_differences )
      ! Lscale
      call check_flipped_results( "Lscale", Lscale, Lscale_flip, nzt, ngrdcol, &
                                  l_differences )
      ! sclr variables
      if ( sclr_dim > 0 ) then
         do sclr = 1, sclr_dim
            call check_flipped_results( "sclrm", sclrm(:,:,sclr), sclrm_flip(:,:,sclr), &
                                        nzt, ngrdcol, &
                                        l_differences )
            call check_flipped_results( "wpsclrp", wpsclrp(:,:,sclr), wpsclrp_flip(:,:,sclr), &
                                        nzm, ngrdcol, &
                                        l_differences )
            call check_flipped_results( "sclrp2", sclrp2(:,:,sclr), sclrp2_flip(:,:,sclr), &
                                        nzm, ngrdcol, &
                                        l_differences )
            call check_flipped_results( "sclrprtp", sclrprtp(:,:,sclr), sclrprtp_flip(:,:,sclr), &
                                        nzm, ngrdcol, &
                                        l_differences )
            call check_flipped_results( "sclrpthlp", sclrpthlp(:,:,sclr), &
                                        sclrpthlp_flip(:,:,sclr), nzm, ngrdcol, &
                                        l_differences )
            call check_flipped_results( "sclrp3", sclrp3(:,:,sclr), sclrp3_flip(:,:,sclr), &
                                        nzt, ngrdcol, &
                                        l_differences )
            call check_flipped_results( "sclrpthvp", sclrpthvp(:,:,sclr), &
                                        sclrpthvp_flip(:,:,sclr), nzm, ngrdcol, &
                                        l_differences )
         enddo ! sclr = 1, sclr_dim
      endif ! sclr_dim > 0
      ! edsclrm
      if ( edsclr_dim > 0 ) then
         do edsclr = 1, edsclr_dim
            call check_flipped_results( "edsclrm", edsclrm(:,:,edsclr), edsclrm_flip(:,:,edsclr), &
                                        nzt, ngrdcol, &
                                        l_differences )
         enddo ! edsclr = 1, edsclr_dim
      endif ! edsclr_dim > 0
      ! pdf_params
      call check_flipped_results( "pdf_params%w_1", pdf_params%w_1, &
                                  pdf_params_flip%w_1, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params%w_2", pdf_params%w_2, &
                                  pdf_params_flip%w_2, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params%varnce_w_1", pdf_params%varnce_w_1, &
                                  pdf_params_flip%varnce_w_1, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params%varnce_w_2", pdf_params%varnce_w_2, &
                                  pdf_params_flip%varnce_w_2, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params%rt_1", pdf_params%rt_1, &
                                  pdf_params_flip%rt_1, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params%rt_2", pdf_params%rt_2, &
                                  pdf_params_flip%rt_2, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params%varnce_rt_1", pdf_params%varnce_rt_1, &
                                  pdf_params_flip%varnce_rt_1, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params%varnce_rt_2", pdf_params%varnce_rt_2, &
                                  pdf_params_flip%varnce_rt_2, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params%thl_1", pdf_params%thl_1, &
                                  pdf_params_flip%thl_1, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params%thl_2", pdf_params%thl_2, &
                                  pdf_params_flip%thl_2, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params%varnce_thl_1", pdf_params%varnce_thl_1, &
                                  pdf_params_flip%varnce_thl_1, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params%varnce_thl_2", pdf_params%varnce_thl_2, &
                                  pdf_params_flip%varnce_thl_2, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params%corr_w_rt_1", pdf_params%corr_w_rt_1, &
                                  pdf_params_flip%corr_w_rt_1, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params%corr_w_rt_2", pdf_params%corr_w_rt_2, &
                                  pdf_params_flip%corr_w_rt_2, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params%corr_w_thl_1", pdf_params%corr_w_thl_1, &
                                  pdf_params_flip%corr_w_thl_1, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params%corr_w_thl_2", pdf_params%corr_w_thl_2, &
                                  pdf_params_flip%corr_w_thl_2, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params%corr_rt_thl_1", pdf_params%corr_rt_thl_1, &
                                  pdf_params_flip%corr_rt_thl_1, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params%corr_rt_thl_2", pdf_params%corr_rt_thl_2, &
                                  pdf_params_flip%corr_rt_thl_2, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params%alpha_thl", pdf_params%alpha_thl, &
                                  pdf_params_flip%alpha_thl, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params%alpha_rt", pdf_params%alpha_rt, &
                                  pdf_params_flip%alpha_rt, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params%crt_1", pdf_params%crt_1, &
                                  pdf_params_flip%crt_1, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params%crt_2", pdf_params%crt_2, &
                                  pdf_params_flip%crt_2, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params%cthl_1", pdf_params%cthl_1, &
                                  pdf_params_flip%cthl_1, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params%cthl_2", pdf_params%cthl_2, &
                                  pdf_params_flip%cthl_2, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params%chi_1", pdf_params%chi_1, &
                                  pdf_params_flip%chi_1, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params%chi_2", pdf_params%chi_2, &
                                  pdf_params_flip%chi_2, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params%stdev_chi_1", pdf_params%stdev_chi_1, &
                                  pdf_params_flip%stdev_chi_1, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params%stdev_chi_2", pdf_params%stdev_chi_2, &
                                  pdf_params_flip%stdev_chi_2, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params%stdev_eta_1", pdf_params%stdev_eta_1, &
                                  pdf_params_flip%stdev_eta_1, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params%stdev_eta_2", pdf_params%stdev_eta_2, &
                                  pdf_params_flip%stdev_eta_2, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params%covar_chi_eta_1", pdf_params%covar_chi_eta_1, &
                                  pdf_params_flip%covar_chi_eta_1, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params%covar_chi_eta_2", pdf_params%covar_chi_eta_2, &
                                  pdf_params_flip%covar_chi_eta_2, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params%corr_w_chi_1", pdf_params%corr_w_chi_1, &
                                  pdf_params_flip%corr_w_chi_1, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params%corr_w_chi_2", pdf_params%corr_w_chi_2, &
                                  pdf_params_flip%corr_w_chi_2, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params%corr_w_eta_1", pdf_params%corr_w_eta_1, &
                                  pdf_params_flip%corr_w_eta_1, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params%corr_w_eta_2", pdf_params%corr_w_eta_2, &
                                  pdf_params_flip%corr_w_eta_2, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params%corr_chi_eta_1", pdf_params%corr_chi_eta_1, &
                                  pdf_params_flip%corr_chi_eta_1, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params%corr_chi_eta_2", pdf_params%corr_chi_eta_2, &
                                  pdf_params_flip%corr_chi_eta_2, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params%rsatl_1", pdf_params%rsatl_1, &
                                  pdf_params_flip%rsatl_1, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params%rsatl_2", pdf_params%rsatl_2, &
                                  pdf_params_flip%rsatl_2, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params%rc_1", pdf_params%rc_1, &
                                  pdf_params_flip%rc_1, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params%rc_2", pdf_params%rc_2, &
                                  pdf_params_flip%rc_2, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params%cloud_frac_1", pdf_params%cloud_frac_1, &
                                  pdf_params_flip%cloud_frac_1, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params%cloud_frac_2", pdf_params%cloud_frac_2, &
                                  pdf_params_flip%cloud_frac_2, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params%mixt_frac", pdf_params%mixt_frac, &
                                  pdf_params_flip%mixt_frac, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params%ice_supersat_frac_1", &
                                  pdf_params%ice_supersat_frac_1, &
                                  pdf_params_flip%ice_supersat_frac_1, nzt, ngrdcol, &
                                  l_differences ) 
      call check_flipped_results( "pdf_params%ice_supersat_frac_2", &
                                  pdf_params%ice_supersat_frac_2, &
                                  pdf_params_flip%ice_supersat_frac_2, nzt, ngrdcol, &
                                  l_differences )
      ! pdf_params_zm
      call check_flipped_results( "pdf_params_zm%w_1", pdf_params_zm%w_1, &
                                  pdf_params_zm_flip%w_1, nzm, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params_zm%w_2", pdf_params_zm%w_2, &
                                  pdf_params_zm_flip%w_2, nzm, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params_zm%varnce_w_1", pdf_params_zm%varnce_w_1, &
                                  pdf_params_zm_flip%varnce_w_1, nzm, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params_zm%varnce_w_2", pdf_params_zm%varnce_w_2, &
                                  pdf_params_zm_flip%varnce_w_2, nzm, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params_zm%rt_1", pdf_params_zm%rt_1, &
                                  pdf_params_zm_flip%rt_1, nzm, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params_zm%rt_2", pdf_params_zm%rt_2, &
                                  pdf_params_zm_flip%rt_2, nzm, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params_zm%varnce_rt_1", pdf_params_zm%varnce_rt_1, &
                                  pdf_params_zm_flip%varnce_rt_1, nzm, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params_zm%varnce_rt_2", pdf_params_zm%varnce_rt_2, &
                                  pdf_params_zm_flip%varnce_rt_2, nzm, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params_zm%thl_1", pdf_params_zm%thl_1, &
                                  pdf_params_zm_flip%thl_1, nzm, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params_zm%thl_2", pdf_params_zm%thl_2, &
                                  pdf_params_zm_flip%thl_2, nzm, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params_zm%varnce_thl_1", pdf_params_zm%varnce_thl_1, &
                                  pdf_params_zm_flip%varnce_thl_1, nzm, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params_zm%varnce_thl_2", pdf_params_zm%varnce_thl_2, &
                                  pdf_params_zm_flip%varnce_thl_2, nzm, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params_zm%corr_w_rt_1", pdf_params_zm%corr_w_rt_1, &
                                  pdf_params_zm_flip%corr_w_rt_1, nzm, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params_zm%corr_w_rt_2", pdf_params_zm%corr_w_rt_2, &
                                  pdf_params_zm_flip%corr_w_rt_2, nzm, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params_zm%corr_w_thl_1", pdf_params_zm%corr_w_thl_1, &
                                  pdf_params_zm_flip%corr_w_thl_1, nzm, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params_zm%corr_w_thl_2", pdf_params_zm%corr_w_thl_2, &
                                  pdf_params_zm_flip%corr_w_thl_2, nzm, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params_zm%corr_rt_thl_1", pdf_params_zm%corr_rt_thl_1, &
                                  pdf_params_zm_flip%corr_rt_thl_1, nzm, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params_zm%corr_rt_thl_2", pdf_params_zm%corr_rt_thl_2, &
                                  pdf_params_zm_flip%corr_rt_thl_2, nzm, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params_zm%alpha_thl", pdf_params_zm%alpha_thl, &
                                  pdf_params_zm_flip%alpha_thl, nzm, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params_zm%alpha_rt", pdf_params_zm%alpha_rt, &
                                  pdf_params_zm_flip%alpha_rt, nzm, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params_zm%crt_1", pdf_params_zm%crt_1, &
                                  pdf_params_zm_flip%crt_1, nzm, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params_zm%crt_2", pdf_params_zm%crt_2, &
                                  pdf_params_zm_flip%crt_2, nzm, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params_zm%cthl_1", pdf_params_zm%cthl_1, &
                                  pdf_params_zm_flip%cthl_1, nzm, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params_zm%cthl_2", pdf_params_zm%cthl_2, &
                                  pdf_params_zm_flip%cthl_2, nzm, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params_zm%chi_1", pdf_params_zm%chi_1, &
                                  pdf_params_zm_flip%chi_1, nzm, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params_zm%chi_2", pdf_params_zm%chi_2, &
                                  pdf_params_zm_flip%chi_2, nzm, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params_zm%stdev_chi_1", pdf_params_zm%stdev_chi_1, &
                                  pdf_params_zm_flip%stdev_chi_1, nzm, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params_zm%stdev_chi_2", pdf_params_zm%stdev_chi_2, &
                                  pdf_params_zm_flip%stdev_chi_2, nzm, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params_zm%stdev_eta_1", pdf_params_zm%stdev_eta_1, &
                                  pdf_params_zm_flip%stdev_eta_1, nzm, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params_zm%stdev_eta_2", pdf_params_zm%stdev_eta_2, &
                                  pdf_params_zm_flip%stdev_eta_2, nzm, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params_zm%covar_chi_eta_1", pdf_params_zm%covar_chi_eta_1, &
                                  pdf_params_zm_flip%covar_chi_eta_1, nzm, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params_zm%covar_chi_eta_2", pdf_params_zm%covar_chi_eta_2, &
                                  pdf_params_zm_flip%covar_chi_eta_2, nzm, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params_zm%corr_w_chi_1", pdf_params_zm%corr_w_chi_1, &
                                  pdf_params_zm_flip%corr_w_chi_1, nzm, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params_zm%corr_w_chi_2", pdf_params_zm%corr_w_chi_2, &
                                  pdf_params_zm_flip%corr_w_chi_2, nzm, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params_zm%corr_w_eta_1", pdf_params_zm%corr_w_eta_1, &
                                  pdf_params_zm_flip%corr_w_eta_1, nzm, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params_zm%corr_w_eta_2", pdf_params_zm%corr_w_eta_2, &
                                  pdf_params_zm_flip%corr_w_eta_2, nzm, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params_zm%corr_chi_eta_1", pdf_params_zm%corr_chi_eta_1, &
                                  pdf_params_zm_flip%corr_chi_eta_1, nzm, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params_zm%corr_chi_eta_2", pdf_params_zm%corr_chi_eta_2, &
                                  pdf_params_zm_flip%corr_chi_eta_2, nzm, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params_zm%rsatl_1", pdf_params_zm%rsatl_1, &
                                  pdf_params_zm_flip%rsatl_1, nzm, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params_zm%rsatl_2", pdf_params_zm%rsatl_2, &
                                  pdf_params_zm_flip%rsatl_2, nzm, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params_zm%rc_1", pdf_params_zm%rc_1, &
                                  pdf_params_zm_flip%rc_1, nzm, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params_zm%rc_2", pdf_params_zm%rc_2, &
                                  pdf_params_zm_flip%rc_2, nzm, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params_zm%cloud_frac_1", pdf_params_zm%cloud_frac_1, &
                                  pdf_params_zm_flip%cloud_frac_1, nzm, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params_zm%cloud_frac_2", pdf_params_zm%cloud_frac_2, &
                                  pdf_params_zm_flip%cloud_frac_2, nzm, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params_zm%mixt_frac", pdf_params_zm%mixt_frac, &
                                  pdf_params_zm_flip%mixt_frac, nzm, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params_zm%ice_supersat_frac_1", &
                                  pdf_params_zm%ice_supersat_frac_1, &
                                  pdf_params_zm_flip%ice_supersat_frac_1, nzm, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_params_zm%ice_supersat_frac_2", &
                                  pdf_params_zm%ice_supersat_frac_2, &
                                  pdf_params_zm_flip%ice_supersat_frac_2, nzm, ngrdcol, &
                                  l_differences )
      ! pdf_implicit_coefs_terms
      call check_flipped_results( "pdf_implicit_coefs_terms%coef_wp4_implicit", &
                                  pdf_implicit_coefs_terms%coef_wp4_implicit, &
                                  pdf_implicit_coefs_terms_flip%coef_wp4_implicit, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_implicit_coefs_terms%coef_wp2rtp_implicit", &
                                  pdf_implicit_coefs_terms%coef_wp2rtp_implicit, &
                                  pdf_implicit_coefs_terms_flip%coef_wp2rtp_implicit, &
                                  nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_implicit_coefs_terms%term_wp2rtp_explicit", &
                                  pdf_implicit_coefs_terms%term_wp2rtp_explicit, &
                                  pdf_implicit_coefs_terms_flip%term_wp2rtp_explicit, &
                                  nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_implicit_coefs_terms%coef_wp2thlp_implicit", &
                                  pdf_implicit_coefs_terms%coef_wp2thlp_implicit, &
                                  pdf_implicit_coefs_terms_flip%coef_wp2thlp_implicit, &
                                  nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_implicit_coefs_terms%term_wp2thlp_explicit", &
                                  pdf_implicit_coefs_terms%term_wp2thlp_explicit, &
                                  pdf_implicit_coefs_terms_flip%term_wp2thlp_explicit, &
                                  nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_implicit_coefs_terms%coef_wp2up_implicit", &
                                  pdf_implicit_coefs_terms%coef_wp2up_implicit, &
                                  pdf_implicit_coefs_terms_flip%coef_wp2up_implicit, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_implicit_coefs_terms%term_wp2up_explicit", &
                                  pdf_implicit_coefs_terms%term_wp2up_explicit, &
                                  pdf_implicit_coefs_terms_flip%term_wp2up_explicit, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_implicit_coefs_terms%coef_wp2vp_implicit", &
                                  pdf_implicit_coefs_terms%coef_wp2vp_implicit, &
                                  pdf_implicit_coefs_terms_flip%coef_wp2vp_implicit, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_implicit_coefs_terms%term_wp2vp_explicit", &
                                  pdf_implicit_coefs_terms%term_wp2vp_explicit, &
                                  pdf_implicit_coefs_terms_flip%term_wp2vp_explicit, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_implicit_coefs_terms%coef_wprtp2_implicit", &
                                  pdf_implicit_coefs_terms%coef_wprtp2_implicit, &
                                  pdf_implicit_coefs_terms_flip%coef_wprtp2_implicit, &
                                  nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_implicit_coefs_terms%term_wprtp2_explicit", &
                                  pdf_implicit_coefs_terms%term_wprtp2_explicit, &
                                  pdf_implicit_coefs_terms_flip%term_wprtp2_explicit, &
                                  nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_implicit_coefs_terms%coef_wpthlp2_implicit", &
                                  pdf_implicit_coefs_terms%coef_wpthlp2_implicit, &
                                  pdf_implicit_coefs_terms_flip%coef_wpthlp2_implicit, &
                                  nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_implicit_coefs_terms%term_wpthlp2_explicit", &
                                  pdf_implicit_coefs_terms%term_wpthlp2_explicit, &
                                  pdf_implicit_coefs_terms_flip%term_wpthlp2_explicit, &
                                  nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_implicit_coefs_terms%coef_wprtpthlp_implicit", &
                                  pdf_implicit_coefs_terms%coef_wprtpthlp_implicit, &
                                  pdf_implicit_coefs_terms_flip%coef_wprtpthlp_implicit, &
                                  nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_implicit_coefs_terms%term_wprtpthlp_explicit", &
                                  pdf_implicit_coefs_terms%term_wprtpthlp_explicit, &
                                  pdf_implicit_coefs_terms_flip%term_wprtpthlp_explicit, &
                                  nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_implicit_coefs_terms%coef_wpup2_implicit", &
                                  pdf_implicit_coefs_terms%coef_wpup2_implicit, &
                                  pdf_implicit_coefs_terms_flip%coef_wpup2_implicit, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_implicit_coefs_terms%term_wpup2_explicit", &
                                  pdf_implicit_coefs_terms%term_wpup2_explicit, &
                                  pdf_implicit_coefs_terms_flip%term_wpup2_explicit, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_implicit_coefs_terms%coef_wpvp2_implicit", &
                                  pdf_implicit_coefs_terms%coef_wpvp2_implicit, &
                                  pdf_implicit_coefs_terms_flip%coef_wpvp2_implicit, nzt, ngrdcol, &
                                  l_differences )
      call check_flipped_results( "pdf_implicit_coefs_terms%term_wpvp2_explicit", &
                                  pdf_implicit_coefs_terms%term_wpvp2_explicit, &
                                  pdf_implicit_coefs_terms_flip%term_wpvp2_explicit, nzt, ngrdcol, &
                                  l_differences )
      if ( sclr_dim > 0 ) then
         do sclr = 1, sclr_dim
            call check_flipped_results( "pdf_implicit_coefs_terms%coef_wp2sclrp_implicit", &
                               pdf_implicit_coefs_terms%coef_wp2sclrp_implicit(:,:,sclr), &
                               pdf_implicit_coefs_terms_flip%coef_wp2sclrp_implicit(:,:,sclr), &
                               nzt, ngrdcol, &
                               l_differences )
            call check_flipped_results( "pdf_implicit_coefs_terms%term_wp2sclrp_explicit", &
                               pdf_implicit_coefs_terms%term_wp2sclrp_explicit(:,:,sclr), &
                               pdf_implicit_coefs_terms_flip%term_wp2sclrp_explicit(:,:,sclr), &
                               nzt, ngrdcol, &
                               l_differences )
            call check_flipped_results( "pdf_implicit_coefs_terms%coef_wpsclrp2_implicit", &
                               pdf_implicit_coefs_terms%coef_wpsclrp2_implicit(:,:,sclr), &
                               pdf_implicit_coefs_terms_flip%coef_wpsclrp2_implicit(:,:,sclr), &
                               nzt, ngrdcol, &
                               l_differences )
            call check_flipped_results( "pdf_implicit_coefs_terms%term_wpsclrp2_explicit", &
                               pdf_implicit_coefs_terms%term_wpsclrp2_explicit(:,:,sclr), &
                               pdf_implicit_coefs_terms_flip%term_wpsclrp2_explicit(:,:,sclr), &
                               nzt, ngrdcol, &
                               l_differences )
            call check_flipped_results( "pdf_implicit_coefs_terms%coef_wprtpsclrp_implicit", &
                               pdf_implicit_coefs_terms%coef_wprtpsclrp_implicit(:,:,sclr), &
                               pdf_implicit_coefs_terms_flip%coef_wprtpsclrp_implicit(:,:,sclr), &
                               nzt, ngrdcol, &
                               l_differences )
            call check_flipped_results( "pdf_implicit_coefs_terms%term_wprtpsclrp_explicit", &
                               pdf_implicit_coefs_terms%term_wprtpsclrp_explicit(:,:,sclr), &
                               pdf_implicit_coefs_terms_flip%term_wprtpsclrp_explicit(:,:,sclr), &
                               nzt, ngrdcol, &
                               l_differences )
            call check_flipped_results( "pdf_implicit_coefs_terms%coef_wpthlpsclrp_implicit", &
                               pdf_implicit_coefs_terms%coef_wpthlpsclrp_implicit(:,:,sclr), &
                               pdf_implicit_coefs_terms_flip%coef_wpthlpsclrp_implicit(:,:,sclr), &
                               nzt, ngrdcol, &
                               l_differences )
            call check_flipped_results( "pdf_implicit_coefs_terms%term_wpthlpsclrp_explicit", &
                               pdf_implicit_coefs_terms%term_wpthlpsclrp_explicit(:,:,sclr), &
                               pdf_implicit_coefs_terms_flip%term_wpthlpsclrp_explicit(:,:,sclr), &
                               nzt, ngrdcol, &
                               l_differences )
         enddo ! sclr = 1, sclr_dim
      endif ! sclr_dim > 0


      ! Print a message and stop the run if there are any discrepanices found
      if ( l_differences ) then
         ! Stop the run and exit
         print *, "##################################################"
         print *, "Discrepancy found in ascending vs. descending grid" &
                  // " direction test. Please see messages listed above."
         err_code = clubb_generalized_grd_test_err
      endif ! l_differences


      return

  end subroutine clubb_generalized_grid_testing

  !-----------------------------------------------------------------------------
  subroutine check_flipped_results( varname, var, var_flip, nz, ngrdcol, &
                                    l_differences )

    use constants_clubb, only: &
        zero

    use grid_class, only: &
        flip

    use clubb_precision, only: &
        core_rknd

    implicit none

    ! Input Variables
    character(len=*), intent(in) :: &
      varname   ! Name of Variable being tested

    integer, intent(in) :: &
      nz,      & ! Number of vertical grid levels
      ngrdcol    ! Number of grid columns

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      var,      & ! Variable (ascending grid)
      var_flip    ! Variable (descending grid; flipped vertical array)

    ! Input/Output Variable
    logical, intent(inout) :: &
      l_differences    ! Flag that gets set to true if any differences are found

    ! Local Variable
    integer :: k, i  ! Loop indices


    ! Variable
    do i = 1, ngrdcol
       if ( any( abs( var(i,:) - flip( var_flip(i,:), nz ) ) > zero ) ) then
          l_differences = .true.
          print *, "**************************************************"
          print *, "Differences found in ", trim( varname )
          do k = 1, nz, 1
             if ( abs( var(i,k) - var_flip(i,nz-k+1) ) > zero ) then
                print *, "array index ascending = ", k
                print *, "array index descending = ", nz-k+1
                print *, trim( varname ), " post-solve ascending = ", var(i,k)
                print *, trim( varname ), " post-solve descending = ", var_flip(i,nz-k+1)
             endif
          enddo
       endif
    enddo


    return

  end subroutine check_flipped_results


end module generalized_grid_test
