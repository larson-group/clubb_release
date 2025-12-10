
!===============================================================================
module generalized_grid_test

  ! Guide:
  !
  ! *** Where to use generalized grid statements:
  !
  ! 1) Where the limits of the loop over the vertical grid are not symmetric
  !
  !    Examples of symmetric grid loops:
  !       do k = 1, nzt
  !       do k = 2, nzt-1
  !       do k = 1, nzm
  !       do k = 2, nzm-1
  !
  !    Examples of grid loops that are not symmetric and need special grid
  !    statements:
  !
  !       do k = 1, nzt-1
  !
  !       should be replaced by:
  !
  !       do k = gr%k_lb_zt, gr%k_ub_zt-grid_dir_indx, gr%grid_dir_indx
  !
  !       For an ascending grid, this results in a loop from 1 (lower boundary)
  !       to nzt-1 with an increment of 1. For a descending grid, this results
  !       in a loop from nzt (lower boundary) to 2 with an increment of -1.
  !
  !       As another example:
  !
  !       do k = 2, nzm
  !
  !       should be replaced by:
  !
  !       do k = gr%k_lb_zm+gr%grid_dir_indx, gr%k_ub_zm, gr%grid_dir_indx
  !
  !       For an ascending grid, this results in a loop from 2 to nzm (upper
  !       boundary) with an increment of 1. For a descending grid, this results
  !       in a loop from nzm-1 to 1 (upper boundary) with an increment of -1.
  !
  ! 2) Statements specific to upper or lower boundary grid levels
  !
  !    For example, setting surface values for some variables:
  !
  !    do i = 1, ngrdcol
  !      wpthlp(i,gr%k_lb_zm) = wpthlp_sfc(i)
  !      wprtp(i,gr%k_lb_zm)  = wprtp_sfc(i)
  !      upwp(i,gr%k_lb_zm)   = upwp_sfc(i)
  !      vpwp(i,gr%k_lb_zm)   = vpwp_sfc(i)
  !    end do
  !
  !    As another example, setting the upper or lower boundary conditons used
  !    in CLUBB's lhs and rhs arrays that are used to advance the predictive
  !    equations:
  !
  !    do i = 1, ngrdcol
  !      rhs(i,gr%k_lb_zm) = xap2(i,gr%k_lb_zm)
  !      ! The value of u'^2 or v'^2 at the upper boundary will be set to the
  !      ! threshold minimum value of w_tol_sqd.
  !      rhs(i,gr%k_ub_zm) = w_tol_sqd
  !    end do
  !
  ! 3) Loops over the vertical grid that need to be handled in a consistent
  !    direction, regardless of what grid direction is used
  !
  !    A great example of this is CLUBB's traditional length scale calculation,
  !    which is found in subroutine compute_mixing_length in mixing_length.F90.
  !    Even though loops are symmetric and the calculation doesn't use
  !    boundary-specific statements, generalized grid statements are still
  !    necessary because the calculation starts at the surface or lower boundary
  !    always integrates upward.
  !
  !    Other examples of the code needing to use generalized grid statements
  !    because calculations are needing to be handled consistently in the same
  !    direction, regardless of grid direction, are found in the monotonic flux
  !    limiter and in the vertical hole filler.
  !
  ! Notes:
  !
  ! All statements where 2 quantities are added together match bit-for-bit when
  ! compiled with -O0 optimization and in debug mode. In other words,
  ! A + B = B + A, and it doesn't matter in which order the addition or
  ! subtraction is handled. However, this isn't necessarily true when 3 or more
  ! quantities are added together. In other words, A + B + C doesn't necessarily
  ! match C + B + A bit-for-bit. In order to promote a bit-for-bit match between
  ! ascending and descending grids and allow for a test to check the grid
  ! direction integrity of the code, grid generalization statements have been
  ! added to some places where 3 or more quantities are added together in a row
  ! to ensure that the addition or subtraction always happens in the same order.
  !
  ! *** Comparing Results:
  !
  ! When comparing arrays in the vertical, follow the method found
  ! in this test, as well as in the G-unit test found in
  ! src/G_unit_test_types/rev_direction_grid_test.F90. The value of a
  ! thermodynamic-level field found at level nzt in the ascending grid needs
  ! to match the value found at level 1 in the descending grid, the value
  ! at level 3 in the ascending grid needs to match the value at level nzt-2
  ! in the descending grid, etc. Additionally, for lhs arrays, the values along
  ! the superdiagonal in the ascending grid need to match (in reverse) the
  ! values along the subdiagonal in the descending grid, etc. Examples of this
  ! comparison can be found in src/G_unit_test_types/rev_direction_grid_test.F90
  ! for some commonly used lhs subroutines.
  !
  ! *** When the generalized vertical grid test fails:
  !
  ! 1) Check which flagset number or numbers fail. The flag changes
  !    associated with each flagset can be found in the run_scripts directory
  !    in the file run_bindiff_w_flags_config_core_flags.json (for the
  !    test clubb_generalized_vertical_grid_test) or in the file
  !    run_bindiff_w_flags_config_host_flags.json (for the
  !    test clubb_generalized_vert_grid_host_flags). Note that the test script
  !    automatically runs the unaltered default flag set as the final flagset.
  !    (The .json file might list 17 flagsets, and then the default
  !    configuration is run as flagset 18.)
  !
  ! 2) Within a test for a flag set, all CLUBB test cases are run. Note which
  !    cases fail and pick a single case to use to debug. When there's a choice,
  !    choose a simpler, shorter case to use for debugging.
  !
  ! 3) Run the grid generalization test manually. The two main necessities are:
  !
  !    a) Setting l_test_grid_generalization to true in the file
  !       src/CLUBB_core/model_flags.F90; and
  !
  !    b) Compiling using the compiler script linux_x86_64_gfortran_debug.bash,
  !       which uses -O0 compiler optimization and debug settings, with the
  !       command (when run from the compile directory):
  !       ./compile.bash -c config/linux_x86_64_gfortran_debug.bash.
  !
  !    Further modifications that you might want or need to make
  !    (especially when running one of the CGILS cases) can be found by
  !    following the steps listed in the Jenkins test recipe in
  !    jenkins_tests/clubb_generalized_vertical_grid_test/Jenkinsfile.
  !
  ! 4) Locally, change the flag settings for the run located in
  !    input/tunable_parameters/configurable_model_flags.in to match those
  !    found in the guilty flag set. Then, run CLUBB for the single case
  !    that you chose to test with. It is best to do multiple runs and change
  !    flag settings one at a time until you find the guilty flag.
  !
  ! 5) The grid generalization test will exit and fail after the first timestep
  !    where output that is not bit-for-bit between ascending and descending
  !    grids is detected. A list of variables (and all grid levels for each
  !    variable) that are not bit-for-bit will be printed to the screen. Figure
  !    out which failing variable was calculated first in the sequence. This
  !    will give a clue to roughly where the error might be located.
  !
  ! 6) Use the information gained from which variable goes wrong first, which
  !    flags are being run (flagset info), and information from the changeset
  !    when the test went wrong to get an idea of where the issue might be
  !    in the code. Also take into account the information listed in the
  !    "Where to use generalized grid statements" section above.
  !
  ! 7) Finally, utilize print statements. Narrow it down and then print out
  !    everything until the culprit is found!
  !
  ! 8) Once the culprit is found, rerun all cases using that flagset to check
  !    that they all pass. Rerun the full Jenkins test after the fix has been
  !    committed so that it can be run in the background.
  !
  ! 9) Congratulations! Now don't break it again!
  !
  ! Brian Griffin; May 23, 2025
 
 
  implicit none

  public :: clubb_generalized_grid_testing, &
            silhs_generalized_grid_testing

  private :: check_flipped_results 

  contains

  !=============================================================================
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
               mixt_frac_max_mag, T0, ts_nudge, &                         ! Intent(in)
               rtm_min, rtm_nudge_max_altitude, &                         ! Intent(in)
               clubb_config_flags, &                                      ! Intent(in)
               stats_metadata, &                                          ! Intent(in)
               stats_zt, stats_zm, stats_sfc, &                           ! intent(inout)
               um, vm, upwp, vpwp, up2, vp2, up3, vp3, &                  ! Intent(inout)
               thlm, rtm, wprtp, wpthlp, &                                ! Intent(inout)
               wp2, wp3, rtp2, rtp3, thlp2, thlp3, rtpthlp, &             ! Intent(inout)
               sclrm, sclrp2, sclrp3, sclrprtp, sclrpthlp, &              ! Intent(inout)
               wpsclrp, edsclrm, err_info, &                              ! Intent(inout)
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

    use constants_clubb, only: &
        fstderr

    use grid_class, only: &
        grid, & ! Type(s)
        flip    ! Procedure(s)

    use clubb_api_module, only: &
        advance_clubb_core_api, &
        init_pdf_implicit_coefs_terms_api

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
        init_pdf_params

    use model_flags, only: &
        clubb_config_flags_type, &
        iiPDF_new, &
        iiPDF_new_hybrid

    use stats_type, only: &
        stats ! Type

    use stats_variables, only: &
        stats_metadata_type

    use error_code, only: &
        clubb_generalized_grd_test_err, &
        clubb_no_error, &
        clubb_fatal_error

    use err_info_type_module, only: &
      err_info_type        ! Type

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

    real( kind = core_rknd ), dimension(ngrdcol,nzm,hydromet_dim), intent(in) :: &
      wphydrometp    ! Covariance of w and a hydrometeor   [(m/s) <hm units>]

    real( kind = core_rknd ), dimension(ngrdcol,nzt,hydromet_dim), intent(in) :: &
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

    type(err_info_type), intent(inout) :: &
      err_info        ! err_info struct containing err_code and err_header

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

    real( kind = core_rknd ), dimension(ngrdcol,nzm,hydromet_dim) :: &
      wphydrometp_flip    ! Covariance of w and a hydrometeor   [(m/s) <hm units>]

    real( kind = core_rknd ), dimension(ngrdcol,nzt,hydromet_dim) :: &
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


      ! Allocate space and initialize flipped pdf parameter terms
      call init_pdf_params( gr%nzt, ngrdcol, pdf_params_flip )
      call init_pdf_params( gr%nzm, ngrdcol, pdf_params_zm_flip )

      call init_pdf_implicit_coefs_terms_api( gr%nzt, ngrdcol, sclr_dim, &   ! Intent(in)
                                              pdf_implicit_coefs_terms_flip ) ! Intent(out)

      ! Set up "flipped" variables for call to descending grid.
      thlm_forcing_flip(:,:) = thlm_forcing(:,nzt:1:-1)
      rtm_forcing_flip(:,:) = rtm_forcing(:,nzt:1:-1)
      um_forcing_flip(:,:) = um_forcing(:,nzt:1:-1)
      vm_forcing_flip(:,:) = vm_forcing(:,nzt:1:-1)
      wprtp_forcing_flip(:,:) = wprtp_forcing(:,nzm:1:-1)
      wpthlp_forcing_flip(:,:) = wpthlp_forcing(:,nzm:1:-1)
      rtp2_forcing_flip(:,:) = rtp2_forcing(:,nzm:1:-1)
      thlp2_forcing_flip(:,:) = thlp2_forcing(:,nzm:1:-1)
      rtpthlp_forcing_flip(:,:) = rtpthlp_forcing(:,nzm:1:-1)
      wm_zm_flip(:,:) = wm_zm(:,nzm:1:-1)
      wm_zt_flip(:,:) = wm_zt(:,nzt:1:-1)
      rtm_ref_flip(:,:) = rtm_ref(:,nzt:1:-1)
      thlm_ref_flip(:,:) = thlm_ref(:,nzt:1:-1)
      um_ref_flip(:,:) = um_ref(:,nzt:1:-1)
      vm_ref_flip(:,:) = vm_ref(:,nzt:1:-1)
      ug_flip(:,:) = ug(:,nzt:1:-1)
      vg_flip(:,:) = vg(:,nzt:1:-1)
      p_in_Pa_flip(:,:) = p_in_Pa(:,nzt:1:-1)
      rho_zm_flip(:,:) = rho_zm(:,nzm:1:-1)
      rho_flip(:,:) = rho(:,nzt:1:-1)
      exner_flip(:,:) = exner(:,nzt:1:-1)
      rho_ds_zm_flip(:,:) = rho_ds_zm(:,nzm:1:-1)
      rho_ds_zt_flip(:,:) = rho_ds_zt(:,nzt:1:-1)
      invrs_rho_ds_zm_flip(:,:) = invrs_rho_ds_zm(:,nzm:1:-1)
      invrs_rho_ds_zt_flip(:,:) = invrs_rho_ds_zt(:,nzt:1:-1)
      thv_ds_zt_flip(:,:) = thv_ds_zt(:,nzt:1:-1)
      thv_ds_zm_flip(:,:) = thv_ds_zm(:,nzm:1:-1)
      rfrzm_flip(:,:) = rfrzm(:,nzt:1:-1)
    
      um_flip(:,:) = um(:,nzt:1:-1)
      vm_flip(:,:) = vm(:,nzt:1:-1)
      upwp_flip(:,:) = upwp(:,nzm:1:-1)
      vpwp_flip(:,:) = vpwp(:,nzm:1:-1)
      up2_flip(:,:) = up2(:,nzm:1:-1)
      vp2_flip(:,:) = vp2(:,nzm:1:-1)
      up3_flip(:,:) = up3(:,nzt:1:-1)
      vp3_flip(:,:) = vp3(:,nzt:1:-1)
      thlm_flip(:,:) = thlm(:,nzt:1:-1)
      rtm_flip(:,:) = rtm(:,nzt:1:-1)
      wpthlp_flip(:,:) = wpthlp(:,nzm:1:-1)
      wprtp_flip(:,:) = wprtp(:,nzm:1:-1)
      wp2_flip(:,:) = wp2(:,nzm:1:-1)
      wp3_flip(:,:) = wp3(:,nzt:1:-1)
      rtp2_flip(:,:) = rtp2(:,nzm:1:-1)
      rtp3_flip(:,:) = rtp3(:,nzt:1:-1)
      thlp2_flip(:,:) = thlp2(:,nzm:1:-1)
      thlp3_flip(:,:) = thlp3(:,nzt:1:-1)
      rtpthlp_flip(:,:) = rtpthlp(:,nzm:1:-1)
      rcm_flip(:,:) = rcm(:,nzt:1:-1)
      cloud_frac_flip(:,:) = cloud_frac(:,nzt:1:-1)
      wpthvp_flip(:,:) = wpthvp(:,nzm:1:-1)
      wp2thvp_flip(:,:) = wp2thvp(:,nzt:1:-1)
      rtpthvp_flip(:,:) = rtpthvp(:,nzm:1:-1)
      thlpthvp_flip(:,:) = thlpthvp(:,nzm:1:-1)
      wp2rtp_flip(:,:) = wp2rtp(:,nzt:1:-1)
      wp2thlp_flip(:,:) = wp2thlp(:,nzt:1:-1)
      uprcp_flip(:,:) = uprcp(:,nzm:1:-1)
      vprcp_flip(:,:) = vprcp(:,nzm:1:-1)
      rc_coef_zm_flip(:,:) = rc_coef_zm(:,nzm:1:-1)
      wp4_flip(:,:) = wp4(:,nzm:1:-1)
      wpup2_flip(:,:) = wpup2(:,nzt:1:-1)
      wpvp2_flip(:,:) = wpvp2(:,nzt:1:-1)
      wp2up2_flip(:,:) = wp2up2(:,nzm:1:-1)
      wp2vp2_flip(:,:) = wp2vp2(:,nzm:1:-1)
      ice_supersat_frac_flip(:,:) = ice_supersat_frac(:,nzt:1:-1)

      if ( sclr_dim > 0 ) then
        sclrm_forcing_flip(:,:,:) = sclrm_forcing(:,nzt:1:-1,:)
        sclrm_flip(:,:,:) = sclrm(:,nzt:1:-1,:)
        wpsclrp_flip(:,:,:) = wpsclrp(:,nzm:1:-1,:)
        sclrprtp_flip(:,:,:) = sclrprtp(:,nzm:1:-1,:)
        sclrpthlp_flip(:,:,:) = sclrpthlp(:,nzm:1:-1,:)
        sclrpthvp_flip(:,:,:) = sclrpthvp(:,nzm:1:-1,:)
        sclrp2_flip(:,:,:) = sclrp2(:,nzm:1:-1,:)
        sclrp3_flip(:,:,:) = sclrp3(:,nzt:1:-1,:)
      end if

      ! edsclr variables
      if ( edsclr_dim > 0 ) then
        edsclrm_forcing_flip(:,:,:) = edsclrm_forcing(:,nzt:1:-1,:)
        edsclrm_flip(:,:,:) = edsclrm(:,nzt:1:-1,:)
      end if

      ! hydrometeor variables
      if ( hydromet_dim > 0 ) then
        wphydrometp_flip(:,:,:) = wphydrometp(:,nzm:1:-1,:)
        wp2hmp_flip(:,:,:) = wp2hmp(:,nzt:1:-1,:)
        rtphmp_zt_flip(:,:,:) = rtphmp_zt(:,nzt:1:-1,:)
        thlphmp_zt_flip(:,:,:) = thlphmp_zt(:,nzt:1:-1,:)
      end if

      ! pdf_params
      pdf_params_flip%w_1(:,:) = pdf_params%w_1(:,nzt:1:-1)
      pdf_params_flip%w_2(:,:) = pdf_params%w_2(:,nzt:1:-1)
      pdf_params_flip%varnce_w_1(:,:) = pdf_params%varnce_w_1(:,nzt:1:-1)
      pdf_params_flip%varnce_w_2(:,:) = pdf_params%varnce_w_2(:,nzt:1:-1)
      pdf_params_flip%rt_1(:,:) = pdf_params%rt_1(:,nzt:1:-1)
      pdf_params_flip%rt_2(:,:) = pdf_params%rt_2(:,nzt:1:-1)
      pdf_params_flip%varnce_rt_1(:,:) = pdf_params%varnce_rt_1(:,nzt:1:-1)
      pdf_params_flip%varnce_rt_2(:,:) = pdf_params%varnce_rt_2(:,nzt:1:-1)
      pdf_params_flip%thl_1(:,:) = pdf_params%thl_1(:,nzt:1:-1)
      pdf_params_flip%thl_2(:,:) = pdf_params%thl_2(:,nzt:1:-1)
      pdf_params_flip%varnce_thl_1(:,:) = pdf_params%varnce_thl_1(:,nzt:1:-1)
      pdf_params_flip%varnce_thl_2(:,:) = pdf_params%varnce_thl_2(:,nzt:1:-1)
      pdf_params_flip%corr_w_rt_1(:,:) = pdf_params%corr_w_rt_1(:,nzt:1:-1)
      pdf_params_flip%corr_w_rt_2(:,:) = pdf_params%corr_w_rt_2(:,nzt:1:-1)
      pdf_params_flip%corr_w_thl_1(:,:) = pdf_params%corr_w_thl_1(:,nzt:1:-1)
      pdf_params_flip%corr_w_thl_2(:,:) = pdf_params%corr_w_thl_2(:,nzt:1:-1)
      pdf_params_flip%corr_rt_thl_1(:,:) = pdf_params%corr_rt_thl_1(:,nzt:1:-1)
      pdf_params_flip%corr_rt_thl_2(:,:) = pdf_params%corr_rt_thl_2(:,nzt:1:-1)
      pdf_params_flip%alpha_thl(:,:) = pdf_params%alpha_thl(:,nzt:1:-1)
      pdf_params_flip%alpha_rt(:,:) = pdf_params%alpha_rt(:,nzt:1:-1)
      pdf_params_flip%crt_1(:,:) = pdf_params%crt_1(:,nzt:1:-1)
      pdf_params_flip%crt_2(:,:) = pdf_params%crt_2(:,nzt:1:-1)
      pdf_params_flip%cthl_1(:,:) = pdf_params%cthl_1(:,nzt:1:-1)
      pdf_params_flip%cthl_2(:,:) = pdf_params%cthl_2(:,nzt:1:-1)
      pdf_params_flip%chi_1(:,:) = pdf_params%chi_1(:,nzt:1:-1)
      pdf_params_flip%chi_2(:,:) = pdf_params%chi_2(:,nzt:1:-1)
      pdf_params_flip%stdev_chi_1(:,:) = pdf_params%stdev_chi_1(:,nzt:1:-1)
      pdf_params_flip%stdev_chi_2(:,:) = pdf_params%stdev_chi_2(:,nzt:1:-1)
      pdf_params_flip%stdev_eta_1(:,:) = pdf_params%stdev_eta_1(:,nzt:1:-1)
      pdf_params_flip%stdev_eta_2(:,:) = pdf_params%stdev_eta_2(:,nzt:1:-1)
      pdf_params_flip%covar_chi_eta_1(:,:) = pdf_params%covar_chi_eta_1(:,nzt:1:-1)
      pdf_params_flip%covar_chi_eta_2(:,:) = pdf_params%covar_chi_eta_2(:,nzt:1:-1)
      pdf_params_flip%corr_w_chi_1(:,:) = pdf_params%corr_w_chi_1(:,nzt:1:-1)
      pdf_params_flip%corr_w_chi_2(:,:) = pdf_params%corr_w_chi_2(:,nzt:1:-1)
      pdf_params_flip%corr_w_eta_1(:,:) = pdf_params%corr_w_eta_1(:,nzt:1:-1)
      pdf_params_flip%corr_w_eta_2(:,:) = pdf_params%corr_w_eta_2(:,nzt:1:-1)
      pdf_params_flip%corr_chi_eta_1(:,:) = pdf_params%corr_chi_eta_1(:,nzt:1:-1)
      pdf_params_flip%corr_chi_eta_2(:,:) = pdf_params%corr_chi_eta_2(:,nzt:1:-1)
      pdf_params_flip%rsatl_1(:,:) = pdf_params%rsatl_1(:,nzt:1:-1)
      pdf_params_flip%rsatl_2(:,:) = pdf_params%rsatl_2(:,nzt:1:-1)
      pdf_params_flip%rc_1(:,:) = pdf_params%rc_1(:,nzt:1:-1)
      pdf_params_flip%rc_2(:,:) = pdf_params%rc_2(:,nzt:1:-1)
      pdf_params_flip%cloud_frac_1(:,:) = pdf_params%cloud_frac_1(:,nzt:1:-1)
      pdf_params_flip%cloud_frac_2(:,:) = pdf_params%cloud_frac_2(:,nzt:1:-1)
      pdf_params_flip%mixt_frac(:,:) = pdf_params%mixt_frac(:,nzt:1:-1)
      pdf_params_flip%ice_supersat_frac_1(:,:) = pdf_params%ice_supersat_frac_1(:,nzt:1:-1)
      pdf_params_flip%ice_supersat_frac_2(:,:) = pdf_params%ice_supersat_frac_2(:,nzt:1:-1)

      ! pdf_params_zm
      pdf_params_zm_flip%w_1(:,:) = pdf_params_zm%w_1(:,nzm:1:-1)
      pdf_params_zm_flip%w_2(:,:) = pdf_params_zm%w_2(:,nzm:1:-1)
      pdf_params_zm_flip%varnce_w_1(:,:) = pdf_params_zm%varnce_w_1(:,nzm:1:-1)
      pdf_params_zm_flip%varnce_w_2(:,:) = pdf_params_zm%varnce_w_2(:,nzm:1:-1)
      pdf_params_zm_flip%rt_1(:,:) = pdf_params_zm%rt_1(:,nzm:1:-1)
      pdf_params_zm_flip%rt_2(:,:) = pdf_params_zm%rt_2(:,nzm:1:-1)
      pdf_params_zm_flip%varnce_rt_1(:,:) = pdf_params_zm%varnce_rt_1(:,nzm:1:-1)
      pdf_params_zm_flip%varnce_rt_2(:,:) = pdf_params_zm%varnce_rt_2(:,nzm:1:-1)
      pdf_params_zm_flip%thl_1(:,:) = pdf_params_zm%thl_1(:,nzm:1:-1)
      pdf_params_zm_flip%thl_2(:,:) = pdf_params_zm%thl_2(:,nzm:1:-1)
      pdf_params_zm_flip%varnce_thl_1(:,:) = pdf_params_zm%varnce_thl_1(:,nzm:1:-1)
      pdf_params_zm_flip%varnce_thl_2(:,:) = pdf_params_zm%varnce_thl_2(:,nzm:1:-1)
      pdf_params_zm_flip%corr_w_rt_1(:,:) = pdf_params_zm%corr_w_rt_1(:,nzm:1:-1)
      pdf_params_zm_flip%corr_w_rt_2(:,:) = pdf_params_zm%corr_w_rt_2(:,nzm:1:-1)
      pdf_params_zm_flip%corr_w_thl_1(:,:) = pdf_params_zm%corr_w_thl_1(:,nzm:1:-1)
      pdf_params_zm_flip%corr_w_thl_2(:,:) = pdf_params_zm%corr_w_thl_2(:,nzm:1:-1)
      pdf_params_zm_flip%corr_rt_thl_1(:,:) = pdf_params_zm%corr_rt_thl_1(:,nzm:1:-1)
      pdf_params_zm_flip%corr_rt_thl_2(:,:) = pdf_params_zm%corr_rt_thl_2(:,nzm:1:-1)
      pdf_params_zm_flip%alpha_thl(:,:) = pdf_params_zm%alpha_thl(:,nzm:1:-1)
      pdf_params_zm_flip%alpha_rt(:,:) = pdf_params_zm%alpha_rt(:,nzm:1:-1)
      pdf_params_zm_flip%crt_1(:,:) = pdf_params_zm%crt_1(:,nzm:1:-1)
      pdf_params_zm_flip%crt_2(:,:) = pdf_params_zm%crt_2(:,nzm:1:-1)
      pdf_params_zm_flip%cthl_1(:,:) = pdf_params_zm%cthl_1(:,nzm:1:-1)
      pdf_params_zm_flip%cthl_2(:,:) = pdf_params_zm%cthl_2(:,nzm:1:-1)
      pdf_params_zm_flip%chi_1(:,:) = pdf_params_zm%chi_1(:,nzm:1:-1)
      pdf_params_zm_flip%chi_2(:,:) = pdf_params_zm%chi_2(:,nzm:1:-1)
      pdf_params_zm_flip%stdev_chi_1(:,:) = pdf_params_zm%stdev_chi_1(:,nzm:1:-1)
      pdf_params_zm_flip%stdev_chi_2(:,:) = pdf_params_zm%stdev_chi_2(:,nzm:1:-1)
      pdf_params_zm_flip%stdev_eta_1(:,:) = pdf_params_zm%stdev_eta_1(:,nzm:1:-1)
      pdf_params_zm_flip%stdev_eta_2(:,:) = pdf_params_zm%stdev_eta_2(:,nzm:1:-1)
      pdf_params_zm_flip%covar_chi_eta_1(:,:) = pdf_params_zm%covar_chi_eta_1(:,nzm:1:-1)
      pdf_params_zm_flip%covar_chi_eta_2(:,:) = pdf_params_zm%covar_chi_eta_2(:,nzm:1:-1)
      pdf_params_zm_flip%corr_w_chi_1(:,:) = pdf_params_zm%corr_w_chi_1(:,nzm:1:-1)
      pdf_params_zm_flip%corr_w_chi_2(:,:) = pdf_params_zm%corr_w_chi_2(:,nzm:1:-1)
      pdf_params_zm_flip%corr_w_eta_1(:,:) = pdf_params_zm%corr_w_eta_1(:,nzm:1:-1)
      pdf_params_zm_flip%corr_w_eta_2(:,:) = pdf_params_zm%corr_w_eta_2(:,nzm:1:-1)
      pdf_params_zm_flip%corr_chi_eta_1(:,:) = pdf_params_zm%corr_chi_eta_1(:,nzm:1:-1)
      pdf_params_zm_flip%corr_chi_eta_2(:,:) = pdf_params_zm%corr_chi_eta_2(:,nzm:1:-1)
      pdf_params_zm_flip%rsatl_1(:,:) = pdf_params_zm%rsatl_1(:,nzm:1:-1)
      pdf_params_zm_flip%rsatl_2(:,:) = pdf_params_zm%rsatl_2(:,nzm:1:-1)
      pdf_params_zm_flip%rc_1(:,:) = pdf_params_zm%rc_1(:,nzm:1:-1)
      pdf_params_zm_flip%rc_2(:,:) = pdf_params_zm%rc_2(:,nzm:1:-1)
      pdf_params_zm_flip%cloud_frac_1(:,:) = pdf_params_zm%cloud_frac_1(:,nzm:1:-1)
      pdf_params_zm_flip%cloud_frac_2(:,:) = pdf_params_zm%cloud_frac_2(:,nzm:1:-1)
      pdf_params_zm_flip%mixt_frac(:,:) = pdf_params_zm%mixt_frac(:,nzm:1:-1)
      pdf_params_zm_flip%ice_supersat_frac_1(:,:) = pdf_params_zm%ice_supersat_frac_1(:,nzm:1:-1)
      pdf_params_zm_flip%ice_supersat_frac_2(:,:) = pdf_params_zm%ice_supersat_frac_2(:,nzm:1:-1)

      if ( clubb_config_flags%iiPDF_type == iiPDF_new .or. clubb_config_flags%iiPDF_type == iiPDF_new_hybrid ) then
        do i = 1, ngrdcol
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
              end do
          end if ! sclr_dim > 0
        end do ! i = 1, ngrdcol
      end if

      Kh_zm_flip(:,:) = Kh_zm(:,nzm:1:-1)
      Kh_zt_flip(:,:) = Kh_zt(:,nzt:1:-1)
      thlprcp_flip(:,:) = thlprcp(:,nzm:1:-1)
      wprcp_flip(:,:) = wprcp(:,nzm:1:-1)
      w_up_in_cloud_flip(:,:) = w_up_in_cloud(:,nzt:1:-1)
      w_down_in_cloud_flip(:,:) = w_down_in_cloud(:,nzt:1:-1)
      cloudy_updraft_frac_flip(:,:) = cloudy_updraft_frac(:,nzt:1:-1)
      cloudy_downdraft_frac_flip(:,:) = cloudy_downdraft_frac(:,nzt:1:-1)
      rcm_in_layer_flip(:,:) = rcm_in_layer(:,nzt:1:-1)
      cloud_cover_flip(:,:) = cloud_cover(:,nzt:1:-1)
      invrs_tau_zm_flip(:,:) = invrs_tau_zm(:,nzm:1:-1)
      Lscale_flip(:,:) = Lscale(:,nzt:1:-1)


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
              mixt_frac_max_mag, T0, ts_nudge, &                   ! Intent(in)
              rtm_min, rtm_nudge_max_altitude, &                   ! Intent(in)
              clubb_config_flags, &                                ! Intent(in)
              stats_metadata, &                                    ! Intent(in)
              stats_zt, stats_zm, stats_sfc, &                     ! intent(inout)
              um, vm, upwp, vpwp, up2, vp2, up3, vp3, &            ! Intent(inout)
              thlm, rtm, wprtp, wpthlp, &                          ! Intent(inout)
              wp2, wp3, rtp2, rtp3, thlp2, thlp3, rtpthlp, &       ! Intent(inout)
              sclrm, sclrp2, sclrp3, sclrprtp, sclrpthlp, &        ! Intent(inout)
              wpsclrp, edsclrm, err_info, &                        ! Intent(inout)
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
      if ( any(err_info%err_code == clubb_fatal_error) ) then
         write(fstderr, *) "Fatal error in advance_clubb_core using ascending grid direction"
         write(fstderr, *) "Generalized grid test continuing anyway"
         ! Reset error code
         err_info%err_code = clubb_no_error
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
              mixt_frac_max_mag, T0, ts_nudge, &                                    ! Intent(in)
              rtm_min, rtm_nudge_max_altitude, &                                    ! Intent(in)
              clubb_config_flags, &                                                 ! Intent(in)
              stats_metadata_flip, &                                                ! Intent(in)
              stats_zt, stats_zm, stats_sfc, &                                      ! intent(inout)
              um_flip, vm_flip, upwp_flip, vpwp_flip, up2_flip, vp2_flip, up3_flip, vp3_flip, & ! Intent(inout)
              thlm_flip, rtm_flip, wprtp_flip, wpthlp_flip, &                       ! Intent(inout)
              wp2_flip, wp3_flip, rtp2_flip, rtp3_flip, thlp2_flip, thlp3_flip, rtpthlp_flip, & ! Intent(inout)
              sclrm_flip, sclrp2_flip, sclrp3_flip, sclrprtp_flip, sclrpthlp_flip, & ! Intent(inout)
              wpsclrp_flip, edsclrm_flip, err_info, &                               ! Intent(inout)
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
      if ( clubb_config_flags%iiPDF_type == iiPDF_new .or. clubb_config_flags%iiPDF_type == iiPDF_new_hybrid ) then
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
      end if


      ! Print a message and stop the run if there are any discrepanices found
      if ( l_differences ) then
         ! Stop the run and exit
         print *, "##################################################"
         print *, "Discrepancy found in ascending vs. descending grid" &
                  // " direction test. Please see messages listed above."
          ! General error -> set all entries to clubb_generalized_grd_test_err
         err_info%err_code = clubb_generalized_grd_test_err
      endif ! l_differences

      return

  end subroutine clubb_generalized_grid_testing

  !=============================================================================
  subroutine silhs_generalized_grid_testing &
             ( gr, gr_desc, ngrdcol, pdf_dim, hydromet_dim,     & ! In
               itime, dt_main, vert_decorr_coef,                & ! In
               Nc_in_cloud, cloud_frac, ice_supersat_frac,      & ! In
               rho_ds_zt, Lscale, Kh_zm, hydromet, wphydrometp, & ! In
               corr_array_n_cloud, corr_array_n_below,          & ! In
               hm_metadata, pdf_params, clubb_params,           & ! In
               clubb_config_flags, silhs_config_flags,          & ! In
               l_rad_itime, stats_metadata,                     & ! In
               stats_zt, stats_zm, stats_sfc,                   & ! In/Out
               stats_lh_zt, stats_lh_sfc, err_info,             & ! In/Out
               time_clubb_pdf, time_stop, time_start,           & ! In/Out
               hydrometp2,                                      & ! Out
               mu_x_1_n, mu_x_2_n,                              & ! Out
               sigma_x_1_n, sigma_x_2_n,                        & ! Out
               corr_array_1_n, corr_array_2_n,                  & ! Out
               corr_cholesky_mtx_1, corr_cholesky_mtx_2,        & ! Out
               precip_fracs,                                    & ! In/Out
               rtphmp_zt, thlphmp_zt, wp2hmp,                   & ! Out
               X_nl_all_levs, X_mixt_comp_all_levs,             & ! Out
               lh_sample_point_weights,                         & ! Out
               lh_rt_clipped, lh_thl_clipped, lh_rc_clipped,    & ! Out
               lh_rv_clipped, lh_Nc_clipped,                    & ! Out
               hydromet_pdf_params )                              ! Optional(out)

    use grid_class, only: &
        grid, & ! Type(s)
        flip    ! Procedure(s)

    use pdf_hydromet_microphys_wrapper, only: &
        pdf_hydromet_microphys_prep    ! Procedure(s)

    use pdf_parameter_module, only: &
        pdf_parameter,   & ! Type(s)
        init_pdf_params    ! Procedure(s)

    use hydromet_pdf_parameter_module, only: &
        hydromet_pdf_parameter,   & ! Type(s)
        precipitation_fractions,  & ! Procedure(s)
        init_hydromet_pdf_params

    use corr_varnce_module, only: &
        hm_metadata_type    ! Type(s)

    use model_flags, only: &
        clubb_config_flags_type, & ! Type(s)
        l_silhs_rad                ! Variable(s)

    use parameters_silhs, only: &
        silhs_config_flags_type    ! Types(s)

    use stats_type, only: &
        stats ! Type(s)

    use stats_variables, only: &
        stats_metadata_type    ! Type(s)

#ifdef SILHS
    use parameters_microphys, only: &
        lh_num_samples    ! Variable(s)
#endif /*SILHS*/

    use parameter_indices, only: &
        nparams    ! Variable(s)

    use parameters_microphys, only: &
        microphys_scheme    ! Variable(s)

#ifdef SILHS
    use parameters_microphys, only: &
        lh_microphys_type,     & ! Variable(s)
        lh_microphys_disabled, &
        lh_num_samples
#endif /*SILHS*/

    use clubb_api_module, only: &
        init_precip_fracs_api    ! Procedure(s)

    use error_code, only: &
        clubb_generalized_grd_test_err, &
        clubb_no_error, &
        clubb_fatal_error

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    use err_info_type_module, only: &
      err_info_type        ! Type

    implicit none

    ! Input Variables
    type (grid), intent(in) :: &
      gr,      & ! Grid variable type for ascending grid
      gr_desc    ! Grid variable type for descending grid

    integer, intent(in) :: &
      ngrdcol,      & ! Number of grid columns
      pdf_dim,      & ! Number of variables in the correlation array
      hydromet_dim, & ! Number of hydrometeor species
      itime

    real( kind = core_rknd ), intent(in) :: &
      dt_main,          & ! Model timestep                              [s]
      vert_decorr_coef    ! Empirically defined de-correlation constant [-]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), intent(in) :: &
      Nc_in_cloud,       & ! Mean (in-cloud) cloud droplet conc.       [num/kg]
      cloud_frac,        & ! Cloud fraction                            [-]
      ice_supersat_frac, & ! Ice supersaturation fraction              [-]
      rho_ds_zt,         & ! Dry, base-state density on thermo. levs.  [kg/m^3]
      Lscale               ! Turbulent Mixing Length                   [m]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzm), intent(in) :: &
      Kh_zm                ! Eddy diffusivity coef. on momentum levels [m^2/s]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt,hydromet_dim), intent(in) :: &
      hydromet       ! Mean of hydrometeor, hm (overall) (t-levs.) [units]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzm,hydromet_dim), intent(in) :: &
      wphydrometp    ! Covariance < w'h_m' > (momentum levels)     [(m/s)units]

    real( kind = core_rknd ), dimension(pdf_dim,pdf_dim), intent(in) :: &
      corr_array_n_cloud, & ! Prescribed normal space corr. array in cloud  [-]
      corr_array_n_below    ! Prescribed normal space corr. array below cl. [-]

    type (hm_metadata_type), intent(in) :: &
      hm_metadata

    type(pdf_parameter), intent(in) :: &
      pdf_params    ! PDF parameters                               [units vary]

    real( kind = core_rknd ), dimension(nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

    type(clubb_config_flags_type), intent(in) :: &
      clubb_config_flags ! Derived type holding all configurable CLUBB flags

    type(silhs_config_flags_type), intent(in) :: &
      silhs_config_flags

    logical, intent(in) :: &
      l_rad_itime

    type (stats_metadata_type), intent(in) :: &
      stats_metadata

    ! Input/Output Variables
    type (stats), dimension(ngrdcol), intent(inout) :: &
      stats_zt, &
      stats_zm, &
      stats_sfc, &
      stats_lh_zt, &
      stats_lh_sfc

    type(err_info_type), intent(inout) :: &
      err_info        ! err_info struct containing err_code and err_header

    real( kind = core_rknd ), intent(inout) :: &
      time_clubb_pdf, & ! time spent in setup_pdf_parameters and hydrometeor_mixed_moments [s]
      time_start,     & ! help variables to measure the time [s]
      time_stop         ! help variables to measure the time [s]

    ! Output Variables
    real( kind = core_rknd ), dimension(ngrdcol,gr%nzm,hydromet_dim), intent(out) :: &
      hydrometp2    ! Variance of a hydrometeor (overall) (m-levs.)   [units^2]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt,pdf_dim), intent(out) :: &
      mu_x_1_n,    & ! Mean array (normal space): PDF vars. (comp. 1) [un. vary]
      mu_x_2_n,    & ! Mean array (normal space): PDF vars. (comp. 2) [un. vary]
      sigma_x_1_n, & ! Std. dev. array (normal space): PDF vars (comp. 1) [u.v.]
      sigma_x_2_n    ! Std. dev. array (normal space): PDF vars (comp. 2) [u.v.]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt,pdf_dim,pdf_dim), intent(out) :: &
      corr_array_1_n, & ! Corr. array (normal space) of PDF vars. (comp. 1)  [-]
      corr_array_2_n    ! Corr. array (normal space) of PDF vars. (comp. 2)  [-]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt,pdf_dim,pdf_dim), intent(out) :: &
      corr_cholesky_mtx_1, & ! Transposed corr. cholesky matrix, 1st comp. [-]
      corr_cholesky_mtx_2    ! Transposed corr. cholesky matrix, 2nd comp. [-]

    ! This is only an output, but it contains allocated arrays, so we need to treat it as inout
    type(precipitation_fractions), intent(inout) :: &
      precip_fracs           ! Precipitation fractions      [-]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt,hydromet_dim), intent(out) :: &
      wp2hmp,     & ! Higher-order mixed moment:  < w'^2 hm' > [(m/s)^2<hm un.>]
      rtphmp_zt,  & ! Covariance of rt and hm (on t-levs.)     [(kg/kg)<hm un.>]
      thlphmp_zt    ! Covariance of thl and hm (on t-levs.)    [K<hm units>]

    real( kind = core_rknd ), dimension(ngrdcol,lh_num_samples,gr%nzt,pdf_dim), intent(out) :: &
      X_nl_all_levs ! Sample that is transformed ultimately to normal-lognormal

    integer, dimension(ngrdcol,lh_num_samples,gr%nzt), intent(out) :: &
      X_mixt_comp_all_levs ! Which mixture component we're in

    real( kind = core_rknd ), dimension(ngrdcol,lh_num_samples,gr%nzt), intent(out) :: &
      lh_sample_point_weights

    real( kind = core_rknd ), dimension(ngrdcol,lh_num_samples,gr%nzt), intent(out) :: &
      lh_rt_clipped,  & ! rt generated from silhs sample points
      lh_thl_clipped, & ! thl generated from silhs sample points
      lh_rc_clipped,  & ! rc generated from silhs sample points
      lh_rv_clipped,  & ! rv generated from silhs sample points
      lh_Nc_clipped     ! Nc generated from silhs sample points

    type(hydromet_pdf_parameter), dimension(ngrdcol,gr%nzt), optional, intent(out) :: &
      hydromet_pdf_params    ! Hydrometeor PDF parameters        [units vary]

    ! Local Variables
    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt) :: &
      Nc_in_cloud_flip,       & ! Mean (in-cloud) cloud droplet conc.       [num/kg]
      cloud_frac_flip,        & ! Cloud fraction                            [-]
      ice_supersat_frac_flip, & ! Ice supersaturation fraction              [-]
      rho_ds_zt_flip,         & ! Dry, base-state density on thermo. levs.  [kg/m^3]
      Lscale_flip               ! Turbulent Mixing Length                   [m]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzm) :: &
      Kh_zm_flip                ! Eddy diffusivity coef. on momentum levels [m^2/s]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt,hydromet_dim) :: &
      hydromet_flip       ! Mean of hydrometeor, hm (overall) (t-levs.) [units]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzm,hydromet_dim) :: &
      wphydrometp_flip    ! Covariance < w'h_m' > (momentum levels)     [(m/s)units]

    type(pdf_parameter) :: &
      pdf_params_flip    ! PDF parameters                               [units vary]

    type (stats_metadata_type) :: &
      stats_metadata_flip

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzm,hydromet_dim) :: &
      hydrometp2_flip    ! Variance of a hydrometeor (overall) (m-levs.)   [units^2]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt,pdf_dim) :: &
      mu_x_1_n_flip,    & ! Mean array (normal space): PDF vars. (comp. 1) [un. vary]
      mu_x_2_n_flip,    & ! Mean array (normal space): PDF vars. (comp. 2) [un. vary]
      sigma_x_1_n_flip, & ! Std. dev. array (normal space): PDF vars (comp. 1) [u.v.]
      sigma_x_2_n_flip    ! Std. dev. array (normal space): PDF vars (comp. 2) [u.v.]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt,pdf_dim,pdf_dim) :: &
      corr_array_1_n_flip, & ! Corr. array (normal space) of PDF vars. (comp. 1)  [-]
      corr_array_2_n_flip    ! Corr. array (normal space) of PDF vars. (comp. 2)  [-]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt,pdf_dim,pdf_dim) :: &
      corr_cholesky_mtx_1_flip, & ! Transposed corr. cholesky matrix, 1st comp. [-]
      corr_cholesky_mtx_2_flip    ! Transposed corr. cholesky matrix, 2nd comp. [-]

    type(precipitation_fractions) :: &
      precip_fracs_flip           ! Precipitation fractions      [-]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt,hydromet_dim) :: &
      wp2hmp_flip,     & ! Higher-order mixed moment:  < w'^2 hm' > [(m/s)^2<hm un.>]
      rtphmp_zt_flip,  & ! Covariance of rt and hm (on t-levs.)     [(kg/kg)<hm un.>]
      thlphmp_zt_flip    ! Covariance of thl and hm (on t-levs.)    [K<hm units>]

    real( kind = core_rknd ), dimension(ngrdcol,lh_num_samples,gr%nzt,pdf_dim) :: &
      X_nl_all_levs_flip ! Sample that is transformed ultimately to normal-lognormal

    integer, dimension(ngrdcol,lh_num_samples,gr%nzt) :: &
      X_mixt_comp_all_levs_flip ! Which mixture component we're in

    real( kind = core_rknd ), dimension(ngrdcol,lh_num_samples,gr%nzt) :: &
      lh_sample_point_weights_flip

    real( kind = core_rknd ), dimension(ngrdcol,lh_num_samples,gr%nzt) :: &
      lh_rt_clipped_flip,  & ! rt generated from silhs sample points
      lh_thl_clipped_flip, & ! thl generated from silhs sample points
      lh_rc_clipped_flip,  & ! rc generated from silhs sample points
      lh_rv_clipped_flip,  & ! rv generated from silhs sample points
      lh_Nc_clipped_flip     ! Nc generated from silhs sample points

    type(hydromet_pdf_parameter), dimension(ngrdcol,gr%nzt) :: &
      hydromet_pdf_params_flip    ! Hydrometeor PDF parameters        [units vary]

    integer :: k, i, hm_idx, hm_indx, pdf_idx, pdf_indx, sample  ! Loop indices

    logical :: l_differences = .false.

    !------------------------- Begin Code -------------------------


      ! Allocate space and initialize flipped pdf parameter terms
      call init_pdf_params( gr%nzt, ngrdcol, pdf_params_flip )

      ! Allocate space and initialize flipped precipitation fraction terms
      call init_precip_fracs_api( gr%nzt, ngrdcol, &
                                  precip_fracs_flip )
                                    
      ! Initialize hydromet_pdf_params_flip to 0.
      do k = 1, gr%nzt
        do i = 1, ngrdcol
          call init_hydromet_pdf_params( hydromet_pdf_params_flip(i,k) )
        end do
      end do


      Nc_in_cloud_flip(:,:) = Nc_in_cloud(:,gr%nzt:1:-1) 
      cloud_frac_flip(:,:) = cloud_frac(:,gr%nzt:1:-1) 
      ice_supersat_frac_flip(:,:) = ice_supersat_frac(:,gr%nzt:1:-1) 
      rho_ds_zt_flip(:,:) = rho_ds_zt(:,gr%nzt:1:-1) 
      Lscale_flip(:,:) = Lscale(:,gr%nzt:1:-1) 
      Kh_zm_flip(:,:) = Kh_zm(:,gr%nzm:1:-1) 

      ! Hydrometeor variables
      if ( hydromet_dim > 0 ) then
        hydromet_flip(:,:,:) = hydromet(:,gr%nzt:1:-1,:)
        wphydrometp_flip(:,:,:) = wphydrometp(:,gr%nzm:1:-1,:)
      end if

      ! pdf_params
      pdf_params_flip%w_1(:,:) = pdf_params%w_1(:,gr%nzt:1:-1)
      pdf_params_flip%w_2(:,:) = pdf_params%w_2(:,gr%nzt:1:-1)
      pdf_params_flip%varnce_w_1(:,:) = pdf_params%varnce_w_1(:,gr%nzt:1:-1)
      pdf_params_flip%varnce_w_2(:,:) = pdf_params%varnce_w_2(:,gr%nzt:1:-1)
      pdf_params_flip%rt_1(:,:) = pdf_params%rt_1(:,gr%nzt:1:-1)
      pdf_params_flip%rt_2(:,:) = pdf_params%rt_2(:,gr%nzt:1:-1)
      pdf_params_flip%varnce_rt_1(:,:) = pdf_params%varnce_rt_1(:,gr%nzt:1:-1)
      pdf_params_flip%varnce_rt_2(:,:) = pdf_params%varnce_rt_2(:,gr%nzt:1:-1)
      pdf_params_flip%thl_1(:,:) = pdf_params%thl_1(:,gr%nzt:1:-1)
      pdf_params_flip%thl_2(:,:) = pdf_params%thl_2(:,gr%nzt:1:-1)
      pdf_params_flip%varnce_thl_1(:,:) = pdf_params%varnce_thl_1(:,gr%nzt:1:-1)
      pdf_params_flip%varnce_thl_2(:,:) = pdf_params%varnce_thl_2(:,gr%nzt:1:-1)
      pdf_params_flip%corr_w_rt_1(:,:) = pdf_params%corr_w_rt_1(:,gr%nzt:1:-1)
      pdf_params_flip%corr_w_rt_2(:,:) = pdf_params%corr_w_rt_2(:,gr%nzt:1:-1)
      pdf_params_flip%corr_w_thl_1(:,:) = pdf_params%corr_w_thl_1(:,gr%nzt:1:-1)
      pdf_params_flip%corr_w_thl_2(:,:) = pdf_params%corr_w_thl_2(:,gr%nzt:1:-1)
      pdf_params_flip%corr_rt_thl_1(:,:) = pdf_params%corr_rt_thl_1(:,gr%nzt:1:-1)
      pdf_params_flip%corr_rt_thl_2(:,:) = pdf_params%corr_rt_thl_2(:,gr%nzt:1:-1)
      pdf_params_flip%alpha_thl(:,:) = pdf_params%alpha_thl(:,gr%nzt:1:-1)
      pdf_params_flip%alpha_rt(:,:) = pdf_params%alpha_rt(:,gr%nzt:1:-1)
      pdf_params_flip%crt_1(:,:) = pdf_params%crt_1(:,gr%nzt:1:-1)
      pdf_params_flip%crt_2(:,:) = pdf_params%crt_2(:,gr%nzt:1:-1)
      pdf_params_flip%cthl_1(:,:) = pdf_params%cthl_1(:,gr%nzt:1:-1)
      pdf_params_flip%cthl_2(:,:) = pdf_params%cthl_2(:,gr%nzt:1:-1)
      pdf_params_flip%chi_1(:,:) = pdf_params%chi_1(:,gr%nzt:1:-1)
      pdf_params_flip%chi_2(:,:) = pdf_params%chi_2(:,gr%nzt:1:-1)
      pdf_params_flip%stdev_chi_1(:,:) = pdf_params%stdev_chi_1(:,gr%nzt:1:-1)
      pdf_params_flip%stdev_chi_2(:,:) = pdf_params%stdev_chi_2(:,gr%nzt:1:-1)
      pdf_params_flip%stdev_eta_1(:,:) = pdf_params%stdev_eta_1(:,gr%nzt:1:-1)
      pdf_params_flip%stdev_eta_2(:,:) = pdf_params%stdev_eta_2(:,gr%nzt:1:-1)
      pdf_params_flip%covar_chi_eta_1(:,:) = pdf_params%covar_chi_eta_1(:,gr%nzt:1:-1)
      pdf_params_flip%covar_chi_eta_2(:,:) = pdf_params%covar_chi_eta_2(:,gr%nzt:1:-1)
      pdf_params_flip%corr_w_chi_1(:,:) = pdf_params%corr_w_chi_1(:,gr%nzt:1:-1)
      pdf_params_flip%corr_w_chi_2(:,:) = pdf_params%corr_w_chi_2(:,gr%nzt:1:-1)
      pdf_params_flip%corr_w_eta_1(:,:) = pdf_params%corr_w_eta_1(:,gr%nzt:1:-1)
      pdf_params_flip%corr_w_eta_2(:,:) = pdf_params%corr_w_eta_2(:,gr%nzt:1:-1)
      pdf_params_flip%corr_chi_eta_1(:,:) = pdf_params%corr_chi_eta_1(:,gr%nzt:1:-1)
      pdf_params_flip%corr_chi_eta_2(:,:) = pdf_params%corr_chi_eta_2(:,gr%nzt:1:-1)
      pdf_params_flip%rsatl_1(:,:) = pdf_params%rsatl_1(:,gr%nzt:1:-1)
      pdf_params_flip%rsatl_2(:,:) = pdf_params%rsatl_2(:,gr%nzt:1:-1)
      pdf_params_flip%rc_1(:,:) = pdf_params%rc_1(:,gr%nzt:1:-1)
      pdf_params_flip%rc_2(:,:) = pdf_params%rc_2(:,gr%nzt:1:-1)
      pdf_params_flip%cloud_frac_1(:,:) = pdf_params%cloud_frac_1(:,gr%nzt:1:-1)
      pdf_params_flip%cloud_frac_2(:,:) = pdf_params%cloud_frac_2(:,gr%nzt:1:-1)
      pdf_params_flip%mixt_frac(:,:) = pdf_params%mixt_frac(:,gr%nzt:1:-1)
      pdf_params_flip%ice_supersat_frac_1(:,:) = pdf_params%ice_supersat_frac_1(:,gr%nzt:1:-1)
      pdf_params_flip%ice_supersat_frac_2(:,:) = pdf_params%ice_supersat_frac_2(:,gr%nzt:1:-1)

      ! Set statistical sampling to false for descending grid to avoid
      ! the error from having too many statistical samples.
      stats_metadata_flip%l_stats_samp = .false.


      ! Call pdf_hydromet_microphys_prep for the ascending grid direction
      call pdf_hydromet_microphys_prep &
           ( gr, ngrdcol, pdf_dim, hydromet_dim,              & ! In
             itime, dt_main, vert_decorr_coef,                & ! In
             Nc_in_cloud, cloud_frac, ice_supersat_frac,      & ! In
             rho_ds_zt, Lscale, Kh_zm, hydromet, wphydrometp, & ! In
             corr_array_n_cloud, corr_array_n_below,          & ! In
             hm_metadata, pdf_params, clubb_params,           & ! In
             clubb_config_flags, silhs_config_flags,          & ! In
             l_rad_itime, stats_metadata,                     & ! In
             stats_zt, stats_zm, stats_sfc,                   & ! In/Out
             stats_lh_zt, stats_lh_sfc, err_info,             & ! In/Out
             time_clubb_pdf, time_stop, time_start,           & ! In/Out
             hydrometp2,                                      & ! Out
             mu_x_1_n, mu_x_2_n,                              & ! Out
             sigma_x_1_n, sigma_x_2_n,                        & ! Out
             corr_array_1_n, corr_array_2_n,                  & ! Out
             corr_cholesky_mtx_1, corr_cholesky_mtx_2,        & ! Out
             precip_fracs,                                    & ! In/Out
             rtphmp_zt, thlphmp_zt, wp2hmp,                   & ! Out
             X_nl_all_levs, X_mixt_comp_all_levs,             & ! Out
             lh_sample_point_weights,                         & ! Out
             lh_rt_clipped, lh_thl_clipped, lh_rc_clipped,    & ! Out
             lh_rv_clipped, lh_Nc_clipped,                    & ! Out
             hydromet_pdf_params )                              ! Optional(out)


      ! In the case of a fatal error during the 1st call (ascending grid) to 
      ! pdf_hydromet_microphys_prep, reset the error code.
      !
      ! Otherwise, a mismatch between the ascending grid and the descending grid
      ! will be caused by the fact that the 2nd call to pdf_hydromet_microphys_prep
      ! will be returned from right away because err_code is already set to 
      ! clubb_fatal_error upon entering pdf_hydromet_microphys_prep. Thus, there
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
      if ( any(err_info%err_code == clubb_fatal_error) ) then
         ! Reset error code
         err_info%err_code = clubb_no_error
      endif


      ! Call pdf_hydromet_microphys_prep for the descending grid direction
      ! All variables with a vertical dimension should be "flip" variables
      ! in this call.
      call pdf_hydromet_microphys_prep &
           ( gr_desc, ngrdcol, pdf_dim, hydromet_dim,         & ! In
             itime, dt_main, vert_decorr_coef,                & ! In
             Nc_in_cloud_flip, cloud_frac_flip, ice_supersat_frac_flip, & ! In
             rho_ds_zt_flip, Lscale_flip, Kh_zm_flip, hydromet_flip, wphydrometp_flip, & ! In
             corr_array_n_cloud, corr_array_n_below,          & ! In
             hm_metadata, pdf_params_flip, clubb_params,      & ! In
             clubb_config_flags, silhs_config_flags,          & ! In
             l_rad_itime, stats_metadata_flip,                & ! In
             stats_zt, stats_zm, stats_sfc,                   & ! In/Out
             stats_lh_zt, stats_lh_sfc, err_info,             & ! In/Out
             time_clubb_pdf, time_stop, time_start,           & ! In/Out
             hydrometp2_flip,                                 & ! Out
             mu_x_1_n_flip, mu_x_2_n_flip,                    & ! Out
             sigma_x_1_n_flip, sigma_x_2_n_flip,              & ! Out
             corr_array_1_n_flip, corr_array_2_n_flip,        & ! Out
             corr_cholesky_mtx_1_flip, corr_cholesky_mtx_2_flip, & ! Out
             precip_fracs_flip,                               & ! In/Out
             rtphmp_zt_flip, thlphmp_zt_flip, wp2hmp_flip,    & ! Out
             X_nl_all_levs_flip, X_mixt_comp_all_levs_flip,   & ! Out
             lh_sample_point_weights_flip,                    & ! Out
             lh_rt_clipped_flip, lh_thl_clipped_flip, lh_rc_clipped_flip,    & ! Out
             lh_rv_clipped_flip, lh_Nc_clipped_flip,          & ! Out
             hydromet_pdf_params_flip )                         ! Optional(out)


      ! Compare the ascending grid variables to the descending grid variables
      if ( .not. trim( microphys_scheme ) == "none" ) then

         ! The variables that are output from the calls to setup_pdf_parameters
         ! and hydrometeor_mixed_moments only have assigned values when a
         ! microphysics scheme is in use. Thus, only run this check when that
         ! is the case.
         do pdf_idx = 1, pdf_dim
            ! mu_x_1_n
            call check_flipped_results( "mu_x_1_n", mu_x_1_n(:,:,pdf_idx), &
                                        mu_x_1_n_flip(:,:,pdf_idx), &
                                        gr%nzt, ngrdcol, &
                                        l_differences )
            ! mu_x_2_n
            call check_flipped_results( "mu_x_2_n", mu_x_2_n(:,:,pdf_idx), &
                                        mu_x_2_n_flip(:,:,pdf_idx), &
                                        gr%nzt, ngrdcol, &
                                        l_differences )
            ! sigma_x_1_n
            call check_flipped_results( "sigma_x_1_n", sigma_x_1_n(:,:,pdf_idx), &
                                        sigma_x_1_n_flip(:,:,pdf_idx), &
                                        gr%nzt, ngrdcol, &
                                        l_differences )
            ! sigma_x_2_n
            call check_flipped_results( "sigma_x_2_n", sigma_x_2_n(:,:,pdf_idx), &
                                        sigma_x_2_n_flip(:,:,pdf_idx), &
                                        gr%nzt, ngrdcol, &
                                        l_differences )
            do pdf_indx = 1, pdf_dim
               ! corr_array_1_n
               call check_flipped_results( "corr_array_1_n", &
                                           corr_array_1_n(:,:,pdf_idx,pdf_indx), &
                                           corr_array_1_n_flip(:,:,pdf_idx,pdf_indx), &
                                           gr%nzt, ngrdcol, &
                                           l_differences )
               ! corr_array_1_n
               call check_flipped_results( "corr_array_2_n", &
                                           corr_array_2_n(:,:,pdf_idx,pdf_indx), &
                                           corr_array_2_n_flip(:,:,pdf_idx,pdf_indx), &
                                           gr%nzt, ngrdcol, &
                                           l_differences )
               ! corr_cholesky_mtx_1
               call check_flipped_results( "corr_cholesky_mtx_1", &
                                           corr_cholesky_mtx_1(:,:,pdf_idx,pdf_indx), &
                                           corr_cholesky_mtx_1_flip(:,:,pdf_idx,pdf_indx), &
                                           gr%nzt, ngrdcol, &
                                           l_differences )
               ! corr_cholesky_mtx_2
               call check_flipped_results( "corr_cholesky_mtx_2", &
                                           corr_cholesky_mtx_2(:,:,pdf_idx,pdf_indx), &
                                           corr_cholesky_mtx_2_flip(:,:,pdf_idx,pdf_indx), &
                                           gr%nzt, ngrdcol, &
                                           l_differences )
            enddo ! pdf_indx = 1, pdf_dim
         enddo ! pdf_idx = 1, pdf_dim

         do hm_idx = 1, hydromet_dim
            ! hydrometp2
            call check_flipped_results( "hydrometp2", hydrometp2(:,:,hm_idx), &
                                        hydrometp2_flip(:,:,hm_idx), &
                                        gr%nzm, ngrdcol, &
                                        l_differences )
         enddo ! hm_idx = 1, hydromet_dim

         ! precip_fracs
         call check_flipped_results( "precip_fracs%precip_frac", precip_fracs%precip_frac, &
                                     precip_fracs_flip%precip_frac, gr%nzt, ngrdcol, &
                                     l_differences )
         call check_flipped_results( "precip_fracs%precip_frac_1", precip_fracs%precip_frac_1, &
                                     precip_fracs_flip%precip_frac_1, gr%nzt, ngrdcol, &
                                     l_differences )
         call check_flipped_results( "precip_fracs%precip_frac_2", precip_fracs%precip_frac_2, &
                                     precip_fracs_flip%precip_frac_2, gr%nzt, ngrdcol, &
                                     l_differences )

         ! hydromet_pdf_params
         do hm_idx = 1, hydromet_dim
            call check_flipped_results( "hydromet_pdf_params%hm_1", &
                                        hydromet_pdf_params%hm_1(hm_idx), &
                                        hydromet_pdf_params_flip%hm_1(hm_idx), &
                                        gr%nzt, ngrdcol, &
                                        l_differences )
            call check_flipped_results( "hydromet_pdf_params%hm_2", &
                                        hydromet_pdf_params%hm_2(hm_idx), &
                                        hydromet_pdf_params_flip%hm_2(hm_idx), &
                                        gr%nzt, ngrdcol, &
                                        l_differences )
            call check_flipped_results( "hydromet_pdf_params%mu_hm_1", &
                                        hydromet_pdf_params%mu_hm_1(hm_idx), &
                                        hydromet_pdf_params_flip%mu_hm_1(hm_idx), &
                                        gr%nzt, ngrdcol, &
                                        l_differences )
            call check_flipped_results( "hydromet_pdf_params%mu_hm_2", &
                                        hydromet_pdf_params%mu_hm_2(hm_idx), &
                                        hydromet_pdf_params_flip%mu_hm_2(hm_idx), &
                                        gr%nzt, ngrdcol, &
                                        l_differences )
            call check_flipped_results( "hydromet_pdf_params%sigma_hm_1", &
                                        hydromet_pdf_params%sigma_hm_1(hm_idx), &
                                        hydromet_pdf_params_flip%sigma_hm_1(hm_idx), &
                                        gr%nzt, ngrdcol, &
                                        l_differences )
            call check_flipped_results( "hydromet_pdf_params%sigma_hm_2", &
                                        hydromet_pdf_params%sigma_hm_2(hm_idx), &
                                        hydromet_pdf_params_flip%sigma_hm_2(hm_idx), &
                                        gr%nzt, ngrdcol, &
                                        l_differences )
            call check_flipped_results( "hydromet_pdf_params%corr_w_hm_1", &
                                        hydromet_pdf_params%corr_w_hm_1(hm_idx), &
                                        hydromet_pdf_params_flip%corr_w_hm_1(hm_idx), &
                                        gr%nzt, ngrdcol, &
                                        l_differences )
            call check_flipped_results( "hydromet_pdf_params%corr_w_hm_2", &
                                        hydromet_pdf_params%corr_w_hm_2(hm_idx), &
                                        hydromet_pdf_params_flip%corr_w_hm_2(hm_idx), &
                                        gr%nzt, ngrdcol, &
                                        l_differences )
            call check_flipped_results( "hydromet_pdf_params%corr_chi_hm_1", &
                                        hydromet_pdf_params%corr_chi_hm_1(hm_idx), &
                                        hydromet_pdf_params_flip%corr_chi_hm_1(hm_idx), &
                                        gr%nzt, ngrdcol, &
                                        l_differences )
            call check_flipped_results( "hydromet_pdf_params%corr_chi_hm_2", &
                                        hydromet_pdf_params%corr_chi_hm_2(hm_idx), &
                                        hydromet_pdf_params_flip%corr_chi_hm_2(hm_idx), &
                                        gr%nzt, ngrdcol, &
                                        l_differences )
            call check_flipped_results( "hydromet_pdf_params%corr_eta_hm_1", &
                                        hydromet_pdf_params%corr_eta_hm_1(hm_idx), &
                                        hydromet_pdf_params_flip%corr_eta_hm_1(hm_idx), &
                                        gr%nzt, ngrdcol, &
                                        l_differences )
            call check_flipped_results( "hydromet_pdf_params%corr_eta_hm_2", &
                                        hydromet_pdf_params%corr_eta_hm_2(hm_idx), &
                                        hydromet_pdf_params_flip%corr_eta_hm_2(hm_idx), &
                                        gr%nzt, ngrdcol, &
                                        l_differences )
            do hm_indx = 1, hydromet_dim
               call check_flipped_results( "hydromet_pdf_params%corr_hmx_hmy_1", &
                                           hydromet_pdf_params%corr_hmx_hmy_1(hm_idx,hm_indx), &
                                           hydromet_pdf_params_flip%corr_hmx_hmy_1(hm_idx,hm_indx),&
                                           gr%nzt, ngrdcol, &
                                           l_differences )
               call check_flipped_results( "hydromet_pdf_params%corr_hmx_hmy_2", &
                                           hydromet_pdf_params%corr_hmx_hmy_2(hm_idx,hm_indx), &
                                           hydromet_pdf_params_flip%corr_hmx_hmy_2(hm_idx,hm_indx),&
                                           gr%nzt, ngrdcol, &
                                           l_differences )
            enddo ! hm_indx = 1, hydromet_dim
         enddo ! hm_idx = 1, hydromet_dim
         call check_flipped_results( "hydromet_pdf_params%mu_Ncn_1", &
                                     hydromet_pdf_params%mu_Ncn_1, &
                                     hydromet_pdf_params_flip%mu_Ncn_1, &
                                     gr%nzt, ngrdcol, &
                                     l_differences )
         call check_flipped_results( "hydromet_pdf_params%mu_Ncn_2", &
                                     hydromet_pdf_params%mu_Ncn_2, &
                                     hydromet_pdf_params_flip%mu_Ncn_2, &
                                     gr%nzt, ngrdcol, &
                                     l_differences )
         call check_flipped_results( "hydromet_pdf_params%sigma_Ncn_1", &
                                     hydromet_pdf_params%sigma_Ncn_1, &
                                     hydromet_pdf_params_flip%sigma_Ncn_1, &
                                     gr%nzt, ngrdcol, &
                                     l_differences )
         call check_flipped_results( "hydromet_pdf_params%sigma_Ncn_2", &
                                     hydromet_pdf_params%sigma_Ncn_2, &
                                     hydromet_pdf_params_flip%sigma_Ncn_2, &
                                     gr%nzt, ngrdcol, &
                                     l_differences )

         do hm_idx = 1, hydromet_dim
            ! rtphmp_zt
            call check_flipped_results( "rtphmp_zt", rtphmp_zt(:,:,hm_idx), &
                                        rtphmp_zt_flip(:,:,hm_idx), &
                                        gr%nzt, ngrdcol, &
                                        l_differences )
            ! thlphmp_zt
            call check_flipped_results( "thlphmp_zt", thlphmp_zt(:,:,hm_idx), &
                                        thlphmp_zt_flip(:,:,hm_idx), &
                                        gr%nzt, ngrdcol, &
                                        l_differences )
            ! wp2hmp
            call check_flipped_results( "wp2hmp", wp2hmp(:,:,hm_idx), &
                                        wp2hmp_flip(:,:,hm_idx), &
                                        gr%nzt, ngrdcol, &
                                        l_differences )
         enddo ! hm_idx = 1, hydromet_dim

      endif ! .not. trim( microphys_scheme ) == "none"

#ifdef SILHS
      if ( lh_microphys_type /= lh_microphys_disabled .or. l_silhs_rad ) then

         ! The variables that are output from the calls to
         ! generate_silhs_sample_api and clip_transform_silhs_output_api only
         ! have assigned values when SILHS is in use. Thus, only run this check
         ! when that is the case.
         do sample = 1, lh_num_samples
            ! X_nl_all_levs
            do pdf_idx = 1, pdf_dim
               call check_flipped_results( "X_nl_all_levs", &
                                           X_nl_all_levs(:,sample,:,pdf_idx), &
                                           X_nl_all_levs_flip(:,sample,:,pdf_idx), &
                                           gr%nzt, ngrdcol, &
                                           l_differences )
            enddo ! pdf_idx = 1, pdf_dim
            ! X_mixt_comp_all_levs
            call check_flipped_results( "X_mixt_comp_all_levs", &
                                        real( X_mixt_comp_all_levs(:,sample,:), &
                                              kind = core_rknd ), &
                                        real( X_mixt_comp_all_levs_flip(:,sample,:), &
                                              kind = core_rknd ), &
                                        gr%nzt, ngrdcol, &
                                        l_differences )
            ! lh_sample_point_weights
            call check_flipped_results( "lh_sample_point_weights", &
                                        lh_sample_point_weights(:,sample,:), &
                                        lh_sample_point_weights_flip(:,sample,:), &
                                        gr%nzt, ngrdcol, &
                                        l_differences )
            ! lh_rt_clipped
            call check_flipped_results( "lh_rt_clipped", lh_rt_clipped(:,sample,:), &
                                        lh_rt_clipped_flip(:,sample,:), &
                                        gr%nzt, ngrdcol, &
                                        l_differences )
            ! lh_thl_clipped
            call check_flipped_results( "lh_thl_clipped", lh_thl_clipped(:,sample,:), &
                                        lh_thl_clipped_flip(:,sample,:), &
                                        gr%nzt, ngrdcol, &
                                        l_differences )
            ! lh_rc_clipped
            call check_flipped_results( "lh_rc_clipped", lh_rc_clipped(:,sample,:), &
                                        lh_rc_clipped_flip(:,sample,:), &
                                        gr%nzt, ngrdcol, &
                                        l_differences )
            ! lh_rv_clipped
            call check_flipped_results( "lh_rv_clipped", lh_rv_clipped(:,sample,:), &
                                        lh_rv_clipped_flip(:,sample,:), &
                                        gr%nzt, ngrdcol, &
                                        l_differences )
            ! lh_Nc_clipped
            call check_flipped_results( "lh_Nc_clipped", lh_Nc_clipped(:,sample,:), &
                                        lh_Nc_clipped_flip(:,sample,:), &
                                        gr%nzt, ngrdcol, &
                                        l_differences )
         enddo ! sample = 1, lh_num_samples

      endif ! lh_microphys_type /= lh_microphys_disabled .or. l_silhs_rad
#endif /*SILHS*/
      

      ! Print a message and stop the run if there are any discrepanices found
      if ( l_differences ) then
         ! Stop the run and exit
         print *, "##################################################"
         print *, "Discrepancy found in SILHS ascending vs. descending grid" &
                  // " direction test. Please see messages listed above."
         err_info%err_code = clubb_generalized_grd_test_err
      endif ! l_differences


      return

  end subroutine silhs_generalized_grid_testing

  !=============================================================================
  subroutine check_flipped_results( varname, var, var_flip, nz, ngrdcol, &
                                    l_differences )

    use constants_clubb, only: &
        zero, &
        eps

    use grid_class, only: &
        flip

    use clubb_precision, only: &
        core_rknd

    use model_flags, only: &
        l_force_descending_solves

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

    ! Local Variables
    integer :: k, i  ! Loop indices

    real( kind = core_rknd ) :: &
      tolerance  ! tolerance used to determine if results match

    !--------------- Begin Code -------------------

    if ( l_force_descending_solves ) then
      ! If we are forcing the matrix solves to be done in descending grid mode, 
      ! then results should be BFB and we can use tolerance 0.0.
      tolerance = zero
    else
      ! If we do not force matrix solves to be done in a specific direction, then
      ! results are not expected to be BFB, so we need a non-zero tolerance, but 
      ! they should be close, so we want a pretty small tolerance
      tolerance = max( 1.e-8_core_rknd, epsilon(tolerance) )    ! max statement for single precision runs
    end if

    if ( any( abs( var - var_flip(:,nz:1:-1) ) > tolerance ) ) then
      l_differences = .true.

      print *, "**************************************************"
      print *, "Differences found in ", trim( varname )

      do i = 1, ngrdcol

        print *, " - column ", i, ": "

        do k = 1, nz, 1
          if ( abs( var(i,k) - var_flip(i,nz-k+1) ) > tolerance ) then
            !print *, " ascending/descending level = ", k, " / ", nz-k+1
            print *, "ascending @ k =", k, ":", var(i,k), "descending @ k =", &
                      nz-k+1, ":", var_flip(i,nz-k+1), " difference = ", &
                      abs(var(i,k)-var_flip(i,nz-k+1))
          endif
        enddo
      enddo

    endif

    return

  end subroutine check_flipped_results

!===============================================================================

end module generalized_grid_test
