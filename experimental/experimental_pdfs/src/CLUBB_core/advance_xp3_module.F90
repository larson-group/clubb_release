! $Id$
!===============================================================================
module advance_xp3_module

  ! Description:
  ! Predicts the value of <x'^3> for <rt'^3>, <thl'^3>, and <sclr'^3>.

  ! References:
  !-------------------------------------------------------------------------

  use stats_netcdf, only: &
      stats_type, &
      stats_update, &
      stats_begin_budget, &
      stats_finalize_budget

  implicit none

  public :: advance_xp3, & ! Procedure(s)
            diagnose_xp3

  private :: advance_xp3_simplified, & ! Procedure(s)
             advance_xp3_trivar, &
             xp3_trivar_lhs, &
             xp3_trivar_solve, &
             term_ta_rhs, &
             term_tp_rhs, &
             term_ac_rhs

  private ! default scope

  integer, parameter, private :: &
    xp3_rtp3 = 1,   & ! Named constant for solving rtp3
    xp3_thlp3 = 2,  & ! Named constant for solving thlp3
    xp3_sclrp3 = 3    ! Named constant for solving sclrp3

  contains

  !=============================================================================
  subroutine advance_xp3( nzm, nzt, ngrdcol, sclr_dim, sclr_tol, gr, dt,   & ! Intent(in)
                          rtm, thlm, rtp2, thlp2, wprtp,                   & ! Intent(in)
                          wpthlp, wprtp2, wpthlp2, rho_ds_zm,              & ! Intent(in)
                          invrs_rho_ds_zt, invrs_tau_zt, tau_max_zt,       & ! Intent(in)
                          sclrm, sclrp2, wpsclrp, wpsclrp2,                & ! Intent(in)
                          wp2, wp3, upwp, vpwp, up2, vp2,                  & ! Intent(in)
                          thvm, wm_zt, Kh_zm, Kh_zt, pdf_params,           & ! Intent(in)
                          clubb_params, nu_vert_res_dep, iiPDF_type,       & ! Intent(in)
                          tridiag_solve_method, l_upwind_xm_ma,             & ! Intent(in)
                          l_implemented, err_info,                          & ! Intent(inout)
                          l_lmm_stepping,                                  & ! Intent(in)
                          stats,                                           & ! Intent(inout)
                          rtp3, thlp3, sclrp3, up3, vp3 )                    ! Intent(inout)

    ! Description:
    ! Advance <rt'^3>, <thl'^3>, and <sclr'^3> one model timestep using a
    ! simplified form of the <x'^3> predictive equation.  The simplified <x'^3>
    ! equation can either be advanced from its previous value or calculated
    ! using a steady-state approximation.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        grid,       & ! Type
        zm2zt_api,  & ! Procedure(s)
        zt2zm_api,  &
        ddzm

    use pdf_parameter_module, only: &
        pdf_parameter

    use parameters_tunable, only: &
        nu_vertical_res_dep

    use err_info_type_module, only: &
        err_info_type

    use model_flags, only: &
        iiPDF_trivariate_transport

    use constants_clubb, only: &
        rt_tol,         & ! Variable(s)
        thl_tol,        &
        w_tol,          &
        w_tol_sqd,      &
        zero_threshold, &
        grav,           &
        zero,           &
        one

    use parameter_indices, only: &
        nparams,          &
        ixp3_coef_base,   &
        ixp3_coef_slope

    use Skx_module, only: &
        Skx_func,           & ! Procedure(s)
        xp3_LG_2005_ansatz

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! --------------------- Input Variables ---------------------
    integer, intent(in) :: &
      nzm,          & ! Number of momentum vertical levels
      nzt,          & ! Number of thermodynamic vertical levels
      ngrdcol,      & ! Number of grid columns
      sclr_dim        ! Number of passive scalars

    real( kind = core_rknd ), intent(in), dimension(sclr_dim) :: &
      sclr_tol          ! Threshold(s) on the passive scalars  [units vary]

    type (grid), intent(in) :: &
      gr

    real( kind = core_rknd ), intent(in) :: &
      dt                 ! Model timestep                            [s]

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      rtm,             & ! Mean (overall) of rt (thermo. levels)                [kg/kg]
      thlm,            & ! Mean (overall) of thl (thermo. levels)               [K]
      wprtp2,          & ! <w'rt'^2> (thermodynamic levels)                     [m/s(kg/kg)^2]
      wpthlp2,         & ! <w'thl'^2> (thermodynamic levels)                    [m/s K^2]
      invrs_rho_ds_zt, & ! Inv. dry, static density at thermo. levels           [m^3/kg]
      invrs_tau_zt,    & ! Inverse time-scale tau on thermodynamic levels       [1/s]
      tau_max_zt         ! Max. allowable eddy dissipation time scale on t-levs [s]

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      wp3,    & ! w'^3 on thermo. grid     [m^3/s^3]
      thvm,     & ! Virtual potential temperature on thermodynamic levels [K]
      wm_zt       ! Mean vertical velocity on thermodynamic levels [m/s]

    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: &
      rtp2,            & ! Variance (overall) of rt (m-levs.)                   [kg^2/kg^2]
      thlp2,           & ! Variance (overall) of thl (m-levs.)                  [K^2]
      wprtp,           & ! Turbulent flux of rt (momentum levs.)                [m/s kg/kg]
      wpthlp,          & ! Turbulent flux of thl (momentum levs.)               [m/s K]
      rho_ds_zm,        & ! Dry, static density on momentum levels               [kg/m^3]
      upwp,             & ! u'w' (momentum levels)                              [m^2/s^2]
      vpwp,             & ! v'w' (momentum levels)                              [m^2/s^2]
      wp2,              & ! w'^2 (momentum levels)                              [m^2/s^2]
      up2,              & ! u'^2 (momentum levels)                              [m^2/s^2]
      vp2                 ! v'^2 (momentum levels)                              [m^2/s^2]

    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: &
      Kh_zm,              & ! Eddy diffusivity on momentum levels [m^2/s]
      Kh_zt                 ! Eddy diffusivity on thermodynamic levels [m^2/s]

    real( kind = core_rknd ), dimension(ngrdcol,nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

    type(pdf_parameter), intent(in) :: &
      pdf_params       ! PDF parameters on thermodynamic levels [units vary]

    type(nu_vertical_res_dep), intent(in) :: &
      nu_vert_res_dep  ! Resolution-dependent background diffusivities [m^2/s]

    integer, intent(in) :: &
      iiPDF_type,              & ! PDF type selector
      tridiag_solve_method        ! Tridiagonal-solver selector

    logical, intent(in) :: &
      l_upwind_xm_ma, & ! Use upwind mean advection for scalar moments
      l_implemented     ! True when CLUBB is embedded in a host model

    type(err_info_type), intent(inout) :: &
      err_info          ! Structured numerical-solver error information

    real( kind = core_rknd ), dimension(ngrdcol,nzt,sclr_dim), intent(in) :: &
      sclrm,    & ! Mean (overall) of sclr (thermo. levels) [sclr units]
      wpsclrp2    ! <w'sclr'^2> (thermodynamic levels)      [m/s(sclr units)^2]

    real( kind = core_rknd ), dimension(ngrdcol,nzm,sclr_dim), intent(in) :: &
      sclrp2,   & ! Variance (overall) of sclr (m-levs.)    [(sclr units)^2]
      wpsclrp     ! Turbulent flux of sclr (momentum levs.) [m/s(sclr units)]

    logical, intent(in) :: &
      l_lmm_stepping    ! Apply Linear Multistep Method (LMM) Stepping

    type(stats_type), intent(inout) :: &
      stats

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(inout) :: &
      rtp3,  & ! <rt'^3> (thermodynamic levels)     [kg^3/kg^3]
      thlp3, & ! <thl'^3> (thermodynamic levels)    [K^3]
      up3,   & ! u'^3 (thermodynamic levels)        [m^3/s^3]
      vp3      ! v'^3 (thermodynamic levels)        [m^3/s^3]

    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
      thvm_zm                           ! Virtual potential temperature on momentum levs. [K]

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      ddzm_thvm_zm,                   & ! d(thvm_zm)/dz [K/m]
      brunt_vaisala_freq_sqd_zt,      & ! Buoyancy frequency squared on t-levs. [s^-2]
      Skw_zt,                         & ! w skewness on thermodynamic levels [-]
      wp2_zt,                         & ! w'^2 on thermo. grid [m^2/s^2]
      thlp2_zt,                       & ! thl'^2 on thermo. grid [K^2]
      wpthlp_zt,                      & ! w'thl' on thermo. grid [m K/s]
      wprtp_zt,                       & ! w'rt' on thermo. grid [m kg/(kg s)]
      rtp2_zt,                        & ! rt'^2 on thermo. grid [(kg/kg)^2]
      upwp_zt,                        & ! u'w' on thermo. grid [m^2/s^2]
      vpwp_zt,                        & ! v'w' on thermo. grid [m^2/s^2]
      up2_zt,                         & ! u'^2 on thermo. grid [m^2/s^2]
      vp2_zt,                         & ! v'^2 on thermo. grid [m^2/s^2]
      xp3_coef_fnc                      ! Coefficient in simple xp3 equation [-]

    real( kind = core_rknd ), dimension(ngrdcol,nzt,sclr_dim), intent(inout) :: &
      sclrp3    ! <sclr'^3> (thermodynamic levels)    [(sclr units)^3]

    ! --------------------- Local Variable ---------------------
    integer :: i, k, sclr    ! Loop index

    !$acc data create( thvm_zm, ddzm_thvm_zm, brunt_vaisala_freq_sqd_zt, &
    !$acc              Skw_zt, wp2_zt, thlp2_zt, wpthlp_zt, wprtp_zt, rtp2_zt, &
    !$acc              upwp_zt, vpwp_zt, up2_zt, vp2_zt, &
    !$acc              xp3_coef_fnc )

    wp2_zt(:,:) = zm2zt_api( nzm, nzt, ngrdcol, gr, wp2(:,:), w_tol_sqd )

    if ( iiPDF_type == iiPDF_trivariate_transport ) then
      ! trivariate transport PDF has an exact two-Gaussian diagnosis for <w'x'^3>.  Use that
      ! conservative flux while treating mean advection, diffusion, and tau
      ! damping implicitly in one thermodynamic-level tridiagonal solve.
      call advance_xp3_trivar( nzm, nzt, ngrdcol, gr, dt,                    & ! Intent(in)
                               rtm, thlm, rtp2, thlp2, wprtp, wpthlp,       & ! Intent(in)
                               wprtp2, wpthlp2, wm_zt, Kh_zm, Kh_zt, rho_ds_zm, & ! Intent(in)
                               invrs_rho_ds_zt, invrs_tau_zt,                & ! Intent(in)
                               pdf_params, clubb_params, nu_vert_res_dep,   & ! Intent(in)
                               tridiag_solve_method, l_upwind_xm_ma,        & ! Intent(in)
                               l_implemented, l_lmm_stepping,               & ! Intent(in)
                               stats, err_info,                              & ! Intent(inout)
                               rtp3, thlp3 )                                   ! Intent(inout)
    else
      ! Preserve the established inexpensive treatment for non-trivariate transport PDF PDFs.
      call advance_xp3_simplified( nzm, nzt, ngrdcol, gr, xp3_rtp3, dt, & ! Intent(in)
                                   rtm, rtp2, wprtp,                    & ! Intent(in)
                                   wprtp2, rho_ds_zm,                   & ! Intent(in)
                                   invrs_rho_ds_zt,                     & ! Intent(in)
                                   invrs_tau_zt, tau_max_zt,            & ! Intent(in)
                                   rt_tol, l_lmm_stepping,              & ! Intent(in)
                                   stats,                                & ! Intent(inout)
                                   rtp3 )                                  ! Intent(inout)

      call advance_xp3_simplified( nzm, nzt, ngrdcol, gr, xp3_thlp3, dt, & ! Intent(in)
                                   thlm, thlp2, wpthlp,                  & ! Intent(in)
                                   wpthlp2, rho_ds_zm,                   & ! Intent(in)
                                   invrs_rho_ds_zt,                      & ! Intent(in)
                                   invrs_tau_zt, tau_max_zt,             & ! Intent(in)
                                   thl_tol, l_lmm_stepping,              & ! Intent(in)
                                   stats,                                 & ! Intent(inout)
                                   thlp3 )                                  ! Intent(inout)
    end if

    ! Advance <sclr'^3> one model timestep or calculate <sclr'^3> using a
    ! steady-state approximation.
    do sclr = 1, sclr_dim, 1
      call advance_xp3_simplified( nzm, nzt, ngrdcol, gr, xp3_sclrp3, dt,                & ! In
                                   sclrm(:,:,sclr), sclrp2(:,:,sclr), wpsclrp(:,:,sclr), & ! In
                                   wpsclrp2(:,:,sclr), rho_ds_zm,                        & ! In
                                   invrs_rho_ds_zt,                                      & ! In
                                   invrs_tau_zt, tau_max_zt,                             & ! In
                                   sclr_tol(sclr), l_lmm_stepping,                       & ! In
                                   stats,                                                & ! In
                                   sclrp3(:,:,sclr) )                                      ! In/Out
      end do ! sclr = 1, sclr_dim

    ! Use a modified form of the Larson and Golaz (2005) ansatz for the
    ! ADG1 PDF to calculate <u'^3> and <v'^3> for another type of PDF.
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

    thvm_zm(:,:)                   = zt2zm_api( nzm, nzt, ngrdcol, gr, thvm(:,:), &
                                                zero_threshold )
    ddzm_thvm_zm(:,:)              = ddzm( nzm, nzt, ngrdcol, gr, thvm_zm(:,:) )
    brunt_vaisala_freq_sqd_zt(:,:) = max( ( grav / thvm(:,:) ) * ddzm_thvm_zm(:,:), zero )

    ! The xp3_coef_fnc provides some extra tunability to the simple xp3
    ! equation.
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

    if ( stats%l_sample ) then
      !$acc update host( thlp2_zt, wpthlp_zt, wprtp_zt, rtp2_zt, &
      !$acc              up2_zt, vp2_zt, upwp_zt, vpwp_zt )

      call stats_update( "thlp2_zt", thlp2_zt, stats )
      call stats_update( "wpthlp_zt", wpthlp_zt, stats )
      call stats_update( "wprtp_zt", wprtp_zt, stats )
      call stats_update( "rtp2_zt", rtp2_zt, stats )
      call stats_update( "up2_zt", up2_zt, stats )
      call stats_update( "vp2_zt", vp2_zt, stats )
      call stats_update( "upwp_zt", upwp_zt, stats )
      call stats_update( "vpwp_zt", vpwp_zt, stats )
    end if

    !$acc end data

    return

  end subroutine advance_xp3

  !=============================================================================
  subroutine diagnose_xp3( nzm, nzt, ngrdcol, sclr_dim, sclr_tol, gr, &
                          iiPDF_type, clubb_params, &
                          wp2, wp3, thvm, &
                          wprtp, wpthlp, rtp2, thlp2, upwp, vpwp, up2, vp2, &
                          sigma_sqd_w, wpsclrp, sclrp2, &
                          stats, &
                          rtp3, thlp3, up3, vp3, &
                          sclrp3 )

    ! Description:
    !   Diagnose third-order scalar and horizontal-velocity moments from the
    !   current second-order moments and PDF width information.
    !
    ! References:
    !   Larson and Golaz (2005) for the diagnostic ansatz.
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd

    use constants_clubb, only: &
        thl_tol,        &
        rt_tol,         &
        w_tol,          &
        w_tol_sqd,      &
        zero_threshold, &
        one,            &
        zero,           &
        grav

    use grid_class, only: &
        grid,       &
        zm2zt_api,  &
        zt2zm_api,  &
        ddzm

    use model_flags, only: &
        iiPDF_ADG1, &
        iiPDF_parcel_ensemble

    use parameter_indices, only: &
        nparams,          &
        ixp3_coef_base,   &
        ixp3_coef_slope

    use Skx_module, only: &
        Skx_func,           &
        xp3_LG_2005_ansatz

    implicit none

    !--------------------------- Input Variables ---------------------------
    integer, intent(in) :: &
      nzm,       & ! Number of momentum levels
      nzt,       & ! Number of thermodynamic levels
      ngrdcol,   & ! Number of grid columns
      sclr_dim,  & ! Number of passive scalars
      iiPDF_type   ! PDF type selector

    real( kind = core_rknd ), dimension(sclr_dim), intent(in) :: &
      sclr_tol   ! Thresholds for passive scalars [units vary]

    type(grid), intent(in) :: &
      gr   ! Grid structure

    real( kind = core_rknd ), dimension(ngrdcol,nparams), intent(in) :: &
      clubb_params   ! Array of CLUBB tunable parameters [units vary]

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      wp3,  & ! Third moment of vertical velocity [m^3/s^3]
      thvm    ! Virtual potential temperature [K]

    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: &
      wp2,         & ! Variance of vertical velocity [m^2/s^2]
      wprtp,       & ! Turbulent flux of rt [m kg/(kg s)]
      wpthlp,      & ! Turbulent flux of thl [m K/s]
      rtp2,        & ! Variance of rt [(kg/kg)^2]
      thlp2,       & ! Variance of thl [K^2]
      upwp,        & ! Turbulent flux of u [m^2/s^2]
      vpwp,        & ! Turbulent flux of v [m^2/s^2]
      up2,         & ! Variance of u wind [m^2/s^2]
      vp2,         & ! Variance of v wind [m^2/s^2]
      sigma_sqd_w    ! PDF width parameter [-]

    real( kind = core_rknd ), dimension(ngrdcol,nzm,sclr_dim), intent(in) :: &
      wpsclrp, & ! Turbulent flux of passive scalars [m sclr/s]
      sclrp2     ! Variance of passive scalars [units vary]

    !----------------------- Input/Output Variables ------------------------
    type(stats_type), intent(inout) :: &
      stats   ! Statistics accumulator

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(inout) :: &
      rtp3,  & ! Third moment of rt [(kg/kg)^3]
      thlp3, & ! Third moment of thl [K^3]
      up3,   & ! Third moment of u wind [m^3/s^3]
      vp3      ! Third moment of v wind [m^3/s^3]

    real( kind = core_rknd ), dimension(ngrdcol,nzt,sclr_dim), intent(inout) :: &
      sclrp3   ! Third moments of passive scalars [units vary]

    !--------------------------- Local Variables ---------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
      thvm_zm                           ! Virtual potential temperature on momentum levs. [K]

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      ddzm_thvm_zm,                   & ! d(thvm_zm)/dz [K/m]
      brunt_vaisala_freq_sqd_zt,      & ! Buoyancy frequency squared on t-levs. [s^-2]
      wp2_zt,                         & ! Variance of vertical velocity on thermodynamic levels [m^2/s^2]
      Skw_zt,                         & ! w skewness on thermodynamic levels [-]
      wpthlp_zt,                      & ! w'thl' on thermo. grid [m K/s]
      wprtp_zt,                       & ! w'rt' on thermo. grid [m kg/(kg s)]
      thlp2_zt,                       & ! thl'^2 on thermo. grid [K^2]
      rtp2_zt,                        & ! rt'^2 on thermo. grid [(kg/kg)^2]
      upwp_zt,                        & ! u'w' on thermo. grid [m^2/s^2]
      vpwp_zt,                        & ! v'w' on thermo. grid [m^2/s^2]
      up2_zt,                         & ! u'^2 on thermo. grid [m^2/s^2]
      vp2_zt,                         & ! v'^2 on thermo. grid [m^2/s^2]
      sigma_sqd_w_zt,                 & ! PDF width parameter (thermodynamic levels) [-]
      xp3_coef_fnc,                   & ! Coefficient in simple xp3 equation [-]
      wpsclrp_zt,                     & ! Scalar flux on thermo. levels [un. vary]
      sclrp2_zt                         ! Scalar variance on thermo. levels [un. vary]

    integer :: &
      i,    & ! Grid-column loop index
      k,    & ! Vertical-level loop index
      sclr    ! Passive scalar loop index

    !----------------------------- Begin Code ------------------------------

    !$acc data create( thvm_zm, ddzm_thvm_zm, brunt_vaisala_freq_sqd_zt, &
    !$acc              wp2_zt, Skw_zt, wpthlp_zt, wprtp_zt, upwp_zt, vpwp_zt, &
    !$acc              thlp2_zt, rtp2_zt, up2_zt, vp2_zt, sigma_sqd_w_zt, &
    !$acc              xp3_coef_fnc, wpsclrp_zt, sclrp2_zt )

    wp2_zt(:,:) = zm2zt_api( nzm, nzt, ngrdcol, gr, wp2(:,:), w_tol_sqd )

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

    if ( iiPDF_type == iiPDF_ADG1 .or. iiPDF_type == iiPDF_parcel_ensemble ) then

      ! Use the Larson and Golaz (2005) ansatz for the ADG1 PDF to
      ! calculate <rt'^3>, <thl'^3>, <u'^3>, <v'^3>, and <sclr'^3>.
      sigma_sqd_w_zt(:,:) = zm2zt_api( nzm, nzt, ngrdcol, gr, sigma_sqd_w(:,:), &
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

        wpsclrp_zt(:,:) = zm2zt_api( nzm, nzt, ngrdcol, gr, wpsclrp(:,:,sclr), &
                                      sclr_tol(sclr)**2 )
        sclrp2_zt(:,:)  = zm2zt_api( nzm, nzt, ngrdcol, gr, sclrp2(:,:,sclr), &
                                      sclr_tol(sclr)**2 )

        call xp3_LG_2005_ansatz( nzt, ngrdcol, Skw_zt, wpsclrp_zt, wp2_zt, &
                                  sclrp2_zt, sigma_sqd_w_zt, &
                                  clubb_params, sclr_tol(sclr), &
                                  sclrp3(:,:,sclr) )

      end do ! sclr = 1, sclr_dim

    else ! iiPDF_type /= iiPDF_ADG1

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
        sclrp2_zt(:,:)  = zm2zt_api( nzm, nzt, ngrdcol, gr, sclrp2(:,:,sclr), &
                                      sclr_tol(sclr)**2 )

        call xp3_LG_2005_ansatz( nzt, ngrdcol, Skw_zt(:,:), wpsclrp_zt(:,:), wp2_zt(:,:), &
                                  sclrp2_zt(:,:), xp3_coef_fnc(:,:), &
                                  clubb_params, sclr_tol(sclr), &
                                  sclrp3(:,:,sclr) )
      end do ! sclr = 1, sclr_dim

    end if ! iiPDF_type == iiPDF_ADG1

    if ( stats%l_sample ) then
      !$acc update host( thlp2_zt, wpthlp_zt, wprtp_zt, rtp2_zt, &
      !$acc              up2_zt, vp2_zt, upwp_zt, vpwp_zt )

      call stats_update( "thlp2_zt", thlp2_zt, stats )
      call stats_update( "wpthlp_zt", wpthlp_zt, stats )
      call stats_update( "wprtp_zt", wprtp_zt, stats )
      call stats_update( "rtp2_zt", rtp2_zt, stats )
      call stats_update( "up2_zt", up2_zt, stats )
      call stats_update( "vp2_zt", vp2_zt, stats )
      call stats_update( "upwp_zt", upwp_zt, stats )
      call stats_update( "vpwp_zt", vpwp_zt, stats )
    end if

    !$acc end data

  end subroutine diagnose_xp3

  !=============================================================================
  subroutine advance_xp3_trivar( nzm, nzt, ngrdcol, gr, dt,                    & ! Intent(in)
                                rtm, thlm, rtp2, thlp2, wprtp, wpthlp,       & ! Intent(in)
                                wprtp2, wpthlp2, wm_zt, Kh_zm, Kh_zt, rho_ds_zm, & ! Intent(in)
                                invrs_rho_ds_zt, invrs_tau_zt,                & ! Intent(in)
                                pdf_params, clubb_params, nu_vert_res_dep,   & ! Intent(in)
                                tridiag_solve_method, l_upwind_xm_ma,        & ! Intent(in)
                                l_implemented, l_lmm_stepping,               & ! Intent(in)
                                stats, err_info,                              & ! Intent(inout)
                                rtp3, thlp3 )                                   ! Intent(inout)

    ! Description:
    !   Advances the trivariate transport PDF scalar third moments <rt'^3> and <thl'^3> with a
    !   conservative, partly implicit discretization.  For either scalar x,
    !
    !   d<x'^3>/dt = - <w> d<x'^3>/dz
    !                - rho_d^-1 d[rho_d <w'x'^3>]/dz
    !                - 3 <w'x'^2> d<x>/dz
    !                + 3 <x'^2> rho_d^-1 d[rho_d <w'x'>]/dz
    !                - <x'^3>/tau
    !                + d[(K_xp3+nu3)d<x'^3>/dz]/dz.
    !
    !   trivariate transport PDF supplies the previously missing flux <w'x'^3> exactly from its
    !   frozen two-Gaussian geometry.  It is diagnosed on thermodynamic levels
    !   from the authoritative PDF parameter state, then interpolated to the
    !   momentum-level flux grid used by the conservative divergence.  This
    !   avoids the optional momentum-level PDF container, which is populated
    !   only when the two-call PDF option is selected.  For a component with
    !   center departures
    !   d_w and d_x, scalar variance V_x, and internal covariance C_wx,
    !
    !   <w'x'^3>_comp = d_w(d_x^3+3d_x V_x) + 3C_wx(d_x^2+V_x).
    !
    !   Its density-weighted divergence, turbulent production, and accumulation
    !   are evaluated explicitly at the old time.  Mean advection, eddy
    !   diffusion K_xp3=c_K3 K_h, and tau damping are assembled in the usual
    !   three-band thermodynamic-level matrix.  Hence each scalar requires one
    !   standard tridiagonal solve, in the same descending ordering used by the
    !   other CLUBB implicit advances.  The top value is fixed to zero, matching
    !   the historical xp3 boundary condition.  The diagnosed turbulent flux is
    !   closed at the lower and upper model boundaries so the new term cannot
    !   inject an unmodelled boundary third-moment flux.
    !
    !   The forcing covariance 3<x'^2 F_x'> is intentionally not included in
    !   this first implementation because microphysical and radiative schemes
    !   do not yet provide that covariance consistently.  This routine therefore
    !   improves resolved vertical transport without claiming a complete scalar
    !   third-moment forcing closure.
    !---------------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd

    use constants_clubb, only: &
        one,      &
        one_half, &
        zero,     &
        rt_tol,   &
        thl_tol,  &
        zero_threshold

    use grid_class, only: &
        grid,       &
        zm2zt_api,  &
        zt2zm_api

    use pdf_parameter_module, only: &
        pdf_parameter

    use pdf_closure_module, only: &
        calc_wpxp3_pdf

    use parameters_tunable, only: &
        nu_vertical_res_dep

    use parameter_indices, only: &
        nparams, &
        ic_K3

    use mean_adv, only: &
        term_ma_zt_lhs

    use diffusion, only: &
        diffusion_zt_lhs

    use err_info_type_module, only: &
        err_info_type

    implicit none

    integer, intent(in) :: &
      nzm, nzt, ngrdcol

    type(grid), intent(in) :: &
      gr

    real(kind=core_rknd), intent(in) :: &
      dt

    real(kind=core_rknd), dimension(ngrdcol,nzt), intent(in) :: &
      rtm, thlm, wprtp2, wpthlp2, wm_zt, invrs_rho_ds_zt, invrs_tau_zt

    real(kind=core_rknd), dimension(ngrdcol,nzm), intent(in) :: &
      rtp2, thlp2, wprtp, wpthlp, Kh_zm, Kh_zt, rho_ds_zm

    type(pdf_parameter), intent(in) :: &
      pdf_params

    real(kind=core_rknd), dimension(ngrdcol,nparams), intent(in) :: &
      clubb_params

    type(nu_vertical_res_dep), intent(in) :: &
      nu_vert_res_dep

    integer, intent(in) :: &
      tridiag_solve_method

    logical, intent(in) :: &
      l_upwind_xm_ma, l_implemented, l_lmm_stepping

    type(stats_type), intent(inout) :: &
      stats

    type(err_info_type), intent(inout) :: &
      err_info

    real(kind=core_rknd), dimension(ngrdcol,nzt), intent(inout) :: &
      rtp3, thlp3

    real(kind=core_rknd), dimension(ngrdcol,nzm) :: &
      rtm_zm, thlm_zm, Kxp3_zm, Kxp3_zt, wprtp3, wpthlp3

    real(kind=core_rknd), dimension(ngrdcol,nzt) :: &
      rtp2_zt, thlp2_zt, wprtp3_zt, wpthlp3_zt,                             &
      term_tp_rt, term_tp_thl, term_ac_rt, term_ac_thl, &
      term_ta_rt, term_ta_thl, term_ma_rt, term_ma_thl, term_df_rt, term_df_thl, &
      term_dp_rt, term_dp_thl, term_bt_rt, term_bt_thl, term_rs_rt, term_rs_thl, &
      rtp3_old, thlp3_old

    real(kind=core_rknd), dimension(3,ngrdcol,nzt) :: &
      lhs_ma, lhs_diff, lhs

    real(kind=core_rknd), dimension(ngrdcol,nzt,1) :: &
      rhs, solution

    integer :: &
      i, k, kp1

    real(kind=core_rknd) :: &
      budget_factor

    ! All new terms are evaluated from the frozen PDF state before either
    ! scalar-third-moment solution is written back.
    rtm_zm = zt2zm_api( nzm, nzt, ngrdcol, gr, rtm, zero_threshold )
    thlm_zm = zt2zm_api( nzm, nzt, ngrdcol, gr, thlm, zero_threshold )
    rtp2_zt = zm2zt_api( nzm, nzt, ngrdcol, gr, rtp2, rt_tol**2 )
    thlp2_zt = zm2zt_api( nzm, nzt, ngrdcol, gr, thlp2, thl_tol**2 )

    call calc_wpxp3_pdf( nzt, ngrdcol, wm=wm_zt, xm=rtm,                    &
                          w_1=pdf_params%w_1, w_2=pdf_params%w_2,           &
                          x_1=pdf_params%rt_1, x_2=pdf_params%rt_2,         &
                          varnce_w_1=pdf_params%varnce_w_1,                 &
                          varnce_w_2=pdf_params%varnce_w_2,                 &
                          varnce_x_1=pdf_params%varnce_rt_1,                &
                          varnce_x_2=pdf_params%varnce_rt_2,                &
                          corr_w_x_1=pdf_params%corr_w_rt_1,                &
                          corr_w_x_2=pdf_params%corr_w_rt_2,                &
                          mixt_frac=pdf_params%mixt_frac, wpxp3=wprtp3_zt )

    call calc_wpxp3_pdf( nzt, ngrdcol, wm=wm_zt, xm=thlm,                   &
                          w_1=pdf_params%w_1, w_2=pdf_params%w_2,           &
                          x_1=pdf_params%thl_1, x_2=pdf_params%thl_2,       &
                          varnce_w_1=pdf_params%varnce_w_1,                 &
                          varnce_w_2=pdf_params%varnce_w_2,                 &
                          varnce_x_1=pdf_params%varnce_thl_1,               &
                          varnce_x_2=pdf_params%varnce_thl_2,               &
                          corr_w_x_1=pdf_params%corr_w_thl_1,               &
                          corr_w_x_2=pdf_params%corr_w_thl_2,               &
                          mixt_frac=pdf_params%mixt_frac, wpxp3=wpthlp3_zt )

    wprtp3 = zt2zm_api( nzm, nzt, ngrdcol, gr, wprtp3_zt, zero )
    wpthlp3 = zt2zm_api( nzm, nzt, ngrdcol, gr, wpthlp3_zt, zero )

    ! No surface or model-top closure exists for <w'x'^3>; use a closed flux
    ! boundary rather than importing a synthetic source from PDF interpolation.
    wprtp3(:,gr%k_lb_zm) = zero
    wprtp3(:,gr%k_ub_zm) = zero
    wpthlp3(:,gr%k_lb_zm) = zero
    wpthlp3(:,gr%k_ub_zm) = zero

    Kxp3_zm = spread(clubb_params(:,ic_K3), dim=2, ncopies=nzm) * Kh_zm
    Kxp3_zt = spread(clubb_params(:,ic_K3), dim=2, ncopies=nzt) * Kh_zt

    call term_ma_zt_lhs( nzm, nzt, ngrdcol, wm_zt, gr%weights_zt2zm, &
                         gr%invrs_dzt, gr%invrs_dzm, l_upwind_xm_ma, &
                         gr%grid_dir, lhs_ma )
    call diffusion_zt_lhs( nzm, nzt, ngrdcol, gr, Kxp3_zm, Kxp3_zt, &
                           nu_vert_res_dep%nu3, invrs_rho_ds_zt, rho_ds_zm, lhs_diff )
    call xp3_trivar_lhs( nzt, ngrdcol, gr, dt, invrs_tau_zt, lhs_ma, lhs_diff, lhs )

    do k = 1, nzt-1
      kp1 = min(k+1,nzm)
      do i = 1, ngrdcol
        term_tp_rt(i,k) = term_tp_rhs( rtp2_zt(i,k), wprtp(i,kp1), wprtp(i,k), &
                                        rho_ds_zm(i,kp1), rho_ds_zm(i,k), &
                                        invrs_rho_ds_zt(i,k), gr%invrs_dzt(i,k) )
        term_tp_thl(i,k) = term_tp_rhs( thlp2_zt(i,k), wpthlp(i,kp1), wpthlp(i,k), &
                                         rho_ds_zm(i,kp1), rho_ds_zm(i,k), &
                                         invrs_rho_ds_zt(i,k), gr%invrs_dzt(i,k) )
        term_ac_rt(i,k) = term_ac_rhs( rtm_zm(i,kp1), rtm_zm(i,k), wprtp2(i,k), gr%invrs_dzt(i,k) )
        term_ac_thl(i,k) = term_ac_rhs( thlm_zm(i,kp1), thlm_zm(i,k), &
                                         wpthlp2(i,k), gr%invrs_dzt(i,k) )
        term_ta_rt(i,k) = term_ta_rhs( wprtp3(i,kp1), wprtp3(i,k), rho_ds_zm(i,kp1), &
                                        rho_ds_zm(i,k), invrs_rho_ds_zt(i,k), gr%invrs_dzt(i,k) )
        term_ta_thl(i,k) = term_ta_rhs( wpthlp3(i,kp1), wpthlp3(i,k), rho_ds_zm(i,kp1), &
                                         rho_ds_zm(i,k), invrs_rho_ds_zt(i,k), gr%invrs_dzt(i,k) )
      end do
    end do
    term_tp_rt(:,nzt) = zero
    term_tp_thl(:,nzt) = zero
    term_ac_rt(:,nzt) = zero
    term_ac_thl(:,nzt) = zero
    term_ta_rt(:,nzt) = zero
    term_ta_thl(:,nzt) = zero

    rtp3_old = rtp3
    thlp3_old = thlp3
    rhs(:,:,1) = rtp3 / dt + term_tp_rt + term_ac_rt + term_ta_rt
    rhs(:,nzt,1) = zero
    call xp3_trivar_solve( nzt, ngrdcol, gr, tridiag_solve_method, l_implemented, &
                          lhs, rhs, err_info, solution )
    rtp3 = solution(:,:,1)

    rhs(:,:,1) = thlp3 / dt + term_tp_thl + term_ac_thl + term_ta_thl
    rhs(:,nzt,1) = zero
    call xp3_trivar_solve( nzt, ngrdcol, gr, tridiag_solve_method, l_implemented, &
                          lhs, rhs, err_info, solution )
    thlp3 = solution(:,:,1)

    if ( stats%l_sample ) then
      ! The LMM option stores the midpoint of the old and fully solved states.
      ! Scale every right- and left-hand-side tendency consistently so that
      ! the printed budget sums to the actual stored-state tendency.
      budget_factor = one
      if ( l_lmm_stepping ) budget_factor = one_half

      call implicit_zt_budget_term( nzt, ngrdcol, gr, lhs_ma, rtp3, term_ma_rt )
      call implicit_zt_budget_term( nzt, ngrdcol, gr, lhs_ma, thlp3, term_ma_thl )
      call implicit_zt_budget_term( nzt, ngrdcol, gr, lhs_diff, rtp3, term_df_rt )
      call implicit_zt_budget_term( nzt, ngrdcol, gr, lhs_diff, thlp3, term_df_thl )

      term_ta_rt = budget_factor * term_ta_rt
      term_ta_thl = budget_factor * term_ta_thl
      term_tp_rt = budget_factor * term_tp_rt
      term_tp_thl = budget_factor * term_tp_thl
      term_ac_rt = budget_factor * term_ac_rt
      term_ac_thl = budget_factor * term_ac_thl
      term_ma_rt = budget_factor * term_ma_rt
      term_ma_thl = budget_factor * term_ma_thl
      term_df_rt = budget_factor * term_df_rt
      term_df_thl = budget_factor * term_df_thl
      term_dp_rt = -budget_factor * invrs_tau_zt * rtp3
      term_dp_thl = -budget_factor * invrs_tau_zt * thlp3
      term_bt_rt = budget_factor * ( rtp3 - rtp3_old ) / dt
      term_bt_thl = budget_factor * ( thlp3 - thlp3_old ) / dt
      term_rs_rt = term_bt_rt - ( term_ta_rt + term_tp_rt + term_ac_rt + term_ma_rt + &
                                   term_df_rt + term_dp_rt )
      term_rs_thl = term_bt_thl - ( term_ta_thl + term_tp_thl + term_ac_thl + term_ma_thl + &
                                    term_df_thl + term_dp_thl )

      ! ``*_tp``/``*_ac`` retain their historic names.  The remaining fields
      ! make the trivariate transport PDF equation observable term by term; ``*_rs`` should be
      ! roundoff-sized except for an intentional future unbudgeted adjustment.
      call stats_update( "rtp3_bt", term_bt_rt, stats )
      call stats_update( "rtp3_ta", term_ta_rt, stats )
      call stats_update( "rtp3_tp", term_tp_rt, stats )
      call stats_update( "rtp3_ac", term_ac_rt, stats )
      call stats_update( "rtp3_ma", term_ma_rt, stats )
      call stats_update( "rtp3_df", term_df_rt, stats )
      call stats_update( "rtp3_dp", term_dp_rt, stats )
      call stats_update( "rtp3_rs", term_rs_rt, stats )
      call stats_update( "thlp3_bt", term_bt_thl, stats )
      call stats_update( "thlp3_ta", term_ta_thl, stats )
      call stats_update( "thlp3_tp", term_tp_thl, stats )
      call stats_update( "thlp3_ac", term_ac_thl, stats )
      call stats_update( "thlp3_ma", term_ma_thl, stats )
      call stats_update( "thlp3_df", term_df_thl, stats )
      call stats_update( "thlp3_dp", term_dp_thl, stats )
      call stats_update( "thlp3_rs", term_rs_thl, stats )
    end if

    if ( l_lmm_stepping ) then
      rtp3 = one_half * (rtp3_old + rtp3)
      thlp3 = one_half * (thlp3_old + thlp3)
      rtp3(:,nzt) = zero
      thlp3(:,nzt) = zero
    end if

  end subroutine advance_xp3_trivar

  !=============================================================================
  subroutine implicit_zt_budget_term( nzt, ngrdcol, gr, lhs_term, field, tendency )

    ! Description:
    !   Reconstructs the tendency represented by a three-band thermodynamic-
    !   level operator on the left-hand side of the trivariate transport PDF scalar-third-moment
    !   solve.  The solver stores the negative of each implicit physical
    !   tendency on its LHS, hence this routine returns ``-L * field``.  It
    !   follows the same boundary/interior stencil used by CLUBB's existing
    !   wp3 and mean-wind budget diagnostics.  The trivariate transport PDF upper boundary is a
    !   fixed zero third moment and has no prognostic tendency.
    !---------------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd

    use constants_clubb, only: &
        zero

    use grid_class, only: &
        grid

    implicit none

    integer, intent(in) :: &
      nzt, ngrdcol

    type(grid), intent(in) :: &
      gr

    real(kind=core_rknd), dimension(3,ngrdcol,nzt), intent(in) :: &
      lhs_term

    real(kind=core_rknd), dimension(ngrdcol,nzt), intent(in) :: &
      field

    real(kind=core_rknd), dimension(ngrdcol,nzt), intent(out) :: &
      tendency

    integer :: &
      i, k, km1, kp1

    tendency = zero
    do i = 1, ngrdcol
      k = gr%k_lb_zt
      kp1 = k + gr%grid_dir_indx
      tendency(i,k) = - lhs_term(2,i,k) * field(i,k) &
                      - lhs_term(2-gr%grid_dir_indx,i,k) * field(i,kp1)

      do k = 2, nzt-1
        km1 = k - gr%grid_dir_indx
        kp1 = k + gr%grid_dir_indx
        tendency(i,k) = - lhs_term(2+gr%grid_dir_indx,i,k) * field(i,km1) &
                        - lhs_term(2,i,k) * field(i,k) &
                        - lhs_term(2-gr%grid_dir_indx,i,k) * field(i,kp1)
      end do
    end do

  end subroutine implicit_zt_budget_term

  !=============================================================================
  subroutine xp3_trivar_lhs( nzt, ngrdcol, gr, dt, invrs_tau_zt, lhs_ma, lhs_diff, lhs )

    ! Description:
    !   Adds the fully implicit scalar-third-moment terms to the standard
    !   thermodynamic-level three-band matrix.  Mean advection and diffusion
    !   enter with the signs already required on the left-hand side; the
    !   backward-Euler time tendency and tau damping add 1/dt+1/tau to its
    !   diagonal.  The upper boundary remains the historical fixed xp3=0.
    !---------------------------------------------------------------------------

    use clubb_precision, only: core_rknd
    use constants_clubb, only: one, zero
    use grid_class, only: grid
    implicit none
    integer, intent(in) :: nzt, ngrdcol
    type(grid), intent(in) :: gr
    real(kind=core_rknd), intent(in) :: dt
    real(kind=core_rknd), dimension(ngrdcol,nzt), intent(in) :: invrs_tau_zt
    real(kind=core_rknd), dimension(3,ngrdcol,nzt), intent(in) :: lhs_ma, lhs_diff
    real(kind=core_rknd), dimension(3,ngrdcol,nzt), intent(out) :: lhs
    integer :: i, k
    do k = 1, nzt-1
      do i = 1, ngrdcol
        lhs(:,i,k) = lhs_ma(:,i,k) + lhs_diff(:,i,k)
        lhs(2,i,k) = lhs(2,i,k) + one / dt + invrs_tau_zt(i,k)
      end do
    end do
    do i = 1, ngrdcol
      lhs(1,i,gr%k_ub_zt) = zero
      lhs(2,i,gr%k_ub_zt) = one
      lhs(3,i,gr%k_ub_zt) = zero
    end do
  end subroutine xp3_trivar_lhs

  !=============================================================================
  subroutine xp3_trivar_solve( nzt, ngrdcol, gr, tridiag_solve_method, l_implemented, &
                               lhs, rhs, err_info, solution )

    ! Description:
    !   Solves one trivariate transport PDF scalar-third-moment tridiagonal system.  Reversing
    !   both level order and matrix bands before the solve matches CLUBB's
    !   descending-order reproducibility convention.
    !---------------------------------------------------------------------------

    use clubb_precision, only: core_rknd
    use grid_class, only: grid
    use matrix_solver_wrapper, only: tridiag_solve
    use model_flags, only: l_force_descending_solves
    use err_info_type_module, only: err_info_type
    implicit none
    integer, intent(in) :: nzt, ngrdcol, tridiag_solve_method
    type(grid), intent(in) :: gr
    logical, intent(in) :: l_implemented
    real(kind=core_rknd), dimension(3,ngrdcol,nzt), intent(inout) :: lhs
    real(kind=core_rknd), dimension(ngrdcol,nzt,1), intent(inout) :: rhs
    type(err_info_type), intent(inout) :: err_info
    real(kind=core_rknd), dimension(ngrdcol,nzt,1), intent(out) :: solution
    if ( l_force_descending_solves .and. gr%grid_dir_indx > 0 ) then
      lhs = lhs(3:1:-1,:,nzt:1:-1)
      rhs = rhs(:,nzt:1:-1,:)
    end if
    call tridiag_solve( "xp3_trivar", tridiag_solve_method, ngrdcol, nzt, 1, &
                        l_implemented, lhs, rhs, err_info, solution )
    if ( l_force_descending_solves .and. gr%grid_dir_indx > 0 ) then
      solution = solution(:,nzt:1:-1,:)
    end if
  end subroutine xp3_trivar_solve

  !=============================================================================
  function term_ta_rhs( wpxpp3, wpxp3, rho_ds_zmp1, rho_ds_zm, &
                        invrs_rho_ds_zt, invrs_dzt ) result(term_ta)

    ! Description:
    !   Conservative explicit turbulent-advection tendency
    !   -rho_d^-1 d[rho_d<w'x'^3>]/dz at one thermodynamic level.
    !---------------------------------------------------------------------------
    use clubb_precision, only: core_rknd
    implicit none
    real(kind=core_rknd), intent(in) :: wpxpp3, wpxp3, rho_ds_zmp1, rho_ds_zm, &
                                        invrs_rho_ds_zt, invrs_dzt
    real(kind=core_rknd) :: term_ta
    term_ta = -invrs_rho_ds_zt * invrs_dzt * (rho_ds_zmp1*wpxpp3 - rho_ds_zm*wpxp3)
  end function term_ta_rhs

  !=============================================================================
  subroutine advance_xp3_simplified( nzm, nzt, ngrdcol, gr, solve_type, dt, & ! Intent(in)
                                     xm, xp2, wpxp,                         & ! Intent(in)
                                     wpxp2, rho_ds_zm,                      & ! Intent(in)
                                     invrs_rho_ds_zt,                       & ! Intent(in)
                                     invrs_tau_zt, tau_max_zt,              & ! Intent(in)
                                     x_tol, l_lmm_stepping,                 & ! Intent(in)
                                     stats,                                 & ! Intent(in)
                                     xp3 )                                    ! Intent(inout)

    ! Description:
    ! Predicts the value of <x'^3> using a simplified form of the <x'^3>
    ! predictive equation.
    !
    ! The full predictive equation for <x'^3>, where <x'^3> can be <rt'^3>,
    ! <thl'^3>, or <sclr'^3>, is:
    !
    ! d<x'^3>/dt = - <w> * d<x'^3>/dz
    !              - (1/rho_ds) * d( rho_ds * <w'x'^3> )/dz
    !              - 3 * <w'x'^2> * d<x>/dz
    !              + 3 * ( <x'^2> / rho_ds ) * d( rho_ds * <w'x'> )/dz
    !              - ( C_xp3_dissipation / tau ) * <x'^3>
    !              + d ( ( K_xp3 + nu_xp3 ) * d<x'^3>/dz )/dz
    !              + 3 * < x'^2 (dx/dt)|_f' >;
    !
    ! where (dx/dt)|_f is the "forcing" term, which may include effects such as
    ! microphysical effects or radiative effects.  The tunable coefficients are
    ! C_xp3_dissipation, K_xp3, and nu_xp3.  The terms are listed as follows:
    !
    ! time tendency: d<x'^3>/dt;
    ! mean advection: - <w> * d<x'^3>/dz;
    ! turbulent advection: - (1/rho_ds) * d( rho_ds * <w'x'^3> )/dz;
    ! accumulation: - 3 * <w'x'^2> * d<x>/dz;
    ! turbulent production: + 3 * ( <x'^2> / rho_ds ) * d( rho_ds * <w'x'> )/dz;
    ! turbulent dissipation: - ( C_xp3_dissipation / tau ) * <x'^3>;
    ! diffusion: + d ( ( K_xp3 + nu_xp3 ) * d<x'^3>/dz )/dz; and
    ! microphysics/other forcing: + 3 * < x'^2 (dx/dt)|_f' >.
    !
    ! The microphysics and turbulent advection terms are both found by
    ! integration over the subgrid PDF.  This requires new integrated terms.
    ! The turbulent advection term may need to be made semi-implicit in order
    ! to aid model stability.  This may be difficult to do for <x'^3>.
    ! Additionally, if it could be made semi-implicit, it involves a derivative
    ! and would require a tridiagonal solver to include contributions from
    ! <x'^3> on three grid levels.  While the microphysics term and turbulent
    ! advection term are important contributors to <x'^3>, they are being
    ! omitted because of the additional complications they bring.
    !
    ! The mean advection and diffusion terms also would require a tridiagonal
    ! solver in order to make the terms implicit because they involve
    ! derivatives and values of <x'^3> on three grid levels.  While tridiagonal
    ! solvers are not very computationally expensive, they are still more
    ! expensive than a simplified one-line equation.  The mean advection and
    ! diffusion terms are also rather small in magnitude, so they are also
    ! being neglected.
    !
    ! This leaves the following equation:
    !
    ! d<x'^3>/dt = - 3 * <w'x'^2> * d<x>/dz
    !              + 3 * ( <x'^2> / rho_ds ) * d( rho_ds * <w'x'> )/dz
    !              - ( C_xp3_dissipation / tau ) * <x'^3>;
    !
    ! which is a balance of time-tendency, accumulation, turbulent production,
    ! and turbulent dissipation.  This equation can be handled semi-implicitly
    ! as:
    !
    ! ( <x'^3>(t+1) - <x'^3>(t) ) / delta_t
    ! = - 3 * <w'x'^2> * d<x>/dz
    !   + 3 * ( <x'^2> / rho_ds ) * d( rho_ds * <w'x'> )/dz
    !   - ( C_xp3_dissipation / tau ) * <x'^3>(t+1);
    !
    ! which can be rewritten as:
    !
    ! ( 1 / delta_t + ( C_xp3_dissipation / tau ) ) * <x'^3>(t+1)
    ! = ( <x'^3>(t) / delta_t )
    !   - 3 * <w'x'^2> * d<x>/dz
    !   + 3 * ( <x'^2> / rho_ds ) * d( rho_ds * <w'x'> )/dz.
    !
    ! The predictive equation can be solved for <x'^3> as:
    !
    ! <x'^3>(t+1)
    ! = ( ( <x'^3>(t) / delta_t )
    !     - 3 * <w'x'^2> * d<x>/dz
    !     + 3 * ( <x'^2> / rho_ds ) * d( rho_ds * <w'x'> )/dz )
    !   / ( 1 / delta_t + ( C_xp3_dissipation / tau ) ).
    !
    ! Alternatively, a steady-state approximation can be used, which
    ! approximates d<x'^3>/dt = 0.  The equation becomes a balance of
    ! accumulation, turbulent production, and turbulent dissipation, and is
    ! written as:
    !
    ! 0 = - 3 * <w'x'^2> * d<x>/dz
    !     + 3 * ( <x'^2> / rho_ds ) * d( rho_ds * <w'x'> )/dz
    !     - ( C_xp3_dissipation / tau ) * <x'^3>.
    !
    ! The equation can be solved for <x'^3> as:
    !
    ! <x'^3>
    ! = ( tau / C_xp3_dissipation )
    !   * ( - 3 * <w'x'^2> * d<x>/dz
    !       + 3 * ( <x'^2> / rho_ds ) * d( rho_ds * <w'x'> )/dz ).
    !
    ! When the flag l_predict_xp3 is enabled, the predictive version of <x'^3>
    ! is used.  When the flag is turned off, the steady-state approximation is
    ! used.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        grid, & ! Type
        zm2zt_api, & ! Procedure(s)
        zt2zm_api

    use constants_clubb, only: &
        one,      & ! Variable(s)
        one_half, &
        zero, &
        zero_threshold

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! ----------------------- Input Variables -----------------------
    integer, intent(in) :: &
      nzm, &
      nzt, &
      ngrdcol

    type (grid), intent(in) :: gr

    integer, intent(in) :: &
      solve_type    ! Flag for solving for rtp3, thlp3, or sclrp3

    real( kind = core_rknd ), intent(in) :: &
      dt                 ! Model timestep                            [s]

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      xm,              & ! Mean (overall) of x (thermo. levels)                 [(x units)]
      wpxp2,           & ! <w'x'^2> (thermodynamic levels)                      [m/s(x units)^2]
      invrs_rho_ds_zt, & ! Inv. dry, static density at thermo. levels           [m^3/kg]
      invrs_tau_zt,    & ! Inverse time-scale tau on thermodynamic levels       [1/s]
      tau_max_zt         ! Max. allowable eddy dissipation time scale on t-levs [s]

    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: &
      xp2,             & ! Variance (overall) of x (m-levs.)       [(x units)^2]
      wpxp,            & ! Turbulent flux of x (momentum levs.)    [m/s(x units)]
      rho_ds_zm          ! Dry, static density on momentum levels  [kg/m^3]

    real( kind = core_rknd ), intent(in) :: &
      x_tol    ! Tolerance value of x                           [(x units)]

    logical, intent(in) :: &
      l_lmm_stepping    ! Apply Linear Multistep Method (LMM) Stepping

    ! --------------------- Input/Output Variables ---------------------
    ! ----------------------- Input/Output Variable -----------------------
    type(stats_type), intent(inout) :: &
      stats

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(inout) :: &
      xp3    ! <x'^3> (thermodynamic levels)    [(x units)^3]

    ! ----------------------- Local Variables -----------------------
    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      xp3_old    ! Saved <x'^3> (thermodynamic levels)    [(x units)^3]

    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
      xm_zm      ! Mean of x interpolated to momentum levels     [(x units)]

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      xp2_zt,  & ! Variance of x interpolated to thermo. levels  [(x units)^2]
      term_tp, & ! <x'^3> turbulent production term              [(x units)^3/s]
      term_ac    ! <x'^3> accumulation term                      [(x units)^3/s]

    integer :: &
      i, k, kp1     ! Grid indices

    character(len=32) :: &
      name_bt, &
      name_tp, &
      name_ac, &
      name_dp

    ! Coefficient in the <x'^3> turbulent dissipation term    [-]
    real( kind = core_rknd ), parameter :: &
      C_xp3_dissipation = 1.0_core_rknd

    ! Flag to either predict <x'^3> or use steady-state approximation.
    logical, parameter :: &
      l_predict_xp3 = .true.

    ! ----------------------- Begin Code -----------------------

    select case ( solve_type )
    case( xp3_rtp3 )
      ! Budget stats for rtp3
      name_bt = "rtp3_bt"
      name_tp = "rtp3_tp"
      name_ac = "rtp3_ac"
      name_dp = "rtp3_dp"
    case( xp3_thlp3 )
      ! Budget stats for thlp3
      name_bt = "thlp3_bt"
      name_tp = "thlp3_tp"
      name_ac = "thlp3_ac"
      name_dp = "thlp3_dp"
    case default
      ! Budgets aren't setup for the passive scalars
      name_bt = ""
      name_tp = ""
      name_ac = ""
      name_dp = ""
    end select ! solve_type

    if ( l_predict_xp3 ) then
      if ( stats%l_sample ) then
        call stats_begin_budget( name_bt, xp3 / dt, stats )
      end if
    end if ! l_predict_xp3

    ! Initialize variables
    term_tp = zero
    term_ac = zero

    ! Interpolate <x> to momentum levels.
    xm_zm = zt2zm_api( nzm, nzt, ngrdcol, gr, xm, zero_threshold )

    ! Interpolate <x'^2> to thermodynamic levels.
    xp2_zt = zm2zt_api( nzm, nzt, ngrdcol, gr, xp2, x_tol**2 )  ! Positive definite quantity

    do k = 1, nzt-1, 1
      do i = 1, ngrdcol

        ! Define the km1 index.
        kp1 = min( k+1, nzt )

        ! Calculate the <x'^3> turbulent production (tp) term.
        term_tp(i,k) = term_tp_rhs( xp2_zt(i,k), wpxp(i,kp1), wpxp(i,k), &
                                    rho_ds_zm(i,kp1), rho_ds_zm(i,k), &
                                    invrs_rho_ds_zt(i,k), &
                                    gr%invrs_dzt(i,k) )

        ! Calculate the <x'^3> accumulation (ac) term.
        term_ac(i,k) = term_ac_rhs( xm_zm(i,kp1), xm_zm(i,k), wpxp2(i,k), &
                                    gr%invrs_dzt(i,k) )

        if ( l_predict_xp3 ) then

           if ( l_lmm_stepping ) then
              xp3_old(i,k) = xp3(i,k)
           endif ! l_lmm_stepping

           ! Advance <x'^3> one time step.
           xp3(i,k) = ( ( xp3(i,k) / dt ) + term_tp(i,k) + term_ac(i,k) ) &
                      / ( ( one / dt ) + ( C_xp3_dissipation * invrs_tau_zt(i,k) ) )

           if ( l_lmm_stepping ) then
              xp3(i,k) = one_half * ( xp3_old(i,k) + xp3(i,k) )
           endif ! l_lmm_stepping

        else

           ! Calculate <x'^3> using the steady-state approximation.
           xp3(i,k) = min( one / invrs_tau_zt(i,k), tau_max_zt(i,k) ) * one / C_xp3_dissipation &
                      * ( term_tp(i,k) + term_ac(i,k) )

        endif ! l_predict_xp3

      end do
    end do ! k = 2, gr%nzt-1, 1

    ! Set Upper Boundary Condition
    xp3(:,nzt) = zero

    if ( stats%l_sample ) then
      call stats_update( name_tp, term_tp, stats )
      call stats_update( name_ac, term_ac, stats )
      call stats_update( name_dp, -(C_xp3_dissipation * invrs_tau_zt) * xp3, stats )
      if ( l_predict_xp3 ) then
        call stats_finalize_budget( name_bt, xp3 / dt, stats )
      end if
    end if

    return

  end subroutine advance_xp3_simplified

  !=============================================================================
  function term_tp_rhs( xp2_zt, wpxpp1, wpxp, &
                        rho_ds_zmp1, rho_ds_zm, &
                        invrs_rho_ds_zt, &
                        invrs_dzt ) &
  result( term_tp )

    ! Description:
    ! Turbulent production of <x'^3>:  explicit portion of the code.
    !
    ! The d<x'^3>/dt equation contains a turbulent production term:
    !
    ! + 3 * ( <x'^2> / rho_ds ) * d( rho_ds * <w'x'> )/dz.
    !
    ! The <x'^3> turbulent production term is completely explicit and is
    ! discretized as follows:
    !
    ! The values of <x'^3> are found on the thermodynamic levels, while the
    ! values of <w'x'> and <x'^2> are found on the momentum levels.
    ! Additionally, the values of rho_ds_zm are found on the momentum levels,
    ! and the values of invrs_rho_ds_zt are found on the thermodynamic levels.
    ! The values of <x'^2> are interpolated to the central thermodynamic level
    ! as <x'^2>|_zt.  On the momentum levels, the values of <w'x'> are
    ! multiplied by rho_ds_zm.  Then, the derivative (d/dz) of
    ! rho_ds_zm * <w'x'> is taken over the central thermodynamic level.  At the
    ! central thermodynamic level, the derivative is multiplied by
    ! invrs_rho_ds_zt, and their product is also multiplied by 3 * <x'^2>|_zt,
    ! yielding the desired results.
    !
    ! =========wpxpp1=========rho_ds_zmp1===========xp2p1================ m(k+1)
    !
    ! --xp3--d( rho_ds_zm * wpxp )/dz--invrs_rho_ds_zt--xp2_zt(interp.)-- t(k)
    !
    ! =========wpxp===========rho_ds_zm=============xp2================== m(k)
    !
    ! The vertical indices m(k+1), t(k), and m(k) correspond with altitudes
    ! zm(k+1), zt(k), and zm(k), respectively.  The letter "t" is used for
    ! thermodynamic levels and the letter "m" is used for momentum levels.
    !
    ! invrs_dzt(k) = 1 / ( zm(k+1) - zm(k) )

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        three    ! Variable(s)

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      xp2_zt,          & ! <x'^2> interp. to thermo. level (k)     [(x units)^2]
      wpxpp1,          & ! <w'x'> at momentum level (k+1)         [m/s(x units)]
      wpxp,            & ! <w'x'> at momentum level (k)           [m/s(x units)]
      rho_ds_zmp1,     & ! Dry, static density on momentum level (k+1)  [kg/m^3]
      rho_ds_zm,       & ! Dry, static density on momentum level (k)    [kg/m^3]
      invrs_rho_ds_zt, & ! Inv. dry, static density at thermo. lev. (k) [m^3/kg]
      invrs_dzt          ! Inverse of grid spacing (k)                     [1/m]

    ! Return Variable
    real( kind = core_rknd ) :: &
      term_tp    ! <x'^3> turbulent production term              [(x units)^3/s]


    ! The <x'^3> turbulent production term.
    term_tp &
    = + three * xp2_zt * invrs_rho_ds_zt &
        * invrs_dzt * ( rho_ds_zmp1 * wpxpp1 - rho_ds_zm * wpxp )


    return

  end function term_tp_rhs

  !=============================================================================
  function term_ac_rhs( xm_zmp1, xm_zm, wpxp2, &
                        invrs_dzt ) &
  result( term_ac )

    ! Description:
    ! Accumulation of <x'^3>:  explicit portion of the code.
    !
    ! The d<x'^3>/dt equation contains an accumulation term:
    !
    ! - 3 * <w'x'^2> * d<x>/dz.
    !
    ! The <x'^3> accumulation term is completely explicit and is discretized as
    ! follows:
    !
    ! The values of <x'^3>, <x>, and <w'x'^2> are found on thermodynamic levels.
    ! The values of <x> are interpolated to the intermediate momentum levels as
    ! <x>|_zm.  Then, the derivative (d/dz) of <x>|_zm is taken over the
    ! central thermodynamic level, where it is multiplied by -3 * <w'x'^2>.
    !
    ! ----------------------xmp1----------------------------------------- t(k+1)
    !
    ! =========================xm_zmp1(interp.)========================== m(k+1)
    !
    ! ----------xp3---------xm---------dxm_zm/dz---------wpxp2----------- t(k)
    !
    ! =========================xm_zm(interp.)============================ m(k)
    !
    ! ----------------------xmm1----------------------------------------- t(k-1)
    !
    ! The vertical indices t(k+1), m(k+1), t(k), m(k), and t(k-1) correspond
    ! with altitudes zt(k+1), zm(k+1), zt(k), zm(k), and zt(k-1), respectively.
    ! The letter "t" is used for thermodynamic levels and the letter "m" is
    ! used for momentum levels.
    !
    ! invrs_dzt(k) = 1 / ( zm(k+1) - zm(k) )

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        three    ! Variable(s)

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      xm_zmp1,   & ! <x> interpolated to momentum level (k+1)  [(x units)]
      xm_zm,     & ! <x> interpolated to momentum level (k)    [(x units)]
      wpxp2,     & ! <w'x'^2> at thermodynamic level (k)       [m/s(x units)^2]
      invrs_dzt    ! Inverse of grid spacing (k)               [1/m]

    ! Return Variable
    real( kind = core_rknd ) :: &
      term_ac    ! <x'^3> accumulation term                    [(x units)^3/s]


    ! The <x'^3> accumulation term.
    term_ac &
    = - three * wpxp2 * invrs_dzt * ( xm_zmp1 - xm_zm )


    return

  end function term_ac_rhs

  !=============================================================================

end module advance_xp3_module
