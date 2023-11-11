! $Id$
!===============================================================================
module advance_xp3_module

  ! Description:
  ! Predicts the value of <x'^3> for <rt'^3>, <thl'^3>, and <sclr'^3>.

  ! References:
  !-------------------------------------------------------------------------

  implicit none

  public :: advance_xp3    ! Procedure(s)

  private :: advance_xp3_simplified, & ! Procedure(s)
             term_tp_rhs, &
             term_ac_rhs

  private ! default scope

  integer, parameter, private :: &
    xp3_rtp3 = 1,   & ! Named constant for solving rtp3
    xp3_thlp3 = 2,  & ! Named constant for solving thlp3
    xp3_sclrp3 = 3    ! Named constant for solving sclrp3

  contains

  !=============================================================================
  subroutine advance_xp3( nz, ngrdcol, gr, dt,                       & ! Intent(in)
                          rtm, thlm, rtp2, thlp2, wprtp,             & ! Intent(in)
                          wpthlp, wprtp2, wpthlp2, rho_ds_zm,        & ! Intent(in)
                          invrs_rho_ds_zt, invrs_tau_zt, tau_max_zt, & ! Intent(in)
                          sclrm, sclrp2, wpsclrp, wpsclrp2,          & ! Intent(in)
                          l_lmm_stepping,                            & ! Intent(in)
                          stats_metadata,                            & ! Intent(in)
                          stats_zt,                                  & ! intent(inout)
                          rtp3, thlp3, sclrp3 )                        ! Intent(inout)

    ! Description:
    ! Advance <rt'^3>, <thl'^3>, and <sclr'^3> one model timestep using a
    ! simplified form of the <x'^3> predictive equation.  The simplified <x'^3>
    ! equation can either be advanced from its previous value or calculated
    ! using a steady-state approximation.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        grid ! Type

    use constants_clubb, only: &
        rt_tol,  & ! Variable(s)
        thl_tol

    use parameters_model, only: &
        sclr_dim, & ! Variable(s)
        sclr_tol

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    use stats_type, only: &
        stats ! Type

    use stats_variables, only: &  
        stats_metadata_type

    implicit none

    ! --------------------- Input Variables ---------------------
    integer, intent(in) :: &
      nz, &
      ngrdcol
    
    type (grid), target, intent(in) :: gr
  
    real( kind = core_rknd ), intent(in) :: &
      dt                 ! Model timestep                            [s]

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      rtm,             & ! Mean (overall) of rt (thermo. levels)  [kg/kg]
      thlm,            & ! Mean (overall) of thl (thermo. levels) [K]
      rtp2,            & ! Variance (overall) of rt (m-levs.)     [kg^2/kg^2]
      thlp2,           & ! Variance (overall) of thl (m-levs.)    [K^2]
      wprtp,           & ! Turbulent flux of rt (momentum levs.)  [m/s kg/kg]
      wpthlp,          & ! Turbulent flux of thl (momentum levs.) [m/s K]
      wprtp2,          & ! <w'rt'^2> (thermodynamic levels)       [m/s(kg/kg)^2]
      wpthlp2,         & ! <w'thl'^2> (thermodynamic levels)      [m/s K^2]
      rho_ds_zm,       & ! Dry, static density on momentum levels      [kg/m^3]
      invrs_rho_ds_zt, & ! Inv. dry, static density at thermo. levels  [m^3/kg]
      invrs_tau_zt,    & ! Inverse time-scale tau on thermodynamic levels [1/s]
      tau_max_zt         ! Max. allowable eddy dissipation time scale on t-levs[s]

    real( kind = core_rknd ), dimension(ngrdcol,nz,sclr_dim), intent(in) :: &
      sclrm,    & ! Mean (overall) of sclr (thermo. levels) [sclr units]
      sclrp2,   & ! Variance (overall) of sclr (m-levs.)    [(sclr units)^2]
      wpsclrp,  & ! Turbulent flux of sclr (momentum levs.) [m/s(sclr units)]
      wpsclrp2    ! <w'sclr'^2> (thermodynamic levels)      [m/s(sclr units)^2]

    logical, intent(in) :: &
      l_lmm_stepping    ! Apply Linear Multistep Method (LMM) Stepping

    type (stats_metadata_type), intent(in) :: &
      stats_metadata

    ! --------------------- Input/Output Variables ---------------------
    type (stats), target, dimension(ngrdcol), intent(inout) :: &
      stats_zt
      
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(inout) :: &
      rtp3,  & ! <rt'^3> (thermodynamic levels)     [kg^3/kg^3]
      thlp3    ! <thl'^3> (thermodynamic levels)    [K^3]

    real( kind = core_rknd ), dimension(ngrdcol,nz,sclr_dim), intent(inout) :: &
      sclrp3    ! <sclr'^3> (thermodynamic levels)    [(sclr units)^3]

    ! --------------------- Local Variable ---------------------
    integer :: i, k, sclr    ! Loop index


    ! Advance <rt'^3> one model timestep or calculate <rt'^3> using a
    ! steady-state approximation.
    call advance_xp3_simplified( nz, ngrdcol, gr, xp3_rtp3, dt, & ! Intent(in)
                                 rtm, rtp2, wprtp,              & ! Intent(in)
                                 wprtp2, rho_ds_zm,             & ! Intent(in)
                                 invrs_rho_ds_zt,               & ! Intent(in)
                                 invrs_tau_zt, tau_max_zt,      & ! Intent(in) 
                                 rt_tol, l_lmm_stepping,        & ! Intent(in)
                                 stats_metadata,                & ! Intent(in)
                                 stats_zt,                      & ! intent(inout)
                                 rtp3 )                           ! Intent(inout)

    ! Advance <thl'^3> one model timestep or calculate <thl'^3> using a
    ! steady-state approximation.
    call advance_xp3_simplified( nz, ngrdcol, gr, xp3_thlp3, dt,  & ! Intent(in)
                                 thlm, thlp2, wpthlp,             & ! Intent(in)
                                 wpthlp2, rho_ds_zm,              & ! Intent(in)
                                 invrs_rho_ds_zt,                 & ! Intent(in)
                                 invrs_tau_zt, tau_max_zt,        & ! Intent(in) 
                                 thl_tol, l_lmm_stepping,         & ! Intent(in)
                                 stats_metadata,                  & ! Intent(in)
                                 stats_zt,                        & ! intent(inout)
                                 thlp3 )                            ! Intent(inout)

    ! Advance <sclr'^3> one model timestep or calculate <sclr'^3> using a
    ! steady-state approximation.
    do sclr = 1, sclr_dim, 1
      call advance_xp3_simplified( nz, ngrdcol, gr, xp3_sclrp3, dt,                     & ! In
                                  sclrm(:,:,sclr), sclrp2(:,:,sclr), wpsclrp(:,:,sclr), & ! In
                                  wpsclrp2(:,:,sclr), rho_ds_zm,                        & ! In
                                  invrs_rho_ds_zt,                                      & ! In
                                  invrs_tau_zt, tau_max_zt,                             & ! In 
                                  sclr_tol(sclr), l_lmm_stepping,                       & ! In
                                  stats_metadata,                                       & ! Intent(in)
                                  stats_zt,                                             & ! In/Out
                                  sclrp3(:,:,sclr) )                                      ! In/Out
    end do ! i = 1, sclr_dim

    return

  end subroutine advance_xp3

  !=============================================================================
  subroutine advance_xp3_simplified( nz, ngrdcol, gr, solve_type, dt, & ! Intent(in)
                                     xm, xp2, wpxp,                   & ! Intent(in)
                                     wpxp2, rho_ds_zm,                & ! Intent(in)
                                     invrs_rho_ds_zt,                 & ! Intent(in)
                                     invrs_tau_zt, tau_max_zt,        & ! Intent(in) 
                                     x_tol, l_lmm_stepping,           & ! Intent(in)
                                     stats_metadata,                  & ! Intent(in)
                                     stats_zt,                        & ! Intent(inout)
                                     xp3 )                              ! Intent(inout)

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
        zm2zt, & ! Procedure(s)
        zt2zm

    use constants_clubb, only: &
        one,      & ! Variable(s)
        one_half, &
        zero

    use stats_type_utilities, only: &
        stat_begin_update, & ! Procedure(s)
        stat_end_update,   &
        stat_update_var

    use stats_variables, only: &
        stats_metadata_type

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    use stats_type, only: stats ! Type

    implicit none
 
    ! ----------------------- Input Variables -----------------------
    integer, intent(in) :: &
      nz, &
      ngrdcol
    
    type (grid), target, intent(in) :: gr
  
    integer, intent(in) :: &
      solve_type    ! Flag for solving for rtp3, thlp3, or sclrp3

    real( kind = core_rknd ), intent(in) :: &
      dt                 ! Model timestep                            [s]

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      xm,              & ! Mean (overall) of x (thermo. levels) [(x units)]
      xp2,             & ! Variance (overall) of x (m-levs.)    [(x units)^2]
      wpxp,            & ! Turbulent flux of x (momentum levs.) [m/s(x units)]
      wpxp2,           & ! <w'x'^2> (thermodynamic levels)      [m/s(x units)^2]
      rho_ds_zm,       & ! Dry, static density on momentum levels      [kg/m^3]
      invrs_rho_ds_zt, & ! Inv. dry, static density at thermo. levels  [m^3/kg]
      invrs_tau_zt,    & ! Inverse time-scale tau on thermodynamic levels  [1/s]
      tau_max_zt         ! Max. allowable eddy dissipation time scale on t-levs[s]

    real( kind = core_rknd ), intent(in) :: &
      x_tol    ! Tolerance value of x                           [(x units)]

    logical, intent(in) :: &
      l_lmm_stepping    ! Apply Linear Multistep Method (LMM) Stepping

    type (stats_metadata_type), intent(in) :: &
      stats_metadata

    ! ----------------------- Input/Output Variable -----------------------
    type (stats), target, dimension(ngrdcol), intent(inout) :: &
      stats_zt
      
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(inout) :: &
      xp3    ! <x'^3> (thermodynamic levels)    [(x units)^3]

    ! ----------------------- Local Variables -----------------------
    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      xp3_old    ! Saved <x'^3> (thermodynamic levels)    [(x units)^3]

    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      xm_zm,   & ! Mean of x interpolated to momentum levels     [(x units)]
      xp2_zt,  & ! Variance of x interpolated to thermo. levels  [(x units)^2]
      term_tp, & ! <x'^3> turbulent production term              [(x units)^3/s]
      term_ac    ! <x'^3> accumulation term                      [(x units)^3/s]

    integer :: &
      i, k, km1     ! Grid indices

    integer :: &
      ixp3_bt, & ! Budget statistics index for <x'^3> time tendency
      ixp3_tp, & ! Budget statistics index for <x'^3> turbulent production
      ixp3_ac, & ! Budget statistics index for <x'^3> accumulation
      ixp3_dp    ! Budget statistics index for <x'^3> dissipation

    ! Coefficient in the <x'^3> turbulent dissipation term    [-]
    real( kind = core_rknd ), parameter :: &
      C_xp3_dissipation = 1.0_core_rknd

    ! Flag to either predict <x'^3> or use steady-state approximation.
    logical, parameter :: &
      l_predict_xp3 = .false.

    ! ----------------------- Begin Code -----------------------

    if ( stats_metadata%l_stats_samp ) then

      select case ( solve_type )
      case( xp3_rtp3 )
        ! Budget stats for rtp3
        ixp3_bt = stats_metadata%irtp3_bt
        ixp3_tp = stats_metadata%irtp3_tp
        ixp3_ac = stats_metadata%irtp3_ac
        ixp3_dp = stats_metadata%irtp3_dp
      case( xp3_thlp3 )
        ! Budget stats for thlp3
        ixp3_bt = stats_metadata%ithlp3_bt
        ixp3_tp = stats_metadata%ithlp3_tp
        ixp3_ac = stats_metadata%ithlp3_ac
        ixp3_dp = stats_metadata%ithlp3_dp
      case default
        ! Budgets aren't setup for the passive scalars
        ixp3_bt = 0
        ixp3_tp = 0
        ixp3_ac = 0
        ixp3_dp = 0
      end select ! solve_type

      if ( l_predict_xp3 ) then
        do i = 1, ngrdcol
          call stat_begin_update( nz, ixp3_bt, xp3(i,:) / dt, & ! Intent(in)
                                  stats_zt(i) )                 ! Intent(inout)
        end do
      end if ! l_predict_xp3

    end if ! stats_metadata%l_stats_samp

    ! Initialize variables
    term_tp = zero
    term_ac = zero

    ! Interpolate <x> to momentum levels.
    xm_zm = zt2zm( nz, ngrdcol, gr, xm )

    ! Interpolate <x'^2> to thermodynamic levels.
    xp2_zt = max( zm2zt( nz, ngrdcol, gr, xp2 ), x_tol**2 )  ! Positive definite quantity

    do k = 2, nz-1, 1
      do i = 1, ngrdcol

        ! Define the km1 index.
        km1 = max( k-1, 1 )

        ! Calculate the <x'^3> turbulent production (tp) term.
        term_tp(i,k) = term_tp_rhs( xp2_zt(i,k), wpxp(i,k), wpxp(i,km1), &
                                    rho_ds_zm(i,k), rho_ds_zm(i,km1), &
                                    invrs_rho_ds_zt(i,k), &
                                    gr%invrs_dzt(i,k) )

        ! Calculate the <x'^3> accumulation (ac) term.
        term_ac(i,k) = term_ac_rhs( xm_zm(i,k), xm_zm(i,km1), wpxp2(i,k), &
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
    end do ! k = 2, gr%nz-1, 1

    ! Set Boundary Conditions
    xp3(:,1) = zero
    xp3(:,nz) = zero

    if ( stats_metadata%l_stats_samp ) then
      do i = 1, ngrdcol
        call stat_update_var( ixp3_tp, term_tp(i,:),  & ! intent(in)
                              stats_zt(i) )             ! intent(inout)
        call stat_update_var( ixp3_ac, term_ac(i,:),  & ! intent(in)
                              stats_zt(i) )             ! intent(inout)
        call stat_update_var( ixp3_dp, -(C_xp3_dissipation * invrs_tau_zt(i,:))*xp3(i,:), & ! intent(in)
                              stats_zt(i) ) ! intent(inout)

        if ( l_predict_xp3 ) then
          call stat_end_update( nz, ixp3_bt, xp3(i,:) / dt, & ! Intent(in)
                                stats_zt(i) )                 ! Intent(inout)
        end if ! l_predict_xp3
      end do
    end if ! stats_metadata%l_stats_samp

    return

  end subroutine advance_xp3_simplified

  !=============================================================================
  function term_tp_rhs( xp2_zt, wpxp, wpxpm1, &
                             rho_ds_zm, rho_ds_zmm1, &
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
    ! =========wpxp===========rho_ds_zm=============xp2================== m(k)
    !
    ! --xp3--d( rho_ds_zm * wpxp )/dz--invrs_rho_ds_zt--xp2_zt(interp.)-- t(k)
    !
    ! =========wpxpm1=========rho_ds_zmm1===========xp2m1================ m(k-1)
    !
    ! The vertical indices m(k), t(k), and m(k-1) correspond with altitudes
    ! zm(k), zt(k), and zm(k-1), respectively.  The letter "t" is used for
    ! thermodynamic levels and the letter "m" is used for momentum levels.
    !
    ! invrs_dzt(k) = 1 / ( zm(k) - zm(k-1) )

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
      wpxp,            & ! <w'x'> at momentum level (k)           [m/s(x units)]
      wpxpm1,          & ! <w'x'> at momentum level (k-1)         [m/s(x units)]
      rho_ds_zm,       & ! Dry, static density on momentum level (k)    [kg/m^3]
      rho_ds_zmm1,     & ! Dry, static density on momentum level (k-1)  [kg/m^3]
      invrs_rho_ds_zt, & ! Inv. dry, static density at thermo. lev. (k) [m^3/kg]
      invrs_dzt          ! Inverse of grid spacing (k)                     [1/m]

    ! Return Variable
    real( kind = core_rknd ) :: &
      term_tp    ! <x'^3> turbulent production term              [(x units)^3/s]


    ! The <x'^3> turbulent production term.
    term_tp &
    = + three * xp2_zt * invrs_rho_ds_zt &
        * invrs_dzt * ( rho_ds_zm * wpxp - rho_ds_zmm1 * wpxpm1 )


    return

  end function term_tp_rhs

  !=============================================================================
  function term_ac_rhs( xm_zm, xm_zmm1, wpxp2, &
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
    ! =========================xm_zm(interp.)============================ m(k)
    !
    ! ----------xp3---------xm---------dxm_zm/dz---------wpxp2----------- t(k)
    !
    ! =========================xm_zmm1(interp.)========================== m(k-1)
    !
    ! ----------------------xmm1----------------------------------------- t(k-1)
    !
    ! The vertical indices t(k+1), m(k), t(k), m(k-1), and t(k-1) correspond
    ! with altitudes zt(k+1), zm(k), zt(k), zm(k-1), and zt(k-1), respectively.
    ! The letter "t" is used for thermodynamic levels and the letter "m" is
    ! used for momentum levels.
    !
    ! invrs_dzt(k) = 1 / ( zm(k) - zm(k-1) )

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        three    ! Variable(s)

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      xm_zm,     & ! <x> interpolated to momentum level (k)    [(x units)]
      xm_zmm1,   & ! <x> interpolated to momentum level (k-1)  [(x units)]
      wpxp2,     & ! <w'x'^2> at thermodynamic level (k)       [m/s(x units)^2]
      invrs_dzt    ! Inverse of grid spacing (k)               [1/m]

    ! Return Variable
    real( kind = core_rknd ) :: &
      term_ac    ! <x'^3> accumulation term                    [(x units)^3/s]


    ! The <x'^3> accumulation term.
    term_ac &
    = - three * wpxp2 * invrs_dzt * ( xm_zm - xm_zmm1 )


    return

  end function term_ac_rhs

  !=============================================================================

end module advance_xp3_module
