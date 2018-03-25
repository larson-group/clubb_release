! $Id$
!===============================================================================
module advance_xp3_module

  ! Description:
  ! Predicts the value of <x'^3> for <rt'^3>, <thl'^3>, and <sclr'^3>.

  ! References:
  !-------------------------------------------------------------------------

  implicit none

  public :: advance_xp3_simplified    ! Procedure(s)

  private :: term_tp_rhs, & ! Procedure(s)
             term_ac_rhs

  private ! default scope

  contains

  !=============================================================================
  subroutine advance_xp3_simplified( dt, xm, xp2,     & ! Intent(in)
                                     wpxp, wpxp2,     & ! Intent(in)
                                     rho_ds_zm,       & ! Intent(in)
                                     invrs_rho_ds_zt, & ! Intent(in)
                                     tau_zt, x_tol,   & ! Intent(in)
                                     xp3              ) ! Intent(inout)

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
        gr,    & ! Variable Type
        zm2zt, & ! Procedure(s)
        zt2zm

    use constants_clubb, only: &
        one,  & ! Variable(s)
        zero

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      dt                 ! Model timestep                            [s]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      xm,              & ! Mean (overall) of x (thermo. levels) [(x units)]
      xp2,             & ! Variance (overall) of x (m-levs.)    [(x units)^2]
      wpxp,            & ! Turbulent flux of x (momentum levs.) [m/s(x units)]
      wpxp2,           & ! <w'x'^2> (thermodynamic levels)      [m/s(x units)^2]
      rho_ds_zm,       & ! Dry, static density on momentum levels      [kg/m^3]
      invrs_rho_ds_zt, & ! Inv. dry, static density at thermo. levels  [m^3/kg]
      tau_zt             ! Time-scale tau on thermodynamic levels      [s]

    real( kind = core_rknd ), intent(in) :: &
      x_tol    ! Tolerance value of x                           [(x units)]

    ! Input/Output Variable
    real( kind = core_rknd ), dimension(gr%nz), intent(inout) :: &
      xp3    ! <x'^3> (thermodynamic levels)    [(x units)^3]

    ! Local Variables
    real( kind = core_rknd ), dimension(gr%nz) :: &
      xm_zm,   & ! Mean of x interpolated to momentum levels     [(x units)]
      xp2_zt,  & ! Variance of x interpolated to thermo. levels  [(x units)^2]
      term_tp, & ! <x'^3> turbulent production term              [(x units)^3/s]
      term_ac    ! <x'^3> accumulation term                      [(x units)^3/s]

    integer :: &
      k, km1     ! Grid indices

    ! Coefficient in the <x'^3> turbulent dissipation term    [-]
    real( kind = core_rknd ), parameter :: &
      C_xp3_dissipation = 1.0_core_rknd

    ! Flag to either predict <x'^3> or use steady-state approximation.
    logical, parameter :: &
      l_predict_xp3 = .false.


    ! Initialize variables
    term_tp = zero
    term_ac = zero

    ! Interpolate <x> to momentum levels.
    xm_zm = zt2zm( xm )

    ! Interpolate <x'^2> to thermodynamic levels.
    xp2_zt = max( zm2zt( xp2 ), x_tol**2 )  ! Positive definite quantity

    do k = 2, gr%nz-1, 1

      ! Define the km1 index.
      km1 = max( k-1, 1 )

      ! Calculate the <x'^3> turbulent production (tp) term.
      term_tp(k) = term_tp_rhs( xp2_zt(k), wpxp(k), wpxp(km1), &
                                rho_ds_zm(k), rho_ds_zm(km1), &
                                invrs_rho_ds_zt(k), &
                                gr%invrs_dzt(k) )

      ! Calculate the <x'^3> accumulation (ac) term.
      term_ac(k) = term_ac_rhs( xm_zm(k), xm_zm(km1), wpxp2(k), &
                                gr%invrs_dzt(k) )

      if ( l_predict_xp3 ) then

         ! Advance <x'^3> one time step.
         xp3(k) = ( ( xp3(k) / dt ) + term_tp(k) + term_ac(k) ) &
                  / ( ( one / dt ) + ( C_xp3_dissipation / tau_zt(k) ) )

      else

         ! Calculate <x'^3> using the steady-state approximation.
         xp3(k) = ( tau_zt(k) / C_xp3_dissipation ) &
                  * ( term_tp(k) + term_ac(k) )

      endif ! l_predict_xp3

    enddo ! k = 2, gr%nz-1, 1

    ! Set Boundary Conditions
    xp3(1) = zero
    xp3(gr%nz) = zero


    return

  end subroutine advance_xp3_simplified

  !=============================================================================
  pure function term_tp_rhs( xp2_zt, wpxp, wpxpm1, &
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
  pure function term_ac_rhs( xm_zm, xm_zmm1, wpxp2, &
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
