!----------------------------------------------------------------------
! $Id$
module astex_a209

  ! Description:
  !   Contains subroutines for the ASTEX KK case.
  ! References:
  !   http://www.euclipse.nl/wp3/ASTEX_Lagrangian/Introduction.shtml
  !----------------------------------------------------------------------

  implicit none

  public :: astex_a209_tndcy, astex_a209_sfclyr

  private ! Default Scope

  contains

  !----------------------------------------------------------------------
  subroutine astex_a209_tndcy( wm_zt, wm_zm,  & 
                               thlm_forcing, rtm_forcing, &
                               sclrm_forcing, edsclrm_forcing )

    ! Description:
    !   Subroutine to set theta and water tendencies for ASTEX KK case
    ! References:
    !   http://www.euclipse.nl/wp3/ASTEX_Lagrangian/Introduction.shtml
    !----------------------------------------------------------------------

    use parameters_model, only: sclr_dim, edsclr_dim ! Variable(s)

    use grid_class, only: gr ! Variable(s)

    use grid_class, only: zt2zm ! Procedure(s)

    use clubb_precision, only: core_rknd ! Variable(s)

    use array_index, only: iisclr_rt, iisclr_thl, iiedsclr_rt, iiedsclr_thl ! Variable(s)

    implicit none

    ! Output Variables
    real( kind = core_rknd ), intent(out), dimension(gr%nz) ::  & 
      wm_zt,         & ! w wind on the thermodynamic grid        [m/s]
      wm_zm,         & ! w wind on the momentum grid             [m/s]
      thlm_forcing,  & ! Liquid potential temperature tendency   [K/s]
      rtm_forcing      ! Total water mixing ratio tendency       [kg/kg/s]

    real( kind = core_rknd ), intent(out), dimension(gr%nz,sclr_dim) ::  & 
      sclrm_forcing ! Passive scalar forcing  [units/s]

    real( kind = core_rknd ), intent(out), dimension(gr%nz,edsclr_dim) ::  & 
      edsclrm_forcing ! Passive scalar forcing  [units/s]

    ! Local variables

    integer :: i

    ! Compute large-scale subsidence

    do i=2,gr%nz

      wm_zt(i) = - 5.e-6_core_rknd * gr%zt(i)

    end do

    ! Lower Boundary condition on zt
    wm_zt(1) = 0.0_core_rknd        ! Below surface

    ! Interpolate to momentum levels
    wm_zm = zt2zm( wm_zt )

    ! Boundary conditions on zm
    wm_zm(1) = 0.0_core_rknd        ! At surface
    wm_zm(gr%nz) = 0.0_core_rknd  ! Model top

    ! Radiative theta-l tendency

    thlm_forcing = 0.0_core_rknd

    ! Large scale advective moisture tendency

    rtm_forcing = 0.0_core_rknd

    ! Test scalars with thetal and rt if desired
    if ( iisclr_thl > 0 ) sclrm_forcing(:,iisclr_thl) = thlm_forcing
    if ( iisclr_rt  > 0 ) sclrm_forcing(:,iisclr_rt)  = rtm_forcing

    if ( iiedsclr_thl > 0 ) edsclrm_forcing(:,iiedsclr_thl) = thlm_forcing
    if ( iiedsclr_rt  > 0 ) edsclrm_forcing(:,iiedsclr_rt)  = rtm_forcing

    return
  end subroutine astex_a209_tndcy

  !----------------------------------------------------------------------
  subroutine astex_a209_sfclyr( time, ubar, rtm, thlm, &
                                lowestlevel, exner_sfc, p_sfc, & 
                                wpthlp_sfc, wprtp_sfc, ustar, T_sfc )

    ! Description:
    !   This subroutine computes surface fluxes of horizontal momentum,
    !   heat and moisture according to ASTEX with Khairoutdinov and Kogan
    !   alteration.

    ! References:
    !   http://www.euclipse.nl/wp3/ASTEX_Lagrangian/Introduction.shtml
    !----------------------------------------------------------------------

    use time_dependent_input, only: T_sfc_given, time_sfc_given, &
                                    time_select ! Variable(s)

    use surface_flux, only: compute_wprtp_sfc, compute_wpthlp_sfc   !Procedure(s)

    use saturation, only: sat_mixrat_liq ! Procedure(s)

    use interpolation, only: linear_interp_factor ! Procedure(s)

    use clubb_precision, only: time_precision, core_rknd ! Variable(s)

    implicit none

    ! Parameter Constants
    integer, parameter :: &
      ntimes = 41

    ! Constants (taken from the rico case)
    real( kind = core_rknd ), parameter :: &
      C_h_20  = 0.001094_core_rknd,  & ! Drag coefficient, defined by RICO 3D specification
      C_q_20  = 0.001133_core_rknd,  & ! Drag coefficient, defined by RICO 3D specification
      z0      = 0.00015_core_rknd      ! Roughness length, defined by ATEX specification

    real( kind = core_rknd ), parameter :: &
      standard_flux_alt = 20._core_rknd ! default height at which the surface flux is computed [m]

    ! Input variables

    real(kind=time_precision), intent(in) :: time     ! Current time

    real( kind = core_rknd ), intent(in) ::  & 
      ubar,        & ! Mean sfc wind speed
      rtm,         & ! This is rt at the lowest above-ground model level.  [kg/kg]
      thlm,        & ! This is theta-l at the lowest above-ground model level.  
                 ! (DOES THIS NEED A CORRECTION FOR THETA-L TO THETA?)  [K]
      lowestlevel, & ! This is z at the lowest above-ground model level.  [m]
      exner_sfc,   & ! This is the surface pressure [Pa].
      p_sfc           ! Sea surface pressure [Pa].
      

    ! Output variables

    real( kind = core_rknd ), intent(out) ::  & 
      wpthlp_sfc,   & ! w'th_l' at (1)   [(m K)/s]  
      wprtp_sfc,    & ! w'r_t'(1) at (1) [(m kg)/(s kg)]
      ustar,        & ! surface friction velocity     [m/s]
      T_sfc           ! Sea surface temperature [K].


    ! Local variables
    integer :: &
      before_time, after_time

    real( kind = core_rknd ) :: &
      Ch,   & ! This is C_h_20 scaled to the height of the lowest model level.
      Cq,   & ! This is C_q_20 scaled to the height of the lowest model level.
      time_frac

    !-----------------BEGIN CODE-------------------------


    ! Compute heat and moisture fluxes

    ! (Stevens, et al. 2000, eq 3)
    ! Modification in case lowest model level isn't at 10 m, from ATEX specification
    Ch   = C_h_20 * ((log(standard_flux_alt/z0))/(log(lowestlevel/z0))) * &
           ((log(standard_flux_alt/z0))/(log(lowestlevel/z0)))
    ! Modification in case lowest model level isn't at 10 m, from ATEX specification
    Cq   = C_q_20 * ((log(standard_flux_alt/z0))/(log(lowestlevel/z0))) * &
           ((log(standard_flux_alt/z0))/(log(lowestlevel/z0)))


    !sensible_heat_flx = 10.0_core_rknd
    !latent_heat_flx = 25.0_core_rknd


    T_sfc = 0.0_core_rknd

    ! We set ustar as it is set in rico
    ustar = 0.155_core_rknd

    ! Use time_select to determine the time indexes before and after time
    ! and to calculate the time fraction necessary for linear_interp_factor
    call time_select(time, ntimes, time_sfc_given, &
                before_time, after_time, time_frac)

    ! Interpolate the value for T_sfc based on time.
    T_sfc = linear_interp_factor( time_frac, T_sfc_given(after_time), &
                                  T_sfc_given(before_time) )
   
    wpthlp_sfc = compute_wpthlp_sfc( Ch, ubar, thlm, T_sfc, exner_sfc )
    wprtp_sfc  = compute_wprtp_sfc( Cq, ubar, rtm, sat_mixrat_liq(p_sfc,T_sfc) )

    !wpthlp_sfc = sensible_heat_flx / ( rho_sfc * Cp )
    !wprtp_sfc  = latent_heat_flx / ( rho_sfc * Lv )

    ! Momentum fluxes are computed elsewhere

    return
  end subroutine astex_a209_sfclyr

end module astex_a209
