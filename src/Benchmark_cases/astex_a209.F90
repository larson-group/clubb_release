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
  subroutine astex_a209_tndcy( sclr_dim, edsclr_dim, sclr_idx, &
                               gr, wm_zt, wm_zm,  & 
                               thlm_forcing, rtm_forcing, &
                               sclrm_forcing, edsclrm_forcing )

    ! Description:
    !   Subroutine to set theta and water tendencies for ASTEX KK case
    ! References:
    !   http://www.euclipse.nl/wp3/ASTEX_Lagrangian/Introduction.shtml
    !----------------------------------------------------------------------

    use grid_class, only: &
      grid ! Type

    use grid_class, only: &
      zt2zm_api ! Procedure(s)

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    use array_index, only: &
      sclr_idx_type

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      sclr_dim, & 
      edsclr_dim

    type (sclr_idx_type), intent(in) :: &
      sclr_idx

    type (grid), intent(in) :: gr

    !--------------------- Output Variables ---------------------
    real( kind = core_rknd ), intent(out), dimension(gr%nzt) ::  & 
      wm_zt,         & ! w wind on the thermodynamic grid        [m/s]
      thlm_forcing,  & ! Liquid potential temperature tendency   [K/s]
      rtm_forcing      ! Total water mixing ratio tendency       [kg/kg/s]

    real( kind = core_rknd ), intent(out), dimension(gr%nzm) ::  & 
      wm_zm            ! w wind on the momentum grid             [m/s]

    real( kind = core_rknd ), intent(out), dimension(gr%nzt,sclr_dim) ::  & 
      sclrm_forcing ! Passive scalar forcing  [units/s]

    real( kind = core_rknd ), intent(out), dimension(gr%nzt,edsclr_dim) ::  & 
      edsclrm_forcing ! Passive scalar forcing  [units/s]

    !--------------------- Local Variables ---------------------

    integer :: i

    !--------------------- Begin Code ---------------------

    ! Compute large-scale subsidence

    do i=1,gr%nzt

      wm_zt(i) = - 5.e-6_core_rknd * gr%zt(1,i)

    end do

    ! Interpolate to momentum levels
    wm_zm = zt2zm_api( gr, wm_zt )

    ! Boundary conditions on zm
    wm_zm(1) = 0.0_core_rknd        ! At surface
    wm_zm(gr%nzm) = 0.0_core_rknd  ! Model top

    ! Radiative theta-l tendency

    thlm_forcing = 0.0_core_rknd

    ! Large scale advective moisture tendency

    rtm_forcing = 0.0_core_rknd

    ! Test scalars with thetal and rt if desired
    if ( sclr_idx%iisclr_thl > 0 ) sclrm_forcing(:,sclr_idx%iisclr_thl) = thlm_forcing
    if ( sclr_idx%iisclr_rt  > 0 ) sclrm_forcing(:,sclr_idx%iisclr_rt)  = rtm_forcing

    if ( sclr_idx%iiedsclr_thl > 0 ) edsclrm_forcing(:,sclr_idx%iiedsclr_thl) = thlm_forcing
    if ( sclr_idx%iiedsclr_rt  > 0 ) edsclrm_forcing(:,sclr_idx%iiedsclr_rt)  = rtm_forcing

    return
  end subroutine astex_a209_tndcy

  !----------------------------------------------------------------------
  subroutine astex_a209_sfclyr( ngrdcol, time, ubar, rtm, thlm, &
                                lowestlevel, exner_sfc, p_sfc, & 
                                saturation_formula, &
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

    use sfc_flux, only: compute_wprtp_sfc, compute_wpthlp_sfc   !Procedure(s)

    use saturation, only: sat_mixrat_liq_api ! Procedure(s)

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

    integer, intent(in) :: &
      ngrdcol

    real(kind=time_precision), intent(in) :: &
      time     ! Current time

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) ::  & 
      ubar,        & ! Mean sfc wind speed
      rtm,         & ! This is rt at the lowest above-ground model level.  [kg/kg]
      thlm,        & ! This is theta-l at the lowest above-ground model level.  
                 ! (DOES THIS NEED A CORRECTION FOR THETA-L TO THETA?)  [K]
      lowestlevel, & ! This is z at the lowest above-ground model level.  [m]
      exner_sfc,   & ! This is the surface pressure [Pa].
      p_sfc           ! Sea surface pressure [Pa].
      
    integer, intent(in) :: &
      saturation_formula ! Integer that stores the saturation formula to be used

    ! Output variables

    real( kind = core_rknd ), dimension(ngrdcol), intent(out) ::  & 
      wpthlp_sfc,   & ! w'th_l' at (1)   [(m K)/s]  
      wprtp_sfc,    & ! w'r_t'(1) at (1) [(m kg)/(s kg)]
      ustar,        & ! surface friction velocity     [m/s]
      T_sfc           ! Sea surface temperature [K].

    ! Local variables
    integer :: &
      before_time, after_time, i

    real( kind = core_rknd ), dimension(ngrdcol) :: &
      rsat, &
      Ch,   & ! This is C_h_20 scaled to the height of the lowest model level.
      Cq      ! This is C_q_20 scaled to the height of the lowest model level.

    real( kind = core_rknd ) :: &
      time_frac, &
      T_sfc_interp

    !-----------------BEGIN CODE-------------------------

    !$acc enter data create( rsat, Ch, Cq )

    !sensible_heat_flx = 10.0_core_rknd
    !latent_heat_flx = 25.0_core_rknd

    ! Use time_select to determine the time indexes before and after time
    ! and to calculate the time fraction necessary for linear_interp_factor
    call time_select(time, ntimes, time_sfc_given, &
                before_time, after_time, time_frac)

    ! Interpolate the value for T_sfc based on time.
    T_sfc_interp = linear_interp_factor( time_frac, T_sfc_given(after_time), &
                                         T_sfc_given(before_time) )

    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      T_sfc(i) = T_sfc_interp

      ! We set ustar as it is set in rico
      ustar(i) = 0.155_core_rknd

      rsat(i) = sat_mixrat_liq_api( p_sfc(i), T_sfc(i), saturation_formula )

      ! (Stevens, et al. 2000, eq 3)
      ! Modification in case lowest model level isn't at 10 m, from ATEX specification
      Ch(i)   = C_h_20 * ((log(standard_flux_alt/z0))/(log(lowestlevel(i)/z0))) * &
            ((log(standard_flux_alt/z0))/(log(lowestlevel(i)/z0)))

      ! Modification in case lowest model level isn't at 10 m, from ATEX specification
      Cq(i)   = C_q_20 * ((log(standard_flux_alt/z0))/(log(lowestlevel(i)/z0))) * &
            ((log(standard_flux_alt/z0))/(log(lowestlevel(i)/z0)))
    end do

    ! Compute heat and moisture fluxes
    call compute_wpthlp_sfc( ngrdcol, Ch, ubar, thlm, T_sfc, exner_sfc, &
                             wpthlp_sfc ) 

    call compute_wprtp_sfc( ngrdcol, Cq, ubar, rtm, rsat, &
                            wprtp_sfc )

    !wpthlp_sfc = sensible_heat_flx / ( rho_sfc * Cp )
    !wprtp_sfc  = latent_heat_flx / ( rho_sfc * Lv )

    ! Momentum fluxes are computed elsewhere

    !$acc exit data delete( rsat, Ch, Cq )

    return

  end subroutine astex_a209_sfclyr

end module astex_a209
