!$Id$
!----------------------------------------------------------------------
module rico
!
! Description:
!   Contains subroutines for the RICO case.
!----------------------------------------------------------------------

  implicit none

  public :: rico_tndcy, rico_sfclyr

  private  ! Default Scope

  contains

!----------------------------------------------------------------------
  subroutine rico_tndcy( exner, &
                         thlm_forcing, rtm_forcing, & 
                         sclrm_forcing, edsclrm_forcing )
!
! Description:
!   Subroutine to apply case-specific forcings to RICO case
!   (Michael Falk, 13 Dec 2006).
!
! References:
!   ATEX: http://www.atmos.ucla.edu/~bstevens/gcss/setup.html
!   RICO: http://www.knmi.nl/samenw/rico/setup3d.html
!-----------------------------------------------------------------------

  use parameters_model, only: sclr_dim, edsclr_dim ! Variable(s)

  use grid_class, only: gr ! Variable(s)

  use grid_class, only: zt2zm ! Procedure(s)

  use array_index, only: iisclr_rt, iisclr_thl, iiedsclr_rt, iiedsclr_thl ! Variable(s)

  use constants_clubb, only: g_per_kg ! Variable(s)

 
!  use stats_variables


  implicit none

  ! Input Variables

  real, dimension(gr%nzmax), intent(in) :: & 
  exner    ! Exner function                         [-]

  ! Output Variables
  real, dimension(gr%nzmax), intent(out) :: & 
    thlm_forcing, & ! Large-scale thlm tendency               [K s^-1]
    rtm_forcing     ! Large-scale rtm tendency                [kg kg^-1 s^-1]

  real, intent(out), dimension(gr%nzmax,sclr_dim) :: & 
    sclrm_forcing ! Passive scalar LS tendency            [units/s]

  real, intent(out), dimension(gr%nzmax,edsclr_dim) :: & 
    edsclrm_forcing ! Passive eddy-scalar LS tendency     [units/s]

  ! Local Variables, general
  integer :: k          ! Loop index
  real    :: t_tendency ! Temperature (not potential temperature) tendency [K s^-1]

  ! Compute large-scale horizontal temperature advection
  ! NEW-- "And Radiation"... 15 Dec 2006, Michael Falk
  ! Equations located in 1D models > Set up short composite run on reference site
  do k=1,gr%nzmax
    if (gr%zt(k) < 4000. ) then
      t_tendency = -2.51 / 86400. + & 
        (-2.18 + 2.51) / (86400.*4000.) * gr%zt(k)  ! Units [K s^-1] - known magic number
    else if (gr%zt(k) < 5000. ) then
      t_tendency = -2.18 / 86400. + & 
        (2.18) / (86400.*(5000.-4000.)) * (gr%zt(k)-4000.)  ! Units [K s^-1] - known magic number
    else
      t_tendency = 0.  ! Units [K s^-1]
    end if
    ! Convert to units of [K s^-1] but potential T instead of T
!          thlm_forcing(k) = (t_tendency * ((p_sfc/p(k)) ** (Rd/Cp)))
    thlm_forcing(k) = (t_tendency / exner(k))
  end do


  ! Compute large-scale horizontal moisture advection [g kg^-1 s^-1]
  ! Equations located in 1D models > Set up short composite run on reference site
  do k=1,gr%nzmax
    if (gr%zt(k) < 3000.) then
      rtm_forcing(k) = - 1.0 / 86400. + & 
        (0.345+1.0) / (86400. * 3000.) * gr%zt(k)  ! Units [g kg^-1 s^-1] - known magic number
    else if (gr%zt(k) < 4000. ) then
      rtm_forcing(k) = 0.345 / 86400.  ! Units [g kg^-1 s^-1] - known magic number
    else if (gr%zt(k) < 5000. ) then
      rtm_forcing(k) = 0.345 / 86400. + & 
       (-0.345) / (86400.*(5000.-4000.)) * (gr%zt(k)-4000.)! Units [g kg^-1 s^-1] known magic number
    else
      rtm_forcing(k) = 0.  ! Units [g kg^-1 s^-1]
    end if
    rtm_forcing(k) = rtm_forcing(k) / g_per_kg  ! Converts [g kg^-1 s^-1] to [kg kg^-1 s^-1]
  end do

  ! Test scalars with thetal and rt if desired
  if ( iisclr_thl > 0 ) sclrm_forcing(:,iisclr_thl) = thlm_forcing
  if ( iisclr_rt  > 0 ) sclrm_forcing(:,iisclr_rt)  = rtm_forcing

  if ( iiedsclr_thl > 0 ) edsclrm_forcing(:,iiedsclr_thl) = thlm_forcing
  if ( iiedsclr_rt  > 0 ) edsclrm_forcing(:,iiedsclr_rt)  = rtm_forcing

  end subroutine rico_tndcy
 !----------------------------------------------------------------------


 !----------------------------------------------------------------------
  subroutine rico_sfclyr( time, um_sfc, vm_sfc, thlm, rtm, &
                          lowestlevel, p_sfc, exner_sfc, & 
                          upwp_sfc, vpwp_sfc, wpthlp_sfc, & 
                          wprtp_sfc, ustar, T_sfc )
  !----------------------------------------------------------------------
  !        Description:
  !          Surface forcing subroutine for RICO case.  Written
  !          December 2006 by Michael Falk.
  !
  !          Updated to use specific formulations for surface fluxes
  !          as specified in the RICO 3D LES specification, in hopes that
  !          they'll be more accurate.
  !
  !        References:
  !          ATEX: http://www.atmos.ucla.edu/~bstevens/gcss/setup.html
  !          RICO: http://www.knmi.nl/samenw/rico/setup3d.html
  !-----------------------------------------------------------------------

  use constants_clubb, only: kappa, p0 ! Variable(s)
  
  use saturation, only: sat_mixrat_liq ! Procedure(s)

  use surface_flux, only: compute_ubar, compute_momentum_flux, &
                          compute_wpthlp_sfc, compute_wprtp_sfc

  use time_dependent_input, only: time_sfc_given, T_sfc_given, &  ! Variable(s)
                                  time_select                     ! Procedure(s)

  use interpolation, only: factor_interp   ! Procedure(s)

  use stats_precision, only: time_precision ! Variable(s)

  implicit none

  intrinsic :: max, log, sqrt

  ! Constants
  real, parameter :: & 
    C_10    = 0.0013,    & ! Drag coefficient, defined by ATEX specification
    C_m_20  = 0.001229,  & ! Drag coefficient, defined by RICO 3D specification
    C_h_20  = 0.001094,  & ! Drag coefficient, defined by RICO 3D specification
    C_q_20  = 0.001133,  & ! Drag coefficient, defined by RICO 3D specification
    z0      = 0.00015      ! Roughness length, defined by ATEX specification

  real, parameter :: &
    standard_flux_alt = 20. ! default height at which surface flux is computed [m]


  ! Internal variables
  real :: & 
    ubar, & ! This is root (u^2 + v^2), per ATEX and RICO spec.
    Cz,   & ! This is C_10 scaled to the height of the lowest model level.
    Cm,   & ! This is C_m_20 scaled to the height of the lowest model level.
    Ch,   & ! This is C_h_20 scaled to the height of the lowest model level.
    Cq,   & ! This is C_q_20 scaled to the height of the lowest model level.
    time_frac ! The time fraction used for interpolation

  integer :: &
      before_time, after_time  ! time indexes used for interpolation

  logical :: & 
    l_use_old_atex  ! if true, use ATEX version; if not, use RICO-specific

  ! Input variables
  real(time_precision), intent(in) :: &
    time ! the current time

  real, intent(in) :: & 
    um_sfc,        & ! This is u at the lowest above-ground model level.  [m/s]
    vm_sfc,        & ! This is v at the lowest above-ground model level.  [m/s]
    thlm,          & ! This is theta-l at the lowest above-ground model level.  
                     ! (DOES THIS NEED A CORRECTION FOR THETA-L TO THETA?)  [K]
    rtm,           & ! This is rt at the lowest above-ground model level.  [kg/kg]
    lowestlevel,   & ! This is z at the lowest above-ground model level.  [m]
    p_sfc,          & ! This is the surface pressure [Pa].
    exner_sfc

  ! Output variables
  real, intent(out) ::  & 
    upwp_sfc,   & ! The upward flux of u-momentum         [(m^2 s^-2]
    vpwp_sfc,   & ! The Upward flux of v-momentum         [(m^2 s^-2]
    wpthlp_sfc, & ! The upward flux of theta-l            [K m s^-1]
    wprtp_sfc,  & ! The upward flux of rtm (total water)  [kg kg^-1 m s^-1]
    ustar,      & ! surface friction velocity             [m/s]
    T_sfc         ! This is the sea surface temperature   [K]

  !--------------------BEGIN CODE----------------------------

  ! interpolate variables from time_dependent_input

  call time_select( time, size(time_sfc_given), time_sfc_given, &
                       before_time, after_time, time_frac )

  T_sfc = factor_interp( time_frac, T_sfc_given(after_time), &
                                       T_sfc_given(before_time) )

  ! Declare the value of ustar.
  ustar = 0.3

  ! Choose which scheme to use
  l_use_old_atex = .FALSE.

  ! Define variable values
  ubar = compute_ubar( um_sfc, vm_sfc )

  ! (Stevens, et al. 2000, eq 3)
  ! Modification in case lowest model level isn't at 10 m, from ATEX specification
  Cz   = C_10 * ((log(10./z0))/(log(lowestlevel/z0))) * & 
         ((log(10./z0))/(log(lowestlevel/z0))) ! Known magic number        
  ! Modification in case lowest model level isn't at 10 m, from ATEX specification
  Cm   = C_m_20 * ((log(standard_flux_alt/z0))/(log(lowestlevel/z0))) * & 
         ((log(standard_flux_alt/z0))/(log(lowestlevel/z0)))
  ! Modification in case lowest model level isn't at 10 m, from ATEX specification
  Ch   = C_h_20 * ((log(standard_flux_alt/z0))/(log(lowestlevel/z0))) * & 
         ((log(standard_flux_alt/z0))/(log(lowestlevel/z0)))
         ! Modification in case lowest model level isn't at 10 m, from ATEX specification
  Cq   = C_q_20 * ((log(standard_flux_alt/z0))/(log(lowestlevel/z0))) * & 
         ((log(standard_flux_alt/z0))/(log(lowestlevel/z0)))

! Compute heat and moisture fluxes
  if (l_use_old_atex) then ! Use ATEX version
    wpthlp_sfc = compute_wpthlp_sfc( Cz, ubar, thlm, T_sfc, exner_sfc )
    wprtp_sfc = compute_wprtp_sfc( Cz, ubar, rtm, sat_mixrat_liq(p_sfc,T_sfc) )
  call compute_momentum_flux( um_sfc, vm_sfc, ubar, ustar, &
                              upwp_sfc, vpwp_sfc )
  else ! Use RICO version
    wpthlp_sfc = compute_wpthlp_sfc( Ch, ubar, thlm, T_sfc, exner_sfc )
!    wprtp_sfc  = -Cz * ubar * ( .01726 - sat_mixrat_liq(p_sfc,T_sfc) ) ! kg kg^-1  m s^-1
!    wprtp_sfc  = -Cz * ubar * ( .01626 - sat_mixrat_liq(p_sfc,T_sfc) ) ! kg kg^-1  m s^-1
    wprtp_sfc  = compute_wprtp_sfc( Cq, ubar, rtm, sat_mixrat_liq(p_sfc,T_sfc) )
    upwp_sfc   = -um_sfc * Cm * ubar  ! m^2 s^-2
    vpwp_sfc   = -vm_sfc * Cm * ubar  ! m^2 s^-2

  end if

  return
  end subroutine rico_sfclyr

!----------------------------------------------------------------------

end module rico
