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
  subroutine rico_tndcy( sclr_dim, edsclr_dim, sclr_idx, &
                         gr, rtm, exner, &
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

  use grid_class, only: &
    grid

  use array_index, only: &
    sclr_idx_type

  use constants_clubb, only: &
    g_per_kg ! Variable(s)

  use spec_hum_to_mixing_ratio, only: &
      force_spec_hum_to_mixing_ratio ! Procedure(s)

  use clubb_precision, only: &
    core_rknd ! Variable(s)


  implicit none

  !--------------------- Input Variables ---------------------
  integer, intent(in) :: &
    sclr_dim, & 
    edsclr_dim

  type (sclr_idx_type), intent(in) :: &
    sclr_idx

  type (grid), target, intent(in) :: &
    gr

  real( kind = core_rknd ), dimension(gr%nzt), intent(in) :: &
    rtm,   & ! Mean total water mixing ratio    [kg/kg]
    exner    ! Exner function                   [-]

  !--------------------- Output Variables ---------------------
  real( kind = core_rknd ), dimension(gr%nzt), intent(out) :: & 
    thlm_forcing, & ! Large-scale thlm tendency               [K s^-1]
    rtm_forcing     ! Large-scale rtm tendency                [kg kg^-1 s^-1]

  real( kind = core_rknd ), intent(out), dimension(gr%nzt,sclr_dim) :: & 
    sclrm_forcing ! Passive scalar LS tendency            [units/s]

  real( kind = core_rknd ), intent(out), dimension(gr%nzt,edsclr_dim) :: & 
    edsclrm_forcing ! Passive eddy-scalar LS tendency     [units/s]

  !--------------------- Local Variables ---------------------
  integer :: k          ! Loop index

  real( kind = core_rknd )    :: &
    t_tendency,  & ! Temperature (not potential temperature) tendency [K s^-1]
    qtm_forcing    ! Large-scale forcing of specific humidity         [kg/kg/s]

  !--------------------- Begin Code ---------------------

  ! Compute large-scale horizontal temperature advection
  ! NEW-- "And Radiation"... 15 Dec 2006, Michael Falk
  ! Equations located in 1D models > Set up short composite run on reference site
  do k=1,gr%nzt
    if (gr%zt(1,k) < 4000._core_rknd ) then
      t_tendency = -2.51_core_rknd / 86400._core_rknd + & 
        (-2.18_core_rknd + 2.51_core_rknd) / (86400._core_rknd*4000._core_rknd) &
        * gr%zt(1,k)  ! Units [K s^-1] - known magic number
    else if (gr%zt(1,k) < 5000._core_rknd ) then
      t_tendency = -2.18_core_rknd / 86400._core_rknd + & 
        (2.18_core_rknd) / (86400._core_rknd*(5000._core_rknd-4000._core_rknd)) &
        * (gr%zt(1,k)-4000._core_rknd)  ! Units [K s^-1] - known magic number
    else
      t_tendency = 0._core_rknd  ! Units [K s^-1]
    end if
    ! Convert to units of [K s^-1] but potential T instead of T
!          thlm_forcing(k) = (t_tendency * ((p_sfc/p(k)) ** (Rd/Cp)))
    thlm_forcing(k) = (t_tendency / exner(k))
  end do


  ! Compute large-scale horizontal moisture advection [g kg^-1 s^-1]
  ! Equations located in 1D models > Set up short composite run on reference site
  do k=1,gr%nzt
    if (gr%zt(1,k) < 3000._core_rknd) then
      qtm_forcing = - 1.0_core_rknd / 86400._core_rknd + & 
        (0.345_core_rknd+1.0_core_rknd) / (86400._core_rknd * 3000._core_rknd) &
        * gr%zt(1,k)  ! Units [g kg^-1 s^-1] - known magic number
    else if (gr%zt(1,k) < 4000._core_rknd ) then
      qtm_forcing = 0.345_core_rknd / 86400._core_rknd  
                ! Units [g kg^-1 s^-1] - known magic number
    else if (gr%zt(1,k) < 5000._core_rknd ) then
      qtm_forcing = 0.345_core_rknd / 86400._core_rknd + & 
       (-0.345_core_rknd) / (86400._core_rknd*(5000._core_rknd-4000._core_rknd)) &
       * (gr%zt(1,k)-4000._core_rknd)! Units [g kg^-1 s^-1] known magic number
    else
      qtm_forcing = 0._core_rknd  ! Units [g kg^-1 s^-1]
    end if
    qtm_forcing = qtm_forcing / g_per_kg  ! Converts [g kg^-1 s^-1] to [kg kg^-1 s^-1]

    ! Convert forcings from terms of total water specific humidity to terms of
    ! total water mixing ratio.
    call force_spec_hum_to_mixing_ratio( rtm(k), qtm_forcing, rtm_forcing(k) )

  end do

  ! Test scalars with thetal and rt if desired
  if ( sclr_idx%iisclr_thl > 0 ) sclrm_forcing(:,sclr_idx%iisclr_thl) = thlm_forcing
  if ( sclr_idx%iisclr_rt  > 0 ) sclrm_forcing(:,sclr_idx%iisclr_rt)  = rtm_forcing

  if ( sclr_idx%iiedsclr_thl > 0 ) edsclrm_forcing(:,sclr_idx%iiedsclr_thl) = thlm_forcing
  if ( sclr_idx%iiedsclr_rt  > 0 ) edsclrm_forcing(:,sclr_idx%iiedsclr_rt)  = rtm_forcing

  end subroutine rico_tndcy
 !----------------------------------------------------------------------


 !----------------------------------------------------------------------
  subroutine rico_sfclyr( time, um_sfc, vm_sfc, thlm, rtm, &
                          lowestlevel, p_sfc, exner_sfc, & 
                          saturation_formula, &
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

  use saturation, only: sat_mixrat_liq ! Procedure(s)

  use sfc_flux, only: compute_ubar, compute_momentum_flux, &
                          compute_wpthlp_sfc, compute_wprtp_sfc

  use time_dependent_input, only: time_sfc_given, T_sfc_given, &  ! Variable(s)
                                  time_select                     ! Procedure(s)

  use interpolation, only: linear_interp_factor   ! Procedure(s)

  use clubb_precision, only: time_precision, core_rknd ! Variable(s)

  implicit none

  intrinsic :: max, log, sqrt

  ! Constants
  real( kind = core_rknd ), parameter :: & 
    C_10    = 0.0013_core_rknd,    & ! Drag coefficient, defined by ATEX specification
    C_m_20  = 0.001229_core_rknd,  & ! Drag coefficient, defined by RICO 3D specification
    C_h_20  = 0.001094_core_rknd,  & ! Drag coefficient, defined by RICO 3D specification
    C_q_20  = 0.001133_core_rknd,  & ! Drag coefficient, defined by RICO 3D specification
    z0      = 0.00015_core_rknd      ! Roughness length, defined by ATEX specification

  real( kind = core_rknd ), parameter :: &
    standard_flux_alt = 20._core_rknd ! default height at which surface flux is computed [m]


  ! Internal variables
  real( kind = core_rknd ) :: & 
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

  real( kind = core_rknd ), intent(in) :: & 
    um_sfc,        & ! This is u at the lowest above-ground model level.  [m/s]
    vm_sfc,        & ! This is v at the lowest above-ground model level.  [m/s]
    thlm,          & ! This is theta-l at the lowest above-ground model level.  
                     ! (DOES THIS NEED A CORRECTION FOR THETA-L TO THETA?)  [K]
    rtm,           & ! This is rt at the lowest above-ground model level.  [kg/kg]
    lowestlevel,   & ! This is z at the lowest above-ground model level.  [m]
    p_sfc,          & ! This is the surface pressure [Pa].
    exner_sfc

  integer, intent(in) :: &
    saturation_formula ! Integer that stores the saturation formula to be used

  ! Output variables
  real( kind = core_rknd ), intent(out) ::  & 
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

  T_sfc = linear_interp_factor( time_frac, T_sfc_given(after_time), &
                                T_sfc_given(before_time) )

  ! Declare the value of ustar.
  ustar = 0.3_core_rknd

  ! Choose which scheme to use
  l_use_old_atex = .FALSE.

  ! Define variable values
  ubar = compute_ubar( um_sfc, vm_sfc )

  ! (Stevens, et al. 2000, eq 3)
  ! Modification in case lowest model level isn't at 10 m, from ATEX specification
  Cz   = C_10 * ((log(10._core_rknd/z0))/(log(lowestlevel/z0))) * & 
         ((log(10._core_rknd/z0))/(log(lowestlevel/z0))) ! Known magic number        
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
    wprtp_sfc = compute_wprtp_sfc( Cz, ubar, rtm, sat_mixrat_liq(p_sfc,T_sfc,saturation_formula) )
  call compute_momentum_flux( um_sfc, vm_sfc, ubar, ustar, &
                              upwp_sfc, vpwp_sfc )
  else ! Use RICO version
    wpthlp_sfc = compute_wpthlp_sfc( Ch, ubar, thlm, T_sfc, exner_sfc )
!    wprtp_sfc  = -Cz * ubar * ( .01726 - sat_mixrat_liq(p_sfc,T_sfc) ) ! kg kg^-1  m s^-1
!    wprtp_sfc  = -Cz * ubar * ( .01626 - sat_mixrat_liq(p_sfc,T_sfc) ) ! kg kg^-1  m s^-1
    wprtp_sfc  = compute_wprtp_sfc( Cq, ubar, rtm, sat_mixrat_liq(p_sfc,T_sfc,saturation_formula) )
    upwp_sfc   = -um_sfc * Cm * ubar  ! m^2 s^-2
    vpwp_sfc   = -vm_sfc * Cm * ubar  ! m^2 s^-2

  end if

  return
  end subroutine rico_sfclyr

!----------------------------------------------------------------------

end module rico
