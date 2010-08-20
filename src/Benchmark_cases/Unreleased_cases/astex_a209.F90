!----------------------------------------------------------------------
! $Id$
module astex_a209

!       Description:
!       Contains subroutines for the ASTEX KK case.
!----------------------------------------------------------------------

implicit none

public :: astex_a209_tndcy, astex_a209_sfclyr

private ! Default Scope

contains

!----------------------------------------------------------------------
subroutine astex_a209_tndcy( wm_zt, wm_zm,  & 
                        thlm_forcing, rtm_forcing, &
                        sclrm_forcing, edsclrm_forcing )

!       Description:
!       Subroutine to set theta and water tendencies for ASTEX KK case
!       References:
!----------------------------------------------------------------------

use parameters_model, only: sclr_dim, edsclr_dim ! Variable(s)

use grid_class, only: gr ! Variable(s)

use grid_class, only: zt2zm ! Procedure(s)

use stats_precision, only: time_precision ! Variable(s)

use array_index, only: iisclr_rt, iisclr_thl, iiedsclr_rt, iiedsclr_thl ! Variable(s)

implicit none

! Output Variables
real, intent(out), dimension(gr%nnzp) ::  & 
  wm_zt,         & ! w wind on the thermodynamic grid        [m/s]
  wm_zm,         & ! w wind on the momentum grid             [m/s]
  thlm_forcing,  & ! Liquid potential temperature tendency   [K/s]
  rtm_forcing      ! Total water mixing ratio tendency       [kg/kg/s]

real, intent(out), dimension(gr%nnzp,sclr_dim) ::  & 
  sclrm_forcing ! Passive scalar forcing  [units/s]

real, intent(out), dimension(gr%nnzp,edsclr_dim) ::  & 
  edsclrm_forcing ! Passive scalar forcing  [units/s]

! Local variables

integer :: i

! Large-scale subsidence

do i=2,gr%nnzp

   wm_zt(i) = - 5.e-6 * gr%zt(i)

end do

! Boundary condition
wm_zt(1) = 0.0        ! Below surface

! Interpolation
wm_zm = zt2zm( wm_zt )

! Boundary condition
wm_zm(1) = 0.0        ! At surface
wm_zm(gr%nnzp) = 0.0  ! Model top

! Radiative theta-l tendency

thlm_forcing = 0.0

! Large scale advective moisture tendency

rtm_forcing = 0.0

! Test scalars with thetal and rt if desired
if ( iisclr_thl > 0 ) sclrm_forcing(:,iisclr_thl) = thlm_forcing
if ( iisclr_rt  > 0 ) sclrm_forcing(:,iisclr_rt)  = rtm_forcing

if ( iiedsclr_thl > 0 ) edsclrm_forcing(:,iiedsclr_thl) = thlm_forcing
if ( iiedsclr_rt  > 0 ) edsclrm_forcing(:,iiedsclr_rt)  = rtm_forcing

return
end subroutine astex_a209_tndcy

!----------------------------------------------------------------------
subroutine astex_a209_sfclyr( time, rho0, ubar, rtm, thlm, &
                         lowestlevel, exner_sfc, psfc, & 
                         wpthlp_sfc, wprtp_sfc, ustar, T_sfc )

!       Description:
!       This subroutine computes surface fluxes of horizontal momentum,
!       heat and moisture according to ASTEX with Khairoutdinov and Kogan
!       alteration.

!       References:
!----------------------------------------------------------------------

use constants_clubb, only: Cp, Lv, fstdout ! Variable(s)

use time_dependent_input, only: T_sfc_given, time_sfc_given ! Variable(s)

use surface_flux, only: compute_wprtp_sfc, compute_wpthlp_sfc   !Procedure(s)

use saturation, only: sat_mixrat_liq ! Procedure(s)

use interpolation, only: factor_interp ! Procedure(s)

use stats_precision, only: time_precision ! Variable(s)

implicit none

! Parameter Constants
integer, parameter :: &
  ntimes = 41

  ! Constants (taken from the rico case)
  real, parameter :: &
    C_h_20  = 0.001094,  & ! Drag coefficient, defined by RICO 3D specification
    C_q_20  = 0.001133,  & ! Drag coefficient, defined by RICO 3D specification
    z0      = 0.00015      ! Roughness length, defined by ATEX specification

    ! Internal variables
  real :: &
    Ch,   & ! This is C_h_20 scaled to the height of the lowest model level.
    Cq      ! This is C_q_20 scaled to the height of the lowest model level.

! Input variables

real(kind=time_precision), intent(in) :: time     ! Current time

real, intent(in) ::  & 
  rho0,        & ! Density at (1)         [kg/m^3]
  ubar,        & ! Mean sfc wind speed
  rtm,         & ! This is rt at the lowest above-ground model level.  [kg/kg]
  thlm,        & ! This is theta-l at the lowest above-ground model level.  
                 ! (DOES THIS NEED A CORRECTION FOR THETA-L TO THETA?)  [K]
  lowestlevel, & ! This is z at the lowest above-ground model level.  [m]
  exner_sfc,   & ! This is the surface pressure [Pa].
  psfc    

! Output variables

real, intent(out) ::  & 
  wpthlp_sfc,   & ! w'th_l' at (1)   [(m K)/s]  
  wprtp_sfc,    & ! w'r_t'(1) at (1) [(m kg)/(s kg)]
  ustar,        & ! surface friction velocity     [m/s]
  T_sfc           ! Sea surface temperature [K].

! Local variables
integer :: &
  i1, i2

real :: &
  true_time, time_frac

!-----------------BEGIN CODE-------------------------


! Compute heat and moisture fluxes

  ! Modification in case lowest model level isn't at 10 m, from ATEX specification
  Ch   = C_h_20 * ((log(20/z0))/(log(lowestlevel/z0))) * &
         ((log(20/z0))/(log(lowestlevel/z0)))
         ! Modification in case lowest model level isn't at 10 m, from ATEX specification
  Cq   = C_q_20 * ((log(20/z0))/(log(lowestlevel/z0))) * &
         ((log(20/z0))/(log(lowestlevel/z0)))


!sensible_heat_flx = 10.0
!latent_heat_flx = 25.0
true_time = real( time )

T_sfc = 0.0

! We set ustar as it is set in rico
ustar = 0.28

if ( true_time <= time_sfc_given(1) ) then
  T_sfc = T_sfc_given(1)

else if (true_time >= time_sfc_given(ntimes) ) then
  T_sfc = T_sfc_given(ntimes)

else ! true_time > time_sfc_given(1) and true_time < time_sfc_given(ntimes)
  i1 = 1

  do while ( i1 <= ntimes-1 )
    i2 = i1 + 1

    if ( true_time >= time_sfc_given(i1) .and. true_time < time_sfc_given(i2) ) then
      time_frac = (true_time-time_sfc_given(i1))/(time_sfc_given(i2)-time_sfc_given(i1))

      T_sfc = factor_interp( time_frac, T_sfc_given(i2), T_sfc_given(i1) )
      i1 = ntimes
    end if ! true_time >= time_sfc_given(i1) & true_time < time_sfc_given(i2)

    i1 = i2
  end do ! while i1 <= ntimes-1

endif ! else

wpthlp_sfc = compute_wpthlp_sfc( Ch, ubar, thlm, T_sfc, exner_sfc )
wprtp_sfc  = compute_wprtp_sfc( Cq, ubar, rtm, sat_mixrat_liq(psfc,T_sfc) )

!wpthlp_sfc = sensible_heat_flx / ( rho0 * Cp )
!wprtp_sfc  = latent_heat_flx / ( rho0 * Lv )

! Compute momentum fluxes

return
end subroutine astex_a209_sfclyr

end module astex_a209
