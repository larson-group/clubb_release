!----------------------------------------------------------------------
! $Id: astex.F90,v 1.6 2008-08-06 13:53:02 faschinj Exp $
module astex

!       Description:
!       Contains subroutines for the ASTEX KK case.
!----------------------------------------------------------------------

implicit none

public :: astex_tndcy, astex_sfclyr

private ! Default Scope

contains

!----------------------------------------------------------------------
subroutine astex_tndcy( wm_zt, wm_zm,  & 
                        thlm_forcing, rtm_forcing, sclrm_forcing )

!       Description:
!       Subroutine to set theta and water tendencies for ASTEX KK case
!       References:
!----------------------------------------------------------------------

use parameters, only: sclr_dim ! Variable(s)

use grid_class, only: gr ! Variable(s)

use grid_class, only: zt2zm ! Procedure(s)

use stats_precision, only: time_precision ! Variable(s)

implicit none

! Output Variables
real, intent(out), dimension(gr%nnzp) ::  & 
  wm_zt,           & ! w wind on the thermodynamic grid        [m/s]
  wm_zm,           & ! w wind on the momentum grid             [m/s]
  thlm_forcing,  & ! Liquid potential temperature tendency   [K/s]
  rtm_forcing   ! Total water mixing ratio tendency       [kg/kg/s]

real, intent(out), dimension(gr%nnzp,sclr_dim) ::  & 
  sclrm_forcing ! Passive scalar forcing  [units/s]

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

! Passive scalar testing

if ( sclr_dim > 0 ) then
  sclrm_forcing(:,:) = 0.0
end if

return
end subroutine astex_tndcy

!----------------------------------------------------------------------
subroutine astex_sfclyr( rho0, & 
                         upwp_sfc, vpwp_sfc, wpthlp_sfc, & 
                         wprtp_sfc, wpsclrp_sfc, wpedsclrp_sfc )

!       Description:
!       This subroutine computes surface fluxes of horizontal momentum,
!       heat and moisture according to ASTEX with Khairoutdinov and Kogan
!       alteration.

!       References:
!----------------------------------------------------------------------

use constants, only: Cp, Lv ! Variable(s)

use parameters, only: sclr_dim ! Variable(s)

use array_index, only: iisclr_rt, iisclr_thl ! Variable(s)

implicit none

! Input variables

real, intent(in) ::  & 
  rho0        ! Density at (1)         [kg/m^3]
! Output variables

real, intent(out) ::  & 
  upwp_sfc,     & ! u'w' at (1)      [m^2/s^2]
  vpwp_sfc,     & ! v'w'at (1)       [m^2/s^2]
  wpthlp_sfc,   & ! w'th_l' at (1)   [(m K)/s]  
  wprtp_sfc    ! w'r_t'(1) at (1) [(m kg)/(s kg)]

real, intent(out), dimension(sclr_dim) ::  & 
  wpsclrp_sfc,   & ! w' scalar at surface [units m/s]
  wpedsclrp_sfc ! w' scalar at surface [units m/s]

! Local variables

real :: sensible_heat_flx,  & ! W/m^2
        latent_heat_flx    ! W/m^2

! Compute heat and moisture fluxes

sensible_heat_flx = 10.0
latent_heat_flx = 25.0

wpthlp_sfc = sensible_heat_flx/( rho0*Cp )
wprtp_sfc  = latent_heat_flx/( rho0*Lv )

! Compute momentum fluxes

upwp_sfc = 0.09
vpwp_sfc = 0.09

! Test scalars
if ( iisclr_rt > 0 ) then
  wpsclrp_sfc(iisclr_rt)    = wprtp_sfc
  wpedsclrp_sfc(iisclr_thl) = wprtp_sfc
end if
if ( iisclr_thl > 0 ) then
  wpsclrp_sfc(iisclr_thl)   = wpthlp_sfc
  wpedsclrp_sfc(iisclr_thl) = wpthlp_sfc
end if

return
end subroutine astex_sfclyr

end module astex
