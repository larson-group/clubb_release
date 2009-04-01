!----------------------------------------------------------------------
! $Id$
module wangara

implicit none

public :: wangara_tndcy, wangara_sfclyr

private ! Default Scope

contains
!----------------------------------------------------------------------
subroutine wangara_tndcy( wm_zt, wm_zm,  & 
                          thlm_forcing, rtm_forcing, & 
                          sclrm_forcing, edsclrm_forcing )
!       Description:
!       Subroutine to set theta and water tendencies for Wangara case
!       References;
!       None
!----------------------------------------------------------------------

use grid_class, only: gr ! Variable(s)

use parameters_model, only: sclr_dim, edsclr_dim ! Variable(s)

use stats_precision, only: time_precision ! Variable(s)

use array_index, only: iisclr_thl, iisclr_rt, iiedsclr_thl, iiedsclr_rt ! Variable(s)

implicit none

! Output Variables
real, intent(out), dimension(gr%nnzp) :: & 
  wm_zt,        & ! w wind on thermodynamic grid                [m/s]
  wm_zm,        & ! w wind on momentum grid                     [m/s]
  thlm_forcing, & ! Liquid water potential temperature tendency [K/s]
  rtm_forcing     ! Total water mixing ratio tendency           [kg/kg/s]

real, intent(out), dimension(gr%nnzp,sclr_dim) :: & 
  sclrm_forcing ! Passive scalar tendency [units/s]

real, intent(out), dimension(gr%nnzp,edsclr_dim) :: & 
  edsclrm_forcing ! Eddy-passive scalar tendency [units/s]

! No large-scale subsidence for now
wm_zt = 0.0
wm_zm = 0.0

! No large-scale water tendency or cooling

rtm_forcing  = 0.0
thlm_forcing = 0.0

! Test scalars with thetal and rt if desired
if ( iisclr_thl > 0 ) sclrm_forcing(:,iisclr_thl) = thlm_forcing
if ( iisclr_rt  > 0 ) sclrm_forcing(:,iisclr_rt)  = rtm_forcing

if ( iiedsclr_thl > 0 ) edsclrm_forcing(:,iiedsclr_thl) = thlm_forcing
if ( iiedsclr_rt  > 0 ) edsclrm_forcing(:,iiedsclr_rt)  = rtm_forcing

return
end subroutine wangara_tndcy

!----------------------------------------------------------------------
subroutine wangara_sfclyr( time, um_sfc, vm_sfc,  & 
                           upwp_sfc, vpwp_sfc,  & 
                           wpthlp_sfc, wprtp_sfc, ustar, & 
                           wpsclrp_sfc, wpedsclrp_sfc )
!       Description:
!       This subroutine computes surface fluxes of horizontal momentum,
!       heat and moisture for Wangara day 33
!----------------------------------------------------------------------

use constants, only: pi, fstderr, sec_per_day ! Variable(s)

use parameters_model, only: sclr_dim ! Variable(s)

use stats_precision, only: time_precision ! Variable(s)

use array_index, only:  & 
    iisclr_thl, iisclr_rt ! Variable(s)

use surface_flux, only: compute_ubar, compute_momentum_flux

implicit none

intrinsic :: mod, max, cos, sqrt, present

! Input variables
real(kind=time_precision), intent(in) ::  & 
  time    ! Current time  [s]

real, intent(in) ::  & 
  um_sfc,  & ! um(2)         [m/s]
  vm_sfc     ! vm(2)         [m/s]

! Output variables
real, intent(out) ::  & 
  upwp_sfc,     & ! u'w' at (1)      [m^2/s^2]
  vpwp_sfc,     & ! v'w'at (1)       [m^2/s^2]
  wpthlp_sfc,   & ! w'th_l' at (1)   [(m K)/s]  
  wprtp_sfc,    & ! w'r_t'(1) at (1) [(m kg)/(s kg)]
  ustar           ! surface friction velocity [m/s]

! Output variables
real, intent(out), dimension(sclr_dim) ::  & 
  wpsclrp_sfc,   & ! Passive scalar surface flux      [units m/s]
  wpedsclrp_sfc    ! Passive eddy-scalar surface flux [units m/s]

! Local variables

real :: ubar 
real(kind=time_precision) :: time_utc, time_est


! Declare the value of ustar.
ustar = 0.13

! Compute UTC time of the day in seconds

time_utc = mod( time, sec_per_day )

! Now convert UTC time to Australia EST (local time)

time_est = mod( time_utc + 36000., sec_per_day )

if ( time_est < 27000 .or. time_est > 63000 ) then
   write(fstderr,*) "wangara_sfclyr: error local time must" & 
     //" be between 730 and 1730."
   write(fstderr,*) 'time_est = ',time_est
   stop
end if

! Compute heat and moisture fluxes

wpthlp_sfc = real(0.18 * cos( (time_est-45000.0)/36000.0 * pi ))
wprtp_sfc  = 1.3e-4 * wpthlp_sfc

! Let the passive scalars be equal to rt and theta_l for now

! Let passive scalars be equal to rt and theta_l for now
if ( iisclr_thl > 0 ) wpsclrp_sfc(iisclr_thl) = wpthlp_sfc
if ( iisclr_rt  > 0 ) wpsclrp_sfc(iisclr_rt)  = wprtp_sfc

if ( iisclr_thl > 0 ) wpedsclrp_sfc(iisclr_thl) = wpthlp_sfc
if ( iisclr_rt  > 0 ) wpedsclrp_sfc(iisclr_rt)  = wprtp_sfc

! Compute momentum fluxes

ubar = compute_ubar( um_sfc, vm_sfc )

call compute_momentum_flux( um_sfc, vm_sfc, ubar, ustar, &
                            upwp_sfc, vpwp_sfc )

return
end subroutine wangara_sfclyr

!----------------------------------------------------------------------

end module wangara
