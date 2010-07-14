!----------------------------------------------------------------------
! $Id$
module cobra
!       Description:
!       Contains subroutines for the COBRA CO2 case.
!----------------------------------------------------------------------

implicit none

private ! Default Scope

public :: & 
  cobra_tndcy, & 
  cobra_sfclyr

contains

!----------------------------------------------------------------------
subroutine cobra_tndcy( thlm_forcing, rtm_forcing, & 
                        sclrm_forcing, edsclrm_forcing )
!       Description:
!       Subroutine to set theta and water tendencies for COBRA CO2 case

!       References:
!       None
!----------------------------------------------------------------------

use grid_class, only: gr ! Variable(s)

use grid_class, only: zt2zm ! Procedure(s)

use constants_clubb, only: fstderr ! Variable(s)

use parameters_model, only: sclr_dim, edsclr_dim ! Variable(s)

use stats_precision, only: time_precision ! Variable(s)

use array_index, only:  & 
    iisclr_thl, iisclr_rt, iisclr_CO2, & ! Variable(s)
    iiedsclr_thl, iiedsclr_rt, iiedsclr_CO2

implicit none

! Output Variables
real, intent(out), dimension(gr%nnzp) :: & 
  thlm_forcing, & ! Liquid water potential temperature tendency  [K/s]
  rtm_forcing     ! Total water mixing ratio tendency            [kg/kg/s]

real, intent(out), dimension(gr%nnzp,sclr_dim) :: & 
  sclrm_forcing ! Passive scalar tendency        [units vary/s]

real, intent(out), dimension(gr%nnzp,edsclr_dim) :: & 
  edsclrm_forcing ! Passive eddy-scalar tendency [units/s]

! Local Variables
integer :: k

DO k = 2, gr%nnzp, 1
   if ( gr%zt(k) < 0. ) then
      write(fstderr,*) "cobra_tndcy:" & 
         //" error in subsidence profile."
      write(fstderr,*) 'Altitude gr%zt = ', gr%zt(k)
      stop
   end if
END DO

! No large-scale water tendency or cooling
rtm_forcing  = 0.0
thlm_forcing = 0.0

! Setup passive scalars, if they're enabled
if ( iisclr_CO2 > 0 ) sclrm_forcing(:,iisclr_CO2) = 0.0
if ( iisclr_thl > 0 ) sclrm_forcing(:,iisclr_thl) = thlm_forcing
if ( iisclr_rt  > 0 ) sclrm_forcing(:,iisclr_rt) = rtm_forcing
if ( iiedsclr_CO2 > 0 ) edsclrm_forcing(:,iiedsclr_CO2) = 0.0
if ( iiedsclr_thl > 0 ) edsclrm_forcing(:,iiedsclr_thl) = thlm_forcing
if ( iiedsclr_rt  > 0 ) edsclrm_forcing(:,iiedsclr_rt) = rtm_forcing

return
end subroutine cobra_tndcy

!-----------------------------------------------------------------------
subroutine cobra_sfclyr( time, z, dn0, thlm_sfc, um_sfc, vm_sfc, & 
                         upwp_sfc, vpwp_sfc,  & 
                         wpthlp_sfc, wprtp_sfc, ustar, & 
                         wpsclrp_sfc, wpedsclrp_sfc )

!       Description:
!       This subroutine computes surface fluxes of horizontal momentum,
!       heat and moisture according to the format used for the GCSS ARM 
!       case.

!       Notes:
!       The data has been altered so it can be used for the COBRA CO2 
!       case.
!-----------------------------------------------------------------------

use constants_clubb, only: Cp, Lv, grav ! Variable(s)

use parameters_model, only: sclr_dim, edsclr_dim ! Variable(s)

use interpolation, only: factor_interp

use stats_precision, only: time_precision ! Variable(s)

use diag_ustar_module, only: diag_ustar ! Variable(s)

use array_index, only: &
  iisclr_rt, iisclr_thl, iisclr_CO2, & ! Variable(s)
  iiedsclr_rt, iiedsclr_thl, iiedsclr_CO2
  
use surface_flux, only: compute_ubar, compute_momentum_flux

use time_dependent_input, only: LH_given, SH_given, time_sfc_given, CO2_sfc_given ! Variable(s)

implicit none

intrinsic :: sqrt, max

! Parameter Constants
real, parameter :: & 
  M_da  = 0.02897  ! Molecular weight of dry air.
integer, parameter :: &
  ntimes = 49

! Input variables
real(kind=time_precision), intent(in) ::  & 
  time      ! Current time                [s]

real, intent(in) :: & 
  z,         & ! Elevation at zt=2           [m]
  dn0,       & ! Air density at surface      [kg/m^3]
  thlm_sfc,  & ! Theta_l at zt(2)            [K]
  um_sfc,    & ! u wind at zt(2)             [m/s]
  vm_sfc       ! v wind at zt(2)             [m/s]

! Output variables
real, intent(out) ::  & 
  upwp_sfc,    & ! u'w' at surface           [m^2/s^2]
  vpwp_sfc,    & ! v'w' at surface           [m^2/s^2]
  wpthlp_sfc,  & ! w'theta_l' surface flux   [(m K)/s]
  wprtp_sfc,   & ! w'rt' surface flux        [(m kg)/(kg s)]
  ustar          ! surface friction velocity [m/s]

! Output variables
real, intent(out), dimension(sclr_dim) ::  & 
  wpsclrp_sfc    ! w'sclr' surface flux          [units m/s]

real, intent(out), dimension(edsclr_dim) ::  & 
  wpedsclrp_sfc  ! w' edsclr' surface flux       [units m/s]

! Local variables
integer :: &
  i1, i2
real ::  & 
  ubar, & 
  true_time, & 
  heat_flx, moisture_flx, &                 ! [W/m^2]
  heat_flx2, moisture_flx2, &
  time_frac, &
  bflx
real :: CO2_flx
real :: CO2_flx2

! COBRA roughness height
! real, parameter :: z0 = 0.035  ! ARM momentum roughness height
real, parameter :: z0 = 1.75   ! momentum roughness height

!-----------------BEGIN CODE-------------------------

! Default Initialization
heat_flx = 0.0
moisture_flx = 0.0
CO2_flx = 0.0
CO2_flx2 = 0.0 ! Default initialization

! Compute heat and moisture fluxes from ARM data in (W/m2)
true_time = real( time )

if ( true_time <= time_sfc_given(1) ) then
   heat_flx     = SH_given(1)
   moisture_flx = LH_given(1)
   CO2_flx      = CO2_sfc_given(1)
   
else if ( true_time >= time_sfc_given(ntimes) ) then
   heat_flx     = SH_given(ntimes)
   moisture_flx = LH_given(ntimes)
   CO2_flx      = CO2_sfc_given(ntimes)
   
else  ! true_time > time_sfc_given(1) and true_time < time_sfc_given(ntimes)
  i1 = 1

  do while ( i1 <= ntimes-1 )
    i2 = i1 + 1

    if ( true_time >= time_sfc_given(i1) .and. true_time < time_sfc_given(i2) ) then
      time_frac = (true_time-time_sfc_given(i1))/(time_sfc_given(i2)-time_sfc_given(i1))

      !heat_flx     = ( 1. - time_frac ) * SH_given(i1) + time_frac * SH_given(i2)
      !moisture_flx = ( 1. - time_frac ) * LH_given(i1) + time_frac * LH_given(i2)
      heat_flx = factor_interp( time_frac, SH_given(i2), SH_given(i1) )
      moisture_flx = factor_interp( time_frac, LH_given(i2), LH_given(i1) )
      !CO2_flx = ( 1. - time_frac ) * CO2_sfc_given(i1) + time_frac * CO2_sfc_given(i2)
      CO2_flx = factor_interp( time_frac, CO2_sfc_given(i2), CO2_sfc_given(i1) )
      i1           = ntimes
    end if ! true_time >= time_sfc_given(i1) & true_time < time_sfc_given(i2)

    i1 = i2
  end do ! while i1 <= ntimes-1

end if ! else 

! Compute momentum fluxes

! Convert heat_flx and moisture_flx to natural units
heat_flx2     = heat_flx / ( Cp * dn0 )    ! (K m/s)
moisture_flx2 = moisture_flx / ( Lv * dn0 )! (m/s)

!       Convert CO2 surface flux to natural units.
!       The CO2 flux has been given in units of:  umol/(m^2 s).
!       umol stands for micromoles.  The CO2 concentration in
!       this code is in units of ppmv, which is also the molar
!       mixing ratio times 10^6.
!       The units are:  10^6 * [ mol (CO2) / mol (dry air) ].
!       w'CO2' = (Flux) * [ M (dry air) / rho (dry air) ];
!       where M is the molecular weight of dry air.
CO2_flx2 = CO2_flx * ( M_da / dn0 )

! Heat flux in units of (m2/s3) (needed by diag_ustar)
bflx = grav/thlm_sfc * heat_flx2

! Surface winds
ubar = compute_ubar( um_sfc, vm_sfc )

! Compute ustar
ustar = diag_ustar( z, bflx, ubar, z0 )

! Assign fluxes

call compute_momentum_flux( um_sfc, vm_sfc, ubar, ustar, &
                            upwp_sfc, vpwp_sfc )

wpthlp_sfc = heat_flx2
wprtp_sfc  = moisture_flx2

if ( iisclr_CO2 > 0 ) wpsclrp_sfc(iisclr_CO2) = CO2_flx2
if ( iisclr_thl > 0 ) wpsclrp_sfc(iisclr_thl) = wpthlp_sfc
if ( iisclr_rt  > 0 ) wpsclrp_sfc(iisclr_rt)  = wprtp_sfc

if ( iiedsclr_CO2 > 0 ) wpedsclrp_sfc(iiedsclr_CO2) = CO2_flx2
if ( iiedsclr_thl > 0 ) wpedsclrp_sfc(iiedsclr_thl) = wpthlp_sfc
if ( iiedsclr_rt  > 0 ) wpedsclrp_sfc(iiedsclr_rt)  = wprtp_sfc

return
end subroutine cobra_sfclyr

end module cobra
