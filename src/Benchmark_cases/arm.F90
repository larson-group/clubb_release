!----------------------------------------------------------------------
! $Id$
module arm

!       Description:
!       Contains subroutines for the GCSS ARM case.
!----------------------------------------------------------------------

implicit none

public :: arm_tndcy, arm_sfclyr

private ! Default Scope

contains

!----------------------------------------------------------------------
subroutine arm_tndcy( time, thlm_forcing, radht, rtm_forcing, &
                      sclrm_forcing, edsclrm_forcing )
!       Description:
!       Subroutine to set theta and water tendencies for ARM case

!       References:
!       None
!----------------------------------------------------------------------

use grid_class, only: gr ! Variable(s)

use parameters_model, only: sclr_dim, edsclr_dim ! Variable(s)

use parameters_radiation, only: rad_scheme ! Variable(s)

use stats_precision, only: time_precision ! Variable(s)

use array_index, only:  & 
    iisclr_thl, iisclr_rt, iiedsclr_thl, iiedsclr_rt ! Variable(s)

use interpolation, only: factor_interp ! Procedure(s) 
 
use stats_type, only: stat_update_var ! Procedure(s)

use stats_variables, only: iradht_LW, zt, l_stats_samp ! Variable(s)

implicit none

! External
intrinsic :: int, min, max

! Constant Parameters
real, parameter, dimension(6) ::  & 
  atheta = (/ 0.000, 0.000,  0.000, -0.080, -0.160, -0.160/), & 
  rtheta = (/-0.125, 0.000,  0.000,  0.000,  0.000, -0.100/), & 
  art    = (/ 0.080, 0.020, -0.040, -0.100, -0.160, -0.300/)

! Input Variables
real(kind=time_precision), intent(in) :: time ! Model time [s]

! Output Variables
real, intent(out), dimension(gr%nnzp) ::  & 
  thlm_forcing,  & ! Liquid water potential temperature tendency [K/s]
  radht,         & ! Radiative heating rate                      [K/s]
  rtm_forcing      ! Total water mixing ratio tendency           [kg/kg/s]

real, intent(out), dimension(gr%nnzp, sclr_dim) :: & 
  sclrm_forcing   ! Passive scalar tendency         [units/s]

real, intent(out), dimension(gr%nnzp, edsclr_dim) :: & 
  edsclrm_forcing ! Eddy-passive scalar tendency    [units/s]

! Local variables
integer :: k, i1, i2 ! Loop indices
real ::  & 
  a, b,       & ! [-]
  true_time,  & ! [s]
  theta_tmp,  & ! [K/s]
  rad_tmp,    & ! [K/s]
  rt_tmp        ! [kg/kg/s]

!-----------------------------------------------------------------------

true_time = real( time )

! Interpolate in time to get theta and rt tendency

i1 = floor( ( true_time - 41400. ) / 10800. ) + 1
i1 = min( max( i1, 1 ), 5 )
i2 = i1 + 1

if (i1 < 5) then
   a = ( true_time - (41400. + 10800. * (i1-1)) ) / 10800.
else
   a = ( true_time - (41400. + 10800. * (i1-1)) ) / 9000.
end if

if ( trim( rad_scheme ) == "simplified" ) then
!  theta_tmp = ( 1. - a ) * ( atheta(i1) ) & 
!            + a * ( atheta(i2) )
  theta_tmp = factor_interp( a, atheta(i2), atheta(i1) )
!  rad_tmp = ( 1. - a ) * ( rtheta(i1) ) & 
!          + a * ( rtheta(i2) )
  rad_tmp = factor_interp( a, rtheta(i2), rtheta(i1) )
else ! Factor in radiation later
!  theta_tmp = ( 1. - a ) * ( atheta(i1) + 0.0 ) & 
!            + a * ( atheta(i2) + 0.0 )
  theta_tmp = factor_interp( a, atheta(i2) + 0.0, atheta(i1) + 0.0 )
  rad_tmp   = 0.0

end if

!rt_tmp = ( 1. - a ) * art(i1) + a * art(i2)
rt_tmp = factor_interp(a, art(i2), art(i1) )

! Convert to the right units

theta_tmp = theta_tmp / 3600.
rad_tmp   = rad_tmp / 3600.
rt_tmp    = rt_tmp / ( 3600. * 1000. )

! Interpolate with respect to height

do k = 2, gr%nnzp
  select case( int( gr%zt(k) ) )
  case ( 0:999 )
    rtm_forcing(k)  = rt_tmp
    thlm_forcing(k) = theta_tmp + rad_tmp
    radht(k)        = rad_tmp

  case ( 1000:2999 )
    b               = 1. - ( gr%zt(k) - 1000. ) / 2000.
    rtm_forcing(k)  = b * rt_tmp
    thlm_forcing(k) = b * ( theta_tmp + rad_tmp )
    radht(k)        = b * rad_tmp

  case default
    rtm_forcing(k)  = 0.0
    thlm_forcing(k) = 0.0
    radht(k)        = 0.0

  end select
end do ! k=2..gr%nnzp

rtm_forcing(1)  = 0.0
thlm_forcing(1) = 0.0
radht(1)        = 0.0

if ( l_stats_samp .and. trim( rad_scheme ) == "simplified" ) then
  call stat_update_var( iradht_LW, radht, zt )
end if

! Test scalars with thetal and rt if desired
if ( iisclr_thl > 0 ) sclrm_forcing(:,iisclr_thl) = thlm_forcing
if ( iisclr_rt  > 0 ) sclrm_forcing(:,iisclr_rt)  = rtm_forcing

if ( iiedsclr_thl > 0 ) edsclrm_forcing(:,iiedsclr_thl) = thlm_forcing
if ( iiedsclr_rt  > 0 ) edsclrm_forcing(:,iiedsclr_rt)  = rtm_forcing

return
end subroutine arm_tndcy
!----------------------------------------------------------------------
subroutine arm_sfclyr( time, z, dn0, thlm_sfc, um_sfc, vm_sfc,  & 
                       upwp_sfc, vpwp_sfc,  & 
                       wpthlp_sfc, wprtp_sfc, ustar, & 
                       wpsclrp_sfc, wpedsclrp_sfc )

!       Description:
!       This subroutine computes surface fluxes of horizontal momentum,
!       heat and moisture according to GCSS ARM specifications
!----------------------------------------------------------------------

use constants, only: Cp, Lv, grav ! Variable(s)

use parameters_model, only: sclr_dim, edsclr_dim ! Variable(s)

use stats_precision, only: time_precision ! Variable(s)

use diag_ustar_mod, only: diag_ustar ! Variable(s)

use array_index, only: iisclr_rt, iisclr_thl, iiedsclr_rt, iiedsclr_thl ! Variable(s)

use surface_flux, only: compute_momentum_flux, compute_ubar

implicit none

intrinsic :: max, sqrt

! ARM roughness height
real, parameter ::  & 
  z0 = 0.035  ! momentum roughness height

! Input variables
real(kind=time_precision), intent(in) ::  & 
  time            ! Current time          [s]

real, intent(in) ::  & 
  z,               & ! Height at zt(2)       [m]
  dn0,             & ! Density at zm(1)      [kg/m^3]
  thlm_sfc,        & ! Theta_l at zt(2)      [K]
  um_sfc,          & ! um at zt(2)           [m/s]
  vm_sfc             ! vm at zt(2)           [m/s]

! Output variables
real, intent(out) ::  & 
  upwp_sfc,    & ! u'w' at surface           [m^2/s^2]
  vpwp_sfc,    & ! v'w' at surface           [m^2/s^2]
  wpthlp_sfc,  & ! w'theta_l' surface flux   [(m K)/s]
  wprtp_sfc,   & ! w'rt' surface flux        [(m kg)/(kg s)]
  ustar          ! surface friction velocity [m/s]

real,  dimension(sclr_dim), intent(out) ::  & 
  wpsclrp_sfc        ! Passive scalar surface flux      [units m/s] 

real,  dimension(edsclr_dim), intent(out) ::  & 
  wpedsclrp_sfc      ! Passive eddy-scalar surface flux [units m/s]

! Internal variables
real ::  & 
  ubar, & 
  true_time, & 
  heat_flx, moisture_flx, & 
  heat_flx2, moisture_flx2, & 
  bflx

! Compute heat and moisture fluxes from ARM data in (W/m2)
true_time = real( time )

call arm_sfcflx( true_time, heat_flx, moisture_flx )

! Compute momentum fluxes

! Convert heat_flx and moisture_flx to natural units
heat_flx2     = heat_flx / ( Cp * dn0 )    ! (K m/s)
moisture_flx2 = moisture_flx / ( Lv * dn0 )! (m/s)

! Heat flux in units of (m2/s3) (needed by diag_ustar)
bflx = grav/thlm_sfc * heat_flx2

! Surface winds

! Compute ubar
ubar = compute_ubar( um_sfc, vm_sfc )
! Compute ustar
ustar = diag_ustar( z, bflx, ubar, z0 )

! Assign fluxes
call compute_momentum_flux( um_sfc, vm_sfc, ubar, ustar, &
                              upwp_sfc, vpwp_sfc )

wpthlp_sfc = heat_flx2
wprtp_sfc  = moisture_flx2

! Let passive scalars be equal to rt and theta_l for now
if ( iisclr_thl > 0 ) wpsclrp_sfc(iisclr_thl) = wpthlp_sfc
if ( iisclr_rt  > 0 ) wpsclrp_sfc(iisclr_rt)  = wprtp_sfc

if ( iiedsclr_thl > 0 ) wpedsclrp_sfc(iiedsclr_thl) = wpthlp_sfc
if ( iiedsclr_rt  > 0 ) wpedsclrp_sfc(iiedsclr_rt)  = wprtp_sfc

return
end subroutine arm_sfclyr

!------------------------------------------------------------------------
subroutine arm_sfcflx( time, heat_flx, moisture_flx )

!       Description:
!       This subroutine computes surface heat and moisture for a specific time
!       according to GCSS ARM specifications. Flux returned are in (W/m2)
!------------------------------------------------------------------------
use interpolation, only: factor_interp

implicit none

! Parameter constants
integer, parameter :: ntimes = 7

real, parameter, dimension(ntimes) ::  & 
  times = (/ 41400., 55800., 64800., 68400., & 
           77400., 86400., 93600. /), & 
  ! H and LE specifications
  H  = (/-30,  90, 140, 140, 100, -10, -10/), & 
  LE = (/  5, 250, 450, 500, 420, 180,   0/)

! Input variable
real, intent(in) :: time !  Current time [s]

! Output variables
real, intent(out) :: heat_flx, moisture_flx

! Local variables
integer :: i1, i2
real :: a

if ( time <= times(1) ) then
   heat_flx     = H(1)
   moisture_flx = LE(1)
else if ( time >= times(ntimes) ) then
   heat_flx     = H(ntimes)
   moisture_flx = LE(ntimes)
else
   i1 = 1
   do while ( i1 <= ntimes-1 )
      i2 = i1 + 1
      if ( time >= times(i1) .and. time < times(i2) ) then
         a            = (time-times(i1))/(times(i2)-times(i1))
!         heat_flx     = ( 1. - a ) * H(i1) + a * H(i2)
         heat_flx = factor_interp( a, H(i2), H(i1) )
!         moisture_flx = ( 1. - a ) * LE(i1) + a * LE(i2)
         moisture_flx = factor_interp( a, LE(i2), LE(i1) )
         i1           = ntimes
      end if
      i1 = i2
   end do
end if ! time <= times(1)

return
end subroutine arm_sfcflx

end module arm
