!----------------------------------------------------------------------
! $Id$
module lba

!       Description:
!       Contains subroutines for the LBA case.
!----------------------------------------------------------------------

implicit none

private ! Default Scope

public :: lba_tndcy, lba_sfclyr, lba_init

private :: zrad, krad, ntimes, nzrad

! Constant Parameters
integer, parameter :: ntimes = 36, nzrad = 33

real, dimension(nzrad) :: & 
  zrad

real, dimension(nzrad, ntimes) :: & 
  krad 

integer, parameter :: per_line = 5

! File path for forcing the forcing files 
character(*), parameter ::  & 
   file_path = '../model/lba_forcings/'


contains

!----------------------------------------------------------------------
subroutine lba_tndcy( time, & 
                      wm_zt, wm_zm, radht, & 
                      thlm_forcing, rtm_forcing, & 
                      sclrm_forcing )
!       Description:
!       Subroutine to set theta and water tendencies for LBA case.

!       References:
!----------------------------------------------------------------------

use grid_class, only: gr !  Variable(s)

use model_flags, only: l_bugsrad ! Variable(s)

use parameters_model, only: sclr_dim ! Variable(s)

use interpolation, only: zlinterp_fnc ! Procedure(s)

use stats_precision, only: time_precision ! Variable(s)

use array_index, only:  & 
    iisclr_thl, iisclr_rt ! Variable(s)

use interpolation, only: factor_interp ! Procedure(s)

implicit none

! Input
real(kind=time_precision), intent(in) :: time ! Model time [s]

! Output Variables
real, intent(out), dimension(gr%nnzp) :: & 
  wm_zt,        & ! w wind on thermodynamic grid                 [m/s]
  wm_zm,        & ! w wind on momentum grid                      [m/s]
  radht,        & ! Radiative heating rate                       [K/s]
  thlm_forcing, & ! Liquid water potential temperature tendency  [K/s]
  rtm_forcing     ! Total water mixing ratio tendency            [kg/kg/s]

! Output Variables (optional)
real, optional, intent(out), dimension(gr%nnzp,sclr_dim) :: & 
  sclrm_forcing ! Passive scalar forcing [units vary]

! Local Variables
real, dimension(nzrad) :: radhtz
real :: a
integer :: i1, i2

! Large scale subsidence
wm_zt(1:gr%nnzp) = 0.0
wm_zm(1:gr%nnzp) = 0.0

if ( .not. l_bugsrad ) then

  ! Calculate radiative heating rate
  if ( time <=  600. ) then
    radhtz = krad(:,1)

  else if ( time >= ntimes * 600. ) then
    radhtz = krad(:,ntimes)

  else
    i1 = 1
    do while ( i1 <= ntimes-1 )
      i2 = i1 + 1
      if ( time >= 600. * i1 .and. time < 600. * i2  ) then
        a  = real(( time - 600. * i1 )/( 600. * i2 - 600. * i1))
        !radhtz(:) = ( 1. - a ) * krad(:,i1) + a * krad(:,i2)
        radhtz(:) = factor_interp( a, krad(:,i2), krad(:,i1) )
        i1     = ntimes
      end if
      i1 = i2
    end do
  end if ! time <= times(1)

  radht = zlinterp_fnc( gr%nnzp, nzrad, gr%zt, zrad, radhtz )

  ! Radiative theta-l tendency

  thlm_forcing = radht

else ! Compute heating rate interactively with BUGSrad

  thlm_forcing = 0.0

end if ! ~l_bugsrad

! Boundary conditions
thlm_forcing(1) = 0.0  ! Below surface

! Large scale advective moisture tendency
rtm_forcing(:) = 0.0

! Test scalars with thetal and rt if desired
if ( iisclr_thl > 0 ) sclrm_forcing(:,iisclr_thl) = thlm_forcing
if ( iisclr_rt  > 0 ) sclrm_forcing(:,iisclr_rt)  = rtm_forcing

return
end subroutine lba_tndcy

!----------------------------------------------------------------------
subroutine lba_sfclyr( time, z, rho0, & 
                       thlm_sfc, um_sfc, vm_sfc,  & 
                       upwp_sfc, vpwp_sfc, & 
                       wpthlp_sfc, wprtp_sfc, ustar, & 
                       wpsclrp_sfc, wpedsclrp_sfc )

!       Description:
!       This subroutine computes surface fluxes of horizontal momentum,
!       heat and moisture according to GCSS BOMEX specifications

!       References:
!       Grabowski, et al. (2005)
!----------------------------------------------------------------------

use constants, only: pi, grav, Lv, Cp ! Variable(s)

use parameters_model, only: sclr_dim ! Variable(s)

use stats_precision, only: time_precision ! Variable(s)

use diag_ustar_mod, only: diag_ustar ! Variable(s)

use array_index, only:  & 
    iisclr_thl, iisclr_rt ! Variable(s)

implicit none

intrinsic :: max, sqrt

! Constant Parameters
real, parameter ::  & 
  ubmin = 0.25,   & ! Minimum value for ubar 
  z0    = 0.035  ! ARM mom. roughness height

! Input Variables
real(kind=time_precision), intent(in) ::  & 
  time      ! Current time        [s] 

real, intent(in) ::  & 
  z,         & ! Height at zt=2      [m] 
  rho0,      & ! Density at zm=1     [kg/m^3] 
  thlm_sfc,  & ! thlm at (2)         [m/s]
  um_sfc,    & ! um at (2)           [m/s]
  vm_sfc       ! vm at (2)           [m/s]

! Output variables
real, intent(out) ::  & 
  upwp_sfc,     & ! u'w' at (1)      [m^2/s^2]
  vpwp_sfc,     & ! v'w'at (1)       [m^2/s^2]
  wpthlp_sfc,   & ! w'th_l' at (1)   [(m K)/s]  
  wprtp_sfc,    & ! w'r_t'(1) at (1) [(m kg)/(s kg)]
  ustar           ! surface friction velocity [m/s]

! Output variables (optional)
real, intent(out), optional, dimension(sclr_dim) ::  & 
  wpsclrp_sfc,   & ! Passive scalar surface flux      [units m/s] 
  wpedsclrp_sfc    ! Passive eddy-scalar surface flux [units m/s]

! Local variables
!        real :: ft, ubar, ustar, bflx
real :: ft, ubar, bflx

! Compute heat and moisture fluxes
! From Table A.1.
ft = real( max( 0._time_precision,  & 
             cos( 0.5 * pi * ( (5.25 - ( time/3600.)) / 5.25 ) ) & 
        ) )

wpthlp_sfc =  ( 270. * ft**1.5 ) / ( rho0 * Cp )
wprtp_sfc  =  ( 554. * ft**1.3 ) / ( rho0 * Lv )


! Let passive scalars be equal to rt and theta_l for now
if ( iisclr_thl > 0 ) wpsclrp_sfc(iisclr_thl) = wpthlp_sfc
if ( iisclr_rt  > 0 ) wpsclrp_sfc(iisclr_rt)  = wprtp_sfc

if ( iisclr_thl > 0 ) wpedsclrp_sfc(iisclr_thl) = wpthlp_sfc
if ( iisclr_rt  > 0 ) wpedsclrp_sfc(iisclr_rt)  = wprtp_sfc

! Compute momentum fluxes using ARM formulae

ubar = max( ubmin, sqrt( um_sfc**2 + vm_sfc**2 ) )

bflx = grav/thlm_sfc * wpthlp_sfc

! Compute ustar
ustar = diag_ustar( z, bflx, ubar, z0 )

upwp_sfc = -um_sfc * ustar**2 / ubar
vpwp_sfc = -vm_sfc * ustar**2 / ubar

return
end subroutine lba_sfclyr

!----------------------------------------------------------------
subroutine lba_init()
!
!       Description:
!       This subroutine initializes the module by reading in forcing
!       data used in the tndcy subroutine.
!----------------------------------------------------------------

  use file_functions, only: file_read_1d, file_read_2d ! Procedure(s)

  implicit none

  call file_read_1d(10, & 
    file_path//'lba_heights.dat', & 
    nzrad, per_line, zrad)

  call file_read_2d(10, & 
    file_path//'lba_rad.dat', & 
    nzrad, ntimes, per_line, krad)
   
end subroutine lba_init


end module lba

