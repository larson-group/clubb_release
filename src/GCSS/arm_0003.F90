!----------------------------------------------------------------------
! $Id$
module arm_0003

!       Description:
!       Contains subroutines for the March 2000 IOP ARM case.
!----------------------------------------------------------------------

implicit none

public :: arm_0003_tndcy, arm_0003_sfclyr, arm_0003_init

private ! Default Scope

! Constant Parameters
integer, parameter :: ntimes = 167, nz = 38, & 
 per_line = 5

! File path for the forcing files
character(*), parameter :: & 
 file_path = '../model/arm_0003_forcings/'

real, dimension(ntimes) :: times       ! Time from day0      [s]
real, dimension(nz) :: z               ! Height              [m]
real, dimension(nz, ntimes) :: thl_ls  ! Potential Temperature
                                       ! Tendency            [K/s]
real, dimension(nz, ntimes) :: rt_ls   ! Water Vapor Advective
                                       ! Tendency            [Kg/Kg/s]
real, dimension(nz, ntimes) :: um_obs  ! Obs. wind u         [m/s]
real, dimension(nz, ntimes) :: vm_obs  ! Obs. wind v         [m/s]
real, dimension(ntimes) :: SE          ! Sensible heat flux  [W/m^2]
real, dimension(ntimes) :: LE          ! Evaporation         [W/m^2] 


contains

!----------------------------------------------------------------------
subroutine arm_0003_tndcy( time, &
                           wm_zt, wm_zm, thlm_forcing, &
                           rtm_forcing, um_hoc_grid, vm_hoc_grid, sclrm_forcing )
!       Description:
!       Subroutine to set thetal and total water tendencies for ARM 0003 case

!       References:
!       None
!----------------------------------------------------------------------

use interpolation, only: factor_interp ! Procedure(s)

use grid_class, only: gr ! Variable(s)

use parameters, only: sclr_dim ! Variable(s)

use interpolation, only: zlinterp_fnc ! Procedure(s)

use stats_precision, only: time_precision ! Variable(s)

use error_code, only: clubb_debug ! Procedure(s)

use array_index, only: iisclr_rt, iisclr_thl

implicit none

! External
intrinsic :: present

! Input Variables
real(kind=time_precision), intent(in) :: time ! Model time [s]

! Output Variables
real, intent(out), dimension(gr%nnzp) ::  & 
  thlm_forcing,  & ! Liquid water potential temperature tendency  [K/s]
  wm_zm,         & ! Vertical velocity on moment. grid            [m/s]
  wm_zt,         & ! Vertical velocity on thermo. grid            [m/s]
  rtm_forcing,   & ! Total water mixing ratio tendency            [kg/kg/s]
  um_hoc_grid,   & ! Observed wind, for nudging                   [m/s]
  vm_hoc_grid      ! Observed wind, for nudging                   [m/s]

! Output (optional) Variables
real, optional, intent(out), dimension(gr%nnzp,sclr_dim) ::  & 
  sclrm_forcing ! Passive scalar tendency [units vary]

real :: time_frac
integer :: i1, i2

real, dimension(nz) :: thlm_t_interp, & 
  rtm_t_interp, um_obs_t_interp, vm_obs_t_interp

!-----------------------------------------------------------------------

! Thetal forcing is equal to the LS tendency given here and the
! interactive calculation done in BUGSrad

time_frac = -1.0 ! Default initialization

if ( time == times(1) ) then
  time_frac = 0.0
  i1 = 1
  i2 = 2
else if ( time >= times(ntimes) ) then
  time_frac = 1.0
  i1 = ntimes-1
  i2 = ntimes
else
  i1 = 1
  do while ( i1 <= ntimes-1 )
    i2 = i1 + 1
    if ( time >= times(i1) .and. time < times(i2) ) then
      time_frac = real((time-times(i1))/(times(i2)-times(i1)))
      exit
    end if
      i1 = i2
  end do
end if ! time <= times(1)

if (time_frac == -1.0) then
        call clubb_debug & 
           (1,"times not sorted in arm_0003_tndcy")
endif

! Interpolate LS thetal tendency to the HOC grid
! Time
!thlm_t_interp = (1.-time_frac) * thl_ls(:,i1) + time_frac * & 
!   thl_ls(:,i2)
thlm_t_interp = factor_interp( time_frac, thl_ls(:,i2), thl_ls(:,i1) )
! Vertical
thlm_forcing(1:gr%nnzp) = zlinterp_fnc( gr%nnzp, nz, gr%zt(:), z, thlm_t_interp )

! Interpolate LS rt tendency to the HOC grid
! Time
!rtm_t_interp = (1.-time_frac) * rt_ls(:,i1) + time_frac * & 
!   rt_ls(:,i2)
rtm_t_interp = factor_interp( time_frac, rt_ls(:,i2), rt_ls(:,i1) )
! Vertical
rtm_forcing(1:gr%nnzp) = zlinterp_fnc & 
         ( gr%nnzp, nz, gr%zt(:), z, rtm_t_interp )
! Modified by Joshua Fasching (October 2007)

! Interpolate um observed to the HOC grid
! Time
!um_obs_t_interp = (1.-time_frac) * um_obs(:,i1) + time_frac *  & 
!   um_obs(:,i2)
um_obs_t_interp = factor_interp( time_frac, um_obs(:,i2), um_obs(:,i1) )
! Vertical
um_hoc_grid(1:gr%nnzp) = zlinterp_fnc & 
         ( gr%nnzp, nz, gr%zt(:), z, um_obs_t_interp )

! Interpolate vm observed to the HOC grid
! Time
!vm_obs_t_interp = (1.-time_frac) * vm_obs(:,i1) + time_frac *  & 
!   vm_obs(:,i2)
vm_obs_t_interp = factor_interp( time_frac, vm_obs(:,i2), vm_obs(:,i1) )
! Vertical
vm_hoc_grid(1:gr%nnzp) = zlinterp_fnc & 
         ( gr%nnzp, nz, gr%zt(:), z, vm_obs_t_interp )

! Added by Joshua Fasching (October 27 2007)
um_hoc_grid (1) = um_hoc_grid(2)
vm_hoc_grid (1) = vm_hoc_grid(2)

! The specification calls for zero imposed subsidence (always)
wm_zm(:) = 0.
wm_zt(:) = 0.

! Test scalars with thetal and rt if desired
if ( iisclr_thl > 0 ) sclrm_forcing(:,iisclr_thl) = thlm_forcing
if ( iisclr_rt  > 0 ) sclrm_forcing(:,iisclr_rt)  = rtm_forcing

return
end subroutine arm_0003_tndcy
!----------------------------------------------------------------------
subroutine arm_0003_sfclyr( time, z, rho0, & 
                            thlm_sfc, um_sfc, vm_sfc,  & 
                            upwp_sfc, vpwp_sfc, & 
                            wpthlp_sfc, wprtp_sfc, ustar, & 
                            wpsclrp_sfc, wpedsclrp_sfc )
!       Description:
!       This subroutine computes surface fluxes of horizontal momentum,
!       heat and moisture according to GCSS ARM specifications
!----------------------------------------------------------------------

use constants, only: Cp, Lv, grav ! Variable(s)

use parameters, only: sclr_dim ! Variable(s)

use stats_precision, only: time_precision ! Variable(s)

use diag_ustar_mod, only: diag_ustar ! Variable(s)

use array_index, only: iisclr_rt, iisclr_thl ! Variable(s)

use interpolation, only: factor_interp ! Procedure(s)

implicit none

intrinsic :: max, sqrt

real, parameter ::  & 
  ubmin = 0.25,    & ! Minimum value for ubar 
  z0    = 0.035   ! ARM Cu mom. roughness height

! Input Variables
real(kind=time_precision), intent(in) ::  & 
  time      ! Current time        [s] 

real, intent(in) ::  & 
  z,         & ! Height at zt=2      [s] 
  rho0,      & ! Density at zm=1     [kg/m^3] 
  um_sfc,    & ! um at (2)           [m/s]
  vm_sfc,    & ! vm at (2)           [m/s]
  thlm_sfc     ! thlm at (2)         [m/s]

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
!        real :: ubar, ustar, bflx, heat_flx, moisture_flx, time_frac
real :: ubar, bflx, heat_flx, moisture_flx, time_frac
integer :: i1, i2
!----------------------------------------------------------------------

! Default Initialization
heat_flx = 0.0
moisture_flx = 0.0

if ( time <= times(1) ) then
  heat_flx     = SE(1)
  moisture_flx = LE(1)
else if ( time >= times(ntimes) ) then
  heat_flx     = SE(ntimes)
  moisture_flx = LE(ntimes)
else
  i1 = 1
  do while ( i1 <= ntimes-1 )
    i2 = i1 + 1
    if ( time >= times(i1) .and. time < times(i2) ) then
      time_frac       = real((time-times(i1))/(times(i2) & 
         - times(i1)))
!      heat_flx     = ( 1. - time_frac ) * SE(i1) +  & 
!         time_frac * SE(i2)
      heat_flx = factor_interp( time_frac, SE(i2), SE(i1) )
!      moisture_flx = ( 1. - time_frac ) * LE(i1) +  & 
!         time_frac * LE(i2)
      moisture_flx = factor_interp( time_frac, SE(i2), SE(i1) )
      i1           = ntimes
    end if
      i1 = i2
   end do
 end if ! time <= times(1)

! Convert W/m^2 into w'thl' w'rt' units
wpthlp_sfc = heat_flx / ( Cp * rho0 )     ! (K m/s)
wprtp_sfc  = moisture_flx / ( Lv * rho0 ) ! (kg m/ kg s)

! Let passive scalars be equal to rt and theta_l for now
if ( iisclr_thl > 0 ) wpsclrp_sfc(iisclr_thl) = wpthlp_sfc
if ( iisclr_rt  > 0 ) wpsclrp_sfc(iisclr_rt)  = wprtp_sfc

if ( iisclr_thl > 0 ) wpedsclrp_sfc(iisclr_thl) = wpthlp_sfc
if ( iisclr_rt  > 0 ) wpedsclrp_sfc(iisclr_rt)  = wprtp_sfc

! Compute momentum fluxes using ARM Cu formulae

ubar = max( ubmin, sqrt( um_sfc**2 + vm_sfc**2 ) )

bflx = grav/thlm_sfc * wpthlp_sfc

! Compute ustar
ustar = diag_ustar( z, bflx, ubar, z0 )

upwp_sfc = -um_sfc * ustar**2 / ubar
vpwp_sfc = -vm_sfc * ustar**2 / ubar

return
end subroutine arm_0003_sfclyr

!----------------------------------------------------------------------
subroutine arm_0003_init()
!
!       Description:
!       This subroutine initializes the module by reading in forcing
!       data used in the tndcy and sfclyr subroutines.
!----------------------------------------------------------------------
   use file_functions, only: file_read_1d, file_read_2d ! Procedure(s)

   implicit none

! Reading in the files located in
! ../model/arm_0003_forcings/ and storing them in arrays
! Joshua Fasching (October 2007)

   call file_read_1d(10, & 
    file_path//'arm_0003_times.dat', & 
    ntimes, per_line, times)

   call file_read_1d(10, & 
    file_path//'arm_0003_heights.dat', & 
    nz, per_line, z)

   call file_read_1d(10, & 
    file_path//'arm_0003_LE.dat', & 
    ntimes, per_line, LE)

   call file_read_1d(10, & 
    file_path//'arm_0003_SE.dat', & 
    ntimes, per_line, SE)

   call file_read_2d(10, & 
    file_path//'arm_0003_dTdt.dat', & 
    nz, ntimes, per_line, thl_ls)

   call file_read_2d(10, & 
     file_path//'arm_0003_dqdt.dat', & 
     nz, ntimes, per_line, rt_ls)

   call file_read_2d(10, & 
     file_path//'arm_0003_um_obs.dat', & 
     nz, ntimes, per_line, um_obs)

   call file_read_2d(10, & 
     file_path//'arm_0003_vm_obs.dat', & 
     nz, ntimes, per_line, vm_obs)

end subroutine arm_0003_init
!----------------------------------------------------------------------
end module arm_0003

