!----------------------------------------------------------------------
! $Id: arm_0003.F90 3363 2009-04-02 21:42:22Z dschanen@uwm.edu $
module twp_ice

!       Description:
!       Contains subroutines for the March 2000 IOP ARM case.
!----------------------------------------------------------------------

implicit none

public :: twp_ice_tndcy, twp_ice_sfclyr, twp_ice_init

private ! Default Scope

! Constant Parameters
integer, parameter :: ntimes = 214, nz = 40, & 
 per_line = 5

! File path for the forcing files
!character(*), parameter :: file_path = '../input/case_setups/twp_ice_forcings/'

real, dimension(ntimes) :: times               ! Time from day0      [s]
real, dimension(nz) :: z                       ! Height              [m]
real, dimension(nz, ntimes) :: thl_ls          ! Potential Temperature
                                               ! Tendency            [K/s]
real, dimension(nz, ntimes) :: rt_ls           ! Water Vapor Advective
                                               ! Tendency            [Kg/Kg/s]
real, dimension(nz, ntimes) :: um_obs          ! Obs. wind u         [m/s]
real, dimension(nz, ntimes) :: vm_obs          ! Obs. wind v         [m/s]
real, dimension(nz, ntimes) :: omega_forcing   ! Vertical velocity forcing [mb/hr]

contains

!----------------------------------------------------------------------
subroutine twp_ice_tndcy( time, p_in_Pa, thvm, &
                           wm_zt, wm_zm, thlm_forcing, &
                           rtm_forcing, um_hoc_grid, vm_hoc_grid, &
                           sclrm_forcing, edsclrm_forcing )
!       Description:
!       Subroutine to set thetal and total water tendencies for ARM 0003 case

!       References:
!       None
!----------------------------------------------------------------------

use constants, only: Rd, grav ! Variable(s)

use interpolation, only: factor_interp ! Procedure(s)

use grid_class, only: gr ! Variable(s)

use grid_class, only: zt2zm ! Procedure(s)

use parameters_model, only: sclr_dim, edsclr_dim ! Variable(s)

use interpolation, only: zlinterp_fnc ! Procedure(s)

use stats_precision, only: time_precision ! Variable(s)

use error_code, only: clubb_debug ! Procedure(s)

use array_index, only: iisclr_rt, iisclr_thl, iiedsclr_rt, iiedsclr_thl

implicit none

! External
intrinsic :: present

! Input Variables
real(kind=time_precision), intent(in) :: time ! Model time [s]

real, dimension(gr%nnzp), intent(in) ::  &
  p_in_Pa,  &    ! Pressure [Pa]
  thvm

! Output Variables
real, intent(out), dimension(gr%nnzp) ::  & 
  thlm_forcing,  & ! Liquid water potential temperature tendency  [K/s]
  wm_zm,         & ! Vertical velocity on moment. grid            [m/s]
  wm_zt,         & ! Vertical velocity on thermo. grid            [m/s]
  rtm_forcing,   & ! Total water mixing ratio tendency            [kg/kg/s]
  um_hoc_grid,   & ! Observed wind, for nudging                   [m/s]
  vm_hoc_grid      ! Observed wind, for nudging                   [m/s]

real, dimension(gr%nnzp) :: omega_hoc_grid

real, intent(out), dimension(gr%nnzp,sclr_dim) ::  & 
  sclrm_forcing ! Passive scalar tendency       [units/s]

real, intent(out), dimension(gr%nnzp,edsclr_dim) ::  & 
  edsclrm_forcing ! Eddy-passive scalar tendency[units/s]

real :: time_frac
integer :: i, i1, i2

real :: velocity_omega

real, dimension(nz) :: thlm_t_interp, & 
  rtm_t_interp, um_obs_t_interp, vm_obs_t_interp, omega_interp

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
           (1,"times not sorted in twp_ice_tndcy")
endif

! Interpolate LS thetal tendency to the HOC grid
! Time
thlm_t_interp = factor_interp( time_frac, thl_ls(:,i2), thl_ls(:,i1) )
! Vertical
thlm_forcing(1:gr%nnzp) = zlinterp_fnc( gr%nnzp, nz, gr%zt(:), z, thlm_t_interp )
!thlm_forcing(1:gr%nnzp) = 0

! Interpolate LS rt tendency to the HOC grid
! Time
rtm_t_interp = factor_interp( time_frac, rt_ls(:,i2), rt_ls(:,i1) )
! Vertical
rtm_forcing(1:gr%nnzp) = zlinterp_fnc & 
         ( gr%nnzp, nz, gr%zt(:), z, rtm_t_interp )
!rtm_forcing(1:gr%nnzp) = 0
! Modified by Joshua Fasching (October 2007)

! Interpolate um observed to the HOC grid
! Time
um_obs_t_interp = factor_interp( time_frac, um_obs(:,i2), um_obs(:,i1) )
! Vertical
um_hoc_grid(1:gr%nnzp) = zlinterp_fnc & 
         ( gr%nnzp, nz, gr%zt(:), z, um_obs_t_interp )

! Interpolate vm observed to the HOC grid
! Time
vm_obs_t_interp = factor_interp( time_frac, vm_obs(:,i2), vm_obs(:,i1) )
! Vertical
vm_hoc_grid(1:gr%nnzp) = zlinterp_fnc & 
         ( gr%nnzp, nz, gr%zt(:), z, vm_obs_t_interp )

omega_interp = factor_interp( time_frac, omega_forcing(:,i2), omega_forcing(:,i1) )
! Vertical
omega_hoc_grid(1:gr%nnzp) = zlinterp_fnc & 
         ( gr%nnzp, nz, gr%zt(:), z, omega_interp )

! Compute vertical motion
do i=2,gr%nnzp
   velocity_omega = omega_hoc_grid(i) * 100 / 3600 ! convering mb/hr to Pa/s
   wm_zt(i) = -velocity_omega * Rd * thvm(i) / p_in_Pa(i) / grav
   !wm_zt(i) = 0.
end do

! Boundary condition
wm_zt(1) = 0.0        ! Below surface

! Interpolation
wm_zm = zt2zm( wm_zt )

! Added by Joshua Fasching (October 27 2007)
um_hoc_grid (1) = um_hoc_grid(2)
vm_hoc_grid (1) = vm_hoc_grid(2)

! Test scalars with thetal and rt if desired
if ( iisclr_thl > 0 ) sclrm_forcing(:,iisclr_thl) = thlm_forcing
if ( iisclr_rt  > 0 ) sclrm_forcing(:,iisclr_rt)  = rtm_forcing

if ( iiedsclr_thl > 0 ) edsclrm_forcing(:,iiedsclr_thl) = thlm_forcing
if ( iiedsclr_rt  > 0 ) edsclrm_forcing(:,iiedsclr_rt)  = rtm_forcing

return
end subroutine twp_ice_tndcy
!----------------------------------------------------------------------
subroutine twp_ice_sfclyr( z, sst, exner_sfc, thlm_sfc, & 
                            um_sfc, vm_sfc, rtm, psfc,  & 
                            upwp_sfc, vpwp_sfc, & 
                            wpthlp_sfc, wprtp_sfc, ustar, & 
                            wpsclrp_sfc, wpedsclrp_sfc )
!       Description:
!       This subroutine computes surface fluxes of horizontal momentum,
!       heat and moisture according to GCSS ARM specifications
!----------------------------------------------------------------------

use constants, only: Cp, Lv, grav ! Variable(s)

use parameters_model, only: sclr_dim, edsclr_dim ! Variable(s)

use saturation, only: sat_mixrat_liq ! Procedure(s)

use stats_precision, only: time_precision ! Variable(s)

use array_index, only: iisclr_rt, iisclr_thl, iiedsclr_rt, iiedsclr_thl ! Variable(s)

use surface_flux, only: compute_ubar, compute_momentum_flux, &
                          compute_wpthlp_sfc, compute_wprtp_sfc
implicit none

intrinsic :: max, sqrt

! Constants
real, parameter :: & 
  C_m_20  = 0.001229,  & ! Drag coefficient, defined by RICO 3D specification
  C_h_20  = 0.001094,  & ! Drag coefficient, defined by RICO 3D specification
  C_q_20  = 0.001133,  & ! Drag coefficient, defined by RICO 3D specification
  z0      = 0.00015      ! Roughness length, defined by ATEX specification

real, intent(in) ::  & 
  z,             & ! Height at zt=2      [s] 
  sst,           & ! Sea surface temp    [K]
  exner_sfc,     & ! Exner function at (2) 
  um_sfc,        & ! um at (2)           [m/s]
  vm_sfc,        & ! vm at (2)           [m/s]
  thlm_sfc,      & ! thlm at (2)         [m/s]
  rtm,           & ! rt at (2)           [kg/kg]
  psfc             ! surface pressure    [Pa]

! Output variables
real, intent(out) ::  & 
  upwp_sfc,     & ! u'w' at (1)      [m^2/s^2]
  vpwp_sfc,     & ! v'w'at (1)       [m^2/s^2]
  wpthlp_sfc,   & ! w'th_l' at (1)   [(m K)/s]  
  wprtp_sfc,    & ! w'r_t'(1) at (1) [(m kg)/(s kg)]
  ustar           ! surface friction velocity [m/s]

real, intent(out), dimension(sclr_dim) ::  & 
  wpsclrp_sfc      ! Passive scalar surface flux      [units m/s] 

real, intent(out), dimension(edsclr_dim) ::  & 
  wpedsclrp_sfc    ! Passive eddy-scalar surface flux [units m/s]

! Internal variables
real :: & 
  ubar, & ! This is root (u^2 + v^2), per ATEX and RICO spec.
  Cm,   & ! This is C_m_20 scaled to the height of the lowest model level.
  Ch,   & ! This is C_h_20 scaled to the height of the lowest model level.
  Cq      ! This is C_q_20 scaled to the height of the lowest model level.
!----------------------------------------------------------------------

! Declare the value of ustar.
ustar = 0.3

! Define variable values
ubar = compute_ubar( um_sfc, vm_sfc )
        
! Modification in case lowest model level isn't at 10 m, from ATEX specification
Cm   = C_m_20 * ((log(20/z0))/(log(z/z0))) * & 
       ((log(20/z0))/(log(z/z0)))             
! Modification in case lowest model level isn't at 10 m, from ATEX specification
Ch   = C_h_20 * ((log(20/z0))/(log(z/z0))) * & 
       ((log(20/z0))/(log(z/z0)))          
! Modification in case lowest model level isn't at 10 m, from ATEX specification
Cq   = C_q_20 * ((log(20/z0))/(log(z/z0))) * & 
       ((log(20/z0))/(log(z/z0)))  

wpthlp_sfc = compute_wpthlp_sfc( Ch, ubar, thlm_sfc, sst, exner_sfc )
wprtp_sfc  = compute_wprtp_sfc( Cq, ubar, rtm, sat_mixrat_liq(psfc,sst) )
upwp_sfc   = -um_sfc * Cm * ubar  ! m^2 s^-2
vpwp_sfc   = -vm_sfc * Cm * ubar  ! m^2 s^-2

! Compute momentum fluxes
call compute_momentum_flux( um_sfc, vm_sfc, ubar, ustar, &
                            upwp_sfc, vpwp_sfc )

! Let passive scalars be equal to rt and theta_l for testing
if ( iisclr_thl > 0 ) wpsclrp_sfc(iisclr_thl) = wpthlp_sfc
if ( iisclr_rt  > 0 ) wpsclrp_sfc(iisclr_rt)  = wprtp_sfc

if ( iiedsclr_thl > 0 ) wpedsclrp_sfc(iiedsclr_thl) = wpthlp_sfc
if ( iiedsclr_rt  > 0 ) wpedsclrp_sfc(iiedsclr_rt)  = wprtp_sfc

return
end subroutine twp_ice_sfclyr

!----------------------------------------------------------------------
subroutine twp_ice_init( iunit, file_path )
!
!       Description:
!       This subroutine initializes the module by reading in forcing
!       data used in the tndcy and sfclyr subroutines.
!----------------------------------------------------------------------
   use file_functions, only: file_read_1d, file_read_2d ! Procedure(s)

   implicit none

   integer, intent(in) :: iunit ! File unit number

   character(len=*), intent(in) :: &
     file_path ! Path to the forcing files

   ! Reading in the files located in
   ! ../model/twp_ice_forcings/ and storing them in arrays
   ! Joshua Fasching (October 2007)

   ! ---- Begin Code ----

   call file_read_1d(iunit, & 
    file_path//'twp_ice_times.dat', & 
    ntimes, per_line, times)

   call file_read_1d(iunit, & 
    file_path//'twp_ice_heights.dat', & 
    nz, per_line, z)

   call file_read_2d(iunit, & 
    file_path//'twp_ice_dTdt.dat', & 
    nz, ntimes, per_line, thl_ls)

   call file_read_2d(iunit, & 
     file_path//'twp_ice_dqdt.dat', & 
     nz, ntimes, per_line, rt_ls)

   call file_read_2d(iunit, & 
     file_path//'twp_ice_um_obs.dat', & 
     nz, ntimes, per_line, um_obs)

   call file_read_2d(iunit, & 
     file_path//'twp_ice_vm_obs.dat', & 
     nz, ntimes, per_line, vm_obs)

   call file_read_2d( iunit, &
     file_path//'twp_ice_omega.dat', &
     nz, ntimes, per_line, omega_forcing)

end subroutine twp_ice_init
!----------------------------------------------------------------------
end module twp_ice

