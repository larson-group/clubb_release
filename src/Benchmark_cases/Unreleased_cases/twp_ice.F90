!----------------------------------------------------------------------
! $Id: arm_0003.F90 3363 2009-04-02 21:42:22Z dschanen@uwm.edu $
module twp_ice

!       Description:
!       Contains subroutines for the March 2000 IOP ARM case.
!----------------------------------------------------------------------

implicit none

public :: twp_ice_tndcy, twp_ice_sfclyr

private ! Default Scope

! File path for the forcing files
!character(*), parameter :: file_path = '../input/case_setups/twp_ice_forcings/'

contains

!----------------------------------------------------------------------
subroutine twp_ice_tndcy( time, rho, &
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

use time_dependant_input, only: &
  time_select,  &
  time_f_given, &
  thlm_f_given, &
  rtm_f_given,  &
  um_given,     &
  vm_given,     &
  wm_given,     &
  l_t_dependant

use stats_precision, only: time_precision ! Variable(s)

use error_code, only: clubb_debug ! Procedure(s)

use array_index, only: iisclr_rt, iisclr_thl, iiedsclr_rt, iiedsclr_thl

implicit none

! External
intrinsic :: present

! Input Variables
real(kind=time_precision), intent(in) :: time ! Model time [s]

real, dimension(gr%nnzp), intent(in) :: rho

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
integer :: i1, i2, i3

real :: velocity_omega

!-----------------------------------------------------------------------

! Thetal forcing is equal to the LS tendency given here and the
! interactive calculation done in BUGSrad

if(l_t_dependant) then

     time_frac = -1.0 ! Default initialization

     call time_select( time,size(time_f_given), time_f_given, i1, i2 )

     time_frac = real((time-time_f_given(i1))/(time_f_given(i2)-time_f_given(i1)))

     if( time_frac == -1.0 ) then
       call clubb_debug(1,"times is not sorted in twp_ice_tndcy")
     endif

   ! Interpolate LS thetal tendency to the HOC grid
   ! Time
   thlm_forcing = factor_interp( time_frac, thlm_f_given(:,i2), thlm_f_given(:,i1) )

   ! Interpolate LS rt tendency to the HOC grid
   ! Time
   rtm_forcing = factor_interp( time_frac, rtm_f_given(:,i2), rtm_f_given(:,i1) )

   ! Interpolate um observed to the HOC grid
   ! Time
   um_hoc_grid = factor_interp( time_frac, um_given(:,i2), um_given(:,i1) )

   ! Interpolate vm observed to the HOC grid
   ! Time
   vm_hoc_grid = factor_interp( time_frac, vm_given(:,i2), vm_given(:,i1) )

   omega_hoc_grid = factor_interp( time_frac, wm_given(:,i2), wm_given(:,i1) )
   ! Compute vertical motion
   do i3=2,gr%nnzp
      velocity_omega = omega_hoc_grid(i3) * 100 / 3600 ! convering mb/hr to Pa/s
      wm_zt(i3) = -velocity_omega / (grav * rho(i3))
   end do

   ! Boundary condition
   wm_zt(1) = 0.0        ! Below surface

   ! Interpolation
   wm_zm = zt2zm( wm_zt )

   um_hoc_grid (1) = um_hoc_grid(2)
   vm_hoc_grid (1) = vm_hoc_grid(2)

   ! Test scalars with thetal and rt if desired
   if ( iisclr_thl > 0 ) sclrm_forcing(:,iisclr_thl) = thlm_forcing
   if ( iisclr_rt  > 0 ) sclrm_forcing(:,iisclr_rt)  = rtm_forcing

   if ( iiedsclr_thl > 0 ) edsclrm_forcing(:,iiedsclr_thl) = thlm_forcing
   if ( iiedsclr_rt  > 0 ) edsclrm_forcing(:,iiedsclr_rt)  = rtm_forcing
endif
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
!Cm   = C_m_20 * ((log(20/z0))/(log(z/z0))) * & 
!       ((log(20/z0))/(log(z/z0)))             
! Modification in case lowest model level isn't at 10 m, from ATEX specification
Ch   = C_h_20 * ((log(20/z0))/(log(z/z0))) * & 
       ((log(20/z0))/(log(z/z0)))          
! Modification in case lowest model level isn't at 10 m, from ATEX specification
Cq   = C_q_20 * ((log(20/z0))/(log(z/z0))) * & 
       ((log(20/z0))/(log(z/z0)))  

wpthlp_sfc = compute_wpthlp_sfc( Ch, ubar, thlm_sfc, sst, exner_sfc )
wprtp_sfc  = compute_wprtp_sfc( Cq, ubar, rtm, sat_mixrat_liq(psfc,sst) )
!upwp_sfc   = -um_sfc * Cm * ubar  ! m^2 s^-2
!vpwp_sfc   = -vm_sfc * Cm * ubar  ! m^2 s^-2

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
end module twp_ice
