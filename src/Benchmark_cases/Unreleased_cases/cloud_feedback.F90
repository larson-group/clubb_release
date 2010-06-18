!----------------------------------------------------------------------
! $Id$
module cloud_feedback

!       Description:
!       Contains subroutines for the Cloud Feedback cases.
!----------------------------------------------------------------------

implicit none

private ! Default Scope

public :: cloud_feedback_sfclyr

contains

!----------------------------------------------------------------------
subroutine cloud_feedback_sfclyr( runtype, sfctype, &
                                  thlm_sfc, rtm_sfc, &
                                  um_sfc, vm_sfc, &
                                  psfc, Tsfc, &
                                  upwp_sfc, vpwp_sfc, &
                                  wpthlp_sfc, wprtp_sfc, ustar, &
                                  wpsclrp_sfc, wpedsclrp_sfc )

!       Description:
!       Sets up surface information for the cloud feedback case

!       References:
!----------------------------------------------------------------------

use constants_clubb, only: pi, grav, Lv, Cp, p0, kappa ! Variable(s)

use parameters_model, only: sclr_dim, edsclr_dim ! Variable(s)

use saturation, only: sat_mixrat_liq ! Variable(s)

use array_index, only: iisclr_rt, iisclr_thl, iiedsclr_rt, iiedsclr_thl ! Variable(s)

use stats_precision, only: time_precision ! Variable(s)

use array_index, only:  & 
    iisclr_thl, iisclr_rt ! Variable(s)

use surface_flux, only: &
    compute_ubar, compute_momentum_flux, compute_wprtp_sfc, compute_wpthlp_sfc

implicit none

intrinsic :: max, sqrt

! Input Variables
character(len=50), intent(in) :: runtype ! The case that is being run

integer, intent(in) :: sfctype

real, intent(in) ::  & 
  thlm_sfc,  & ! thlm at (2)         [m/s]
  rtm_sfc,   & ! rtm at (2)          [kg/kg]
  Tsfc,      & ! Temperature         [K]
  psfc,      & ! Surface pressure    [Pa]
  um_sfc,    & ! um at (2)           [m/s]
  vm_sfc       ! vm at (2)           [m/s]

! Output variables
real, intent(out) ::  & 
  upwp_sfc,     & ! u'w' at (1)      [m^2/s^2]
  vpwp_sfc,     & ! v'w'at (1)       [m^2/s^2]
  wpthlp_sfc,   & ! w'th_l' at (1)   [(m K)/s]  
  wprtp_sfc,    & ! w'r_t'(1) at (1) [(m kg)/(s kg)]
  ustar           ! surface friction velocity [m/s]

real, intent(out), dimension(sclr_dim) ::  & 
  wpsclrp_sfc     ! Passive scalar surface flux      [units m/s] 

real, intent(out), dimension(edsclr_dim) ::  & 
  wpedsclrp_sfc   ! Passive eddy-scalar surface flux [units m/s]

! Constants
real, parameter :: & 
!  rho_sfc_flux = 1.0, &
  C_10    = 0.0013      ! Drag coefficient, defined by ATEX specification

! Internal variables
real :: &
  exner_sfc ! Value of exner at the surface [-]

real :: & 
  ubar 


! Calculate exner_sfc based on psfc.
exner_sfc = ( psfc / p0 )**kappa

ubar = compute_ubar( um_sfc, vm_sfc )

! Just set ustar = 0.3
ustar = 0.3

! Get rid of a compiler warning
if( runtype == "anything" .or. thlm_sfc == 1 .or. rtm_sfc == 1 .or. &
        exner_sfc == 1 .or. psfc == 1 .or. Tsfc == 1) then
    ustar = 0.3
end if

!--------------------------------------------------------------------------------
! Email way
! wprtp = value_from_forcings_file_in_W_m**2 / ( rho_sfc_flux * Lv )
! wpthlp = value_from_forcings_file_in_W_m**2 / ( rho_sfc_flux * Cp )

!lhflx(1) = 0.001 * ubar * rho_sfc_flux * Lv * ( sat_mixrat_liq( psfc, Tsfc ) - & 
!                                                sat_mixrat_liq( psfc, T_in_K ) * 0.8 )
!shflx(1) = 0.001 * ubar * rho_sfc_flux * Cp * ( Tsfc - T_in_K )

! If this is the S6 case, fudge the values of the fluxes using values from the forcings
!if ( runtype == "cloud_feedback_s6" .or. runtype == "cloud_feedback_s6_p2k" ) then
!    wprtp_sfc = lhflx(1) / ( 1.0 * Lv )
!    wpthlp_sfc = shflx(1) / ( 1.0 * Cp )
!else
!    wprtp_sfc = compute_wprtp_sfc( C_10, ubar, rtm_sfc, sat_mixrat_liq( psfc, Tsfc ) )
!    wpthlp_sfc = compute_wpthlp_sfc( C_10, ubar, thlm_sfc, & 
!                                     Tsfc, exner_sfc )
!end if

call compute_momentum_flux( um_sfc, vm_sfc, ubar, ustar, &
                            upwp_sfc, vpwp_sfc )

! 
if ( sfctype == 1 ) then
  wprtp_sfc = compute_wprtp_sfc( C_10, ubar, rtm_sfc, sat_mixrat_liq( psfc, Tsfc ) )
  wpthlp_sfc = compute_wpthlp_sfc( C_10, ubar, thlm_sfc, & 
                                   Tsfc, exner_sfc )

end if

wpsclrp_sfc(:)   = 0.0 ! Initialize flux to 0 
wpedsclrp_sfc(:) = 0.0 ! Initialize flux to 0 

! Set passive scalars be equal to rt and theta_l when checking the scalar code
if ( iisclr_thl > 0 ) wpsclrp_sfc(iisclr_thl) = wpthlp_sfc
if ( iisclr_rt  > 0 ) wpsclrp_sfc(iisclr_rt)  = wprtp_sfc

if ( iiedsclr_thl > 0 ) wpedsclrp_sfc(iiedsclr_thl) = wpthlp_sfc
if ( iiedsclr_rt  > 0 ) wpedsclrp_sfc(iiedsclr_rt)  = wprtp_sfc


return
end subroutine cloud_feedback_sfclyr

end module cloud_feedback
