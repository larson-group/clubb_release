!----------------------------------------------------------------------
! $Id$
module bomex

!       Description:
!       Contains subroutines for the GCSS BOMEX case.
!----------------------------------------------------------------------

implicit none

public :: bomex_tndcy, bomex_sfclyr

private ! Default Scope

contains

!----------------------------------------------------------------------
subroutine bomex_tndcy( rtm, radht, & 
                        thlm_forcing, rtm_forcing, & 
                        sclrm_forcing, edsclrm_forcing )
!       Description:
!       Subroutine to set theta and water tendencies for BOMEX case

!       References:
!       <http://www.knmi.nl/~siebesma/gcss/bomexcomp.init.html>
!----------------------------------------------------------------------

use grid_class, only: gr ! Variable(s)

use grid_class, only: zt2zm ! Procedure(s)

use spec_hum_to_mixing_ratio, only: &
    force_spec_hum_to_mixing_ratio ! Procedure(s)

use parameters_radiation, only: rad_scheme ! Variable(s)

use parameters_model, only: sclr_dim, edsclr_dim ! Variable(s)

use array_index, only: iisclr_rt, iisclr_thl, iiedsclr_rt, iiedsclr_thl ! Variable(s)

implicit none

! Input Variable
real, intent(in), dimension(gr%nnzp) :: &
  rtm    ! Total water mixing ratio (thermodynamic levels)        [kg/kg]

! Output Variables
real, intent(out), dimension(gr%nnzp) :: & 
  radht,         & ! Radiative heating rate                       [K/s]
  thlm_forcing,  & ! Liquid water potential temperature tendency  [K/s]
  rtm_forcing      ! Total water mixing ratio tendency            [kg/kg/s]

real, intent(out), dimension(gr%nnzp,sclr_dim) :: & 
  sclrm_forcing ! Passive scalar forcing        [units vary/s]

real, intent(out), dimension(gr%nnzp,edsclr_dim) :: & 
  edsclrm_forcing ! Eddy-passive scalar forcing [units vary/s]

! Local Variables
real, dimension(gr%nnzp) :: &
  qtm_forcing  ! Specified total water spec. humidity tendency    [kg/kg/s]

integer :: k

if ( trim( rad_scheme ) == "simplified" ) then

! Radiative theta-l tendency
  do k = 2, gr%nnzp

    if ( gr%zt(k) >= 0. .and. gr%zt(k) < 1500. ) then
      radht(k) = -2.315e-5
    else if ( gr%zt(k) >= 1500. .and. gr%zt(k) < 2500. ) then
      radht(k) & 
        = - 2.315e-5  & 
          + 2.315e-5  & 
            * ( gr%zt(k) - 1500. ) / ( 2500. - 1500. )
    else
      radht(k) = 0.
    end if

  end do ! k=2..gr%nnzp

  ! Boundary condition
  radht(1) = 0.0

  thlm_forcing = radht
else ! Compute radht interactively with BUGSrad

  thlm_forcing = 0.0

end if ! simplified

! Large scale advective moisture tendency
! The BOMEX specifications give large-scale advective moisture tendency in
! terms of total water specific humidity.
do k = 2, gr%nnzp

   if ( gr%zt(k) >= 0. .and. gr%zt(k) < 300. ) then
      qtm_forcing(k) = -1.2e-8
   else if ( gr%zt(k) >= 300. .and. gr%zt(k) < 500. ) then
      qtm_forcing(k)  & 
        = - 1.2e-8  & 
            * ( 1. - ( gr%zt(k) - 300. )/( 500. - 300. ) )
   else
      qtm_forcing(k) = 0.
   end if

   ! Convert forcings from terms of total water specific humidity to terms of
   ! total water mixing ratio.
   call force_spec_hum_to_mixing_ratio( rtm(k), qtm_forcing(k), rtm_forcing(k) )

end do


! Boundary conditions
thlm_forcing(1) = 0.0  ! Below surface
rtm_forcing(1)  = 0.0  ! Below surface

! Test scalars with thetal and rt if desired
if ( iisclr_thl > 0 ) sclrm_forcing(:,iisclr_thl) = thlm_forcing
if ( iisclr_rt  > 0 ) sclrm_forcing(:,iisclr_rt)  = rtm_forcing

if ( iiedsclr_thl > 0 ) edsclrm_forcing(:,iiedsclr_thl) = thlm_forcing
if ( iiedsclr_rt  > 0 ) edsclrm_forcing(:,iiedsclr_rt)  = rtm_forcing

return
end subroutine bomex_tndcy

!----------------------------------------------------------------------
subroutine bomex_sfclyr( um_sfc, vm_sfc, rtm_sfc, & 
                         upwp_sfc, vpwp_sfc, & 
                         wpthlp_sfc, wprtp_sfc, ustar, & 
                         wpsclrp_sfc, wpedsclrp_sfc )

!       Description:
!       This subroutine computes surface fluxes of horizontal momentum,
!       heat and moisture according to GCSS BOMEX specifications

!       References:
!----------------------------------------------------------------------

use parameters_model, only: sclr_dim, edsclr_dim  ! Variable(s)

use array_index, only: iisclr_rt, iisclr_thl, iiedsclr_rt, iiedsclr_thl

use surface_flux, only: compute_momentum_flux, compute_ubar

use spec_hum_to_mixing_ratio, only: &
    flux_spec_hum_to_mixing_ratio ! Procedure(s)

implicit none


! Input Variables
real, intent(in) ::  & 
  um_sfc,  & ! um(2)  [m/s]
  vm_sfc,  & ! vm(2)  [m/s]
  rtm_sfc    ! rtm(2) [kg/kg]

! Output variables
real, intent(out) ::  & 
  upwp_sfc,     & ! u'w' at (1)      [m^2/s^2]
  vpwp_sfc,     & ! v'w'at (1)       [m^2/s^2]
  wpthlp_sfc,   & ! w'th_l' at (1)   [(m K)/s]  
  wprtp_sfc,    & ! w'r_t' at (1)    [(m kg)/(s kg)]
  ustar           ! surface friction velocity [m/s]

real,  dimension(sclr_dim), intent(out) ::  & 
  wpsclrp_sfc        ! Passive scalar surface flux      [units m/s] 

real,  dimension(edsclr_dim), intent(out) ::  & 
  wpedsclrp_sfc      ! Passive eddy-scalar surface flux [units m/s]

! Local variables
real :: wpqtp_sfc  ! w'q_t' at (1)         [(m kg)/(s kg)]
real :: ubar       ! mean sfc wind speed   [m/s]

! Declare the value of ustar.
ustar = 0.28

! Compute heat and moisture fluxes

wpthlp_sfc = 8.e-3
! The BOMEX specifications give surface moisture flux in terms of total water
! specific humidity.
wpqtp_sfc  = 5.2e-5

! Convert flux from terms of total water specific humidity to terms of total
! water mixing ratio.
call flux_spec_hum_to_mixing_ratio( rtm_sfc, wpqtp_sfc, wprtp_sfc )

! Compute momentum fluxes

ubar = compute_ubar( um_sfc, vm_sfc )

call compute_momentum_flux( um_sfc, vm_sfc, ubar, ustar, &
                            upwp_sfc, vpwp_sfc )

! Let passive scalars be equal to rt and theta_l for now
if ( iisclr_thl > 0 ) wpsclrp_sfc(iisclr_thl) = wpthlp_sfc
if ( iisclr_rt  > 0 ) wpsclrp_sfc(iisclr_rt)  = wprtp_sfc

if ( iiedsclr_thl > 0 ) wpedsclrp_sfc(iiedsclr_thl) = wpthlp_sfc
if ( iiedsclr_rt  > 0 ) wpedsclrp_sfc(iiedsclr_rt)  = wprtp_sfc

return
end subroutine bomex_sfclyr

end module bomex
