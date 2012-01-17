!----------------------------------------------------------------------
!$Id$
module cloud_sed_module

implicit none

public :: cloud_drop_sed

private ! Default Scope

contains
!----------------------------------------------------------------------
subroutine cloud_drop_sed( rcm, Ncm, rho_zm, rho, exner, sigma_g, &
                           rcm_mc, thlm_mc )

!       Description:
!       Account for cloud droplet sedimentation in cases like DYCOMS II RF 02

!       References:
!       http://journals.ametsoc.org/doi/abs/10.1175/2008MWR2582.1
!----------------------------------------------------------------------

use grid_class, only: gr ! Variable(s)

use grid_class, only: zt2zm, ddzm ! Procedure(s)

use constants_clubb, only: rho_lw, pi, Cp, Lv ! Variable(s)

 
use stats_type, only: stat_update_var ! Procedure(s)

use stats_variables, only:  & 
    ised_rcm, iFcsed, zt, zm, l_stats_samp ! Variable(s)

use clubb_precision, only: core_rknd ! Variable(s)

implicit none

! External
intrinsic :: exp, log

! Input Variables
real( kind = core_rknd ), intent(in), dimension(gr%nz) :: & 
  rcm,    & ! Liquid water mixing ratio.   [kg/kg]
  rho_zm, & ! Density on moment. grid      [kg/m^3]
  rho,    & ! Density on thermo. grid      [kg/m^3]
  exner,  & ! Exner function               [-]
  Ncm       ! Cloud droplet number conc.   [num/kg]

real( kind = core_rknd ), intent(in) :: &
  sigma_g   ! Geometric std. dev. of cloud droplets falling in a stokes regime.

! real, parameter :: sigma_g = 1.5_core_rknd 
! The DYCOMS2 RF02 intercomparison specified value is 1.5.

! sigma_g is now called into this subroutine set as 1.2 for the "astex_a209" case, and
! as 1.5 for all others. 

! Input/Output Variables
real( kind = core_rknd ), intent(inout), dimension(gr%nz) ::  & 
  rcm_mc, & ! r_t change due to microphysics     [kg/kg)/s] 
  thlm_mc   ! thlm change due to microphysics    [K/s] 

! Local Variables
! Addition for DYCOMS_2
real( kind = core_rknd ), dimension(gr%nz) ::  & 
  Fcsed, & ! Cloud water sedimentation flux       [kg/(m^2 s)]
  sed_rcm  ! d(rcm)/dt due to cloud sedimentation [kg/(m^2 s)]

integer :: k


!=====================================================================
! NOTE:  ADDITION OF RAIN EFFECTS AND CLOUD WATER SEDIMENTATION ON
!        RTM AND THLM.
!=====================================================================
!
! Equations:  rtm = rvm + rcm
!            thlm = thm - ( Lv / (Cp*exner) ) * rcm
!
! When water condenses, latent heat is given off and theta (thm)
! increases by a factor of ( Lv / (Cp*exner) ) * rcm(condensed).
! The opposite effect occurs with evaporation.
!
!=====================================================================
!||     Effect     |  rvm  |  rcm  |  rtm  |  rrainm  |  thm  |  thlm  ||
!|===================================================================|
!|| Sedimentation  |       |       |       |       |       |        ||
!|| Effects of     | stays | incr. | incr. | stays | stays | decr.  ||
!|| Cloud Water.   | same  |       |       | same  | same  |        ||
!|| sed_rcm > 0    |       |       |       |       |       |        ||
!|===================================================================|
!|| Evaporation    |       |       |       |       |       |        ||
!|| of rain to     | incr. | stays | incr. | decr. | decr. | decr.  ||
!|| water vapor.   |       | same  |       |       |       |        ||
!|| cond_rrainm < 0   |       |       |       |       |       |        ||
!|-------------------------------------------------------------------|
!|| Autoconversion |       |       |       |       |       |        ||
!|| of cloud water | stays | decr. | decr. | incr. | stays | incr.  ||
!|| to rain water. | same  |       |       |       | same  |        ||
!|| auto_rrainm > 0   |       |       |       |       |       |        ||
!|-------------------------------------------------------------------|
!|| Accretion of   |       |       |       |       |       |        ||
!|| cloud water by | stays | decr. | decr. | incr. | stays | incr.  ||
!|| rain water.    | same  |       |       |       | same  |        ||
!|| accr_rrainm > 0   |       |       |       |       |       |        ||
!=====================================================================
!
! Note: In HOC, cond_rrainm will always be either negative or zero.
!
! Overall effects of rain and cloud water sedimentation:
!
! (drtm/dt)t  = (drtm/dt)0 
!                       + sed_rcm - cond_rrainm - auto_rrainm - accr_rrainm
! (dthlm/dt)t = (dthlm/dt)0  -  ( Lv / (Cp*exner) ) 
!                       * ( sed_rcm - cond_rrainm - auto_rrainm - accr_rrainm )
!
! Note by Brian Griffin.
!=====================================================================

! Code addition by Brian for cloud water sedimentation.
!
! Sedimentation flux of cloud droplets should be treated by assuming
! a log-normal size distribution of droplets falling in a Stoke's
! regime, in which the sedimentation flux (Fcsed) is given by:
!
! Sedimentation Flux = constant
!                     *[(3/(4*pi*rho_lw*NcV))^(2/3)]
!                     *[(rho*rc)^(5/3)]
!                     *EXP[5*((LOG(sigma_g))^2)]
! constant = 1.19 x 10^8 (m^-1 s^-1)
! sigma_g:  geometric standard deviation
!
! When written for a mass-dependent cloud droplet concentration, Nc:
!
! Sedimentation Flux = constant
!                     *[(3/(4*pi*rho_lw*Nc*rho))^(2/3)]
!                     *[(rho*rc)^(5/3)]
!                     *EXP[5*((LOG(sigma_g))^2)]
!
! According to the above equation, sedimentation flux
! is defined positive downwards.  Therefore, 
!
! (drc/dt)Fcsed = (1.0/rho) * d(Fcsed)/dz


! Define cloud water sedimentation flux on momentum levels.

DO k = 2, gr%nz-1, 1

  IF ( zt2zm(rcm,k) > 0.0_core_rknd .AND. zt2zm( Ncm, k ) > 0.0_core_rknd ) THEN
    Fcsed(k) = 1.19E8_core_rknd & 
                * (   (  3._core_rknd / ( 4.0_core_rknd*pi*rho_lw & 
                                 *zt2zm( Ncm, k )*rho_zm(k) )  ) & 
                   **(2.0_core_rknd/3.0_core_rknd)   ) & 
                * ( ( rho_zm(k)*zt2zm( rcm, k ) )**(5.0_core_rknd/3.0_core_rknd) ) & 
                * EXP( 5.0_core_rknd*( (LOG( sigma_g ))**(2) ) ) ! Known magic number
                                 ! See Ackerman - eq. no. 7
  ELSE
    Fcsed(k) = 0.0_core_rknd
  END IF

END DO ! k=2..gr%nz-1

! Boundary conditions.
Fcsed(1)       = 0.0_core_rknd
Fcsed(gr%nz) = 0.0_core_rknd

! Find drc/dt due to cloud water sedimentation flux.
! This value is defined on thermodynamic levels.
! Fcsed units:  kg (liquid) / [ m^2 * s ]
! Multiply by Lv for units of W / m^2.
! sed_rcm units:  [ kg (liquid) / kg (air) ] / s
sed_rcm = (1.0_core_rknd/rho) * ddzm( Fcsed )

if ( l_stats_samp ) then
 
  call stat_update_var( ised_rcm, sed_rcm, zt )

  call stat_update_var( iFcsed, Fcsed, zm ) 
end if


! + thlm/rtm_microphysics -- cloud water sedimentation.
! Code addition by Brian for cloud water sedimentation.

rcm_mc  = rcm_mc + sed_rcm
thlm_mc = thlm_mc - ( Lv / (Cp*exner) ) * sed_rcm

return
end subroutine cloud_drop_sed

end module cloud_sed_module
