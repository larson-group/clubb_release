!----------------------------------------------------------------------
! $Id$
module neutral_case

!       Description:
!       Contains subroutines for the Moeng and Sullivan (1994) shear-driven case.
!       Modeled after George Bryan's shear-driven case in CM1.
!
!       References:
!----------------------------------------------------------------------

  implicit none

  public :: neutral_case_sfclyr

  private ! Default Scope

  contains

!----------------------------------------------------------------------
  subroutine neutral_case_sfclyr( time, &
                                  ! z, thlm_sfc, &
                                  um_sfc, vm_sfc, ubar, &
                                  upwp_sfc, vpwp_sfc, &
                                  wpthlp_sfc, wprtp_sfc, ustar )

!       Description:
!       This subroutine computes surface fluxes of horizontal momentum,
!       heat and moisture for a shear-driven, neutral case. 

!       References:
!       Moeng and Sullivan (1994). 
!----------------------------------------------------------------------

    use sfc_flux, only: compute_momentum_flux

    !use constants_clubb, only: grav ! Variable(s)

    use clubb_precision, only: time_precision, core_rknd ! Variable(s)

    use diag_ustar_module, only: diag_ustar ! Variable(s)

    implicit none

    ! Input Variables
    real(time_precision), intent(in) ::  & 
      time    ! the current time [s]

    real( kind = core_rknd ), intent(in) ::  &
      ! z,               & ! Height at zt(2)       [m]
      ! thlm_sfc,        & ! Theta_l at zt(2)      [K]
      um_sfc,          & ! um at zt(2)           [m/s]
      vm_sfc,          & ! vm at zt(2)           [m/s]
      ubar

    ! Output variables
    real( kind = core_rknd ), intent(out) ::  &
      upwp_sfc,     &
      vpwp_sfc,     & 
      wpthlp_sfc,   & ! w'th_l' at (1)   [(m K)/s]  
      wprtp_sfc,    & ! w'r_t' at (1)    [(m kg)/(s kg)]
      ustar           ! surface friction velocity [m/s]

    ! real( kind = core_rknd ), parameter ::  &
    !  z0 = 0.16_core_rknd ! Roughness length  [m]

    ! real( kind = core_rknd ) :: bflx 

!---- Begin code

    ! Compute heat and moisture fluxes --- turn off heat flux after 3000 s
    if (time > 80880_time_precision) then
      wpthlp_sfc = 0.0_core_rknd
    else
      wpthlp_sfc = 0.05_core_rknd 
    end if

    ! No moisture in this case.
    wprtp_sfc  = 0.0_core_rknd 

    ! Heat flux in units of (m2/s3) (needed by diag_ustar)
    !bflx = grav/thlm_sfc * wpthlp_sfc ! grav variable is commented out

    ! Compute ustar
    ustar = 0.5_core_rknd !diag_ustar( z, bflx, ubar, z0 )

    call compute_momentum_flux( um_sfc, vm_sfc, ubar, ustar, &
                              upwp_sfc, vpwp_sfc )

    return
  end subroutine neutral_case_sfclyr

end module neutral_case
