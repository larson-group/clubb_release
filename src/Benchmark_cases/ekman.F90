!----------------------------------------------------------------------
! $Id$
module ekman

!       Description:
!       Contains subroutines for Hing Ong's Ekman layer case.
!
!       References:
!----------------------------------------------------------------------

  implicit none

  public :: ekman_sfclyr

  private ! Default Scope

  contains

!----------------------------------------------------------------------
  subroutine ekman_sfclyr( ngrdcol, z, &
                           um_sfc, vm_sfc, ubar, &
                           upwp_sfc, vpwp_sfc, &
                           wpthlp_sfc, wprtp_sfc, ustar )

!       Description:
!       This subroutine computes surface fluxes of horizontal momentum,
!       heat and moisture for a neutral Ekman layer case. 

!       References:
!----------------------------------------------------------------------

    use sfc_flux, only: compute_momentum_flux

    use diag_ustar_module, only: diag_ustar ! Variable(s)

    use clubb_precision, only: core_rknd ! Variable(s)

    implicit none

    ! Constant
    real( kind = core_rknd ), parameter ::  &
      z0 = 0.0001_core_rknd ! Momentum roughness height

    ! Input Variables
    integer, intent(in) :: &
      ngrdcol

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) ::  &
      z,               & ! Height at zt(2)       [m]
      um_sfc,          & ! um at zt(2)           [m/s]
      vm_sfc,          & ! vm at zt(2)           [m/s]
      ubar

    ! Output variables
    real( kind = core_rknd ), dimension(ngrdcol), intent(out) ::  &
      upwp_sfc,     & ! u'w' at (1)      [(m/s)^2]
      vpwp_sfc,     & ! v'w' at (1)      [(m/s)^2]
      wpthlp_sfc,   & ! w'th_l' at (1)   [(m K)/s]  
      wprtp_sfc,    & ! w'r_t' at (1)    [(m kg)/(s kg)]
      ustar           ! surface friction velocity [m/s]

    ! Local variable
    integer :: i

    !---- Begin code

    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol

      ! No heat and moisture fluxes
      wpthlp_sfc(i) = 0.0_core_rknd
      wprtp_sfc(i)  = 0.0_core_rknd 

      ! Set ustar
      ustar(i) = 0.3_core_rknd ! diag_ustar( z(i), 0.0_core_rknd, ubar(i), z0 )

    end do

    call compute_momentum_flux( ngrdcol, um_sfc, vm_sfc, ubar, ustar, &
                                upwp_sfc, vpwp_sfc )

    return
    
  end subroutine ekman_sfclyr

end module ekman
