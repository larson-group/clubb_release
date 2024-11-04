!-----------------------------------------------------------------------
! $Id$
!===============================================================================
module spec_hum_to_mixing_ratio

  ! Description:
  ! This module contains subroutines that take various moisture variables that
  ! were given in terms of specific humidity and converts them to terms of
  ! mixing ratio.

  implicit none

  private ! Default Scope

  public :: flux_spec_hum_to_mixing_ratio, &
            force_spec_hum_to_mixing_ratio

  contains

  !=============================================================================
  subroutine flux_spec_hum_to_mixing_ratio( ngrdcol, rtm_zm, wpqtp, wprtp )

    ! Description:
    ! This subroutine takes a flux given in terms of specific humidity, w'q_t',
    ! and converts it to a flux given in terms of mixing ratio, w'r_t'.
    !
    ! This relationship is given by the equation:
    !
    ! w'r_t' = ( 1 + r_tm )^2 * w'q_t';
    !
    ! where r_tm is the mean total water mixing ratio.  Higher-order terms,
    ! which were comparatively small in magnitude, were neglected for purposes
    ! of linearizing the equation.
    !
    ! This subroutine has been written very generally, in order to allow for
    ! this conversion at any vertical level.  However, the most likely
    ! application of this is at the surface, taking a specified surface flux in
    ! terms of specific humidity, w'q_t'|_sfc, and converting it to a surface
    ! flux in terms of mixing ratio, w'r_t'|_sfc.
    
    !-------------------------------------------------------------------

    use clubb_precision, only: core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      ngrdcol

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) ::  &
      rtm_zm, & ! r_tm interpolated to momentum level (k)   [kg/kg]
      wpqtp     ! w'q_t' flux at momentum level (k)         [(kg/kg) m/s]

    ! Output Variable
    real( kind = core_rknd ), dimension(ngrdcol), intent(out) ::  &
      wprtp     ! w'r_t' flux at momentum level (k)         [(kg/kg) m/s]

    integer :: i

    ! Solve for flux in terms of total water mixing ratio.
    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      wprtp(i) = ( 1.0_core_rknd + rtm_zm(i) )**2 * wpqtp(i)
    end do

    return
  end subroutine flux_spec_hum_to_mixing_ratio

  !=============================================================================
  subroutine force_spec_hum_to_mixing_ratio( ngrdcol, nzt, rtm, qtm_forcing, rtm_forcing )

    ! Description:
    ! This subroutine takes a forcing given in terms of specific humidity,
    ! d(q_tm)/dt|_f, and converts it to a forcing given in terms of mixing
    ! ratio, d(r_tm)/dt|_f.
    !
    ! This relationship is given by the equation:
    !
    ! d(r_tm)/dt|_f = ( 1 + r_tm )^2 * d(q_tm)/dt|_f;
    !
    ! where r_tm is the mean total water mixing ratio.

    !-------------------------------------------------------------------

    use clubb_precision, only: core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      ngrdcol, &
      nzt

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) ::  &
      rtm,        & ! r_tm at thermodynamic level (k)            [kg/kg]
      qtm_forcing   ! d(q_tm)/dt|_f at thermodynamic level (k)   [(kg/kg)/s]

    ! Output Variable
    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) ::  &
      rtm_forcing   ! d(r_tm)/dt|_f at thermodynamic level (k)   [(kg/kg)/s]

    integer :: i, k

    ! Solve for forcing in terms of total water mixing ratio.
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzt
      do i = 1, ngrdcol
        rtm_forcing(i,k) = ( 1.0_core_rknd + rtm(i,k) )**2 * qtm_forcing(i,k)
      end do
    end do

    return
  end subroutine force_spec_hum_to_mixing_ratio

!===============================================================================

end module spec_hum_to_mixing_ratio
