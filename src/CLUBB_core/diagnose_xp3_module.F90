!-------------------------------------------------------------------------
!$Id:$
!===============================================================================
module diagnose_xp3_module

  implicit none

  public :: diagnose_xp3

  private

  contains

!-------------------------------------------------------------------------------
  subroutine diagnose_xp3( xm, xp2, wp2, wpxp, tau_zt, small_m_w, big_m_x, x_tol, &
                           xp3 )

! Description:
!   Diagnose the third order moment of x

! References:
!   None
!-------------------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    use grid_class, only: &
        gr,    & ! Variable(s)
        zm2zt    ! Procedure(s)

    use constants_clubb, only: &
        w_tol_sqd,   & ! Constant(s
        one_half,    &
        one,         &
        eps

    implicit none

    ! External
    intrinsic :: min, max

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      xm,          &
      xp2,         &
      wp2,         &
      wpxp,        &
      tau_zt,      &
      small_m_w,   &
      big_m_x

    real( kind = core_rknd), intent(in) :: &
      x_tol

    ! Output Variable
    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      xp3

    ! Local Variables
    integer :: &
      k, &
      kp1

    real( kind = core_rknd), dimension(gr%nz) :: &
      dxm_dz,  &
      dxp2_dz, &
      wp2_over_xp2

    ! ---- Begin Code ----

      do k = 1, gr%nz, 1

        kp1 = min( k+1, gr%nz )

        ! Compute the vertical derivative of xm and xp2
        dxm_dz(k) = gr%invrs_dzm(k) * ( xm(kp1)  - xm(k) )
        dxp2_dz(k) = gr%invrs_dzm(k) * ( xp2(kp1)  - xp2(k) )
        wp2_over_xp2(k) =  wp2(k) / max( xp2(k), x_tol**2 )

        xp3(k) = one

      end do

  end subroutine diagnose_xp3
!-----------------------------------------------------------------------

end module diagnose_xp3_module
