module w_up_in_cloud_tests

  implicit none

  private

  public :: w_up_in_cloud_tests_driver

  private :: w_up_in_cloud_setup_tests

contains

  !=============================================================================
  function w_up_in_cloud_tests_driver()

    ! Description:
    !
    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd

    implicit none

    ! Output var
    integer :: w_up_in_cloud_tests_driver

    ! Local vars
    real(kind = core_rknd) :: tol = 1.0e-6_core_rknd

    print *, "=================================================="
    print *, " "
    print *, "Performing 'mean vertical velocity in clouds' tests"
    print *, " "


  end function w_up_in_cloud_tests_driver


  !=============================================================================
  subroutine w_up_in_cloud_setup_tests()
    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one, &  ! Constant(s)
        zero

    use clubb_precision, only: &
        core_rknd

    implicit none    

  end subroutine w_up_in_cloud_setup_tests


end module w_up_in_cloud_tests
