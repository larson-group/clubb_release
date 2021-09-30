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
    integer :: total_mismatches = 0
    real(kind = core_rknd) :: rel_tol, abs_tol
    
    rel_tol = 1.0e-6_core_rknd ! Acceptable percent difference between results
    abs_tol = 1.0e-10_core_rknd ! Acceptable abs difference between results

    print *, "=================================================="
    print *, " "
    print *, "Performing 'mean vertical velocity in clouds' tests"
    print *, " "
    
    call w_up_in_cloud_setup_tests()
    
    print *, "=================================================="
    print *, " "
    
    if ( total_mismatches == 0 ) then
      print *, "Success!"
      w_up_in_cloud_tests_driver = 0 ! Exit Code = 0, Success!
    else  ! total_mismatches > 0
      print *, "There were", total_mismatches, "total mismatches found."
      w_up_in_cloud_tests_driver = 1 ! Exit Code = 1, Fail

endif


return


  end function w_up_in_cloud_tests_driver


  !=============================================================================
  subroutine w_up_in_cloud_setup_tests( &
                            abs_tol, rel_tol, total_mismatches, &
                            results &
                            )
    
    ! Description:
    !
    ! References:
    !-----------------------------------------------------------------------
    use grid_class, only: &
        grid ! Type
        
    use constants_clubb, only: &
      one, &  ! Constant(s)
      zero

    use clubb_precision, only: &
      core_rknd

    implicit none
    
    ! Output
    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      results ! mean updraft over clouds                               [m/s]
      
      
    !----------- Local Variables -----------
    type (grid), target, intent(in) :: gr
    
    real(kind = core_rknd), dimension(gr%nz) &
      mixt_frac, &      ! mixture fraction                             [-]
      cloud_frac_1, &   ! cloud fraction (1st PDF component)           [-]
      cloud_frac_2, &   ! cloud fraction (2nd PDF component)           [-]
      w_1, &            ! upward velocity (1st PDF component)          [m/s]
      w_2, &            ! upward velocity (2nd PDF component)          [m/s]
      varnce_w_1, &     ! standard deviation of w (1st PDF component)  [m^2/s^2]
      varnce_w_2        ! standard deviation of w (2nd PDF component)  [m^2/s^2]

  end subroutine w_up_in_cloud_setup_tests


end module w_up_in_cloud_tests
