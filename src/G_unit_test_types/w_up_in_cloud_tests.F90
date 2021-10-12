module w_up_in_cloud_tests

  implicit none

  private

  public :: w_up_in_cloud_tests_driver

  private :: w_up_in_cloud_setup_tests

contains

  !=============================================================================
  function w_up_in_cloud_tests_driver(gr)

    ! Description:
    !   Returns the driver with exit code 0 if tests were successful or exit 
    !   code 1 if any test failed.
    !   Tests check whether computation of w_up_in_cloud within pdf_closure_module
    !   satisfy some given boundary cases where an exact solution is known.
    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd
        
    use grid_class, only: &
        grid ! Type

    implicit none
    
    ! Input var
    type (grid), target, intent(in) :: gr

    ! Output var
    integer :: w_up_in_cloud_tests_driver ! Returns the output of the test

    ! Local vars
    integer :: total_mismatches

    print *, "=================================================="
    print *, " "
    print *, "Performing 'mean vertical velocity in clouds' tests"
    print *, " "
    
    call w_up_in_cloud_setup_tests(gr, total_mismatches)
    
    print *, "=================================================="
    print *, " "
    
    if ( total_mismatches .eq. 0 ) then
      print *, "Success!"
      w_up_in_cloud_tests_driver = 0 ! Exit Code = 0, Success!
    else  ! total_mismatches > 0
      print *, "The computed result was wrong!"
      w_up_in_cloud_tests_driver = 1 ! Exit Code = 1, Fail
    endif


return


  end function w_up_in_cloud_tests_driver


  !=============================================================================
  subroutine w_up_in_cloud_setup_tests(gr, total_mismatches)
    
    ! Description:
    !   Tests the subroutine calc_w_up_in_cloud. Test cases for both w_1 and w_2 are:
    !   Case 1: Cloud_frac = 0
    !     if Cloud_frac = 0, then we would divide by 0 in the heaviside function.
    !     Obviously this makes little sense and what we would expect as output
    !     is 0 since if there are no clouds there will be no updraft inside clouds.
    !   
    !   Case 2: sigma = 0
    !     if the variation in a component is zero, then the vertical velocity
    !     is constant in that area. Inside the heaviside function, there are two
    !     terms where we divide by the standard deviation. If it is 0, those
    !     will go to infinity (which is not a problem in exact arithmetic since 
    !     erf(infty) = 1 and exp(-infty) = 0). Thus, for sigma = 0 we expect 
    !     simply max(0, w) to be output from the product of w and the heaviside function
    !
    !   Case 3: w strictly non-Positive WITH small sigma
    !     if w is small enough so that considering the standard deviation, it
    !     will never go above 0, then we expect the average updraft to be 0.
    !
    !  Test cases for testing w_1 and w_2 simultaneously are:
    !  Case 1: C_1=C_2=1, a=0.5, sigma_w_1=sigma_w_2=0
    !     Output is the average of w_1 and w_2
    ! References:
    !-----------------------------------------------------------------------
    use grid_class, only: &
        grid ! Type
        
    use constants_clubb, only: &
      one,      &  ! Constant(s)
      two,      &
      three,    &
      zero,     &
      one_half, &
      eps

    use clubb_precision, only: &
      core_rknd
      
    use pdf_closure_module, only: &
      calc_w_up_in_cloud

    implicit none
    
    !----------- Input Variables -----------
    type (grid), target, intent(in) :: gr
    
    !----------- Local Variables -----------
    real( kind = core_rknd ), dimension(gr%nz) :: result
    
    !----------- Output Variables -----------
    integer, intent(out) :: total_mismatches
    
    real(kind = core_rknd), dimension(gr%nz) :: &
      mixt_frac, &      ! mixture fraction                             [-]
      cloud_frac_1, &   ! cloud fraction (1st PDF component)           [-]
      cloud_frac_2, &   ! cloud fraction (2nd PDF component)           [-]
      w_1, &            ! upward velocity (1st PDF component)          [m/s]
      w_2, &            ! upward velocity (2nd PDF component)          [m/s]
      varnce_w_1, &     ! standard deviation of w (1st PDF component)  [m^2/s^2]
      varnce_w_2, &     ! standard deviation of w (2nd PDF component)  [m^2/s^2]
      result_cmp        ! correct result for w_up_in_cloud

      total_mismatches = 0

      ! For w_1:
        ! Case 1
      mixt_frac = one
      cloud_frac_1 = zero
      cloud_frac_2 = zero
      w_1 = one
      varnce_w_1 = one
      w_2 = zero
      varnce_w_2 = zero
      result_cmp = zero
      call calc_w_up_in_cloud(gr, mixt_frac, &
                              cloud_frac_1, cloud_frac_2, &
                              w_1, w_2, &
                              varnce_w_1, varnce_w_2, &
                              result)
      total_mismatches = total_mismatches + COUNT(abs(result - result_cmp) .ge. eps)
      
        ! Case 2
      cloud_frac_1 = one
      varnce_w_1 = zero
      result_cmp = one
      call calc_w_up_in_cloud(gr, mixt_frac, &
                              cloud_frac_1, cloud_frac_2, &
                              w_1, w_2, &
                              varnce_w_1, varnce_w_2, &
                              result)
      total_mismatches = total_mismatches + COUNT(abs(result - result_cmp) .ge. eps)
      
        ! Case 3
      w_1 = -one
      result_cmp = zero
      call calc_w_up_in_cloud(gr, mixt_frac, &
                              cloud_frac_1, cloud_frac_2, &
                              w_1, w_2, &
                              varnce_w_1, varnce_w_2, &
                              result)
      total_mismatches = total_mismatches + COUNT(abs(result - result_cmp) .ge. eps)
      
      
      ! For w_2:
        ! Case 1
      mixt_frac = zero
      cloud_frac_1 = zero
      w_1 = zero
      varnce_w_1 = zero
      cloud_frac_2 = zero
      w_2 = one
      varnce_w_2 = one
      result_cmp = zero
      call calc_w_up_in_cloud(gr, mixt_frac, &
                              cloud_frac_1, cloud_frac_2, &
                              w_1, w_2, &
                              varnce_w_1, varnce_w_2, &
                              result)
      total_mismatches = total_mismatches + COUNT(abs(result - result_cmp) .ge. eps)
      
        ! Case 2
      cloud_frac_2 = one
      varnce_w_2 = zero
      result_cmp = one
      call calc_w_up_in_cloud(gr, mixt_frac, &
                              cloud_frac_1, cloud_frac_2, &
                              w_1, w_2, &
                              varnce_w_1, varnce_w_2, &
                              result)
      total_mismatches = total_mismatches + COUNT(abs(result - result_cmp) .ge. eps)
      
        ! Case 3
      w_2 = -one
      result_cmp = zero
      call calc_w_up_in_cloud(gr, mixt_frac, &
                              cloud_frac_1, cloud_frac_2, &
                              w_1, w_2, &
                              varnce_w_1, varnce_w_2, &
                              result)
      total_mismatches = total_mismatches + COUNT(abs(result - result_cmp) .ge. eps)
      
      ! Test w_1 and w_2 together
        ! Case 1
      mixt_frac = one_half
      cloud_frac_1 = one
      w_1 = one
      varnce_w_1 = zero
      cloud_frac_2 = one
      w_2 = three
      varnce_w_2 = zero
      result_cmp = two
      call calc_w_up_in_cloud(gr, mixt_frac, &
                              cloud_frac_1, cloud_frac_2, &
                              w_1, w_2, &
                              varnce_w_1, varnce_w_2, &
                              result)
      total_mismatches = total_mismatches + COUNT(abs(result - result_cmp) .ge. eps)
      

  end subroutine w_up_in_cloud_setup_tests


end module w_up_in_cloud_tests
