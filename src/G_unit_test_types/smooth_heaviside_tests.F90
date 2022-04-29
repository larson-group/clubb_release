module smooth_heaviside_tests
  
  implicit none
  
  private
  
  public :: smooth_heaviside_tests_driver
  
  private :: smooth_heaviside_setup_tests
  
contains
  
  !=============================================================================
  function smooth_heaviside_tests_driver()
    
    ! Description:
    !   Returns the driver with exit code 0 if tests were successful or exit 
    !   code 1 if any test failed.
    !   Tests check whether output of smooth_heaviside_peskin within 
    !   advance_helper_module satisfies some cases where a solution has
    !   previously been computed.
    !
    ! References:
    !-----------------------------------------------------------------------

    implicit none

    ! Output var
    integer :: smooth_heaviside_tests_driver ! Returns the output of the test

    ! Local vars
    integer :: total_mismatches

    print *, "=================================================="
    print *, " "
    print *, "Performing 'Smooth Heaviside Peskin' tests"
    print *, " "
    
    call smooth_heaviside_setup_tests(total_mismatches)
    
    if ( total_mismatches == 0 ) then
      print *, "Success!"
      smooth_heaviside_tests_driver = 0 ! Exit Code = 0, Success!
    else  ! total_mismatches > 0
      print *, "The computed result was wrong!"
      smooth_heaviside_tests_driver = 1 ! Exit Code = 1, Fail
    endif

    print *, " "
    print *, "=================================================="
    print *, " "
    
    return
  end function smooth_heaviside_tests_driver
  
  !=============================================================================
  ! I am not completely happy with giving gr as an input argument, but
  ! unfortunately smooth_heaviside_peskin requires the grid as an argument
  !   ~ Jan Gruenenwald
  subroutine smooth_heaviside_setup_tests(total_mismatches)
    
    ! Description:
    !   Tests the function smooth_heaviside_peskin. Test cases are:
    !   Case 1: heaviside_smth_range=zero, output ought to be the same as Heaviside
    !           step function
    !   Case 2: heaviside_smth_range=one (arbitrary value), output values have been
    !           calculated by hand
    !
    ! References:
    !-----------------------------------------------------------------------
        
    use constants_clubb, only: &
      one,      &  ! Constant(s)
      two,      &
      zero,     &
      eps

    use clubb_precision, only: &
      core_rknd
      
    use advance_helper_module, only: &
      smooth_heaviside_peskin

    implicit none
    
    !----------- Output Variables -----------
    integer, intent(out) :: &
      total_mismatches
    
    !----------- Local Variables -----------
    real( kind = core_rknd ), dimension(4, 1) :: &
      result_cmp, result, input
    
    real( kind = core_rknd ) :: &
      heaviside_smth_range
      
    integer :: &
      i
    
    !=========== Code ==================
    
    ! Case 1: Heaviside step function
    print *, "Case 1: heaviside_smth_range = 0, smooth Heaviside collapses to Heaviside"
    heaviside_smth_range = zero
    
    input(1, 1) = -one
    input(2, 1) = -eps
    input(3, 1) = eps
    input(4, 1) = one
    
    result_cmp(1, 1) = zero
    result_cmp(2, 1) = zero
    result_cmp(3, 1) = one
    result_cmp(4, 1) = one
      
    result = smooth_heaviside_peskin(input, heaviside_smth_range)
    print *, "Input: ", input
    print *, "Expected outcome: ", result_cmp
    print *, "True outcome:     ", result
    total_mismatches = total_mismatches + COUNT(abs(result - result_cmp) >= eps)
                          
    ! Case 2: Testing arbitrary smoothing range for precomputed values
    heaviside_smth_range = two
    print *, NEW_LINE('a'), "Case 2: heaviside_smth_range = ", heaviside_smth_range
    
    input(1, 1) = -two - eps
    input(2, 1) = -one
    input(3, 1) = one
    input(4, 1) = two + eps
    
    result_cmp(1, 1) = zero
    result_cmp(2, 1) = 0.0908450569081_core_rknd
    result_cmp(3, 1) = 0.9091549430900_core_rknd
    result_cmp(4, 1) = one
    
    result = smooth_heaviside_peskin(input, heaviside_smth_range)
    print *, "Input: ", input
    print *, "Expected outcome: ", result_cmp
    print *, "True outcome:     ", result
    total_mismatches = total_mismatches + COUNT(abs(result - result_cmp) >= eps)
    
  end subroutine smooth_heaviside_setup_tests
end module smooth_heaviside_tests
