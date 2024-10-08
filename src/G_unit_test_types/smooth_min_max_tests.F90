module smooth_min_max_tests
  
  implicit none
  
  private
  
  public :: smooth_min_max_tests_driver
  
  private :: smooth_min_max_setup_tests
  
contains
  
  !=============================================================================
  function smooth_min_max_tests_driver()
    
    ! Description:
    !   Returns the driver with exit code 0 if tests were successful or exit 
    !   code 1 if any test failed.
    !   Tests check whether output of smooth_min and smooth_max functions 
    !   within advance_helper_module satisfies some cases where a solution 
    !   has previously been computed.
    !
    ! References:
    !-----------------------------------------------------------------------

    implicit none

    ! Output var
    integer :: smooth_min_max_tests_driver ! Returns the output of the test

    ! Local vars
    integer :: total_mismatches

    print *, "=================================================="
    print *, " "
    print *, "Performing 'Smooth Min / Max' tests"
    print *, " "
    
    call smooth_min_max_setup_tests(total_mismatches)
    
    if ( total_mismatches == 0 ) then
      print *, "Success!"
      smooth_min_max_tests_driver = 0 ! Exit Code = 0, Success!
    else  ! total_mismatches > 0
      print *, "The computed result was wrong!"
      smooth_min_max_tests_driver = 1 ! Exit Code = 1, Fail
    endif

    print *, " "
    print *, "=================================================="
    print *, " "
    
    return
  end function smooth_min_max_tests_driver
  
  !=============================================================================
  subroutine smooth_min_max_setup_tests(total_mismatches)
    
    ! Description:
    !   Tests the interfaces smooth_min and smooth_max. Test cases are:
    !   Case 1: min_max_smth_mag=zero. The smooth min and max functions ought to
    !           collapse to their standard min and max variants.
    !   Case 2: min_max_smth_mag=one (arbitrary value), output values have been
    !           calculated by hand.
    !   Case 3: It is checked whether the smooth_min is consistently smaller than the hard min, 
    !           and vice versa for the smooth_max and hard max.
    !   Case 4: We check whether for input_var1 << input_var2, the smooth_min returns results
    !           that lie sufficiently close to input_var1. Similarly for input_var1 >> input_var2
    !           and the smooth_max.
    !   Case 5: Test the scalar variants of smooth_min and smooth_max.
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
      smooth_min, smooth_max

    implicit none
    
    !----------- Output Variables -----------
    integer, intent(out) :: &
      total_mismatches
    
    !----------- Local Variables -----------
    real( kind = core_rknd ), dimension(1, 4) :: &
      result_cmp, result, input
      
    real( kind = core_rknd ), dimension(1, 1000) :: &
      result_cmp_pt2, result_pt2, input_pt2
    
    real( kind = core_rknd ) :: &
      min_max_smth_mag, result_cmp_scalar, result_scalar
      
    integer :: &
      i
    
    total_mismatches = 0
    
    ! Case 1: Simple min / max function
    min_max_smth_mag = zero
    
    input(1, 1) = -one
    input(1, 2) = -eps
    input(1, 3) = eps
    input(1, 4) = one
    
    ! part a) min function 
    result_cmp(1, 1) = -one
    result_cmp(1, 2) = -eps
    result_cmp(1, 3) = zero
    result_cmp(1, 4) = zero
    result = smooth_min(4, 1, input, zero, min_max_smth_mag) ! Order of nz and ngrdcol is 
                                                             ! opposite to order of indexing! 
                                                             ! This is very unintuitive!
    print *, "Input: ", input
    print *, "(Simple) min:"
    print *, "Expected outcome: ", result_cmp
    print *, "True outcome:     ", result, NEW_LINE('A')
    total_mismatches = total_mismatches + COUNT(abs(result - result_cmp) >= eps)

    ! part b) max function
    result_cmp(1, 1) = zero
    result_cmp(1, 2) = zero
    result_cmp(1, 3) = eps
    result_cmp(1, 4) = one
    result = smooth_max(4, 1, input, zero, min_max_smth_mag)
    print *, "(Simple) max:"
    print *, "Expected outcome: ", result_cmp
    print *, "True outcome:     ", result, NEW_LINE('A'), NEW_LINE('A')
    total_mismatches = total_mismatches + COUNT(abs(result - result_cmp) >= eps)
                          
    ! Case 2: Testing arbitrary smoothing range for precomputed values
    min_max_smth_mag = sqrt(two)
    
    input(1, 1) = -two - eps
    input(1, 2) = -one
    input(1, 3) = one
    input(1, 4) = two + eps
    
    ! Since we are taking the max or min with 0, the smooth functions collapse
    ! to the simple Heaviside function and should resemble our previous
    ! results in test case 2 of smooth_heaviside_tests.
    
    ! part a) min function
    result_cmp(1, 1) = -2.22474487148241410_core_rknd
    result_cmp(1, 2) = -1.36602540378443880_core_rknd
    result_cmp(1, 3) = -0.36602540378443871_core_rknd
    result_cmp(1, 4) = -0.22474487138241406_core_rknd
    result = smooth_min(4, 1, input, zero, min_max_smth_mag)
    print *, "Input: ", input
    print *, "Smooth min:"
    print *, "Expected outcome: ", result_cmp
    print *, "True outcome:     ", result
    total_mismatches = total_mismatches + COUNT(abs(result - result_cmp) >= eps)
                          
    ! part b) max function
    result_cmp(1, 1) = 0.22474487138241406_core_rknd
    result_cmp(1, 2) = 0.36602540378443871_core_rknd
    result_cmp(1, 3) = 1.36602540378443880_core_rknd
    result_cmp(1, 4) = 2.22474487148241410_core_rknd
    result = smooth_max(4, 1, input, zero, min_max_smth_mag)
    print *, "Smooth max:"
    print *, "Expected outcome: ", result_cmp
    print *, "True outcome:     ", result, NEW_LINE('A'), NEW_LINE('A')
    total_mismatches = total_mismatches + COUNT(abs(result - result_cmp) >= eps)
    
    ! Case 3: Make sure that on a large number of arbitrarily chosen points with small
    !         min_max_smth_mag, smooth_min < min and smooth_max > max
    print *, "Testing smooth min and max with min_max_smth_mag=1e-7, ", &
             "input_var1 = equidistant grid on [-1, 1), input_var2 = 0", &
             NEW_LINE('A')
    min_max_smth_mag = 1.e-7_core_rknd
    do i=1,1000
      input_pt2(1, i) = -one + i / 500.0_core_rknd  ! input goes from -1 to 1 in even increments
    end do
    result_pt2     = smooth_min(1000, 1, input_pt2, zero, min_max_smth_mag)
    result_cmp_pt2 = smooth_min(1000, 1, input_pt2, zero, zero)
    total_mismatches = total_mismatches + COUNT(result_pt2 > result_cmp_pt2)
    
    result_pt2     = smooth_max(1000, 1, input_pt2, zero, min_max_smth_mag)
    result_cmp_pt2 = smooth_max(1000, 1, input_pt2, zero, zero)
    total_mismatches = total_mismatches + COUNT(result_pt2 < result_cmp_pt2)
    
    ! Case 4: input_var1 << input_var2 or input_var1 >> input_var2
    min_max_smth_mag = one
    
    input(1, 1) = one
    input(1, 2) = 1.0e1_core_rknd * one
    input(1, 3) = 1.0e2_core_rknd * one
    input(1, 4) = 1.0e3_core_rknd * one
    
    ! part a) min function 
    result_cmp(1, 1) = -1.1180339887498949_core_rknd
    result_cmp(1, 2) = -1.0226805085936306_core_rknd
    result_cmp(1, 3) = -1.0024751868658441_core_rknd
    result_cmp(1, 4) = -1.0002497501874359_core_rknd
    result = smooth_min(4, 1, input, -one, min_max_smth_mag)
    print *, "Input: ", input
    print *, "Smooth min:"
    print *, "Expected outcome: ", result_cmp
    print *, "True outcome:     ", result, NEW_LINE('A')
    total_mismatches = total_mismatches + COUNT(abs(result - result_cmp) >= eps)

    ! part b) max function
    result = smooth_max(4, 1, -input, one, min_max_smth_mag)
    print *, "Smooth max:"
    print *, "Expected outcome: ", -result_cmp
    print *, "True outcome:     ", result, NEW_LINE('A'), NEW_LINE('A')
    total_mismatches = total_mismatches + COUNT(abs(result - (-result_cmp)) >= eps)
    
    
    ! Testing smooth_min and smooth_max for large values
    result_cmp(1, 1) = -1.0000000000002501e6_core_rknd
    result_cmp(1, 2) = -1.0000000000000026e7_core_rknd
    result_cmp(1, 3) = -1.0000000000000000e8_core_rknd
    result_cmp(1, 4) = -1.0000000000000000e9_core_rknd
    
    result = smooth_min(4, 1, -1.0e6_core_rknd * input, input, min_max_smth_mag)
    print *, "Input_var1: ", -1.0e6_core_rknd * input
    print *, "Input_var2: ", input
    print *, "Smooth min:"
    print *, "Expected outcome: ", result_cmp
    print *, "True outcome:     ", result, NEW_LINE('A')
    total_mismatches = total_mismatches + COUNT(abs(result - result_cmp) >= eps)
    
    result_cmp(1, 1) = 1.0000002500019036e0_core_rknd
    result_cmp(1, 2) = 1.0000000025145710e1_core_rknd
    result_cmp(1, 3) = 1.0000000000000000e2_core_rknd
    result_cmp(1, 4) = 1.0000000000000000e3_core_rknd
    result = smooth_max(4, 1, -1.0e6_core_rknd * input, input, min_max_smth_mag)
    print *, "Smooth max:"
    print *, "Expected outcome: ", result_cmp
    print *, "True outcome:     ", result, NEW_LINE('A')
    total_mismatches = total_mismatches + COUNT(abs(result - result_cmp) >= eps)
    
    ! Testing scalar versions of smooth_min and smooth_max
    result_cmp_scalar = -1.0000000000002501e6_core_rknd
    result_scalar = smooth_min(-1.0e6_core_rknd, one, min_max_smth_mag)
    print *, "Input_var1: ", -1.0e6_core_rknd
    print *, "Input_var2: ", one
    print *, "Smooth min:"
    print *, "Expected outcome: ", result_cmp_scalar
    print *, "True outcome:     ", result_scalar, NEW_LINE('A')
    if (abs(result_scalar - result_cmp_scalar) >= eps) total_mismatches = total_mismatches + 1
    
    result_cmp_scalar = 1.0000002500019036e0_core_rknd
    result_scalar = smooth_max(-1.0e6_core_rknd, one, min_max_smth_mag)
    print *, "Smooth max:"
    print *, "Expected outcome: ", result_cmp_scalar
    print *, "True outcome:     ", result_scalar, NEW_LINE('A')
    if (abs(result_scalar - result_cmp_scalar) >= eps) total_mismatches = total_mismatches + 1
    
  end subroutine smooth_min_max_setup_tests
end module smooth_min_max_tests
