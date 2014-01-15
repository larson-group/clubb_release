!$Id$
module hole_filling_tests

  implicit none

  private

  public :: hole_filling_tests_driver

  private :: hole_filling_one_lev_tests

  contains

  !=============================================================================
  subroutine hole_filling_tests_driver

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd

    implicit none

    ! Local Variables
    integer :: total_errors = 0

    real( kind = core_rknd ) :: &
      tol = 1.0e-6_core_rknd ! Acceptable percent difference between the results    [-]

    print *, "=================================================="
    print *, " "
    print *, "Performing hole_filling_one_lev tests"
    print *, " "

    call hole_filling_one_lev_tests( tol, total_errors )

    print *, "=================================================="
    print *, " "

    if ( total_errors == 0 ) then

       print *, "Success!"

    else  ! total_mismatches > 0

       print *, "There were", total_errors, "total errors found."

    endif


    return

  end subroutine hole_filling_tests_driver

  !=============================================================================
  subroutine hole_filling_one_lev_tests( tol, total_errors )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one, &  ! Constant(s)
        zero

    use clubb_precision, only: &
        core_rknd

    use fill_holes, only: &
        hole_filling_one_lev

    implicit none

    ! Input Variables
    real( kind = core_rknd ) :: tol

    ! Input/Output Variable
    integer, intent(inout) :: total_errors

    ! Output Variables

    ! Local Variables
    integer :: num_hm_fill

    real( kind = core_rknd ), dimension(4) :: &
      testset1_in, &
      testset1_out, &
      testset1_comp

    real( kind = core_rknd ), dimension(4) :: &
      testset2_in, &
      testset2_out

    real( kind = core_rknd ), dimension(4) :: &
      testset3_in, &
      testset3_out

    ! ---- Begin Code ----
    ! Testset 1
    num_hm_fill = 4
    testset1_in(1) = -2._core_rknd
    testset1_in(2) = 5._core_rknd
    testset1_in(3) = 6._core_rknd
    testset1_in(4) = -3._core_rknd

    testset1_comp(1) = 0._core_rknd
    testset1_comp(2) = 5._core_rknd-25./11.
    testset1_comp(3) = 6._core_rknd-30./11.
    testset1_comp(4) = 0._core_rknd

    call hole_filling_one_lev( num_hm_fill, testset1_in, testset1_out )

    if( any( abs(testset1_out - testset1_comp) > tol ) ) then
      total_errors = total_errors + 1
      print *, "---- Mismatch: ----"
      print *, "testset1_out = ", testset1_out
      print *, "testset1_comp = ", testset1_comp
    endif

    call check_results( num_hm_fill, testset1_in, testset1_out, total_errors )

    ! Testset 2 -- should cause an error
    num_hm_fill = 4
    testset2_in(1) = 2._core_rknd
    testset2_in(2) = 1._core_rknd
    testset2_in(3) = -5._core_rknd
    testset2_in(4) = 1._core_rknd

    call hole_filling_one_lev( num_hm_fill, testset2_in, testset2_out )

    call check_results( num_hm_fill, testset2_in, testset2_out, total_errors )

    ! Testset 3
    num_hm_fill = 4
    testset3_in(1) = -2._core_rknd
    testset3_in(2) = 6._core_rknd
    testset3_in(3) = -2._core_rknd
    testset3_in(4) = -2._core_rknd

    call hole_filling_one_lev( num_hm_fill, testset3_in, testset3_out )

    call check_results( num_hm_fill, testset3_in, testset3_out, total_errors )

    return

  end subroutine hole_filling_one_lev_tests
  !=============================================================================

  subroutine check_results(num_hm_fill, testset_in, testset_result, total_errors)
    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        zero ! Constant(s)

    use clubb_precision, only: &
        core_rknd

    use fill_holes, only: &
        hole_filling_one_lev

    implicit none

    intrinsic :: sum

    ! Input Variables
    integer, intent(in) :: num_hm_fill

    real( kind = core_rknd ), dimension(num_hm_fill), intent(in) :: &
      testset_result, &
      testset_in

    ! Input/Output Variable
    integer, intent(inout) :: total_errors

    ! ---- Begin Code ----

    ! Check if the hole filling was conservative
    if ( sum(testset_in) /= sum(testset_result) ) then
       total_errors = total_errors + 1
       return
    endif

    ! Check if all holes are filled
    if ( any( testset_result < zero ) ) then
       total_errors = total_errors + 1
       return
    endif

    return

  end subroutine check_results

end module hole_filling_tests
