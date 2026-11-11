module cubic_solve_tests

  implicit none

  private

  public :: cubic_solve_tests_driver

  private :: cubic_solve_setup_tests

contains

  !=============================================================================
  function cubic_solve_tests_driver()

    ! Description:
    !   Returns the driver with exit code 0 if tests were successful or exit
    !   code 1 if any test failed.
    !   Tests check whether the roots returned by cubic_solve within calc_roots
    !   match the known roots of cubic equations that were built from them.
    !
    ! References:
    !-----------------------------------------------------------------------

    implicit none

    ! Output var
    integer :: cubic_solve_tests_driver ! Returns the output of the test

    ! Local vars
    integer :: total_mismatches

    print *, "=================================================="
    print *, " "
    print *, "Performing 'Cubic Solve' tests"
    print *, " "

    total_mismatches = 0

    call cubic_solve_setup_tests(total_mismatches)

    if ( total_mismatches == 0 ) then
      print *, "Success!"
      cubic_solve_tests_driver = 0 ! Exit Code = 0, Success!
    else  ! total_mismatches > 0
      print *, "The computed result was wrong!"
      cubic_solve_tests_driver = 1 ! Exit Code = 1, Fail
    endif

    print *, " "
    print *, "=================================================="
    print *, " "

    return
  end function cubic_solve_tests_driver

  !=============================================================================
  subroutine cubic_solve_setup_tests(total_mismatches)

    ! Description:
    !   Tests the function cubic_solve.  Each column is a cubic equation that
    !   was built by multiplying out known roots, and all columns are solved
    !   in a single call.  Test cases are:
    !   Case 1: (x-1)(x^2+4)       = x^3 - x^2 + 4x - 4,  D > 0, R > 0,
    !           R - sqrt(D) < 0;  roots 1, +2i, -2i
    !   Case 2: (x+2)(x^2-2x+2)    = x^3 - 2x + 4,        D > 0, R < 0;
    !           roots -2, 1+i, 1-i
    !   Case 3: (x-1)^2(x+2)       = x^3 - 3x + 2,        D = 0, R < 0;
    !           roots -2, 1, 1
    !   Case 4: (x-1)(x-2)(x-3)    = x^3 - 6x^2 + 11x - 6, D < 0, R = 0;
    !           roots 1, 2, 3
    !   Case 5: x(x-4)(x-5)        = x^3 - 9x^2 + 20x,    D < 0, R < 0;
    !           roots 0, 4, 5
    !   Cases 1-3 have a negative real Cardano argument R +/- sqrt(D), where
    !   the principal complex cube root is not the real cube root.
    !
    !   When D > 0, roots(1) must be the real root and roots(2) the root with a
    !   positive imaginary part, since adg1_adg2_3d_luhar_pdf relies on
    !   roots(1) being real.  When D <= 0, all three roots must be real; their
    !   order is not checked.
    !
    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
      one,      &  ! Constant(s)
      two,      &
      zero

    use clubb_precision, only: &
      core_rknd

    use calc_roots, only: &
      cubic_solve

    implicit none

    !----------- Input/Output Variables -----------
    integer, intent(inout) :: &
      total_mismatches

    !----------- Local Constants -----------
    integer, parameter :: &
      nz = 5    ! Number of test cases (columns)

    !----------- Local Variables -----------
    real( kind = core_rknd ), dimension(nz) :: &
      a_coef, b_coef, c_coef, d_coef

    complex( kind = core_rknd ), dimension(nz,3) :: &
      roots, roots_cmp

    real( kind = core_rknd ), dimension(3) :: &
      roots_sorted, roots_sorted_cmp

    real( kind = core_rknd ) :: &
      tol    ! Tolerance for a root to match the expected root

    integer :: k

    !=========== Code ==================

    tol = 1000.0_core_rknd * epsilon( one )

    ! Case 1: x^3 - x^2 + 4x - 4
    a_coef(1) = one
    b_coef(1) = -one
    c_coef(1) = 4.0_core_rknd
    d_coef(1) = -4.0_core_rknd
    roots_cmp(1,1) = cmplx( one, zero, kind = core_rknd )
    roots_cmp(1,2) = cmplx( zero, two, kind = core_rknd )
    roots_cmp(1,3) = cmplx( zero, -two, kind = core_rknd )

    ! Case 2: x^3 - 2x + 4
    a_coef(2) = one
    b_coef(2) = zero
    c_coef(2) = -two
    d_coef(2) = 4.0_core_rknd
    roots_cmp(2,1) = cmplx( -two, zero, kind = core_rknd )
    roots_cmp(2,2) = cmplx( one, one, kind = core_rknd )
    roots_cmp(2,3) = cmplx( one, -one, kind = core_rknd )

    ! Case 3: x^3 - 3x + 2
    a_coef(3) = one
    b_coef(3) = zero
    c_coef(3) = -3.0_core_rknd
    d_coef(3) = two
    roots_cmp(3,1) = cmplx( -two, zero, kind = core_rknd )
    roots_cmp(3,2) = cmplx( one, zero, kind = core_rknd )
    roots_cmp(3,3) = cmplx( one, zero, kind = core_rknd )

    ! Case 4: x^3 - 6x^2 + 11x - 6
    a_coef(4) = one
    b_coef(4) = -6.0_core_rknd
    c_coef(4) = 11.0_core_rknd
    d_coef(4) = -6.0_core_rknd
    roots_cmp(4,1) = cmplx( one, zero, kind = core_rknd )
    roots_cmp(4,2) = cmplx( two, zero, kind = core_rknd )
    roots_cmp(4,3) = cmplx( 3.0_core_rknd, zero, kind = core_rknd )

    ! Case 5: x^3 - 9x^2 + 20x
    a_coef(5) = one
    b_coef(5) = -9.0_core_rknd
    c_coef(5) = 20.0_core_rknd
    d_coef(5) = zero
    roots_cmp(5,1) = cmplx( zero, zero, kind = core_rknd )
    roots_cmp(5,2) = cmplx( 4.0_core_rknd, zero, kind = core_rknd )
    roots_cmp(5,3) = cmplx( 5.0_core_rknd, zero, kind = core_rknd )

    roots = cubic_solve( nz, a_coef, b_coef, c_coef, d_coef )

    do k = 1, nz

      print *, NEW_LINE('a'), "Case ", k, ": a, b, c, d = ", &
               a_coef(k), b_coef(k), c_coef(k), d_coef(k)
      print *, "Expected roots: ", roots_cmp(k,:)
      print *, "True roots:     ", roots(k,:)

      if ( any( aimag( roots_cmp(k,:) ) /= zero ) ) then

        ! D > 0: the order of the roots is fixed.
        total_mismatches = total_mismatches &
                           + count( abs( roots(k,:) - roots_cmp(k,:) ) &
                                    >= tol * max( one, abs( roots_cmp(k,:) ) ) )

      else

        ! D <= 0: all roots are real, compare them independent of order.
        roots_sorted = sort_three( real( roots(k,:), kind = core_rknd ) )
        roots_sorted_cmp = sort_three( real( roots_cmp(k,:), kind = core_rknd ) )

        total_mismatches = total_mismatches &
                           + count( abs( roots_sorted - roots_sorted_cmp ) &
                                    >= tol * max( one, abs( roots_sorted_cmp ) ) ) &
                           + count( abs( aimag( roots(k,:) ) ) &
                                    >= tol * max( one, maxval( abs( roots_sorted_cmp ) ) ) )

      endif

    enddo

  end subroutine cubic_solve_setup_tests

  !=============================================================================
  pure function sort_three( x ) result( x_sorted )

    ! Description:
    !   Sorts three values from smallest to largest.
    !
    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd

    implicit none

    real( kind = core_rknd ), dimension(3), intent(in) :: &
      x

    real( kind = core_rknd ), dimension(3) :: &
      x_sorted

    x_sorted(1) = minval( x )
    x_sorted(3) = maxval( x )
    x_sorted(2) = sum( x ) - x_sorted(1) - x_sorted(3)

  end function sort_three

end module cubic_solve_tests
