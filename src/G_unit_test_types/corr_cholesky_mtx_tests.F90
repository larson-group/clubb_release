! $Id$
module corr_cholesky_mtx_tests

  implicit none

  private

  public :: corr_cholesky_mtx_tests_driver

  private :: corr_cholesky_setup_tests, corr_mtx_approx_tests, &
             print_matrix, percent_difference

  contains

  !=============================================================================
  function corr_cholesky_mtx_tests_driver()

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd

    implicit none

    ! Output Vars
    integer :: corr_cholesky_mtx_tests_driver ! Returns the exit code of the tes

    ! Local Variables
    integer :: &
      total_mismatches   ! Total number of mismatches

    real( kind = core_rknd ) :: &
      tol      ! Acceptable percent difference between the results    [-]

    real( kind = core_rknd ), dimension(3,3) :: &
      corr_cholesky_array_t_1,     &
      corr_cholesky_array_t_2,     &
      corr_cholesky_array_t_3


    ! Declare the tolerance value, which is the maximum acceptable percent
    ! difference between the results.
    tol = 1.0e-6_core_rknd

    ! Initialize total number of mismatches.
    total_mismatches = 0

    print *, "=================================================="
    print *, " "
    print *, "Performing correlation cholesky matrix setup tests"
    print *, " "

    call corr_cholesky_setup_tests( tol, total_mismatches, &
                                    corr_cholesky_array_t_1, &
                                    corr_cholesky_array_t_2, &
                                    corr_cholesky_array_t_3)

    print *, "=================================================="
    print *, " "
    print *, "Performing correlation matrix approximation tests"
    print *, " "

    call corr_mtx_approx_tests( corr_cholesky_array_t_1, corr_cholesky_array_t_2, &
                                            corr_cholesky_array_t_3, tol, &
                                            total_mismatches )

    print *, "=================================================="
    print *, " "

    if ( total_mismatches == 0 ) then

       print *, "Success!"
       corr_cholesky_mtx_tests_driver = 0 ! Exit Code = 0, Success!

    else  ! total_mismatches > 0

       print *, "There were", total_mismatches, "total mismatches found."
       corr_cholesky_mtx_tests_driver = 1 ! Exit Code = 1, Fail

    endif


    return

  end function corr_cholesky_mtx_tests_driver

  !=============================================================================
  subroutine corr_cholesky_setup_tests( tol, total_mismatches, &
                                        corr_cholesky_array_t_1, &
                                        corr_cholesky_array_t_2, &
                                        corr_cholesky_array_t_3)

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one, &  ! Constant(s)
        zero

    use clubb_precision, only: &
        core_rknd

    use diagnose_correlations_module, only: &
        setup_corr_cholesky_mtx

    implicit none

    intrinsic :: reshape, sqrt

    ! Input Variables

    real( kind = core_rknd ), intent(in) :: &
      tol      ! Acceptable percent difference between the results          [-]

    ! Input/Output Variable
    integer, intent(inout) :: &
      total_mismatches   ! Total number of mismatches

    ! Output Variables
    real( kind = core_rknd ), dimension(3,3), intent(out) :: &
      corr_cholesky_array_t_1,     &
      corr_cholesky_array_t_2,     &
      corr_cholesky_array_t_3

    ! Local Variables
    real( kind = core_rknd ), dimension(3,3) :: &
      corr_array_1, &
      corr_array_2, &
      corr_array_3

    real( kind = core_rknd ), dimension(3,3) :: &
      corr_cholesky_array_t_1_cmp, &
      corr_cholesky_array_t_2_cmp, &
      corr_cholesky_array_t_3_cmp

    integer :: &
      count_incr    ! Increment total number of mismatches

    ! ---- Begin Code ----

    ! Initialize the correlation arrays for the input
    corr_array_1 = reshape( (/one, zero, zero, zero, one, zero, zero, zero, one/), (/3, 3/) )
    corr_array_2 = reshape( (/one, one, one, zero, one, one, zero, zero, one/), (/3, 3/) )
    corr_array_3 = reshape( (/one         , sqrt(0.75_core_rknd), sqrt(0.75_core_rknd), &
                              zero         ,                 one, sqrt(0.75_core_rknd), &
                              zero         ,                zero,                one/), &
                            (/3, 3/) )

    ! Initialize the corresponding solutions for the correlation cholesky matrices
    corr_cholesky_array_t_1_cmp = reshape( (/one, zero, zero, zero, one, zero, zero, zero, one/), &
                                           (/3, 3/) )
    corr_cholesky_array_t_2_cmp = reshape( (/one, one, one, zero, zero, zero, zero, zero, zero/), &
                                           (/3, 3/) )
    corr_cholesky_array_t_3_cmp = reshape( &
                                 (/one     , sqrt(0.75_core_rknd), sqrt(0.75_core_rknd), &
                                   zero    ,       0.5_core_rknd, 0.433012702_core_rknd, &
                                   zero    ,             zero,         0.25_core_rknd/), &
                                 (/3, 3/) )

    print *, "correlation matrix :"
    call print_matrix(3, corr_array_1)

    call setup_corr_cholesky_mtx( 3, corr_array_1, &        ! intent(in)
                                  corr_cholesky_array_t_1 ) ! intent(out)

    call percent_difference( 3, corr_cholesky_array_t_1,   &     ! intent(in)
                             corr_cholesky_array_t_1_cmp, tol, & ! intent(in)
                             count_incr )                        ! intent(inout)

    print *, "correlation matrix :"
    call print_matrix(3, corr_array_2)

    call setup_corr_cholesky_mtx( 3, corr_array_2, &        ! intent(in)
                                  corr_cholesky_array_t_2 ) ! intent(out)

    call percent_difference( 3, corr_cholesky_array_t_2,   &     ! intent(in)
                             corr_cholesky_array_t_2_cmp, tol, & ! intent(in)
                             count_incr )                        ! intent(inout)

    print *, "correlation matrix :"
    call print_matrix(3, corr_array_3)

    call setup_corr_cholesky_mtx( 3, corr_array_3, &        ! intent(in)
                                  corr_cholesky_array_t_3 ) ! intent(out)

    call percent_difference( 3, corr_cholesky_array_t_3,   &     ! intent(in)
                             corr_cholesky_array_t_3_cmp, tol, & ! intent(in)
                             count_incr )                        ! intent(inout)

    total_mismatches = total_mismatches + count_incr

    return

  end subroutine corr_cholesky_setup_tests
  !=============================================================================

  !=============================================================================
  subroutine corr_mtx_approx_tests( corr_cholesky_array_t_1, corr_cholesky_array_t_2, &
                                    corr_cholesky_array_t_3, tol, &
                                    total_mismatches )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one, &  ! Constant(s)
        zero

    use clubb_precision, only: &
        core_rknd

    use diagnose_correlations_module, only: &
        cholesky_to_corr_mtx_approx

    implicit none

    intrinsic :: reshape, sqrt

    ! Input Variables
    real( kind = core_rknd ), dimension(3,3), intent(in) :: &
      corr_cholesky_array_t_1,     &
      corr_cholesky_array_t_2,     &
      corr_cholesky_array_t_3

    real( kind = core_rknd ), intent(in) :: &
      tol      ! Acceptable percent difference between the results          [-]

    ! Input/Output Variable
    integer, intent(inout) :: &
      total_mismatches   ! Total number of mismatches

    ! Local Variables
    real( kind = core_rknd ), dimension(3,3) :: &
      corr_array_approx_1, &
      corr_array_approx_2, &
      corr_array_approx_3

    real( kind = core_rknd ), dimension(3,3) :: &
      corr_array_approx_1_cmp, &
      corr_array_approx_2_cmp, &
      corr_array_approx_3_cmp

    integer :: &
      count_incr    ! Increment total number of mismatches

    ! ---- Begin Code ----

    ! Set the corresponding solutions for the approximated correlation matrices
    corr_array_approx_1_cmp = reshape( (/one, zero, zero, zero, one, zero, zero, zero, one/), &
                                       (/3, 3/) )
    corr_array_approx_2_cmp = reshape( (/one, one, one, one, one, one, one, one, one/), &
                                       (/3, 3/) )
    corr_array_approx_3_cmp = reshape( &
                              (/one        ,  sqrt(0.75_core_rknd),  sqrt(0.75_core_rknd), &
                                sqrt(0.75_core_rknd),          one, 0.966506351_core_rknd, &
                                sqrt(0.75_core_rknd), 0.966506351_core_rknd,        one/), &
                              (/3, 3/) )

    print *, "correlation cholesky matrix :"
    call print_matrix(3, corr_cholesky_array_t_1)

    call cholesky_to_corr_mtx_approx( 3, corr_cholesky_array_t_1, & ! intent(in)
                                      corr_array_approx_1 )         ! intent(out)

    call percent_difference( 3, corr_array_approx_1,   &     ! intent(in)
                             corr_array_approx_1_cmp, tol, & ! intent(in)
                             count_incr )                    ! intent(inout)

    print *, "correlation cholesky matrix :"
    call print_matrix(3, corr_cholesky_array_t_2)

    call cholesky_to_corr_mtx_approx( 3, corr_cholesky_array_t_2, & ! intent(in)
                                      corr_array_approx_2 )         ! intent(out)

    call percent_difference( 3, corr_array_approx_2,   &     ! intent(in)
                             corr_array_approx_2_cmp, tol, & ! intent(in)
                             count_incr )                    ! intent(inout)

    print *, "correlation cholesky matrix :"
    call print_matrix(3, corr_cholesky_array_t_3)

    call cholesky_to_corr_mtx_approx( 3, corr_cholesky_array_t_3, & ! intent(in)
                                      corr_array_approx_3 )         ! intent(out)

    call percent_difference( 3, corr_array_approx_3,   &     ! intent(in)
                             corr_array_approx_3_cmp, tol, & ! intent(in)
                             count_incr )                    ! intent(inout)

    total_mismatches = total_mismatches + count_incr

    return

  end subroutine corr_mtx_approx_tests
  !=============================================================================

  !=============================================================================
  subroutine percent_difference( n_var, obtained_result, comparison_result, &
                                 tol, count_incr )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        eps    ! Constant(s)
        

    use clubb_precision, only: &
        core_rknd

    implicit none

    intrinsic :: abs, epsilon, maxval

    ! Input Variables
    integer, intent(in) :: n_var ! number of elements in the correlation matrix

    real( kind = core_rknd ), dimension(n_var, n_var), intent(in) :: &
      obtained_result,   & ! Result obtained for the integral in CLUBB   [-]
      comparison_result    ! Result calculated "by hand"                 [-]

    real( kind = core_rknd ), intent(in) :: &
      tol                  ! Acceptable percent diff. betw. the results  [-]

    ! Output Variable
    integer, intent(out) :: &
      count_incr    ! Increment total number of mismatches.

    ! Local Variable
    real( kind = core_rknd ) :: &
      percent_diff    ! Percent diff. betw. the CLUBB and MATLAB results    [-]

    real( kind = core_rknd ), dimension(n_var, n_var) :: &
      diff_matrix ! Difference matrix

    integer :: i, j ! Loop iterators

    ! ---- Begin Code ----

    ! Percent difference between the result obtained for the integral using
    ! the CLUBB code and the result obtained by hand
    do i = 1,n_var
       do j = 1, n_var
          if ( abs(comparison_result(i, j)) > eps ) then
             diff_matrix(i,j) = abs( (comparison_result(i,j)-obtained_result(i,j)) &
                                     / comparison_result(i,j) )
          else
             diff_matrix(i,j) = abs( (comparison_result(i,j)-obtained_result(i,j)) &
                                     / epsilon(comparison_result(i,j)) )
          endif
       end do
    end do

    percent_diff = maxval(diff_matrix)


    if ( percent_diff <= tol ) then

       ! The percent difference between the obtained result and the MATLAB
       ! result is within an acceptable tolerance.
       print *, "Agreement"
       print *, "Result:  "
       call print_matrix(3, obtained_result)
       print *, "--"
       print *, " "

       count_incr = 0

    else ! percent_diff > tol

       ! The percent difference between the obtained result and the MATLAB
       ! result is beyond an acceptable tolerance.  Print an error message
       ! telling the user to please check for any changes made to the
       ! relevant portion(s) of the CLUBB model code.
       print *, "Mismatch"
       print *, "Obtained result:  "
       call print_matrix(3, obtained_result)
       print *, "Comparison result:  "
       call print_matrix(3, comparison_result)
       print *, "Percent difference:  ", percent_diff*100._core_rknd, "%"
       print *, "--"
       print *, " "

       count_incr = 1

    endif


    return

    end subroutine percent_difference

!===============================================================================

  !=============================================================================
  subroutine print_matrix( n_var, matrix )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd

    implicit none

    ! Input Variables
    integer, intent(in):: n_var

    real( kind = core_rknd ), dimension( n_var, n_var ), intent(in) :: &
      matrix

    ! Local Variables
    integer :: i  ! Loop iterator

    ! ---- Begin Code ----

    do i = 1, n_var
       print *, matrix(:,i)
    end do

    print *, " "

    return

    end subroutine print_matrix

!===============================================================================


end module corr_cholesky_mtx_tests
