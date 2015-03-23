!-----------------------------------------------------------------------
! $Id$
!-----------------------------------------------------------------------
module silhs_category_test

  implicit none

  private

  public :: silhs_category_test_driver

  contains

  !=============================================================================
  function silhs_category_test_driver()

    ! Description:
    !   This function tests the SILHS category picker algorithm.

    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd

    use silhs_importance_sample_module, only: &
      generate_strat_uniform_variate, & ! Procedure(s)
      pick_sample_categories, &
      num_importance_categories ! Constant

    use constants_clubb, only: &
      fstderr

    implicit none

    ! Local Constants
    integer, parameter :: &
      num_samples = 100       ! Number of samples for this test

    ! Output Variable
    integer :: silhs_category_test_driver ! Returns the exit code of the test

    ! Local Variables
    integer, dimension(num_importance_categories) :: &
      n_sample_points_per_category

    integer, dimension(num_samples) :: &
      int_sample_category

    real( kind = core_rknd ), dimension(num_samples) :: &
      rand_vect

    real( kind = core_rknd ), dimension(num_importance_categories) :: &
      category_prescribed_probs

    integer :: isample, icategory

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    silhs_category_test_driver = 0

    print *, "=================================================="
    print *, ""
    print *, "Performing SILHS category test"
    print *, ""
    print *, "=================================================="

    ! Assign desired numbers for each category. Note: these should add up to
    ! num_sample_points, otherwise the test will fail!!
    n_sample_points_per_category(1) = 13
    n_sample_points_per_category(2) = 12
    n_sample_points_per_category(3) = 26
    n_sample_points_per_category(4) = 31
    n_sample_points_per_category(5) = 5
    n_sample_points_per_category(6) = 8
    n_sample_points_per_category(7) = 3
    n_sample_points_per_category(8) = 2

    category_prescribed_probs(:) = real( n_sample_points_per_category(:), kind=core_rknd ) / &
                                   real( num_samples, kind=core_rknd )

    rand_vect = generate_strat_uniform_variate( num_samples )

    int_sample_category = pick_sample_categories( num_samples, category_prescribed_probs, &
                                                  rand_vect )

    ! Here's a sly trick to check that the generated sample points per category equals
    ! what we prescribed above.
    do isample=1, num_samples
      n_sample_points_per_category(int_sample_category(isample)) = &
        n_sample_points_per_category(int_sample_category(isample)) - 1
    end do ! isample=1, num_samples

    ! They should all be zero.
    do icategory=1, num_importance_categories
      if ( n_sample_points_per_category(icategory) /= 0 ) then
        silhs_category_test_driver = 1 ! Failure code
      end if
    end do

    if ( silhs_category_test_driver /= 0 ) then
      write(fstderr,*) 'Error: the picked number of samples in a category was not equal to &
                       &what was prescribed.'
      write(fstderr,*) 'n_sample_points_per_category(:) = ', n_sample_points_per_category(:)

    else ! silhs_category_test_driver == 1

      print *, 'Success!'

    end if

    return
  end function silhs_category_test_driver

end module silhs_category_test
