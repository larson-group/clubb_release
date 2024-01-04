!$Id$
module hole_filling_tests

  implicit none

  private

  public :: hole_filling_tests_driver

  private :: hole_filling_hm_one_lev_tests, fill_holes_hydromet_tests

  contains

  !=============================================================================
  function hole_filling_tests_driver()

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd

    implicit none

    ! Output Vars
    integer :: hole_filling_tests_driver ! Returns the exit code of the test

    ! Local Variables
    integer :: total_errors = 0

    real( kind = core_rknd ) :: &
      tol = 1.0e-10_core_rknd ! Acceptable absolute diff. between results    [-]

    print *, "=================================================="
    print *, " "
    print *, "Performing hole_filling_hm_one_lev_tests"
    print *, " "

    call hole_filling_hm_one_lev_tests( tol, total_errors )

    print *, "=================================================="
    print *, " "

    if ( total_errors == 0 ) then

       print *, "Success!"
       hole_filling_tests_driver = 0 ! Exit Code = 0, Success!

    else  ! total_mismatches > 0

       print *, "There were", total_errors, "total errors found."
       hole_filling_tests_driver = 1 ! Exit Code = 1, Fail

    endif

    total_errors = 0

    print *, "=================================================="
    print *, " "
    print *, "Performing fill_holes_hydromet_tests"
    print *, " "

    call fill_holes_hydromet_tests( tol, total_errors )

    print *, "=================================================="
    print *, " "

    if ( total_errors == 0 ) then

       print *, "Success!"

    else  ! total_mismatches > 0

       print *, "There were", total_errors, "total errors found."
       hole_filling_tests_driver = 1 ! Exit Code = 1, Fail

    endif


    return

  end function hole_filling_tests_driver

  !=============================================================================
  subroutine hole_filling_hm_one_lev_tests( tol, total_errors )

    ! Description:
    ! Tests the subroutine hole_filling_hm_one_lev for conservation and
    ! non-negativity.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        five,  & ! Constant(s)
        three, &
        two,   &
        one,   &
        zero

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    use fill_holes, only: &
        hole_filling_hm_one_lev  ! Procedure(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      tol    ! Acceptable absolute difference between results      [-]

    ! Input/Output Variable
    integer, intent(inout) :: total_errors

    ! Output Variables

    ! Local Variables
    integer :: num_hm_fill

    real( kind = core_rknd ), dimension(4) :: &
      testset1_in,   &
      testset1_out,  &
      testset1_comp

    real( kind = core_rknd ), dimension(4) :: &
      testset2_in,   &
      testset2_out,  &
      testset2_comp

    real( kind = core_rknd ), dimension(4) :: &
      testset3_in,   &
      testset3_out,  &
      testset3_comp

    logical :: &
      l_expect_not_filled  ! Flag for expected the hole not to be filled

    ! ---- Begin Code ----

    ! Testset 1 -- hole filled (total mass > total hole)
    num_hm_fill = 4
    testset1_in(1) = -two
    testset1_in(2) = five
    testset1_in(3) = 6.0_core_rknd
    testset1_in(4) = -three

    testset1_comp(1) = zero
    testset1_comp(2) = five-25.0_core_rknd/11.0_core_rknd
    testset1_comp(3) = 6.0_core_rknd-30.0_core_rknd/11.0_core_rknd
    testset1_comp(4) = zero

    l_expect_not_filled = .false.

    call hole_filling_hm_one_lev( num_hm_fill, testset1_in, testset1_out )

    if( any( abs(testset1_out - testset1_comp) > tol ) ) then
      total_errors = total_errors + 1
      print *, "---- Mismatch: ----"
      print *, "testset1_out = ", testset1_out
      print *, "testset1_comp = ", testset1_comp
    endif

    call check_results_one_lev( num_hm_fill, testset1_in, testset1_out, &
                                tol, l_expect_not_filled, &
                                total_errors )

    ! Testset 2 -- does not fill hole completely (total mass < total hole)
    num_hm_fill = 4
    testset2_in(1) = two
    testset2_in(2) = one
    testset2_in(3) = -five
    testset2_in(4) = one

    testset2_comp(1) = zero
    testset2_comp(2) = zero
    testset2_comp(3) = -one
    testset2_comp(4) = zero

    l_expect_not_filled = .true.

    call hole_filling_hm_one_lev( num_hm_fill, testset2_in, testset2_out )

    if( any( abs(testset2_out - testset2_comp) > tol ) ) then
      total_errors = total_errors + 1
      print *, "---- Mismatch: ----"
      print *, "testset2_out = ", testset2_out
      print *, "testset2_comp = ", testset2_comp
    endif

    call check_results_one_lev( num_hm_fill, testset2_in, testset2_out, &
                                tol, l_expect_not_filled, &
                                total_errors )

    ! Testset 3 -- total hole equals total mass
    num_hm_fill = 4
    testset3_in(1) = -two
    testset3_in(2) = 6.0_core_rknd
    testset3_in(3) = -two
    testset3_in(4) = -two

    testset3_comp(1) = zero
    testset3_comp(2) = zero
    testset3_comp(3) = zero
    testset3_comp(4) = zero

    l_expect_not_filled = .false.

    call hole_filling_hm_one_lev( num_hm_fill, testset3_in, testset3_out )

    if( any( abs(testset3_out - testset3_comp) > tol ) ) then
      total_errors = total_errors + 1
      print *, "---- Mismatch: ----"
      print *, "testset3_out = ", testset3_out
      print *, "testset3_comp = ", testset3_comp
    endif

    call check_results_one_lev( num_hm_fill, testset3_in, testset3_out, &
                                tol, l_expect_not_filled, &
                                total_errors )


    return

  end subroutine hole_filling_hm_one_lev_tests
  !=============================================================================

  subroutine fill_holes_hydromet_tests( tol, total_errors )

    ! Description:
    ! Tests the subroutine fill_holes_hydromet on water substance conservation
    ! and non-negativity.
    !
    ! Expected number of errors: 2
    !
    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd

    use fill_holes, only: &
        fill_holes_hydromet

    implicit none

    intrinsic :: reshape, transpose

    integer, parameter :: c = core_rknd

    integer, parameter :: &
      hm_dim = 7, &
      hgt_dim = 6

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      tol    ! Acceptable absolute difference between results      [-]

    ! Input/Output Variable
    integer, intent(inout) :: total_errors

    ! Output Variables

    ! Local Variables
    integer :: i

    real( kind = core_rknd ), dimension(hgt_dim,hm_dim) :: &
      testset_in,  & ! Testset input
      testset_out    ! Testset after hole filling

    logical, dimension(hgt_dim) :: &
      l_expect_not_filled_lev  ! Flag for expected the hole not to be filled

    logical, dimension(hm_dim) :: &
      l_frozen_hm, & ! if true, then the hydrometeor is frozen; otherwise liquid
      l_mix_rat_hm   ! if true, then the quantity is a hydrometeor mixing ratio

    ! ---- Begin Code ----

    l_frozen_hm  = (/ .true., .false., .false., .true., .true., .true., .true. /)
    l_mix_rat_hm = (/ .false., .true., .false., .true., .true., .true., .true. /)

!    Testset 1
!    F   F   F   T   T   T   T  -->  Frozen mixratio (T) or not (F)
!    1,  6,  3,  1,  1,  2,  1
!    2,  5,  3,  3, -2, -1,  4
!    3,  4,  3,  2, -4, -1,  2
!    4,  3,  4,  1, -3,  1,  3
!    5,  2,  4, -1, -1, -1, -1
!    6,  1,  4,  0,  0,  0,  0

   testset_in = reshape( (/ 1._c,  2._c,  3._c,  4._c,  5._c,  6._c, &
                            6._c,  5._c,  4._c,  3._c,  2._c,  1._c, &
                            3._c,  3._c,  3._c,  4._c,  4._c,  4._c, &
                            1._c,  3._c,  2._c,  1._c, -1._c,  0._c, &
                            1._c, -2._c, -4._c, -3._c, -1._c,  0._c, &
                            2._c, -1._c, -1._c,  1._c, -1._c,  0._c, &
                            1._c,  4._c,  2._c,  3._c, -1._c,  0._c /), &
                            (/ hgt_dim, hm_dim /) )

    l_expect_not_filled_lev &
    = (/ .false., .false., .true., .false., .true., .false. /)

    do i = 1, hm_dim
       print *, "testset_in(:,",i,") = ", testset_in(:,i)
    enddo

    call fill_holes_hydromet( hgt_dim, hm_dim, testset_in, &
                              l_frozen_hm, l_mix_rat_hm, &
                              testset_out )

    call check_results_hm_filling( hgt_dim, hm_dim, & ! Intent(in)
                                   l_frozen_hm, l_mix_rat_hm, & ! Intent(in)
                                   testset_in, testset_out, & ! Intent(in)
                                   tol, l_expect_not_filled_lev, & ! Intent(in)
                                   total_errors ) ! Intent(inout)

    return

  end subroutine fill_holes_hydromet_tests
  !=============================================================================

  subroutine check_results_hm_filling( hgt_dim, hm_dim, & ! Intent(in)
                                       l_frozen_hm, l_mix_rat_hm, & ! Intent(in)
                                       testset_in, testset_result, &
                                       tol, l_expect_not_filled_lev, & ! Intent(in)
                                       total_errors) ! Intent(inout)
    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd

    implicit none

    intrinsic :: sum

    ! Input Variables
    integer, intent(in) :: hgt_dim, hm_dim

    real( kind = core_rknd ), dimension(hgt_dim, hm_dim), intent(in) :: &
      testset_result, &
      testset_in

    real( kind = core_rknd ), intent(in) :: &
      tol    ! Acceptable absolute difference between results      [-]

    logical, dimension(hgt_dim), intent(in) :: &
      l_expect_not_filled_lev  ! Flag for expected the hole not to be filled

    logical, dimension(hm_dim) :: &
      l_frozen_hm, & ! if true, then the hydrometeor is frozen; otherwise liquid
      l_mix_rat_hm   ! if true, then the quantity is a hydrometeor mixing ratio

    ! Input/Output Variables
    integer, intent(inout) :: total_errors

    ! Local Variables
    integer :: i,j ! Loop iterators

    integer :: num_frozen_hm = 0

    real( kind = core_rknd ), dimension( :, :), allocatable :: &
      testset_frozen_res, &
      testset_frozen_in

    ! ---- Begin Code ----

    ! Determine the number of frozen hydrometeor mixing ratios
    do i = 1, hm_dim
       if ( l_frozen_hm(i) .and. l_mix_rat_hm(i) ) then
          num_frozen_hm = num_frozen_hm + 1
       endif
    enddo

    ! Allocation
    allocate( testset_frozen_in(hgt_dim,num_frozen_hm) )
    allocate( testset_frozen_res(hgt_dim,num_frozen_hm) )

    ! Determine frozen hydrometeor mixing ratios
    j = 1
    do i = 1, hm_dim
       if ( l_frozen_hm(i) .and. l_mix_rat_hm(i) ) then
          testset_frozen_res(:,j) = testset_result(:,i)
          testset_frozen_in(:,j) = testset_in(:,i)
          j = j+1
       endif
    enddo

    ! Check if the hole filling was conservative
    do i = 1, hgt_dim
       call check_results_one_lev( num_frozen_hm, testset_frozen_in(i,:), & ! In
                                   testset_frozen_res(i,:), & ! In
                                   tol, l_expect_not_filled_lev(i), & ! In
                                   total_errors ) ! Out
    enddo

    return

  end subroutine check_results_hm_filling

  !=============================================================================
  subroutine check_results_one_lev( num_hm_fill, testset_in, testset_result, & ! Intent(in)
                                    tol, l_expect_not_filled, & ! Intent(in)
                                    total_errors ) ! Intent(inout)
    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        zero,    & ! Constant(s)
        fstdout

    use clubb_precision, only: &
        core_rknd

    implicit none

    intrinsic :: sum

    ! Input Variables
    integer, intent(in) :: num_hm_fill

    real( kind = core_rknd ), dimension(num_hm_fill), intent(in) :: &
      testset_result, &
      testset_in

    real( kind = core_rknd ), intent(in) :: &
      tol    ! Acceptable absolute difference between results      [-]

    logical, intent(in) :: &
      l_expect_not_filled  ! Flag for expected the hole not to be filled

    ! Input/Output Variable
    integer, intent(inout) :: total_errors

    ! ---- Begin Code ----

    ! Check if the hole filling was conservative
    if ( abs( sum( testset_in ) - sum( testset_result ) ) <= tol ) then
       write (fstdout,*) "Hole-filling was conservative"
    else ! abs( sum( testset_in ) - sum( testset_result ) ) > tol
       write (fstdout,*) "Hole-filling was not conservative"
       total_errors = total_errors + 1
    endif

    ! Check if all holes are filled
    if ( any( testset_result < zero ) ) then
       if ( l_expect_not_filled ) then
          write (fstdout,*) "Hole not entirely filled, as expected"
       else ! not l_expect_not_filled
          write (fstdout,*) "Error: hole not entirely filled"
          total_errors = total_errors + 1
       endif ! l_expect_not_filled
    endif ! any( testset_result ) < 0


    return

  end subroutine check_results_one_lev

end module hole_filling_tests
