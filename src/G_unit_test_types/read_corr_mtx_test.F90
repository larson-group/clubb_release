! $Id$
module read_corr_mtx_test
  implicit none

  private

  public :: read_corr_mtx_unit_test

  contains

  !-----------------------------------------------------------------------
  subroutine read_corr_mtx_unit_test(show_corr_arrays_input)
  ! Description:
  ! Tests the read_correlation_matrix subroutine in corr_matrix_module.F90
  !-----------------------------------------------------------------------

  ! Included Modules
  use corr_matrix_module, &
    only: read_correlation_matrix, & !Subroutine(s)
      setup_pdf_indices, &
      print_corr_matrix

  use clubb_precision, only: core_rknd

  use matrix_operations, only: print_lower_triangular_matrix ! Procedure

  use constants_clubb, only: &
        one,  & ! Constant(s)
        zero

  implicit none

  ! Input var
  logical, optional :: &
    show_corr_arrays_input ! If true, the correlation arrays will be printed when the test is ran

  ! Local Constants
  ! Passed in to read_correlation_matrix
  integer, parameter :: &
    iunit = 5, &    !input unit
    d_variables = 12  !number of variables

  character(LEN=*), parameter :: &
    input_file = "../input_misc/corr_array.in" ! Path to the file

  ! tolerance used for real precision testing
  real( kind = core_rknd ), parameter :: tol = 1.0e-6_core_rknd

  ! This "renaming" is used to shorten the matrix declarations below.
  integer, parameter :: c = core_rknd

  ! Local Variables
  real( kind = core_rknd ), dimension(:,:), allocatable :: &
      corr_array, & ! Retrieved correlation variance array
      test_corr_array! Expected correlation variance array

  integer:: n_row, n_col, & !indeces
    errors !Number of errors found between the two arrays

  logical :: show_corr_arrays ! the member variable of show_corr_arrays_input

  !-----------------------------------------------------------------------
    !----- Begin Code -----
    print *, "Performing corr_mtx_test"
    print *, ""
    print *, "=================================================="

    ! Handle the optional argument
    if (present(show_corr_arrays_input)) then
      show_corr_arrays = show_corr_arrays_input
    else
      show_corr_arrays = .false.
    end if

    ! Allocate Arrays
    allocate( test_corr_array(d_variables,d_variables) )
    allocate( corr_array(d_variables,d_variables) )

    ! Initialize all values to 0.
    test_corr_array = zero
    corr_array = zero

    !setup the test array
    test_corr_array = reshape( (/ &
      1.0_c, 2.0_c, 3.0_c, 4.0_c,  5.0_c,  6.0_c,  &
      7.0_c,     8.0_c,     9.0_c,      1.0_c,     1.1_c,       1.2_c,     &! s

      0.0_c, 1.0_c, 0.2_c, 0.3_c,  0.4_c,  0.5_c,  &
      0.6_c,     0.7_c,     0.8_c,      0.9_c,     0.10_c,      0.11_c,    &! t

      0.0_c, 0.0_c, 1.0_c, 0.02_c, 0.03_c, 0.04_c, &
      0.05_c,    0.06_c,    0.07_c,     0.08_c,    0.09_c,      0.010_c,   &! w

      0.0_c, 0.0_c, 0.0_c, 1.0_c,  -1.0_c, -2.0_c, &
      -3.0_c,    -4.0_c,    -5.0_c,     -6.0_c,    -7.0_c,      -8.0_c,    &! Ncn

      0.0_c, 0.0_c, 0.0_c, 0.0_c,  1.0_c,  -0.1_c, &
      -0.2_c,    -0.3_c,    -0.4_c,     -0.5_c,    -0.6_c,      -0.7_c,    &! rrain

      0.0_c, 0.0_c, 0.0_c, 0.0_c,  0.0_c,  1.0_c,  &
      0.00001_c, 0.00002_c, 0.000003_c, 0.00004_c, 0.0000005_c, 0.000006_c,&! Nr

      0.0_c, 0.0_c, 0.0_c, 0.0_c,  0.0_c,  0.0_c,  &
      1.0_c,     -0.0001_c, -0.00002_c, -0.0003_c, -0.000004_c, -0.00005_c,&! rice

      0.0_c, 0.0_c, 0.0_c, 0.0_c,  0.0_c,  0.0_c,  &
      0.0_c,     1.0_c,     0.0_c,      0.0_c,     0.0_c,       1.0_c,     &! Ni

      0.0_c, 0.0_c, 0.0_c, 0.0_c,  0.0_c,  0.0_c,  &
      0.0_c,     0.0_c,     1.0_c,      0.0_c,     0.0_c,       2.0_c,     &! rsnow

      0.0_c, 0.0_c, 0.0_c, 0.0_c,  0.0_c,  0.0_c,  &
      0.0_c,     0.0_c,     0.0_c,      1.0_c,     0.0_c,       3.0_c,     &! Nsnow

      0.0_c, 0.0_c, 0.0_c, 0.0_c,  0.0_c,  0.0_c,  &
      0.0_c,     0.0_c,     0.0_c,      0.0_c,     1.0_c,       4.0_c,     &! rqraupel

      0.0_c, 0.0_c, 0.0_c, 0.0_c,  0.0_c,  0.0_c,  &
      0.0_c,     0.0_c,     0.0_c,      0.0_c,     0.0_c,       1.0_c      &! Ng
    /), shape(test_corr_array) )

    ! setup the system under test
    call setup_pdf_indices(12, 5, 6, &
                           7, 8, 9, 10, &
                           11, 12)

    ! Read the correlation aray in the given test file
    call read_correlation_matrix(iunit, input_file, &
      d_variables, corr_array)

    ! Print out the array if required
    if (show_corr_arrays) then
      print *, "Expected:"
      call print_corr_matrix(d_variables, test_corr_array)
      print *, "Actual:"
      call print_corr_matrix(d_variables, corr_array)
    end if

    errors = 0

    ! Check for differences (between the original array and the one
    ! read in) greater than the acceptable tolerance.
    do n_col = 1, d_variables
      do n_row = 1, d_variables
        if ( abs(test_corr_array(n_row, n_col)-corr_array(n_row, n_col)) > tol ) then
              errors = errors + 1
        end if
      end do
    end do

    if (errors > 0) then
      print *, "There were ", errors, " errors found"
    else
      print *, "Success!"
    end if

    print *, "=================================================="

    return
  end subroutine read_corr_mtx_unit_test
  !-----------------------------------------------------------------------

end module read_corr_mtx_test