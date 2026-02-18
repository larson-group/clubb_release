! $Id$
module read_corr_mtx_test
  implicit none

  private

  public :: read_corr_mtx_unit_test

  contains

  !-----------------------------------------------------------------------
  function read_corr_mtx_unit_test(show_corr_arrays_input)
  ! Description:
  ! Tests the read_correlation_matrix subroutine in corr_varnce_module.F90
  ! (current implementation is in corr_varnce_input_reader.F90).
  ! Tests the read_correlation_matrix subroutine in corr_varnce_input_reader.F90
  !-----------------------------------------------------------------------

  ! Included Modules
  use corr_varnce_module, &
    only: init_pdf_hydromet_arrays_api, &
      print_corr_matrix, &
      hm_metadata_type

  use corr_varnce_input_reader, only: &
    read_correlation_matrix

  use clubb_precision, only: core_rknd

  use matrix_operations, only: print_lower_triangular_matrix ! Procedure

  use constants_clubb, only: &
        zero, &
        one

  implicit none

  ! Local Constants
  ! Passed in to read_correlation_matrix
  integer, parameter :: &
    iunit = 5  ! input unit

  character(LEN=*), parameter :: &
    input_file = "../input_misc/corr_array.in" ! Path to the file

  ! tolerance used for real precision testing
  real( kind = core_rknd ), parameter :: tol = 1.0e-6_core_rknd

  ! This "renaming" is used to shorten the matrix declarations below.
  integer, parameter :: c = core_rknd

  ! Input vars
  logical, optional :: &
    show_corr_arrays_input ! If true, the correlation arrays will be printed when the test is ran

  ! Output Vars
  integer :: read_corr_mtx_unit_test ! Returns the exit code of the test

  ! Local Variables
  real( kind = core_rknd ), dimension(:,:), allocatable :: &
      corr_array, & ! Retrieved correlation variance array
      test_corr_array! Expected correlation variance array

  integer :: n_row, n_col, & !indeces
    errors !Number of errors found between the two arrays

  logical :: show_corr_arrays ! the member variable of show_corr_arrays_input

  type (hm_metadata_type) :: &
      hm_metadata

  integer :: &
    pdf_dim

  !-----------------------------------------------------------------------
    !----- Begin Code -----
    print *, ""
    print *, "Performing corr_mtx_test"
    print *, ""
    print *, "=================================================="

    ! Handle the optional argument
    if (present(show_corr_arrays_input)) then
      show_corr_arrays = show_corr_arrays_input
    else
      show_corr_arrays = .false.
    end if

    ! setup the system under test
    call init_pdf_hydromet_arrays_api( one, one, 12,              & ! intent(in)
                                       5, 6, 7, 8,                & ! intent(in)
                                       9, 10, 11, 12,             & ! intent(in)
                                       one,                       & ! intent(in)
                                       hm_metadata, pdf_dim )       ! intent(out)

    if ( pdf_dim /= 12 ) then
      print *, "Error calling init_pdf_hydromet_arrays_api"
      print *, "--- pdf_dim expected: 12"
      print *, "--- pdf_dim actual:", pdf_dim
      read_corr_mtx_unit_test = 1 ! Exit Code = 1, Fail
    end if


    ! Allocate Arrays
    allocate( test_corr_array(pdf_dim,pdf_dim) )
    allocate( corr_array(pdf_dim,pdf_dim) )

    ! Initialize all values to 0.
    test_corr_array = zero
    corr_array = zero

    !setup the test array
    test_corr_array = reshape( (/ &
      1.0_c, 2.0_c, 3.0_c, 4.0_c,  5.0_c,  6.0_c,  &
      7.0_c,     8.0_c,     9.0_c,      1.0_c,     1.1_c,       1.2_c,     &! chi

      0.0_c, 1.0_c, 0.2_c, 0.3_c,  0.4_c,  0.5_c,  &
      0.6_c,     0.7_c,     0.8_c,      0.9_c,     0.10_c,      0.11_c,    &! eta

      0.0_c, 0.0_c, 1.0_c, 0.02_c, 0.03_c, 0.04_c, &
      0.05_c,    0.06_c,    0.07_c,     0.08_c,    0.09_c,      0.010_c,   &! w

      0.0_c, 0.0_c, 0.0_c, 1.0_c,  -1.0_c, -2.0_c, &
      -3.0_c,    -4.0_c,    -5.0_c,     -6.0_c,    -7.0_c,      -8.0_c,    &! Ncn

      0.0_c, 0.0_c, 0.0_c, 0.0_c,  1.0_c,  -0.1_c, &
      -0.2_c,    -0.3_c,    -0.4_c,     -0.5_c,    -0.6_c,      -0.7_c,    &! rr

      0.0_c, 0.0_c, 0.0_c, 0.0_c,  0.0_c,  1.0_c,  &
      0.00001_c, 0.00002_c, 0.000003_c, 0.00004_c, 0.0000005_c, 0.000006_c,&! Nr

      0.0_c, 0.0_c, 0.0_c, 0.0_c,  0.0_c,  0.0_c,  &
      1.0_c,     -0.0001_c, -0.00002_c, -0.0003_c, -0.000004_c, -0.00005_c,&! ri

      0.0_c, 0.0_c, 0.0_c, 0.0_c,  0.0_c,  0.0_c,  &
      0.0_c,     1.0_c,     0.0_c,      0.0_c,     0.0_c,       1.0_c,     &! Ni

      0.0_c, 0.0_c, 0.0_c, 0.0_c,  0.0_c,  0.0_c,  &
      0.0_c,     0.0_c,     1.0_c,      0.0_c,     0.0_c,       2.0_c,     &! rs

      0.0_c, 0.0_c, 0.0_c, 0.0_c,  0.0_c,  0.0_c,  &
      0.0_c,     0.0_c,     0.0_c,      1.0_c,     0.0_c,       3.0_c,     &! Ns

      0.0_c, 0.0_c, 0.0_c, 0.0_c,  0.0_c,  0.0_c,  &
      0.0_c,     0.0_c,     0.0_c,      0.0_c,     1.0_c,       4.0_c,     &! rg

      0.0_c, 0.0_c, 0.0_c, 0.0_c,  0.0_c,  0.0_c,  &
      0.0_c,     0.0_c,     0.0_c,      0.0_c,     0.0_c,       1.0_c      &! Ng
    /), shape(test_corr_array) )

    ! Read the correlation aray in the given test file
    call read_correlation_matrix( iunit, input_file,        & ! In
                                  pdf_dim, hm_metadata,  & ! In
                                  corr_array )

    ! Print out the array if required
    if (show_corr_arrays) then
      print *, "Expected:"
      call print_corr_matrix(pdf_dim, test_corr_array)
      print *, "Actual:"
      call print_corr_matrix(pdf_dim, corr_array)
    end if

    errors = 0

    ! Check for differences (between the original array and the one
    ! read in) greater than the acceptable tolerance.
    do n_col = 1, pdf_dim
      do n_row = 1, pdf_dim
        if ( abs(test_corr_array(n_row, n_col)-corr_array(n_row, n_col)) > tol ) then
              errors = errors + 1
        end if
      end do
    end do

    if (errors > 0) then
      print *, "There were ", errors, " errors found"
      read_corr_mtx_unit_test = 1 ! Exit Code = 1, Fail
    else
      print *, "Success!"
      read_corr_mtx_unit_test = 0 ! Exit Code = 0, Success!
    end if

    print *, "=================================================="

    return
  end function read_corr_mtx_unit_test
  !-----------------------------------------------------------------------

end module read_corr_mtx_test
