!---------------------------------------------------------------------------
! $Id$
!===============================================================================
program G_unit_tests

  ! Description:

  ! References:
  !-------------------------------------------------------------------------

  use constants_clubb, only: &
      fstdout  ! Constant(s)

  use KK_integrals_tests, only: &
      KK_integrals_tests_driver  ! Procedure(s)

  use corr_cholesky_mtx_tests, only: &
      corr_cholesky_mtx_tests_driver  ! Procedure(s)

  use hole_filling_tests, only: &
      hole_filling_tests_driver ! Procedure(s)

  use Nc_Ncn_test, only: &
      Nc_Ncn_unit_test ! Procedure(s)

  use read_corr_mtx_test, only: &
      read_corr_mtx_unit_test ! Procedure(s)

  implicit none

  ! Local Constants
  integer, parameter :: iunit = 25

  ! Local Variables
  integer :: exit_code = 0 ! Assume all of the tests worked.

  ! Initialize flags
  logical :: &
    l_KK_unit_tests = .true.,           & ! Flag for KK integrals tests
    l_corr_cholesky_mtx_tests = .true., & ! Flag for corr_cholesky_mtx_tests
    l_hole_filling_tests = .true.,      & ! Flag for hole filling tests
    l_Nc_Ncn_test = .true.,             & ! Flag for Nc-Ncn Equations tests
    l_read_corr_mtx_test = .true.,      & ! Flag for corr matrix read test
      show_read_test_arrays = .true.      ! If true, the arrays used in
                                          !   read_corr_mtx_test will be shown

  ! Definition of namelist
  namelist /G_unit_namelist/ &
    l_KK_unit_tests, l_corr_cholesky_mtx_tests, l_hole_filling_tests, &
    l_Nc_Ncn_test, l_read_corr_mtx_test


  ! Read namelist file
  open(unit=iunit, file="G_unit_tests.in", status='old')
  read(unit=iunit, nml=G_unit_namelist)
  close(unit=iunit)


  write(fstdout,'(A)') "Running G_unit_tests"
  write(fstdout,'(A)') " "


  if ( l_KK_unit_tests ) then
     if (KK_integrals_tests_driver() /= 0) then
       exit_code = 1
     end if
  end if

  if ( l_corr_cholesky_mtx_tests ) then
     if (corr_cholesky_mtx_tests_driver() /= 0) then
       exit_code = 1
     end if
  end if

  if ( l_hole_filling_tests ) then
     if (hole_filling_tests_driver() /= 0) then
       exit_code = 1
     end if
  end if

  if ( l_Nc_Ncn_test ) then
     if (Nc_Ncn_unit_test() /= 0) then
       exit_code = 1
     end if
  end if

  if ( l_read_corr_mtx_test ) then
     if (read_corr_mtx_unit_test(show_read_test_arrays) /= 0) then
       exit_code = 1
     end if
  end if

  ! Stop with exit code if error found
  if (exit_code /= 0) then
    call exit(1)
  end if

!===============================================================================

end program G_unit_tests

