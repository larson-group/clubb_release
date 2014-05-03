!---------------------------------------------------------------------------
! $Id$
!===============================================================================
program G_unit_tests

  ! Description:

  ! References:
  !-------------------------------------------------------------------------

  use KK_integrals_tests, only: &
      KK_integrals_tests_driver  ! Procedure(s)

  use corr_cholesky_mtx_tests, only: &
      corr_cholesky_mtx_tests_driver  ! Procedure(s)

  use hole_filling_tests, only: &
      hole_filling_tests_driver ! Procedure(s)

  use Nc_Ncn_test, only: &
      Nc_Ncn_unit_test ! Procedure(s)

  implicit none

  ! Variables
  logical, parameter :: &
    l_KK_unit_tests = .true.,         &
    l_corr_cholesky_mtx_tests = .true., &
    l_hole_filling_tests = .true., &
    l_Nc_Ncn_test = .true.


  if ( l_KK_unit_tests ) then
     call KK_integrals_tests_driver
  endif

  if ( l_corr_cholesky_mtx_tests ) then
     call corr_cholesky_mtx_tests_driver
  endif

  if ( l_hole_filling_tests ) then
     call hole_filling_tests_driver
  endif

  if ( l_Nc_Ncn_test ) then
     call Nc_Ncn_unit_test
  endif


!===============================================================================

end program G_unit_tests

