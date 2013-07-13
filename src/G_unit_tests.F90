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

  implicit none

  ! Variables
  logical, parameter :: &
    l_KK_unit_tests = .false.,         &
    l_corr_cholesky_mtx_tests = .true.


  if ( l_KK_unit_tests ) then
     call KK_integrals_tests_driver
  endif


  if ( l_corr_cholesky_mtx_tests ) then
     call corr_cholesky_mtx_tests_driver
  endif

!===============================================================================

end program G_unit_tests

