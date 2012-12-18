! $Id$
!===============================================================================
program G_unit_tests

  ! Description:

  ! References:
  !-------------------------------------------------------------------------

  use KK_integrals_tests, only: &
      KK_integrals_tests_driver  ! Procedure(s)

  implicit none

  ! Variables
  logical, parameter :: l_KK_unit_tests = .true.


  if ( l_KK_unit_tests ) then
     call KK_integrals_tests_driver
  endif


!===============================================================================

end program G_unit_tests

