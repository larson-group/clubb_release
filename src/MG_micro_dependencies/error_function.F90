! $Id$
module error_function
!
! Dummy module for importing variables into morrison-gettelman microphysics
!---------------------------------------------------------------------------------------------------

  implicit none

  private

  public :: erf,erfc
  
  ! These functions are not used anywhere in MG, they are just imported. Because of this we
  ! are setting them to dummy values
  integer :: erf = 0, erfc = 0
  
end module error_function
