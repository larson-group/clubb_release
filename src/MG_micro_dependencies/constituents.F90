! $Id$
module constituents
!
! Dummy module for importing variables into morrison-gettelman microphysics
!---------------------------------------------------------------------------------------------------

  implicit none

  private

  public :: pcnst

  ! This variable is not used anywhere in MG, it is just imported. Because of this we
  ! are setting it to a dummy value so MG will compile.
  integer :: &
    pcnst = 0

end module constituents
