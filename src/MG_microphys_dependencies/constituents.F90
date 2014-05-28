! $Id$
module constituents
!
! Dummy module for importing variables into morrison_gettelman microphysics
!---------------------------------------------------------------------------------------------------

  implicit none

  private

  public :: pcnst, qmin

  ! This variable is not used anywhere in MG, it is just imported. Because of this we
  ! are setting it to a dummy value so MG will compile.
  integer :: &
    pcnst = 0

  ! Should not be used, only here so unused subroutines will compile.
  real(8) :: &
    qmin(3) = -9999.99

end module constituents
